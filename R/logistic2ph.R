#' Sieve maximum likelihood estimator (SMLE) for two-phase logistic regression problems
#'
#' This function returns the sieve maximum likelihood estimators (SMLE) for the logistic regression model from Lotspeich et al. (2021). See pacakge vigenette for code examples.
#'
#' @param y_unval Column name of the error-prone or unvalidated binary outcome. This argument is optional. If \code{y_unval = NULL} (the default), \code{y} is treated as error-free.
#' @param y Column name that stores the validated value of \code{y_unval} in the second phase. Subjects with missing values of \code{y} are considered as those not selected in the second phase. This argument is required.
#' @param x_unval Specifies the columns of the error-prone covariates. This argument is required.
#' @param x Specifies the columns that store the validated values of \code{x_unval} in the second phase. Subjects with missing values of \code{x} are considered as those not selected in the second phase. This argument is required.
#' @param z Specifies the columns of the accurately measured covariates. This argument is optional.
#' @param data Specifies the name of the dataset. This argument is required.
#' @param hn_scale Specifies the scale of the perturbation constant in the variance estimation. For example, if \code{hn_scale = 0.5}, then the perturbation constant is \eqn{0.5n^{-1/2}}, where \eqn{n} is the first-phase sample size. The default value is \code{1}. This argument is optional.
#' @param max_iter Maximum number of iterations in the EM algorithm. The default number is \code{1000}. This argument is optional.
#' @param tol Specifies the convergence criterion in the EM algorithm. The default value is \code{1E-4}. This argument is optional.
#' @param se If \code{FALSE}, then the variances of the parameter estimators will not be estimated. The default value is \code{TRUE}. This argument is optional.
#' @param verbose If \code{TRUE}, then show details of the analysis. The default value is \code{FALSE}.
#'
#' @details
#' Models for \code{logistic2ph()} are specified through the arguments. The dataset input should at least contain columns for unvalidated error-prone outcome, validated error-prone outcome,
#' unvalidated error-prone covariate(s), validated error-prone covariate(s), and B-spline basis. B-spline basis can be generated from \code{splines::bs()} function, with argument \code{x}
#' being the unvalidated error-prone covariate(s). See vignette for options in tuning the B-spline basis.
#'
#' @return
#' `logistic2ph()` returns an object of class `"logistic2ph"`. The function `coef()` is used to obtain the coefficients of the fitted model. The function `summary()` is used to obtain and print a summary of results.
#'
#' An object of class `"logistic2ph"` is a list containing at least the following components:
#' \item{call}{the matched call.}
#' \item{coefficients}{A named vector of the logistic regression coefficient estimates.}
#' \item{covariance}{The covariance matrix of the logistic regression coefficient estimates.}
#' \item{converge}{In parameter estimation, if the EM algorithm converges, then \code{converge = TRUE}. Otherwise, \code{converge = FALSE}.}
#' \item{converge_cov}{In variance estimation, if the EM algorithm converges, then \code{converge_cov = TRUE}. Otherwise, \code{converge_cov = FALSE}.}
#'
#' @references
#' Lotspeich, S. C., Shepherd, B. E., Amorim, G. G. C., Shaw, P. A., & Tao, R. (2021). Efficient odds ratio estimation under two-phase sampling using error-prone data from a multi-national HIV research cohort. *Biometrics, biom.13512.* https://doi.org/10.1111/biom.13512
#'
#' @importFrom stats as.formula
#' @importFrom stats glm
#'
#' @export

logistic2ph <- function(y_unval = NULL, y = NULL, x_unval = NULL, x = NULL, z = NULL, data = NULL, hn_scale = 1, se = TRUE, tol = 1E-4, max_iter = 1000, verbose = FALSE) {
  # Store the function call
  model_call <- match.call()

  # variable name change
  Y_unval = y_unval; Y = y ; X_unval = x_unval; X = x; Z = z
  Bspline = attr(data, "bs_name"); noSE = !se; TOL = tol; MAX_ITER = max_iter

  if (missing(data)) {
    stop("No dataset is provided!")
  }

  if (missing(Bspline)) {
    stop("The B-spline basis is not specified!")
  }

  if (missing(Y)) {
    stop("The accurately measured response Y is not specified!")
  }

  if (xor(missing(X), missing(X_unval))) {
    stop("If X_unval and X are NULL, all predictors are assumed to be error-free. You must define both variables or neither!")
  }

  if (length(data[,X_unval]) != length(data[,X])) {
    stop("The number of columns in X_unval and X is different!")
  }

  # Convert data to data.frame to avoid conflicts with tbl
  data <- data.frame(data)

  # Save sample size
  N <- nrow(data)

  # Calculate the validated subjects
  Validated <- logical(N) # initialize a logical vector of length N
  for (i in 1:N) {
    Validated[i] <- !(is.na(data[i,X]) || is.na(data[i,Y]))
  }

  # n is how many validated subjects there are in data
  n <- sum(Validated)

  # Reorder so that the n validated subjects are first ------------
  data <- data[order(as.numeric(Validated), decreasing = TRUE), ]

  # Create vector of predictors for each model ----------------------
  theta_pred <- c(X, Z)
  gamma_pred <- c(X_unval, Y, X, Z)
  pred <- unique(c(theta_pred, gamma_pred))
  # ----------------------  Create vector of predictors for each model

  # Determine error setting -----------------------------------------
  ## If unvalidated outcome was left blank, assume error-free -------
  errorsY <- !is.null(Y_unval)
  if (!errorsY & any(is.na(data[, Y]))) {
    stop("If Y_unval is NULL, the outcome is assumed to be error-free, but Y has missing values.")
  }
  # ----------------------------------------- Determine error setting

  # Add the B spline basis ------------------------------------------
  sn <- ncol(data[, Bspline])
  if(0 %in% colSums(data[c(1:n), Bspline])) {
    warning("Empty sieve in validated data. Reconstruct B-spline basis and try again.", call. = FALSE)

    res_coefficients <- data.frame(Estimate = rep(NA, length(theta_pred) + 1),
                                   SE = NA,
                                   Statistic = NA,
                                   pvalue = NA)
    colnames(res_coefficients) <- c("Estimate", "SE", "Statistic", "p-value")

    res_final = list(coefficients = res_coefficients,
                     covariance = NA,
                     converge = FALSE,
                     converge_cov = NA)

    return(res_final)
  }
  # ------------------------------------------ Add the B spline basis

  # Create "complete" data ------------------------------------------
  ## Save distinct X ------------------------------------------------
  x_obs <- data.frame(unique(data[1:n, c(X)]))
  x_obs <- data.frame(x_obs[order(x_obs[, 1]), ])
  m <- nrow(x_obs)
  x_obs_stacked <- do.call(rbind, replicate(n = (N - n), expr = x_obs, simplify = FALSE))
  x_obs_stacked <- data.frame(x_obs_stacked[order(x_obs_stacked[, 1]), ])
  colnames(x_obs) <- colnames(x_obs_stacked) <- c(X)
  ## ------------------------------------------------ Save distinct X
  ## Save static (X*,Y*,X,Y,Z) since they don't change --------------
  comp_dat_val <- data[c(1:n), c(Y_unval, X_unval, Z, Bspline, X, Y)]
  comp_dat_val <- merge(x = comp_dat_val, y = data.frame(x_obs, k = 1:m), all.x = TRUE)
  comp_dat_val <- comp_dat_val[, c(Y_unval, pred, Bspline, "k")]
  comp_dat_val <- data.matrix(comp_dat_val)
  ## -------------- Save static (X*,Y*,X,Y,Z) since they don't change
  # 2 (m x n)xd matrices (y=0/y=1) of each (one column per person, --
  # one row per x) --------------------------------------------------
  comp_dat_unval <- cbind(data[rep(seq(from = (n+1), N), times = m), c(Y_unval, setdiff(x = pred, y = X), Bspline)],
                          x_obs_stacked,
                          k = rep(seq(1, m), each = (N - n)),
                          row.names = NULL)
  if (errorsY) {
    ### Create two versions of the complete data
    comp_dat_y0 <- comp_dat_y1 <- comp_dat_unval
    comp_dat_y0[, Y] <- 0 ### One with Y = 0
    comp_dat_y1[, Y] <- 1 ### And one with Y = 1
    comp_dat_unval <- data.matrix(rbind(comp_dat_y0, comp_dat_y1))
  }
  comp_dat_unval <- comp_dat_unval[, c(Y_unval, pred, Bspline, "k")]
  comp_dat_all <- rbind(comp_dat_val,
                        comp_dat_unval)

  # Initialize B-spline coefficients {p_kj}  -----------------------
  ## Numerators sum B(Xi*) over k = 1,...,m ------------------------
  ## Save as p_val_num for updates (contributions don't change) ----
  p_val_num <- rowsum(x = comp_dat_val[, Bspline],
                      group = comp_dat_val[, "k"],
                      reorder = TRUE)
  prev_p <- p0 <-  t(t(p_val_num) / colSums(p_val_num))

  # Create model formulas
  theta_formula <- as.formula(paste0(Y, "~", paste(theta_pred, collapse = "+")))
  theta_design_mat <- data.matrix(cbind(int = 1, comp_dat_all[, theta_pred]))
  gamma_formula <- as.formula(paste0(Y_unval, "~", paste(gamma_pred, collapse = "+")))
  gamma_design_mat <- data.matrix(cbind(int = 1, comp_dat_all[, gamma_pred]))

  # Initialize parameter values -------------------------------------
  ## theta, gamma ---------------------------------------------------
  prev_theta <- theta0 <- matrix(0, nrow = ncol(theta_design_mat), ncol = 1)
  prev_gamma <- gamma0 <- matrix(0, nrow = ncol(gamma_design_mat), ncol = 1)

  CONVERGED <- FALSE
  CONVERGED_MSG <- "Unknown"
  it <- 1

  # pre-allocate memory for loop variables
  mus_theta <- vector("numeric", nrow(theta_design_mat) * ncol(prev_theta))
  mus_gamma <- vector("numeric", nrow(gamma_design_mat) * ncol(prev_gamma))

  # Estimate theta using EM -------------------------------------------
  while(it <= MAX_ITER & !CONVERGED) {
    # E Step ----------------------------------------------------------
    ## Update the psi_kyji for unvalidated subjects -------------------
    ### P(Y|X) --------------------------------------------------------
    mu_theta <- as.numeric(theta_design_mat[-c(1:n), ] %*% prev_theta)
    pY_X <- 1 / (1 + exp(- mu_theta))
    I_y0 <- comp_dat_unval[, Y] == 0
    pY_X[I_y0] <- 1 - pY_X[I_y0]
    ### -------------------------------------------------------- P(Y|X)
    ###################################################################
    ### P(Y*|X*,Y,X) --------------------------------------------------
    if (errorsY) {
      pYstar <- 1 / (1 + exp(- as.numeric((gamma_design_mat[-c(1:n), ] %*% prev_gamma))))
      I_ystar0 <- comp_dat_unval[, Y_unval] == 0
      pYstar[I_ystar0] <- 1 - pYstar[I_ystar0]
    } else {
      pYstar <- 1
    }
    ### -------------------------------------------------- P(Y*|X*,Y,X)
    ###################################################################
    ### P(X|X*) -------------------------------------------------------
    ### p_kj ------------------------------------------------------
    ### need to reorder pX so that it's x1, ..., x1, ...., xm, ..., xm-
    ### multiply by the B-spline terms
    if (errorsY) {
      pX <- prev_p[rep(rep(seq(1, m), each = (N - n)), times = 2), ] * comp_dat_unval[, Bspline]
    } else {
      pX <- prev_p[rep(seq(1, m), each = (N - n)), ] * comp_dat_unval[, Bspline]
    }
    ### ------------------------------------------------------- P(X|X*)
    ###################################################################
    ### Estimate conditional expectations -----------------------------
    ### P(Y|X,Z)P(Y*|X*,Y,X,Z)p_kjB(X*) -----------------------------
    psi_num <- c(pY_X * pYstar) * pX
    ### Update denominator ------------------------------------------
    #### Sum up all rows per id (e.g. sum over xk/y) ----------------
    if (errorsY) {
      psi_denom <- rowsum(x = psi_num,
                          group = rep(seq(1, (N - n)),
                                      times = 2 * m))
    } else {
      psi_denom <- rowsum(x = psi_num,
                          group = rep(seq(1, (N - n)),
                                      times = m))
    }
    #### Then sum over the sn splines -------------------------------
    psi_denom <- rowSums(psi_denom)
    #### Avoid NaN resulting from dividing by 0 ---------------------
    psi_denom[psi_denom == 0] <- 1
    ### And divide them! --------------------------------------------
    psi_t <- psi_num / psi_denom
    ### Update the w_kyi for unvalidated subjects -------------------
    ### by summing across the splines/ columns of psi_t -------------
    w_t <- rowSums(psi_t)
    ### Update the u_kji for unvalidated subjects ------------------
    ### by summing over Y = 0/1 w/i each i, k ----------------------
    ### add top half of psi_t (y = 0) to bottom half (y = 1) -------
    u_t <- psi_t[c(1:(m * (N - n))), ]
    if (errorsY) {
      u_t <- u_t + psi_t[- c(1:(m * (N - n))), ]
    }
    # ---------------------------------------------------------- E Step
    ###################################################################

    ###################################################################
    # M Step ----------------------------------------------------------
    ###################################################################
    ## Update theta using weighted logistic regression ----------------
    w_t <- .lengthenWT(w_t, n)

    # calculateMu returns exp(-mu) / (1 + exp(-mu))
    muVector <- .calculateMu(theta_design_mat, prev_theta)
    gradient_theta <- .calculateGradient(w_t, n, theta_design_mat, comp_dat_all[, Y], muVector)
    hessian_theta <- .calculateHessian(theta_design_mat, w_t, muVector, n, mus_theta);
    new_theta <- tryCatch(expr = prev_theta - (solve(hessian_theta) %*% gradient_theta),
      error = function(err) {
        matrix(NA, nrow = nrow(prev_theta))
        })
    if (any(is.na(new_theta))) {
      suppressWarnings(new_theta <- matrix(glm(formula = theta_formula, family = "binomial", data = data.frame(comp_dat_all), weights = w_t)$coefficients, ncol = 1))
      # browser()
    }

    ### Check for convergence -----------------------------------------
    theta_conv <- abs(new_theta - prev_theta) < TOL

    ## --------------------------------------------------- Update theta
    ###################################################################
    # w_t is already the proper size, no need to run .lengthenWT again

    if (errorsY) {
      ## Update gamma using weighted logistic regression ----------------
      muVector <- .calculateMu(gamma_design_mat, prev_gamma)
      gradient_gamma <- .calculateGradient(w_t, n, gamma_design_mat, comp_dat_all[, c(Y_unval)], muVector)
      hessian_gamma <- .calculateHessian(gamma_design_mat, w_t, muVector, n, mus_gamma)
      new_gamma <- tryCatch(expr = prev_gamma - (solve(hessian_gamma) %*% gradient_gamma),
                            error = function(err) {
                              matrix(NA, nrow = nrow(prev_gamma))
                            })
      if (any(is.na(new_gamma))) {
        suppressWarnings(new_gamma <- matrix(glm(formula = gamma_formula, family = "binomial", data = data.frame(comp_dat_all), weights = w_t)$coefficients, ncol = 1))
      }

      # Check for convergence -----------------------------------------
      gamma_conv <- abs(new_gamma - prev_gamma) < TOL
      ## ---------------- Update gamma using weighted logistic regression
    } else {
      new_gamma <- NULL
      gamma_conv <- TRUE
    }
    ###################################################################
    ## Update {p_kj} --------------------------------------------------
    ### Update numerators by summing u_t over i = 1, ..., N ---------
    new_p_num <- p_val_num +
      rowsum(x = u_t,
             group = rep(seq(1, m), each = (N - n)),
             reorder = TRUE)
    new_p <- t(t(new_p_num) / colSums(new_p_num))
    ### Check for convergence ---------------------------------------
    p_conv <- abs(new_p - prev_p) < TOL
    ## -------------------------------------------------- Update {p_kj}

    ###################################################################
    # M Step ----------------------------------------------------------
    ###################################################################
    all_conv <- c(theta_conv, gamma_conv, p_conv)
    if (mean(all_conv) == 1) { CONVERGED <- TRUE }

    # Update values for next iteration  -------------------------------
    it <- it + 1
    prev_theta <- new_theta
    prev_gamma <- new_gamma
    prev_p <- new_p
    #  ------------------------------- Update values for next iteration
  }

  rownames(new_theta) <- c("Intercept", theta_pred)
  if (errorsY) {
    rownames(new_gamma) <- c("Intercept", gamma_pred)
  }

  if(!CONVERGED) {
    if(it > MAX_ITER) {
      CONVERGED_MSG = "MAX_ITER reached"
    }

    res_coefficients <- data.frame(Estimate = rep(NA, length(new_theta)),
                                   SE = NA,
                                   Statistic = NA,
                                   pvalue = NA)
    colnames(res_coefficients) <- c("Estimate", "SE", "Statistic", "p-value")

    if (errorsY) {
      res_final = list(coefficients = res_coefficients,
                       covariance = NA,
                       converge = FALSE,
                       converge_cov = NA)
    } else {
      res_final = list(coefficients = res_coefficients,
                       covariance = NA,
                       converge = FALSE,
                       converge_cov = NA)
    }
    return(res_final)
  }

  if(CONVERGED) { CONVERGED_MSG <- "Converged" }

  # ---------------------------------------------- Estimate theta using EM
  if(noSE) {
    res_coefficients <- data.frame(Estimate = new_theta,
                                   SE = NA,
                                   Statistic = NA,
                                   pvalue = NA)
    colnames(res_coefficients) <- c("Estimate", "SE", "Statistic", "p-value")

    if (errorsY) {
      res_final = list(coefficients = res_coefficients,
                       covariance = NA,
                       converge = CONVERGED,
                       converge_cov = NA)
    } else {
      res_final = list(coefficients = res_coefficients,
                       covariance = NA,
                       converge = CONVERGED,
                       converge_cov = NA)
    }
    return(res_final)
  } else {
    # Estimate Cov(theta) using profile likelihood -------------------------
    h_N <- hn_scale * N ^ ( - 1 / 2) # perturbation ----------------------------

    ## Calculate pl(theta) -------------------------------------------------
    od_loglik_theta <- observed_data_loglik(N = N,
      n = n,
      Y_unval = Y_unval,
      Y = Y,
      X_unval = X_unval,
      X = X,
      Z = Z,
      Bspline = Bspline,
      comp_dat_all = comp_dat_all,
      theta_pred = theta_pred,
      gamma_pred = gamma_pred,
      theta = new_theta,
      gamma = new_gamma,
      p = new_p)

    I_theta <- matrix(data = od_loglik_theta,
                      nrow = nrow(new_theta),
                      ncol = nrow(new_theta))

    single_pert_theta <- sapply(X = seq(1, ncol(I_theta)),
      FUN = pl_theta,
      theta = new_theta,
      h_N = h_N,
      n = n,
      N = N,
      Y_unval = Y_unval,
      Y = Y,
      X_unval = X_unval,
      X,
      Z = Z,
      Bspline = Bspline,
      comp_dat_all = comp_dat_all,
      theta_pred = theta_pred,
      gamma_pred = gamma_pred,
      gamma0 = new_gamma,
      p0 = new_p,
      p_val_num = p_val_num,
      TOL = TOL,
      MAX_ITER = MAX_ITER)

    if (any(is.na(single_pert_theta))) {
      I_theta <- matrix(data = NA,
                        nrow = nrow(new_theta),
                        ncol = nrow(new_theta))
      SE_CONVERGED <- FALSE
    } else {
      spt_wide <- matrix(data = rep(c(single_pert_theta),
                             times = ncol(I_theta)),
                         ncol = ncol(I_theta),
                         byrow = FALSE)
      #for the each kth row of single_pert_theta add to the kth row / kth column of I_theta
      I_theta <- I_theta - spt_wide - t(spt_wide)
      SE_CONVERGED <- TRUE
    }

    for (c in 1:ncol(I_theta)) {
      pert_theta <- new_theta
      pert_theta[c] <- pert_theta[c] + h_N
      double_pert_theta <- sapply(X = seq(c, ncol(I_theta)),
        FUN = pl_theta,
        theta = pert_theta,
        h_N = h_N,
        n = n,
        N = N,
        Y_unval = Y_unval,
        Y = Y,
        X_unval = X_unval,
        X,
        Z = Z,
        Bspline = Bspline,
        comp_dat_all = comp_dat_all,
        theta_pred = theta_pred,
        gamma_pred = gamma_pred,
        gamma0 = new_gamma,
        p0 = new_p,
        p_val_num = p_val_num,
        MAX_ITER = MAX_ITER,
        TOL = TOL)
      dpt <- matrix(data = 0,
                    nrow = nrow(I_theta),
                    ncol = ncol(I_theta))
      dpt[c,c] <- double_pert_theta[1] #Put double on the diagonal
      if(c < ncol(I_theta)) {
        ## And fill the others in on the cth row/ column
        dpt[c, -(1:c)] <- dpt[-(1:c), c] <- double_pert_theta[-1]
      }

      I_theta <- I_theta + dpt
    }
    I_theta <- h_N ^ (- 2) * I_theta
    cov_theta <- tryCatch(expr = - solve(I_theta),
      error = function(err) { matrix(NA, nrow = nrow(I_theta), ncol = ncol(I_theta)) }
      )
    # ------------------------- Estimate Cov(theta) using profile likelihood

    se_theta <- tryCatch(expr = sqrt(diag(cov_theta)),
      warning = function(w) { matrix(NA, nrow = nrow(prev_theta)) })
    if (any(is.na(se_theta))) {
      SE_CONVERGED <- FALSE
    } else {
      SE_CONVERGED <- TRUE
    }

    if (verbose) {
      message(CONVERGED_MSG)
    }

    res_coefficients <- data.frame(Estimate = new_theta, SE = se_theta)
    res_coefficients$Statistic <- res_coefficients$Estimate / res_coefficients$SE
    res_coefficients$pvalue <- 1 - pchisq(res_coefficients$Statistic ^ 2, df = 1)
    colnames(res_coefficients) <- c("Estimate", "SE", "Statistic", "p-value")

    coef_return <- as.numeric(new_theta)
    names(coef_return) <- rownames(res_coefficients)

    res_final = list(
                    call = model_call,  # Store the call in the object
                    coefficients = coef_return,
                    covariance = cov_theta,
                    converge = CONVERGED,
                    converge_cov = SE_CONVERGED)

    res_final <- logistic2ph_class(res_final)

    return(res_final)
  }
}
