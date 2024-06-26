#' Profiles out nuisance parameters from the observed-data log-likelihood for a given value of theta
#'
#' For a given vector `theta` to parameterize P(Y|X,Z), this function repeats the EM algorithm to find
#' the values of `gamma` and `p` at convergence. The resulting parameters are used to find the profile
#' log-likelihood for `theta` by plugging them into the observed-data log-likelihood.
#' This function is used by `pl_theta()`.
#
#' @param theta Parameters for the analysis model (a column vector)
#' @param N Phase I sample size
#' @param n Phase II sample size
#' @param Y_unval Column with the unvalidated outcome (can be name or numeric index)
#' @param Y Column with the validated outcome (can be name or numeric index)
#' @param X_unval Column(s) with the unvalidated predictors (can be name or numeric index)
#' @param X Column(s) with the validated predictors (can be name or numeric index)
#' @param Z (Optional) Column(s) with additional error-free covariates (can be name or numeric index)
#' @param Bspline Vector of columns containing the B-spline basis functions (can be name or numeric index)
#' @param comp_dat_all Augmented dataset containing rows for each combination of unvalidated subjects' data with values from Phase II (a matrix)
#' @param theta_pred Vector of columns in \code{data} that pertain to the predictors in the analysis model.
#' @param gamma_pred Vector of columns in \code{data} that pertain to the predictors in the outcome error model.
#' @param gamma0 Starting values for `gamma`, the parameters for the outcome error model (a column vector)
#' @param p0 Starting values for `p`, the B-spline coefficients for the approximated covariate error model (a matrix)
#' @param p_val_num Contributions of validated subjects to the numerator for `p`, which are fixed (a matrix)
#' @param TOL Tolerance between iterations in the EM algorithm used to define convergence.
#' @param MAX_ITER Maximum number of iterations allowed in the EM algorithm.
#'
#' @return Profile likelihood for `theta`: the value of the observed-data log-likelihood after profiling out other parameters.
#'
#' @importFrom stats as.formula
#' @importFrom stats glm
#'
#' @noRd

profile_out <- function(theta, n, N, Y_unval = NULL, Y = NULL, X_unval = NULL, X = NULL, Z = NULL, Bspline = NULL,
                           comp_dat_all, theta_pred, gamma_pred, gamma0, p0, p_val_num, TOL, MAX_ITER) {
  # Determine error setting -----------------------------------------
  ## If unvalidated variable was left blank, assume error-free ------
  errorsY <- errorsX <- TRUE
  if (is.null(Y_unval)) {errorsY <- FALSE}
  if (is.null(X_unval)) {errorsX <- FALSE}
  # ----------------------------------------- Determine error setting

  sn <- ncol(p0)
  m <- nrow(p0)

  prev_gamma <- gamma0
  prev_p <- p0

  # Convert to matrices
  theta_design_mat <- as.matrix(cbind(int = 1,
                                      comp_dat_all[-c(1:n), theta_pred]))
  comp_dat_all <- as.matrix(comp_dat_all)

  # For the E-step, save static P(Y|X) for unvalidated --------------
  # browser()
  pY_X <- .pYstarCalc(theta_design_mat, 1L, theta, comp_dat_all, match(Y, colnames(comp_dat_all))-1)
  if (errorsY) {
    gamma_formula <- as.formula(paste0(Y_unval, "~", paste(gamma_pred, collapse = "+")))
    gamma_design_mat <- as.matrix(cbind(int = 1,
                                        comp_dat_all[, gamma_pred]))
  }

  CONVERGED <- FALSE
  CONVERGED_MSG <- "Unknown"
  it <- 1

  # pre-allocate memory for our loop variables
  # improves performance
  if (errorsY) {
    pYstar <- vector(mode="numeric", length = nrow(gamma_design_mat) - n)
    mu_gamma <- vector(mode="numeric", length = length(pYstar))
    mus_gamma <- vector("numeric", nrow(gamma_design_mat) * ncol(prev_gamma))
    pX <- matrix(,nrow = m * (N-n) * 2, ncol = length(Bspline))
  } else {
    pX <- matrix(,nrow = m * (N-n), ncol = length(Bspline))
  }
  psi_num <- matrix(,nrow = nrow(comp_dat_all)-n, ncol=length(Bspline))
  psi_t <- matrix(,nrow = nrow(psi_num), ncol = ncol(psi_num))
  w_t <- vector("numeric", length = nrow(psi_t))
  u_t <- matrix(,nrow = m * (N-n), ncol = ncol(psi_t))

  # Estimate gamma/p using EM -----------------------------------------
  #browser()
  while(it <= MAX_ITER & !CONVERGED) {
    # E Step ----------------------------------------------------------
    # P(Y*|X*,Y,X) ---------------------------------------------------
    if (errorsY) {
      pYstar <- .pYstarCalc(gamma_design_mat, n + 1, prev_gamma, comp_dat_all, match(Y_unval, colnames(comp_dat_all))-1)
    } else {
      pYstar <- 1
    }
    ### -------------------------------------------------- P(Y*|X*,Y,X)
    ###################################################################
    ### P(X|X*) -------------------------------------------------------

    # these are the slowest lines in this function (13280 ms), but I can't seem to get any faster with C++
    # pX <- pXCalc(n, comp_dat_all[-c(1:n), Bspline], errorsX, errorsY, pX, prev_p[rep(seq(1, m), each = (N - n)), ])
    ### p_kj ------------------------------------------------------
    ### need to reorder pX so that it's x1, ..., x1, ...., xm, ..., xm-
    ### multiply by the B-spline terms
    if (errorsY) {
      pX <- prev_p[rep(rep(seq(1, m), each = (N - n)), times = 2), ] * comp_dat_all[-c(1:n), Bspline]
    } else {
      pX <- prev_p[rep(seq(1, m), each = (N - n)), ] * comp_dat_all[-c(1:n), Bspline]
    }
    ### ---------------------------------------------------------- p_kj
    ### ------------------------------------------------------- P(X|X*)
    ###################################################################
    ### Estimate conditional expectations -----------------------------


    # this function modifies w_t, u_t, psi_num, and psi_t internally
    # condExp <- conditionalExpectations(errorsX, errorsY, pX, pY_X, pYstar, N - n, m, w_t, u_t, psi_num, psi_t)
    # w_t <- condExp[["w_t"]]
    # u_t <- condExp[["u_t"]]
    # psi_t <- condExp[["psi_t"]]

    ### P(Y|X,Z)P(Y*|X*,Y,X,Z)p_kjB(X*) -----------------------------
    psi_num <- c(pY_X * pYstar) * pX
    ### Update denominator ------------------------------------------
    #### Sum up all rows per id (e.g. sum over xk/y) ----------------
    if (errorsY) {
      psi_denom <- rowsum(psi_num, group = rep(seq(1, (N - n)), times = 2 * m))
    } else {
      psi_denom <- rowsum(psi_num, group = rep(seq(1, (N - n)), times = m))
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
    ### ----------------------------- Estimate conditional expectations
    # ---------------------------------------------------------- E Step
    ###################################################################

    ###################################################################
    # M Step ----------------------------------------------------------
    ###################################################################
    ## Update gamma using weighted logistic regression ----------------
    if (errorsY) {
      w_t <- .lengthenWT(w_t, n)
      muVector <- .calculateMu(gamma_design_mat, prev_gamma)
      gradient_gamma <- .calculateGradient(w_t, n, gamma_design_mat, comp_dat_all[,c(Y_unval)], muVector)
      hessian_gamma <- .calculateHessian(gamma_design_mat, w_t, muVector, n, mus_gamma)

      new_gamma <- tryCatch(expr = prev_gamma - (solve(hessian_gamma) %*% gradient_gamma),
                            error = function(err) {
                              matrix(NA, nrow = nrow(prev_gamma))
                            })
      if (any(is.na(new_gamma))) {
        suppressWarnings(new_gamma <- matrix(glm(formula = gamma_formula,
                                                 family = "binomial",
                                                 data = data.frame(comp_dat_all),
                                                 weights = w_t)$coefficients,
                                             ncol = 1))
        # browser()
      }

      # Check for convergence -----------------------------------------
      gamma_conv <- abs(new_gamma - prev_gamma) < TOL
    } else {
      new_gamma <- NULL
      gamma_conv <- TRUE
    }
    ## ---------------- Update gamma using weighted logistic regression
    ###################################################################
    ## Update {p_kj} --------------------------------------------------
    ### Update numerators by summing u_t over i = 1, ..., N ---------
    new_p_num <- p_val_num +
      rowsum(u_t, group = rep(seq(1, m), each = (N - n)), reorder = TRUE)
    new_p <- t(t(new_p_num) / colSums(new_p_num))
    ### Check for convergence ---------------------------------------
    p_conv <- abs(new_p - prev_p) < TOL
    ## -------------------------------------------------- Update {p_kj}
    ###################################################################
    # M Step ----------------------------------------------------------
    ###################################################################

    all_conv <- c(gamma_conv, p_conv)
    if (mean(all_conv) == 1) { CONVERGED <- TRUE }

    # Update values for next iteration  -------------------------------
    it <- it + 1
    prev_gamma <- new_gamma
    prev_p <- new_p
    #  ------------------------------- Update values for next iteration

  }

  if(it == MAX_ITER & !CONVERGED) {
    CONVERGED_MSG <- "MAX_ITER reached"
    if (errorsY) {
      new_gamma <- matrix(data = NA,
                          nrow = nrow(gamma0),
                          ncol = 1)
    } else {
      new_gamma <- NULL
    }
    new_p <- matrix(data = NA,
                    nrow = nrow(p0),
                    ncol = ncol(p0))
  }
  if(CONVERGED) CONVERGED_MSG <- "converged"
  # ---------------------------------------------- Estimate theta using EM
  return(list("psi_at_conv" = psi_t,
              "gamma_at_conv" = new_gamma,
              "p_at_conv" = new_p,
              "converged" = CONVERGED,
              "converged_msg" = CONVERGED_MSG))
}
