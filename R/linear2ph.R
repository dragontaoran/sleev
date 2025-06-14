#' Sieve maximum likelihood estimator (SMLE) for two-phase linear regression problems
#'
#' Performs efficient semiparametric estimation for general two-phase measurement error models when there are errors in both the outcome and covariates. See pacakge vigenette for code examples.
#'
#' @param y_unval Column name of the error-prone or unvalidated continuous outcome. Subjects with missing values of \code{y_unval} are omitted from the analysis. This argument is required.
#' @param y Column name that stores the validated value of \code{y_unval} in the second phase. Subjects with missing values of \code{y} are considered as those not selected in the second phase. This argument is required.
#' @param x_unval Specifies the columns of the error-prone covariates. Subjects with missing values of \code{x_unval} are omitted from the analysis. This argument is required.
#' @param x Specifies the columns that store the validated values of \code{x_unval} in the second phase. Subjects with missing values of \code{x} are considered as those not selected in the second phase. This argument is required.
#' @param z Specifies the columns of the accurately measured covariates. Subjects with missing values of \code{z} are omitted from the analysis. This argument is optional.
#' @param data Specifies the name of the dataset. This argument is required.
#' @param hn_scale Specifies the scale of the perturbation constant in the variance estimation. For example, if \code{hn_scale = 0.5}, then the perturbation constant is \eqn{0.5n^{-1/2}}, where \eqn{n} is the first-phase sample size. The default value is \code{1}. This argument is optional.
#' @param max_iter Maximum number of iterations in the EM algorithm. The default number is \code{1000}. This argument is optional.
#' @param tol Specifies the convergence criterion in the EM algorithm. The default value is \code{1E-4}. This argument is optional.
#' @param se If \code{FALSE}, then the variances of the parameter estimators will not be estimated. The default value is \code{TRUE}. This argument is optional.
#' @param verbose If \code{TRUE}, then show details of the analysis. The default value is \code{FALSE}.
#'
#' @details
#' Models for \code{linear2ph()} are specified through the arguments. The dataset input should at least contain columns for unvalidated error-prone outcome, validated error-prone outcome,
#' unvalidated error-prone covariate(s), validated error-prone covariate(s), and B-spline basis. B-spline basis can be generated from \code{splines::bs()} function, with argument \code{x}
#' being the unvalidated error-prone covariate(s). See vignette for options in tuning the B-spline basis.
#'
#' @return
#' `linear2ph()` returns an object of class `"linear2ph"`. The function `coef()` is used to obtain the coefficients of the fitted model. The function `summary()` is used to obtain and print a summary of results.
#'
#' An object of class `"linear2ph"` is a list containing at least the following components:
#' \item{call}{the matched call.}
#' \item{coefficients}{A named vector of the linear regression coefficient estimates.}
#' \item{sigma}{The residual standard error.}
#' \item{covariance}{The covariance matrix of the linear regression coefficient estimates.}
#' \item{converge}{In parameter estimation, if the EM algorithm converges, then \code{converge = TRUE}. Otherwise, \code{converge = FALSE}.}
#' \item{converge_cov}{In variance estimation, if the EM algorithm converges, then \code{converge_cov = TRUE}. Otherwise, \code{converge_cov = FALSE}.}
#'
#' @importFrom Rcpp evalCpp
#' @importFrom stats pchisq
#'
#' @references
#' Tao, R., Mercaldo, N. D., Haneuse, S., Maronge, J. M., Rathouz, P. J., Heagerty, P. J., & Schildcrout, J. S. (2021). Two-wave two-phase outcome-dependent sampling designs, with applications to longitudinal binary data. *Statistics in Medicine, 40*(8), 1863–1876. https://doi.org/10.1002/sim.8876
#'
#' @seealso [cv_linear2ph()] to calculate the average predicted log likelihood of this function.
#' @export
#'
linear2ph <- function (y_unval=NULL, y=NULL, x_unval=NULL, x=NULL, z=NULL, data=NULL, hn_scale=1, se=TRUE, tol=1E-4, max_iter=1000, verbose=FALSE)
{
  # Store the function call
  model_call <- match.call()

  # variable name change
  Y_unval = y_unval; Y = y ; X_unval = x_unval; X = x; Z = z
  Bspline =  attr(data, "bs_name"); noSE = !se; TOL = tol; MAX_ITER = max_iter
### linear2ph
    ###############################################################################################################
    #### check data ###############################################################################################
    storage.mode(MAX_ITER) = "integer"
	storage.mode(TOL) = "double"
	storage.mode(noSE) = "integer"

	if (missing(data)) {
	    stop("No dataset is provided!")
	}

	if (missing(Y_unval)) {
		stop("The error-prone response Y_unval is not specified!")
	} else {
		vars_ph1 = Y_unval
	}

	if (missing(X_unval)) {
		stop("The error-prone covariates X_unval is not specified!")
	} else {
		vars_ph1 = c(vars_ph1, X_unval)
	}

	if (missing(Bspline)) {
	    stop("The B-spline basis is not specified!")
	} else {
	    vars_ph1 = c(vars_ph1, Bspline)
	}

	if (missing(Y)) {
		stop("The accurately measured response Y is not specified!")
	}

	if (missing(X)) {
		stop("The validated covariates in the second-phase are not specified!")
	}

	if (length(X_unval) != length(X)) {
	    stop("The number of columns in X_unval and X is different!")
	}

	if (!missing(Z)) {
		vars_ph1 = c(vars_ph1, Z)
	}

    id_exclude = c()
    for (var in vars_ph1) {
        id_exclude = union(id_exclude, which(is.na(data[,var])))
    }

	if (verbose) {
    	message("There are ", nrow(data), " observations in the dataset.")
    	message(length(id_exclude), " observations are excluded due to missing Y_unval, X_unval, or Z.")
	}
	if (length(id_exclude) > 0) {
		data = data[-id_exclude,]
	}

    n = nrow(data)
	if (verbose) {
    	message("There are ", n, " observations in the analysis.")
	}

    id_phase1 = which(is.na(data[,Y]))
    for (var in X) {
        id_phase1 = union(id_phase1, which(is.na(data[,var])))
    }
	if (verbose) {
		message("There are ", n-length(id_phase1), " observations validated in the second phase.")
	}
    #### check data ###############################################################################################
	###############################################################################################################



	###############################################################################################################
	#### prepare analysis #########################################################################################
    Y_unval_vec = c(as.vector(data[-id_phase1,Y_unval]), as.vector(data[id_phase1,Y_unval]))
	storage.mode(Y_unval_vec) = "double"

	X_tilde_mat = rbind(as.matrix(data[-id_phase1,X_unval]), as.matrix(data[id_phase1,X_unval]))
	storage.mode(X_tilde_mat) = "double"

	Bspline_mat = rbind(as.matrix(data[-id_phase1,Bspline]), as.matrix(data[id_phase1,Bspline]))
	storage.mode(Bspline_mat) = "double"

	Y_vec = as.vector(data[-id_phase1,Y])
	storage.mode(Y_vec) = "double"

    X_mat = as.matrix(data[-id_phase1,X])
	storage.mode(X_mat) = "double"

    if (!is.null(Z)) {
        Z_mat = rbind(as.matrix(data[-id_phase1,Z]), as.matrix(data[id_phase1,Z]))
		storage.mode(Z_mat) = "double"
    }

    cov_names = c("Intercept", X)
	if (!is.null(Z)) {
		cov_names = c(cov_names, Z)
	}

	ncov = length(cov_names)
	X_nc = length(X)
	rowmap = rep(NA, ncov)
	res_coefficients = matrix(NA, nrow=ncov, ncol=4)
	colnames(res_coefficients) = c("Estimate", "SE", "Statistic", "p-value")
	rownames(res_coefficients) = cov_names
	res_cov = matrix(NA, nrow=ncov, ncol=ncov)
	colnames(res_cov) = cov_names
	rownames(res_cov) = cov_names

	if (is.null(Z)) {
		Z_mat = as.matrix(rep(1., n))
		rowmap[1] = ncov
		rowmap[2:ncov] = 1:X_nc
	} else {
		Z_mat = cbind(1, Z_mat)
		rowmap[1] = X_nc+1
		rowmap[2:(X_nc+1)] = 1:X_nc
		rowmap[(X_nc+2):ncov] = (X_nc+2):ncov
	}

	hn = hn_scale/sqrt(n)
	#### prepare analysis #########################################################################################
	###############################################################################################################

	if (verbose)
	{
		message("Calling C++ function TwoPhase_MLE0_MEXY")
	}

	## Ensure every variable is the correct type
	if (!is.vector(Y_unval_vec))
	{
		warning("Y_unval_vec is not a vector!")
	}
	if (!is.matrix(X_tilde_mat))
	{
		warning("X_tilde_mat is not a matrix!")
	}
	if (!is.vector(Y_vec))
	{
		warning("Y_vec is not a vector!")
	}
	if (!is.matrix(X_mat))
	{
		warning("X_mat is not a matrix!")
	}
	if (!is.matrix(Z_mat))
	{
		warning("Z_mat is not a matrix!")
	}
	if (!is.matrix(Bspline_mat))
	{
		warning("Bspline_mat is not a matrix!")
	}

	###############################################################################################################
	#### analysis #################################################################################################
	res = .TwoPhase_MLE0_MEXY(Y_unval_vec, X_tilde_mat, Y_vec, X_mat, Z_mat, Bspline_mat, hn, MAX_ITER, TOL, noSE)
    #### analysis #################################################################################################
	###############################################################################################################



    ###############################################################################################################
    #### return results ###########################################################################################
 	res_coefficients[,1] = res$theta[rowmap]
	res_coefficients[which(res_coefficients[,1] == -999),1] = NA
	res_cov = res$cov_theta[rowmap, rowmap]
	res_cov[which(res_cov == -999)] = NA
	diag(res_cov)[which(diag(res_cov) < 0)] = NA

	res_coefficients[,2] = diag(res_cov)
	res_coefficients[which(res_coefficients[,2] > 0),2] = sqrt(res_coefficients[which(res_coefficients[,2] > 0),2])

	id_NA = which(is.na(res_coefficients[,1]) | is.na(res_coefficients[,2]))
	if (length(id_NA) > 0)
	{
	    res_coefficients[-id_NA,3] = res_coefficients[-id_NA,1]/res_coefficients[-id_NA,2]
	    res_coefficients[-id_NA,4] = 1-pchisq(res_coefficients[-id_NA,3]^2, df=1)
	}
	else
	{
	    res_coefficients[,3] = res_coefficients[,1]/res_coefficients[,2]
	    res_coefficients[,4] = 1-pchisq(res_coefficients[,3]^2, df=1)
	}

	res_final = list(
	  call = model_call,  # Store the call in the object
	  coefficients=res_coefficients[,1],
		sigma=sqrt(res$sigma_sq),
		covariance=res_cov,
		converge=!res$flag_nonconvergence,
		converge_cov=!res$flag_nonconvergence_cov)


	res_final <- linear2ph_class(res_final)

	return(res_final)
    #### return results ###########################################################################################
    ###############################################################################################################
}
