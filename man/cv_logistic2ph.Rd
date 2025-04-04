% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv_logistic2ph.R
\name{cv_logistic2ph}
\alias{cv_logistic2ph}
\title{Cross-validation log-likelihood prediction for \code{logistic2ph}}
\usage{
cv_logistic2ph(
  y_unval = NULL,
  y = NULL,
  x_unval = NULL,
  x = NULL,
  z = NULL,
  b_spline = NULL,
  data,
  nfolds = 5,
  tol = 1e-04,
  max_iter = 1000,
  verbose = FALSE
)
}
\arguments{
\item{y_unval}{Column name of the error-prone or unvalidated binary outcome. This argument is optional. If \code{y_unval = NULL} (the default), \code{y} is treated as error-free.}

\item{y}{Column name that stores the validated value of \code{y_unval} in the second phase. Subjects with missing values of \code{y} are considered as those not selected in the second phase. This argument is required.}

\item{x_unval}{Specifies the columns of the error-prone covariates. This argument is required.}

\item{x}{Specifies the columns that store the validated values of \code{x_unval} in the second phase. Subjects with missing values of \code{x} are considered as those not selected in the second phase. This argument is required.}

\item{z}{Specifies the columns of the accurately measured covariates. Subjects with missing values of \code{z} are omitted from the analysis. This argument is optional.}

\item{b_spline}{Specifies the columns of the B-spline basis. Subjects with missing values of \code{b_spline} are omitted from the analysis. This argument is required.}

\item{data}{Specifies the name of the dataset. This argument is required.}

\item{nfolds}{Specifies the number of cross-validation folds. The default value is \code{5}. Although \code{nfolds} can be as large as the sample size (leave-one-out cross-validation), it is not recommended for large datasets. The smallest value allowable is \code{3}.}

\item{tol}{Specifies the convergence criterion in the EM algorithm. The default value is \code{1E-4}. This argument is optional.}

\item{max_iter}{Specifies the maximum number of iterations in the EM algorithm. The default number is \code{2000}. This argument is optional.}

\item{verbose}{If \code{TRUE}, then show details of the analysis. The default value is \code{FALSE}.}
}
\value{
\code{cv_logistic2ph()} returns a list that includes the following components:
\item{avg_pred_loglike}{Stores the average predicted log likelihood.}
\item{pred_loglike}{Stores the predicted log likelihood in each fold.}
\item{converge}{Stores the convergence status of the EM algorithm in each run.}
}
\description{
Performs cross-validation to calculate the average predicted log likelihood for the \code{logistic2ph} function. This function can be used to select the B-spline basis that yields the largest average predicted log likelihood.
}
\details{
\code{cv_logistic2ph} gives log-likelihood prediction for models and data like those in \code{logistic2ph}. Therefore, the arguments of \code{cv_logistic2ph} is analogous to that of \code{logistic2ph}.
}
\examples{
set.seed(1)
# different B-spline sizes
sns <- c(15, 20, 25, 30, 35, 40)
# vector to hold mean log-likelihood and run time for each sn
pred_loglike.1 <- run.time.secs <- rep(NA, length(sns))
# get number of rows of the dataset
n <- nrow(mock.vccc)
# specify number of folds in the cross validation
k <- 5
# calculate proportion of female patients in the dataset
sex_ratio <- sum(mock.vccc$Sex==0)/n
for (i in 1:length(sns)) {
  # constructing B-spline basis using the same process as in Section 4.3.1
  print(i)
  sn <- sns[i]
  art_ratio <- sum(mock.vccc$prior_ART==0)/n
  sn_0 <- round(sn*art_ratio)
  sn_1 <- sn-sn_0
  Bspline_0 <- splines::bs(x=mock.vccc$VL_unval_l10[mock.vccc$prior_ART==0],
                           df=sn_0, degree=3, intercept=TRUE)
  Bspline_1 <- splines::bs(x=mock.vccc$VL_unval_l10[mock.vccc$prior_ART==1],
                           df=sn_1, degree=3, intercept=TRUE)
  Bspline <- matrix(0, n, sn)
  Bspline[mock.vccc$prior_ART==0,1:sn_0] <- Bspline_0
  Bspline[mock.vccc$prior_ART==1,(sn_0+1):sn] <- Bspline_1
  colnames(Bspline) <- paste("bs", 1:sn, sep="")
  data.sieve <- data.frame(cbind(mock.vccc, Bspline))
  # cross validation, produce mean log-likelihood
  start.time <- Sys.time()
  res.1 <- cv_logistic2ph(y="ADE_val", y_unval="ADE_unval",
                          x="VL_val_l10", x_unval="VL_unval_l10",
                          z="prior_ART", b_spline=colnames(Bspline),
                          data=data.sieve,
                          nfolds=5, max_iter = 2000, tol = 1e-04)
  # save run time
  run.time.secs[i] <- difftime(Sys.time(), start.time, units = "secs")
  # save mean log-likelihood result
  pred_loglike.1[i] <- res.1$avg_pred_loglik
}

}
