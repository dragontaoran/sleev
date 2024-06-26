# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @noRd
.TwoPhase_MLE0_MEXY_CV_loglik <- function(Y_tilde, X_tilde, Y, X, Z, Bspline, MAX_ITER, TOL, Train) {
    .Call(`_sleev_TwoPhase_MLE0_MEXY_CV_loglik`, Y_tilde, X_tilde, Y, X, Z, Bspline, MAX_ITER, TOL, Train)
}

#' Two Phase MLE0 MEXY
#' 
#' TODO
#' 
#' @param Y_unval Unvalidated Y variables
#' @param X_unval Unvalidated X variables
#' @param Y Validated Y variables
#' @param X Validated X variables
#' @param Z True covariates
#' @param Bspline Matrix of B splines
#' @param hn Scaling of hn
#' @param MAX_ITER Max iterations to perform when calculating convergence
#' @param TOL Maximum difference between iteration that satisfies convergence requirements
#' @param noSE Skips general spline profiling if converged
#' @noRd
.TwoPhase_MLE0_MEXY <- function(Y_unval, X_unval, Y, X, Z, Bspline, hn, MAX_ITER, TOL, noSE) {
    .Call(`_sleev_TwoPhase_MLE0_MEXY`, Y_unval, X_unval, Y, X, Z, Bspline, hn, MAX_ITER, TOL, noSE)
}

#' Multiply matrix by vector
#'
#' Multiplies each column of a matrix by a vector
#'
#' @param mat The matrix
#' @param v The vector
#' @return mat * v
#' @noRd
NULL

#' Divide matrix by vector
#'
#' Divides each column of a matrix by a vector
#'
#' @param mat The matrix
#' @param v The vector divisor
#' @return mat / v, or mat * v^-1
#' @noRd
NULL

#' Prepend ones to a w_t
#'
#' Lengthens a vector by prepending n ones
#'
#' @param w_t_original The original vector
#' @param n The number of ones to add to the front of the vector
#' @param modifyW_T If false, instantly returns w_t_original
#' @noRd
.lengthenWT <- function(w_t_original, n, modifyW_T = TRUE) {
    .Call(`_sleev_lengthenWT`, w_t_original, n, modifyW_T)
}

#' Calculate Mu
#'
#' Calculates the value of mu according to two variables
#' Small helper function
#'
#' @param design_mat The design matrix
#' @param prev The previous iteration of the design matrix
#' @noRd
.calculateMu <- function(design_mat, prev) {
    .Call(`_sleev_calculateMu`, design_mat, prev)
}

#' Calculate gradient
#'
#' Calculates a gradient given w_t and a design matrix
#' TODO
#'
#' @param w_t A vector indicating ??
#' @param n The number of ones to prepend to w_t
#' @param design_mat The design matrix
#' @param Y_col The column of validated Y values from the complete data matrix
#' @param muVector The vector calculated by calculateMu
#' @param modifyW_T Whether to add ones to the beginning of w_t
#' @noRd
.calculateGradient <- function(w_t, n, design_mat, Y_col, muVector, modifyW_T = FALSE) {
    .Call(`_sleev_calculateGradient`, w_t, n, design_mat, Y_col, muVector, modifyW_T)
}

#' Calculate Hessian Matrix
#'
#' Calculates the Hessian Matrix and lengthens w_t by n
#'
#' @param design_mat The design matrix
#' @param w_t The vector ??
#' @param muVector The vector returned by calculateMu
#' @param n The number of ones to prepend to w_t
#' @param mus An empty, pre-allocated vector of the same length as muVector, pre-allocated memory saves time
#' @param modifyW_T Whether to add ones to the beginning of w_t
#' @noRd
.calculateHessian <- function(design_mat, w_t, muVector, n, mus, modifyW_T = FALSE) {
    .Call(`_sleev_calculateHessian`, design_mat, w_t, muVector, n, mus, modifyW_T)
}

#' Calculate pYstar
#'
#' TODO
#'
#' @param gamma_design_mat The gamma design matrix
#' @param startRow The number of the rows to start within gamma_design_mat
#' @param prev_gamma The previous iteration of gamma_design_mat
#' @param comp_dat_all The complete dataset
#' @param Y_unval_index Which column of comp_dat_all houses the unvalidated Y variable
#' @noRd
.pYstarCalc <- function(gamma_design_mat, startRow, prev_gamma, comp_dat_all, Y_unval_index) {
    .Call(`_sleev_pYstarCalc`, gamma_design_mat, startRow, prev_gamma, comp_dat_all, Y_unval_index)
}

