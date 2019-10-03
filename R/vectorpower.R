#------------------------------------------------
#' @title Power Calculation for Cluster Randomised Vector Control Trials
#'
#' @description Allows the user to estimate the sample size required for a given
#'   intervention efficacny in order to detect a significant effect, using data
#'   simulated from a cluster randomised trial in a malaria transmission model.
#'
#' @docType package
#' @name vectorpower
NULL

#------------------------------------------------
#' @useDynLib vectorpower, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#------------------------------------------------
# unload DLL when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("vectorpower", libpath)
}
NULL
