#' HSA (Hydrological Signature Analysis)
#' 
#' This package is a set of function to compute hydrological signatures
#' from streamflow, precipitation and air temperature time series.
#' 
#' @docType package
#' @author Ivan Horner
#' @import dplyr, segmented, RcppRoll, Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib HSA
#' @name HSA
#' @export
NULL

# instruction for preparing and installing package:
# delete NAMESPACE file
# library(roxygen2)
# roxygenize()
# edit NAMESPACE file: remove 'export HSA' line and remove commas after external dependency names
# library(devtools)
# C:\Users\ihach\Documents\Programming\R\hsa-R-package
# remove.packages("HSA")
# install_local(file.path("C:","Users", "ihach", "Documents", "Programming", "R", "hsa-R-package"))