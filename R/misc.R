# TODO: implement an option to verify the 'validity' of the interpolation
#       For example, it should not be longer than 'n' time step or I could
#       provide a value indicating how much the interpolated value weight
#       on the resulting cumulative curve
#------------------------------------------------------------------------------
#' Wrapper around the \code{cumsum} function where missing values can be 
#' interpolated.
#' 
#' @param x vector to use to compute cumulated values
#' @param na.action action to undertake when missing values are found (default:
#' 'interpolate', any other value leads to the original \code{cumsum} function)
#' @return the cumulative sum of x.
#' @export
cumsum_interpolate <- function(x, na.action = "interpolate") {
  isna <- is.na(x)
  if (all(isna)) stop("'x' cannot be only missing values!")
  if (any(isna)) {
    if (na.action == "interpolate") {
      y <- 1:length(x)
      x[isna] <- approx(x = y[!isna], y = x[!isna], xout = y[isna], method = "linear")$y
    } else if (is.numeric(na.action)) {
      x[isna] <- na.action
    }
    isna <- is.na(x)
    if (any(isna)) {
      warning("'x' has missing values at the start/end! They were treated as zeros.")
      x[isna] <- 0
      y <- cumsum(x)
      y[isna] <- NA
      return(y)
    }
  }
  cumsum(x)
}



# This function is here only to avoid code repetition (BFI, RC, etc.)
#------------------------------------------------------------------------------
#' Compute the ratio of two sums.
#' 
#' Compute the ratio between the sum of 'top' and the sum of 'bottom'.
#' 
#' @param top numeric vector. Its sum will be the nominator of the ratio
#' @param bottom numeric vector. Its sum will be the denominator of the ratio
#' @param na.rm logical. Should missing values be omited?
#' @return A single numeric value corresponding to the ratio of the sum of
#' 'top' and the sum of 'bottom'
#' @export
sum_ratio <- function(top, bottom, na.rm = TRUE) sum(top, na.rm = na.rm) / sum(bottom, na.rm = na.rm)
