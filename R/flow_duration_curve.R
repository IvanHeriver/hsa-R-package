#------------------------------------------------------------------------------
#' Compute the slope of the mid-segment slope of the flow duration curve
#' 
#' Given two probability of streamflow exceedance (probs) which define the
#' start and end of the mid-segment of the flow duration curve, compute the 
#' slope of the mid-segment slope of the flow duration curve.
#' 
#' The computation of the slope is done with the following function:
#' \code{- (log10(Q_high) - log10(Q_low)) / diff(probs)}
#' where \code{log10(Q_low)} (resp. \code{log10(Q_high)}) is the log transformed 
#' values of streamflow corresponding to the low (resp. high) end of the 
#' mid-segment of the flow duration curve
#' 
#' @param Q numeric vector. Streamflow values.
#' @param probs numeric vector of length 2 containing the exceedance
#' probabilities that define the start and end of the mid-segment of the 
#' flow duration curve.
#' @return A single numeric value corresponding to the slope of the mid-segment
#' of the flow duration curve.
#' @seealso \code{\link{fdc_percentiles}} and \code{\link{fdc_slope_percentiles}}
#' @export
fdc_slope <- function(Q, probs = c(0.33, 0.66)) {
  Qp <- quantile(Q, probs = probs, na.rm = TRUE, names = FALSE)
  - (log10(Qp[1L]) - log10(Qp[2L])) / diff(probs)
}


#------------------------------------------------------------------------------
#' Compute any flow duration curve percentiles
#' 
#' Given streamflow exceedance probability values (\code{probs}) compute the
#' corresponding streamflow values/percentiles.
#' 
#' @param Q numeric vector. Streamflow values.
#' @param probs numeric vector containing the exceedance probabilities for the
#' streamflow percentiles to compute.
#' @return A numeric vector of the same length as \code{probs} containing the
#' streamflow percentiles corresponding to the exceedance probabilities
#' specified with argument \code{probs}.
#' @seealso \code{\link{fdc_slope}} and \code{\link{fdc_slope_percentiles}}
#' @export
fdc_percentiles <- function(Q, probs = c(0.1, 0.9)) {
  out <- quantile(Q, probs = 1 - probs, na.rm = TRUE, names = FALSE)
  names(out) <- paste0('Q', probs) # This might be useless and slow down execution
  out
}

# This functions is only intended to speed up code execution by avoiding
# calling the quantile function twice to compute the slope and the percentils.
# However, if the exact same vector is fed the quantile function,
# R will remember it allready did the computation and return the results.
# Therefore, I recommend using the individual functions instead.
#------------------------------------------------------------------------------
#' Compute flow duration curve based hydrological signatures
#' 
#' This function compute both (1) the mid-segment flow duration curve slope 
#' (see \code{\link{fdc_slope}}) and (2) any streamflow percentiles (see 
#' \code{\link{fdc_percentiles}}).
#' 
#' @param Q numeric vector. Streamflow values.
#' @param slope_probs numeric vector of length 2 containing the exceedance
#' probabilities that define the start and end of the mid-segment of the 
#' flow duration curve.
#' @param percentile_probs numeric vector containing the exceedance
#' probabilities for the streamflow percentiles to compute
#' @return A numeric vector of length n+1 where n is the length of
#' \code{percentile_probs} with the first element being the slope of the
#' flow duration curve mid-segment.
#' @seealso \code{\link{fdc_slope}} and \code{\link{fdc_percentiles}}
#' @export 
fdc_slope_percentiles <- function(Q, slope_probs = c(0.33, 0.66), percentile_probs = c(0.1, 0.9)) {
  Qp <- quantile(Q, probs = c(slope_probs, 1 - percentile_probs), na.rm = TRUE, names = FALSE)
  out <- c(- (log10(Qp[1L]) - log10(Qp[2L])) / diff(slope_probs), Qp[-c(1L, 2L)])
  names(out) <- c("slope", paste0('Q', percentile_probs)) # I think this is useless and might slow down execution
  out
}

#------------------------------------------------------------------------------
#' Compute the flow duration curve (FDC)
#' 
#' Given a vector of streamflow values \code{Q}, this function computes a data.frame
#' with two columns:
#' * \code{p}: probabilities of exceedance
#' * \code{Q}: corresponding streamflow values
# Two methods can be used: simply sorting the data (not recommended)
# or using the \code{quantile} function.
#' 
#' @param Q numeric vector. Streamflow values.
#' @param n number of rows in the resulting data.frame (should be smaller than
#'  the length of \code{Q}).
#' @param sort logical. Whether the fdc should be computed by sorting the 
#' streamflow values (TRUE) or by using the \code{quantile} function (FALSE, 
#' default).
#' @param na.rm logical. Should the missing values be ignored? (must be TRUE if
#' the \code{quantile} function is used and \code{Q} contains missing values)
#' @return a data.frame with two columns:
#' * \code{p}: probabilities of exceedance
#' * \code{Q}: corresponding streamflow values
#' @seealso \code{\link{fdc_slope}}, \code{\link{fdc_percentiles}} and \code{\link{fdc_slope_percentiles}}
#' @export 
fdc_values <- function(Q, n = 1000, sort = FALSE, na.rm = TRUE) {
  if (na.rm) Q <- Q[!is.na(Q)]
  if (sort) {
    m <- length(Q)
    pfdc <- 1-1:m/m
    Qfdc <- sort(Q, na.last = ifelse(na.rm, NA, FALSE))
  } else {
    if (n > length(Q)) warning("'n' is larger than the number of values in 'Q'!")
    pfdc <- seq(0, 1, length.out = n)
    Qfdc <- quantile(Q, probs = 1 - pfdc)
  }
  return(data.frame(p = pfdc, Q = Qfdc))
}
