#------------------------------------------------------------------------------
#' Compute the runoff coefficient
#' 
#' compute the runoff coefficient or ratio: the ratio between the total 
#' strealfow volume and the total precipitaion volume. 
#' 
#' Warning: no checks are done on inputs! You should have selected the proper
#' period beforehand: both vector should be of same length and span over the
#' same period. In addition the runoff coefficient should be compute from the
#' start of a hydrological year to the end of a hydrological year.
#' 
#' @param Q numeric vector. Streamflow values.
#' @param P numeric vector. Precipitation values.
#' @param na.rm logical. Should missing values be omited?
#' @return A single numeric value corresponding to the runoff coefficient
#' @export
runoff_coefficient <- function(Q, P, na.rm = TRUE) {
  if (length(Q) != length(P)) warning("'Q' and 'P' don't have the same length!")
  hsaSumRatio(Q, P, na.rm = na.rm)
}
