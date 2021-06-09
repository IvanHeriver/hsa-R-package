

#------------------------------------------------------------------------------
#' Compute the inflexion point of the streamflow mass curve
#' 
#' From a streamflow regime vector \code{Qr} (i.e. a inter-annual average of
#' streamflow for each calendar day), it searches for the inflexion points in
#' a given time window \code{p}. This functions (1) computes the derivative of
#' \code{Qr}, (2) and smooth it using a rolling mean with a window \code{n}.
#' 
#' @param Qcum numeric vector. Cumulative streamflow regime vector.
#' @return If \code{return_derivative} is TRUE, a list with 3 components is
#' returned:
#' * \code{data}: a data.frame containing two column: the derivative of \code{Qr}
#' (dQ) and its smoothed version (dQs).
#'  \code{inflexion_points} the indices of the two inflexions points found
#' If \code{return_derivative} is FALSE, only the two inflexions points are 
#' returned in a length 2 numeric vector.
#' @seealso [mct_snowstorage]
#' @export
mct_inflexionpoints <- function(Qcum, periods = c(100, 300), n = 30L, return_derivative = FALSE) {
  dQ   <- c(diff(Qcum), NA)
  dQs  <- roll_mean(dQ, n=n[1L], fill=NA)
  p <- periods[1L]:periods[2L]
  inflex_points <- c(which.max(dQs[p]), which.min(dQs[p])) + p[1L] - 1L
  if (return_derivative) return(list(data = data.frame(dQ = dQ, dQs = dQs), inflexion_points = inflex_points))
  return(inflex_points)
}


# TODO: check dims of inputs!
#------------------------------------------------------------------------------
#' Compute snow storage estimate using the MCT technique
#' 
#' Given cumulative regime curve of precipitation \code{Pcum} and streamflow
#' \code{Qcum}, the inflexions points of the streamflow regime curve returned
#' by the function \code{mct_inflexionpoints()}, this function compute the 
#' one or two snow storage estimates depending on the method used (see details).
#' 
#' If \code{use_tangentes} is FALSE, the snow storage estimate is simply the 
#' difference \code{Pcum[inflexion_points[1L]] - Qcum[inflexion_points[1L]]}.
#' If \code{use_tangentes} is TRUE, a linear regression is performed on 
#' Pcum to retrieve the average slope in the period defined by \code{period}.
#' If \code{fixed_intercept} is TRUE, a fixed intercept at (0, 0) is used 
#' in the linear regression. If \code{use_tangentes} is TRUE, two snow estimates
#' are computed: (1) from the accumulation period (Sa) and (2) from the 
#' metling period (Sm). For details see Schaefli (2016)
#' 
#' @param Pcum numeric vector. Cumulative precipitation regime curve.
#' @param Qcum numeric vector. Cumulative streamflow regime curve.
#' @param inflexion_points numeric vector. A length 2 numeric vector containing
#' the two inflexion points of the cumulative streamflow regime curve (see
#' \code{mct_inflexionpoints()})
#' @param use_tangentes logical. Should the methodology of Schaefli (2016) be used
#' to compute two snow estimates (see details).
#' @param period. numeric vector. A length 2 numeric vector defining the period
#' to use to perform the linear regression on \code{Pcum} (see details).
#' @param fixed_intercept logical. Should a fixed intercept at (0, 0) be used
#' in the linear regression (see details).
#' @param return_tangents logical. Should the tangents values be returned as well.
#' @return There are three cases: 
#'  1. \code{use_tangentes} is FALSE: a single numeric value corresponding to one 
#' storage estiamte.
#'  2. \code{use_tangentes} is TRUE and \code{return_tangents} is FALSE: a length 2
#' numeric vector is returned containing 2 snow storage estimates: (1) from the
#' accumulation period (Sa) and (2) from the  metling period (Sm).
#'  3. \code{use_tangentes} is TRUE and \code{return_tangents} is TRUE: a list is 
#' returned containing the the two snow storage estimates (S), the slope of resulting
#' from the linear regression over \code{Pcum} (P_slp), the intercept and slope of the
#' the two tangents to streamflow at inflexion point 1 (Qa_slp) and inflexion point 2
#' (Qm_slp). These last two components are length 2 numeric vectors.
#' @seealso [mct_snowstorage]
#' @references
#' B. Schaefli, “Snow hydrology signatures for model identification within a limits-of-acceptability approach”, 
#' Hydrological Processes, vol. 30, no. 22, Art. no. 22, Aug. 2016, doi: 10.1002/hyp.10972.
#' @export
mct_snowstorage <- function(Pcum, Qcum, inflexion_points, 
                            use_tangentes = FALSE, 
                            period = c(1, length(Pcum)),
                            fixed_intercept = FALSE,
                            return_tangents = FALSE) {
  if (use_tangentes) {
    x <- period[1L]:period[2L]
    if (fixed_intercept) {
      Pslp <- c(0, lm(Pcum[x]~x+0)$coefficients) # FIXME: this is slow, find a way using .lm.fit() function
    } else {
      Pslp <- .fast_lm(x = x, y = Pcum[x])$coefficients
    }
    b <- Qcum[inflexion_points] - (Pslp[2L] * inflexion_points)
    if (return_tangents) {
      return(list(S = c(Sa = Pslp[1L] - b[1L],
                        Sm = b[2L] - b[1L]),
                  P_slp = Pslp,
                  Qa_slp = c(b[1L], Pslp[2L]),
                  Qm_slp = c(b[2L], Pslp[2L])))
    } else {
      return(c(Sa = Pslp[1L] - b[1L],
               Sm = b[2L] - b[1L]))
    }
  } else {
    return(Pcum[inflexion_points[1L]] - Qcum[inflexion_points[1L]])
  }
}

