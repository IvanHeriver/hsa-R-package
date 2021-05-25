# TODO: check dims of inputs!
# FIXME: this function might be obsolete: use hsaRegimeCurve() and hsaCumSum() 
# with apply()
# Note: In this version the inter-annual calendar day mean are computed
# first and then the cumulative curves are computed
# Might be obsolete
#------------------------------------------------------------------------------
#' Computes the cumulative inter-annual daily average of P and Q as
#' well as the difference P-Q.
#' 
#' @param Q numeric vector. Streamflow values.
#' @param P numeric vector. Precipitation values.
#' @param hdays numeric vector. Days of the (Hydrological) years vector
#' @return A 3-columns data.frame is returned (Q, P and P_Q) containing the
#' cumulative curve of streamflow (Q), precipitation (P) and the differences 
#' (P-Q)
#' @export
pq_cumsum <- function(Q, P, hdays) {
  data.frame(Q, P, hdays) %>% 
    group_by(hdays) %>%
    summarise_all(mean, na.rm = TRUE) %>%
    mutate_at(vars(Q, P), cumsum_interpolate) %>%
    mutate(P_Q = P - Q) %>%
    data.frame()
}

#------------------------------------------------------------------------------
#' P-Q signatures computations using P-Q cumulative curve
#' 
#' Compute the seasonal response change signatures: breackpoint date,
#' first period slope (dry), second period slope (wet) and the intercepts 
#' associated with these two slopes (only there for plotting purposes).
#' 
#' This function uses a segmented regression to find the dry and wet slopes.
#' It also returns the intercepts of the two lines (which are not signatures)
#' that can be used to plot the lines. 
#' 
#' @param PQ numeric vector. The cumulative precipitation minus 
#' cumulative streamflow (P-Q) vector.
#' @param start numeric. The day number of the start of the period to 
#' search for a change of trend in P-Q.
#' @param end numeric. The day number of the end of the period to 
#' search for a change of trend in P-Q.
#' @param bp numeric value. The day number of the initial guess of threshold
#' date
#' @param intercept logical. Should the intercept be estimated (default: TRUE)
#' or fixed to c(0, 0) (FALSE)?
#' @return A named vector of length 6 containing the breakpoint date ('bp'),
#' breakpoint strength signature ('bp_strength'), the slope of the first (dry)
#' period ('slp_dry') and its corresponding intercept ('b_dry'), the slope of
#' the second (wet) period ('slp_wet') and its corresponding intercept 
#' ('b_wet')
#' @export
pq_curve_slopes <- function(PQ, start = 15, end = 183, bp = mean(c(start, end)), intercept = TRUE) {
  x <- start:end
  y <- PQ[x]
  if (intercept) {
    reg <- segmented(lm(y ~ x), psi = bp)
    coefs <- reg$coefficients
  } else {
    reg <- segmented(lm(y ~ x + 0), psi = bp)
    coefs <- c(0, reg$coefficients)
  }
  coefs[3L] <- sum(coefs[2:3])
  bp <- reg$psi[, 2L]
  out <- c(bp, 
           1 - coefs[3L] / coefs[2L],
           coefs[2L],
           coefs[1L],
           coefs[3L],
           coefs[2L] * bp + coefs[1L] - coefs[3L] * bp)
  names(out) <- c("bp", "bp_strength", "slp_dry", "b_dry", "slp_wet", "b_wet")
  out
}

# TODO: check dims of inputs!
#------------------------------------------------------------------------------
#' P-Q signatures computations using Q, P and hdays vectors
#' 
#' It computes the cumulative inter-annual daily average of P and Q as
#' well as the difference P-Q. Then, it computes the seasonal response change
#' signatures following the so-called P-Q approach: breackpoint date,
#' first period slope (dry), second period slope (wet) and the intercepts 
#' associated with these two slopes (only there for plotting purposes).
#' 
#' This function uses a segmented regression to find the dry and wet slopes.
#' It also returns the intercepts of the two lines (which are not signatures)
#' that can be used to plot the lines. 
#' 
#' @param Q numeric vector. Streamflow values.
#' @param P numeric vector. Precipitation values.
#' @param hdays numeric vector. Days of the (Hydrological) years vector
#' @param start numeric. The day number of the start of the period to 
#' search for a change of trend in P-Q.
#' @param end numeric. The day number of the end of the period to 
#' search for a change of trend in P-Q.
#' @param bp numeric value. The day number of the initial guess of threshold
#' date
#' @param intercept logical. Should the intercept be estimated (default: TRUE)
#' or fixed to c(0, 0) (FALSE)?
#' @return Returned
#' @export
pq_slopes <- function(Q, P, hdays, start = 15, end = 183, bp = mean(c(start, end)), intercept = TRUE) {
  PQ <- data.frame(Q, P, hdays) %>% 
    group_by(hdays) %>%
    summarise_all(mean, na.rm = TRUE) %>%
    mutate_at(vars(Q, P), hsaCumsum) %>%
    mutate(P_Q = P - Q)
  PQ <- PQ$P_Q
  x <- start:end
  y <- PQ[x]
  if (intercept) {
    reg <- segmented(lm(y ~ x), psi = bp)
    coefs <- reg$coefficients
  } else {
    reg <- segmented(lm(y ~ x + 0), psi = bp)
    coefs <- c(0, reg$coefficients)
  }
  coefs[3L] <- sum(coefs[2:3])
  bp <- reg$psi[, 2L]
  
  out <- c(bp, 1 - coefs[3L] / coefs[2L], coefs[2L], coefs[1L], coefs[3L], coefs[2L] * bp + coefs[1L] - coefs[3L] * bp)
  names(out) <- c("bp", "bp_strength", "slp_dry", "b_dry", "slp_wet", "b_wet")
  out
}

# TODO: check dims of inputs!
#------------------------------------------------------------------------------
#' Compute P-Q snow storage estimate
#' 
#' Given the cumulative curve of precipitation (Pcum) and streamflow (Qcum), it
#' computes the differences (Pcum - Qcum) and computes a snow storage estimate
#' corresponding to the maximum of the Pcum - Qcum curve.
#' 
#' @param Pcum Numeric vector. Cumulative precipitation over one hydrological 
#' year (or inter-annual average hydrological year).
#' @param Qcum Numeric vector. Cumulative streamflow over one hydrological 
#' year (or inter-annual average hydrological year).
#' @param period NULL or a numeric vector of length 2 that defines the end and
#' start (number of days) of the period where the maximum of the P-Q curve 
#' is to be found.
#' @param return_data logical. Should the intermediate results (P-Q curve and
#' number of the day corresponding to the maximum of the P-Q curve) be also
#' returned? Default is FALSE. and a single value corresponding to the snow 
#' storage estimate is returned. If TRUE, a list with the snow storage estimate
#' ('S'), and the P-Q curve ('PQ') and the date of the maximum of the P-Q curve
#' ('tQcumMax') is returned.
#' @return If return_data is FALSE, a single value corresponding to the snow 
#' storage estimate is returned. Otherwise, a list with the snow storage
#'  estimate ('S'), and the P-Q curve ('PQ') and the date of the maximum of
#'  the P-Q curve ('tQcumMax') is returned.
#' @export
pq_snowstorage <- function(Pcum, Qcum, period=NULL, return_data = FALSE) {
  PQ <- Pcum - Qcum
  if (is.null(period)) {
    iQmax <- which.max(PQ) 
  } else {
    x <- period[1L]:period[2L]
    iQmax <- which.max(PQ[x]) + period[1L] - 1
  }
  if (return_data) return(list(S = PQ[iQmax], PQ = PQ, tPQcumMax = iQmax))  else return(PQ[iQmax])
}