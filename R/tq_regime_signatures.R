

# TODO: check dims of inputs!
# -------------------------------------------------------------------
# Compute list containing the indices of each periods necessary to
# compute the slopes from the T-Q regime cycle:
# > T+Q+ (TpQp) period
# > T+Q- (TpQm) period
# > T-Q-1 (TmQm1) period
# > T-Q-2 (TmQm2) period
# -------------------------------------------------------------------
# T:            air temperature regime (365 or 366 values)
# Q:            streamflow regime (365 or 366 values)
# TpQp_start, TpQp_end, TpQm_start, TpQm_end, TmQm1_start, TmQm1_end,
# TmQm2_start, TmQm2_end:
#               indices (relative to a hydrological year) of the start
#               and end of each period.
#               Note that:
#                if TpQp_end, TpQm_start are set to NA the timing of 
#                  the maximum of Q (streamflow regime) is used
#                if TpQm_end, TmQm1_start are set to NA the timing of
#                  the maximum of  T (air temperature regime) is used
# periods:      a vector of length 3 defining (1) the start of the 
#               T+Q+ period, (2) the start (resp. end) of the
#               T-Q-2 (resp. T-Q-1) period and (3) the end of the 
#               T-Q-2 period.
#               Note: the other period boundaries are computed using
#               the air temperature and streamflow regime maxima timings.
#------------------------------------------------------------------------------
#' Compute temperature-streamflow regime cycle periods
#' 
#' Compute list containing the indices of each periods necessary to compute 
#' the slopes from the T-Q regime cycle:
#' * T+Q+ (TpQp) period
#' * T+Q- (TpQm) period
#' * T-Q-1 (TmQm1) period
#' * T-Q-2 (TmQm2) period
#'
#' @param T numeric vector. air temperature regime (365 or 366 values)
#' @param Q numeric vector. streamflow regime (365 or 366 values)
#' @param peridos numeric vector. a vector of length 3 defining (1) the start
#' of the T+Q+ period, (2) the start (resp. end) of the T-Q-2 (resp. T-Q-1) 
#' period and (3) the end of the T-Q-2 period. Note: the other period 
#' boundaries are computed using the air temperature and streamflow regime
#' maxima timings.
#' @references
#' I. Horner, F. Branger, H. McMillan, O. Vannier, and I. Braud,
#' “Information content of snow hydrological signatures based on streamflow, 
#' precipitation and air temperature,” 
#' Hydrological Processes, vol. 34, no. 12, Art. no. 12, 2020, doi: 10.1002/hyp.13762.
#' @return A list containing the indices of the different periods in 4 
#' components: "TpQp", "TpQm", "TmQm1", "TmQm2"
#' @export
tq_periods <- function(T, Q, periods = c(183, 336, 62)) {
  i_Qmax <- which.max(Q)[1L]
  i_Tmax <- which.max(T)[1L]
  TQ_indices <- list()
  TQ_indices[[1L]] <- periods[1L]:i_Qmax   # TpQp
  TQ_indices[[2L]] <- i_Qmax:i_Tmax        # TpQm
  TQ_indices[[3L]] <- i_Tmax:periods[2L]   # TmQm1
  TQ_indices[[4L]] <- c(periods[2L]:length(Q), 1:periods[3L]) # TmQm2
  names(TQ_indices) <- c("TpQp", "TpQm", "TmQm1", "TmQm2")
  TQ_indices
}

# TODO: check dims of inputs!
#------------------------------------------------------------------------------
#' Compute the temperature-streamflow regime cycle slopes
#' 
#' Compute slopes (and intercepts) of each periods from the air 
#' temperature - streamflow regime cycle using the periods defined in a list of
#' indices vectors (see \code{tq_periods()})
#' 
#' @param T numeric vector. air temperature regime (365 or 366 values)
#' @param Q numeric vector. streamflow regime (365 or 366 values)
#' @param indices list. list of  length 'n'  containing the vectors of indices
#'  defining the different period to use to compute the slopes/intercept in
#' `T = b + slp * Q`
#' @references
#' I. Horner, F. Branger, H. McMillan, O. Vannier, and I. Braud,
#' “Information content of snow hydrological signatures based on streamflow, 
#' precipitation and air temperature,” 
#' Hydrological Processes, vol. 34, no. 12, Art. no. 12, 2020, doi: 10.1002/hyp.13762.
#' @seealso [tq_periods]
#' @return A matrix with 'n' rows (the different periods) and 2 columns: 'intercept'
#' and 'slope'
#' @export
tq_slopes <- function(T, Q, indices) {
  coefs <- vapply(indices, function(e) .fast_lm(Q[e],T[e])$coefficients, numeric(2L))
  rownames(coefs) <- c("intercept", "slope")
  coefs
}

#------------------------------------------------------------------------------
#' Compute streamflow regime maximum timing
#' 
#' Compute the timing of the maximum of the streamflow regime. This function 
#' is only a wrapper arount which.max() doing ... nothing.
#' 
#' @param Q numeric vector. streamflow regime (365 or 366 values)
#' @return the index of the maximum of input 'Q'
#' @export
q_timing <- function(Q) which.max(Q)

# TODO: check dims of inputs!
# -------------------------------------------------------------------
# Compute signatures based on the air temperature (T) and streamflow (Q)
# regime: slopes of 4 periods (T+Q+, T+Q-, T-Q-1, T-Q-2  periods) of 
# the T-Q cycle and timing of the Q maximum
# -------------------------------------------------------------------
# hdays:        vector of 'julian' date of length m which should be
#               computed according to hydrological year
#               (see function 'hsaHydroYearSeasons()')
# T:            air temperature vector (of length m)
# Q:            streamflow vector (of length m)
# n:            smoothing window
# remove_last:  should the last day of the year (366) be removed
#               recommended
# periods:      a vector of length 3 defining (1) the start of the 
#               T+Q+ period, (2) the start (resp. end) of the
#               T-Q-2 (resp. T-Q-1) period and (3) the end of the 
#               T-Q-2 period.
#               Note: the other period boundaries are computed using
#               the air temperature and streamflow regime maxima timings.
#------------------------------------------------------------------------------
#' Compute signatures based on the air temperature (T) and streamflow (Q) regime
#' 
#' This function computes signatures based on the air temperature (T) and 
#' streamflow (Q) regime cycle: the slopes of 4 periods (T+Q+, T+Q-, T-Q-1, 
#' T-Q-2  periods) of the T-Q cycle and timing of the Q maximum. 
#' 
#' @param hdays numeric vector. vector of 'julian' date of length m which 
#' should be computed according to hydrological year (see function 
#' \code{get_hydro_years_seasons()})
#' @param T numeric vector. air temperature regime (365 or 366 values)
#' @param Q numeric vector. streamflow regime (365 or 366 values)
#' @param n integer. smoothing window
#' @param remove_last logical. should the last day of the year (366) be 
#' @param peridos numeric vector. a vector of length 3 defining (1) the start
#' of the T+Q+ period, (2) the start (resp. end) of the T-Q-2 (resp. T-Q-1) 
#' period and (3) the end of the T-Q-2 period. Note: the other period 
#' boundaries are computed using the air temperature and streamflow regime
#' maxima timings.
#' @references
#' I. Horner, F. Branger, H. McMillan, O. Vannier, and I. Braud,
#' “Information content of snow hydrological signatures based on streamflow, 
#' precipitation and air temperature,” 
#' Hydrological Processes, vol. 34, no. 12, Art. no. 12, 2020, doi: 10.1002/hyp.13762.
#' @seealso [tq_periods][regime_curve][q_timing][tq_slopes]
#' @return A vector of length 5 with the following named components:
#' "TQslp_TpQp", "TQslp_TpQm", "TQslp_TmQm1", "TQslp_TmQm2", "tQmax"
#' @export
tq_signatures <- function(hdays, T, Q, n = 30, remove_last = TRUE, periods = c(183, 336, 62)) {
  TQ_regime <- regime_curve(hdays = hdays, x = data.frame(T, Q), n = n, remove_last = remove_last)
  TQ_indices <- tq_periods(T = TQ_regime$T, Q = TQ_regime$Q, periods = periods)
  TQ_sigs <- c(vapply(TQ_indices, function(e) .fast_lm(TQ_regime$Q[e], TQ_regime$T[e])$coefficients[2L], numeric(1L)), periods[[2L]][1L])
  names(TQ_sigs) <- c("TQslp_TpQp", "TQslp_TpQm", "TQslp_TmQm1", "TQslp_TmQm2", "tQmax")
  TQ_sigs
}

