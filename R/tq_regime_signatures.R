

# FIXME: these functions could also be used for: PQ slopes, BFRmag!
# -------------------------------------------------------------------
# Compute regime curves from the data provided in 'x'.
# It first apply smoothing (rolling mean) with a wondow 'n' of the 
# data provided in 'x'. Then, it compute the inter-annual mean of 
# each calendar day (provided in vector 'hdays').
# It recommended to ignore the last day using 'remove_last = TRUE'.
# -------------------------------------------------------------------
# hdays:        vector of 'julian' date of length m which should be
#               computed according to hydrological year
#               (see function 'hsaHydroYearSeasons()')
# x:            a vector (of length m) or a matrix (with m rows)
# n:            smoothing window
# remove_last:  should the last day of the year (366) be removed
#               (recommended)
hsaRegimeCurve <- function(hdays, x, n = 1L, na.rm = TRUE, remove_last = TRUE) {
  if (is.null(dim(x))) x <- matrix(x, dimnames = list(NULL, as.character(substitute(x))))
  if (n > 1L) x <- apply(x, 2L, roll_mean, n = n, fill = NA)
  y <- data.frame(hdays, x)  %>% 
    group_by(hdays) %>% 
    summarise_all(mean, na.rm = na.rm) %>% 
    data.frame()
  if (remove_last) return(y[-nrow(y), ]) else return(y)
}

# -------------------------------------------------------------------
# Given a vector or a matrix/data.frame 'x', compute the corresponding
# cumulative curve.
# This function is only a wrapper around hsaCumsum() with an apply() 
# call.
# -------------------------------------------------------------------
# x:            a vector (of length m) or a matrix/data.frame (with m rows)
# na.action:    a character string indicating whether missing values
#               should be interpolated or not (see hsaCumsum() function)
hsaCumulativeCurve <- function(x, na.action = "interpolate") {
  if (is.null(dim(x))) x <- matrix(x)
  y <- apply(x, 2L, hsaCumsum, na.action = na.action)
  as.data.frame(y)
}

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
# -------------------------------------------------------------------
# Compute slopes (and intercepts) of each periods from the air
# temperature - streamflow regime cycle using the periods defined
# in a list of indices vectors
# -------------------------------------------------------------------
# T:            air temperature regime (365 or 366 values)
# Q:            streamflow regime (365 or 366 values)
# indices:      list of vectors of indices defining the different
#               period to use to compute the slopes/intercept in
#               T = b + slp * Q
tq_slopes <- function(T, Q, indices) {
  coefs <- vapply(indices, function(e) .hsaFastLm(Q[e],T[e])$coefficients, numeric(2L))
  rownames(coefs) <- c("intercept", "slope")
  coefs
}

# -------------------------------------------------------------------
# Compute the timing of the maximum of the streamflow regime. This
# function is only a wrapper arount which.max() doing ... nothing.
# -------------------------------------------------------------------
# Q:            streamflow regime (365 or 366 values)
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
tq_signatures <- function(hdays, T, Q, n = 30, remove_last = TRUE, periods = c(183, 336, 62)) {
  TQ_regime <- hsaRegimeCurve(hdays = hdays, x = data.frame(T, Q), n = n, remove_last = remove_last)
  TQ_indices <- tq_periods(T = TQ_regime$T, Q = TQ_regime$Q, periods = periods)
  TQ_sigs <- c(vapply(TQ_indices, function(e) .hsaFastLm(TQ_regime$Q[e], TQ_regime$T[e])$coefficients[2L], numeric(1L)), periods[[2L]][1L])
  names(TQ_sigs) <- c("TQslp_TpQp", "TQslp_TpQm", "TQslp_TmQm1", "TQslp_TmQm2", "tQmax")
  TQ_sigs
}

