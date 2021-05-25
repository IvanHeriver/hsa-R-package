
# -------------------------------------------------------------------
# Most of the baseflow extraction algotrithm are coded in C++ and 
# can be called in a similar way as the baseflow_Gustard() function
# given below (which also depends on C++ code).
# The following algorithms are available:
# > baseflow_LyneHollick(Q, k)
# > baseflow_ChapmanMaxwell(Q, k)
# > baseflow_Boughton(Q, k, C)
# > baseflow_Eckhardt(Q, a, BFImax)
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Institute of Hydrology algorithm for the computation of baseflow
# -------------------------------------------------------------------
# Q:            Streamflow vector
# d:            Window size (in days) to look for local minima
# k:            Factors to apply to local minima when looking for
#               pivot points (i.e. local local minima)
#------------------------------------------------------------------------------
#' compute the baseflow time series using the Gustard algorithm.
#' 
#' This function compute baseflow time series using the algorithm of Gustard
#' et al. (XXXX) (a.k.a. method of the Institude of hydrology).
#' 
#' @param Q numeric vector. Streamflow vector.
#' @param d integer. Window size (in days) to look for local minima
#' (Default: 5).
#' @param k float. Factor to apply to local minima when looking for
#' pivot points (i.e. local local minima)
#' @return Returns a baseflow time series.
#' @export
baseflow_Gustard <- function(Q, d = 5, k = 0.9) {
  n <- length(Q)
  isna <- is.na(Q)
  Q[isna] <- 99999999
  i <- c(1, ioh_min_pivots(Q = Q, d = d, k = k) + 1, n) # C++ code for faster execution
  Qbf <- approx(i, Q[i], xout = 1:n)$y 
  a <- Qbf > Q
  Qbf[a] <- Q[a]
  Qbf[isna] <- NA
  Qbf
}

# -------------------------------------------------------------------
# compute the baseflow index: the ratio between the total baseflow
# volume and the total streamflow volume. 
# Warning:
# No checks are done on inputs!
# You should have selected the proper period before hand: both vector
# should be of same length and span over the same period.
# -------------------------------------------------------------------
# Q:            Streamflow vector
# Qbf:          Baseflow vector
# na.rm:        Should missing values be omited?
baseflow_index <- function(Q, Qbf, na.rm = TRUE) {
  if (length(Q) != length(Qbf)) warning("'Q' and 'Qbf' don't have the same length!")
  hsaSumRatio(Qbf, Q, na.rm = na.rm)
} 

# TODO: check dims of inputs!
# FIXME: do not call roll_mean when n == 1
# -------------------------------------------------------------------
# Compute the baseflow magnitude: it is the relative difference
# between the maxima and minima of the 'n'-days smoothed inter-annual
# daily average of baseflow time series
# Note in this second version, the min and the max are also returned
# and the default smoothing is 1 (i.e. no smoothing)
# -------------------------------------------------------------------
# Qbf:          Baseflow vector
# days:         Vector containting the days of the (hydrological) year 
# n:            Window size for the rolling mean function
baseflow_regime_magnitude <- function(Qbf, hdays, n = 1) {
  x <- data.frame(Qbf, hdays)
  y <- x %>% group_by(hdays) %>% summarise(Qbf_dailymean = mean(Qbf, na.rm = TRUE))
  y <- roll_mean(y[["Qbf_dailymean"]], n = n, fill = NA, align = "center")
  max_y <- max(y, na.rm = TRUE)
  min_y <- min(y, na.rm = TRUE)
  c(mag = (max_y - min_y) / max_y, min = min_y, max = max_y)
}

# FIXME: write documentation
# TODO: check dims of inputs!
# FIXME: do not call roll_mean when n == 1
baseflow_regime <- function(Qbf, hdays, n = 1) {
  x <- data.frame(Qbf, hdays)
  y <- x %>% group_by(hdays) %>% summarise(Qbf_dailymean = mean(Qbf, na.rm = TRUE))
  roll_mean(y[["Qbf_dailymean"]], n = n, fill = NA, align = "center")
}

# FIXME: write documentation
baseflow_magnitude <- function(Qbfr) {
  max_y <- max(Qbfr, na.rm = TRUE)
  min_y <- min(Qbfr, na.rm = TRUE)
  c(mag = (max_y - min_y) / max_y, min = min_y, max = max_y)
}

