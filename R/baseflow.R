
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
  # C++ code for faster execution
  i <- c(1, ioh_min_pivots(Q = Q, d = d, k = k) + 1, n)
  Qbf <- approx(i, Q[i], xout = 1:n)$y 
  a <- Qbf > Q
  Qbf[a] <- Q[a]
  Qbf[isna] <- NA
  Qbf
}

#------------------------------------------------------------------------------
#' Compute the baseflow index
#' 
#' This function compute the baseflow index: the ratio between the total
#' baseflow volume and the total streamflow volume.
#' 
#' Warning: No checks are done on inputs except the time series length
#' of Q and Qbf that should match. You should have selected the proper
#' period before hand: both vector should be of same length and span over
#' the same period.
#' 
#' @param Q numeric vector. Streamflow vector.
#' @param Qbf numeric vector. Baseflow vector.
#' @param na.rm logical. Should missing values be omited?
#' @return A single value corresponding to the baseflow index
#' @export
baseflow_index <- function(Q, Qbf, na.rm = TRUE) {
  if (length(Q) != length(Qbf)) warning("'Q' and 'Qbf' don't have the same length!")
  sum_ratio(Qbf, Q, na.rm = na.rm)
} 

# TODO: check dims of inputs!
# FIXME: do not call roll_mean when n == 1
#------------------------------------------------------------------------------
#' Compute the baseflow regime magnitude
#' 
#' Compute the baseflow magnitude i.e. the relative difference between the 
#' maxima and minima of the 'n'-days smoothed inter-annual daily average of
#' the baseflow time series
#' 
#' A length 3-vector is return with the min and the max in second and third
#' position respectively
#' 
#' @param Qbf numeric vector. Baseflow vector.
#' @param hdays numeric vector. Days of the (hydrological) year. see
#' get_hydro_years_seasons() function.
#' @param n integer. Window size for the rolling mean function
#' @return A length 3-vector is return with the baseflow regime magnitude, 
#' the min and the max of the baseflow regime
#' @seealso [baseflow_regime][baseflow_magnitude]
#' @export
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
#------------------------------------------------------------------------------
#' Compute the baseflow regime 
#' 
#' Compute the baseflow regime vector i.e. the 'n'-days smoothed inter-annual 
#' daily average of the baseflow time series
#' 
#' @param Qbf numeric vector. Baseflow vector.
#' @param hdays numeric vector. Days of the (hydrological) year. see
#' get_hydro_years_seasons() function.
#' @param n integer. Window size for the rolling mean function
#' @return A 366-long (or as many different hydrological days found in hdays 
#' vector) numeric vector containing the baseflow regime values for each 
#' calendar day.
#' @seealso [baseflow_regime_magnitude][baseflow_magnitude]
#' @export
baseflow_regime <- function(Qbf, hdays, n = 1) {
  x <- data.frame(Qbf, hdays)
  y <- x %>% group_by(hdays) %>% summarise(Qbf_dailymean = mean(Qbf, na.rm = TRUE))
  roll_mean(y[["Qbf_dailymean"]], n = n, fill = NA, align = "center")
}

#------------------------------------------------------------------------------
#' Compute the baseflow regime magniture
#' 
#' Compute the baseflow regime vector using the already computed baseflow 
#' regime resulting from \code{baseflow_regime()} function
#' 
#' @param Qbfr numeric vector. Baseflow regime vector resulting from the 
#' \code{baseflow_regime()} function.
#' @return A length 3-vector is return with the baseflow regime magnitude, 
#' the min and the max of the baseflow regime
#' @seealso [baseflow_regime_magnitude][baseflow_regime]
#' @export
baseflow_magnitude <- function(Qbfr) {
  max_y <- max(Qbfr, na.rm = TRUE)
  min_y <- min(Qbfr, na.rm = TRUE)
  c(mag = (max_y - min_y) / max_y, min = min_y, max = max_y)
}
