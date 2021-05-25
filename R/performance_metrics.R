
# FIXME: add reference 
#------------------------------------------------------------------------------
#' Computes the Nash-Sutcfliff efficiency coefficient
#' 
#' Given two streamflow vectors of same length (observer 'obs' and simulated
#' 'sim' streamflow), this functions computes the corresponding Nash-Sutcliff
#' efficiency. 
#' 
#' @param obs numeric vector. Observed streamflow vector.
#' @param sim numeric vector. Simulated streamflow vector.
#' @param na.rm logical. Whether missing values should be ignored.
#' @return Returns a float value between -Inf and 1 corresponding to the 
#' Nash-Sutcfliff efficiency coefficient.
#' @export
nse <- function(obs, sim, na.rm = TRUE) {
  if (na.rm) {
    isna <- is.na(obs) | is.na(sim)
    obs <- obs[!isna]
    sim <- sim[!isna]
  }
  1 - sum((sim - obs) ^ 2) / sum((obs - mean(obs)) ^ 2)
}


# FIXME: add the actual references
#------------------------------------------------------------------------------
#' Computes log transformed streamflow vector.
#' 
#' compute the log transformed value of streamflow using either the approach 
#' proposed by Pushpalatha et al. (2012) (i.e. a 100th of the mean value is
#' added to each streamflow value), or by adding any given value to all 
#' streamflow values or simply by setting to NA any non finite value resulting
#' from the application of the log() function.
#' 
#' @param x numeric vector. Streamflow vector to log transform.
#' @param method Method used to deal with missing values: 'Pushpalatha2012'
#' a numeric value (to be added to 'x') or 'inf.na'. If a charcter string 
#' equal to "Pushpalatha2012", the method of Pushpalatha et al. (2012) is used.
#' If a character string equal to "inf.na", any non-finite value are set to NA.
#' If a flow value, it is added to each streamflow value before transofmation.
#' @return Returns a numeric vector of the same length than x containing the log
#' transformed streamflow values.
#' @export
log_transform <- function(x, method = "Pushpalatha2012") {
  if (method == "Pushpalatha2012") x <- x + mean(x, na.rm = TRUE) / 100
  else if (is.numeric(method))     x <- x + method
  if (method == "inf.na")  {
    y <- suppressWarnings(log(x))
    y[!is.finite(y)] <- NA
  } else {
    y <- log(x)
  }
  y
}

#------------------------------------------------------------------------------
#' Computes the Nash-Sutcfliff efficiency coefficient on the log transformed
#' streamflow values
#' 
#' This function computes the Nash-Sutcliff efficiency coefficient using the 
#' nse() function on the streamflow values which are first log transformed 
#' using the log_transform()vfunction. 
#' 
#' @param obs numeric vector. Observed streamflow vector.
#' @param sim numeric vector. Simulated streamflow vector.
#' @param na.rm logical. Whether missing values should be ignored.
#' @param log_method character string or numeric. See log_transform().
#' @return Returns a float value between -Inf and 1 corresponding to the 
#' Nash-Sutcfliff efficiency coefficient computed on the log transofmed 
#' streamflow values.
#' @export
nselog <- function(obs, sim, na.rm = TRUE, log_method = "inf.na") {
  nse(obs = log_transform(obs, method = log_method),
      sim = log_transform(sim, method = log_method), na.rm = na.rm)
}

# FIXME: add the actual references...
#------------------------------------------------------------------------------
#' Computes the Kling-Gupta efficiency coefficient.
#' 
#' This function computes the Kling-Gupta efficiency coefficient. It islargely 
#' inspired by the similar function in the HydroGOF package. Here, only a 
#' simplified version is provided with no checks on inputs and very little 
#' formatting.
#' 
#' @param obs numeric vector. Observed streamflow vector.
#' @param sim numeric vector. Simulated streamflow vector.
#' @param na.rm logical. Whether missing values should be ignored.
#' @param method two methods are implemented (see hydroGOF::KGE): one from 
#' Gupta et al. (2009) and the other one from Kling et al. (2012)
#' @return Returns a float value between -Inf and 1 corresponding to the 
#' Kling-Gupta efficiency coefficient.
#' @export
kge <- function(obs, sim, na.rm = TRUE, method = 1) {
  if (na.rm) {
    isna <- is.na(obs) | is.na(sim)
    obs <- obs[!isna]
    sim <- sim[!isna]
  }
  mobs <- mean(obs)
  msim <- mean(sim)
  sobs <- sd(obs)
  ssim <- sd(sim)
  rso  <- cor(sim, obs)
  ALPHA <- ssim / sobs
  BETA  <- msim / mobs
  if (method == 1) {
    .kge(R = rso, AG = ALPHA, BETA = BETA)
  } else if (method == 2){
    cvobs <- sobs / mobs
    cvsim <- ssim / msim
    GAMMA <- cvsim / cvobs
    .kge(R = rso, AG = GAMMA, BETA = BETA)
  }  else {
    warning("Unknown method, only 1 and 2 are supported. Default method 1 used.")
    .kge(R = rso, AG = ALPHA, BETA = BETA)
  }
}
# an intermediate function in the computation of KGE.
# Here only to avoid code repetition
.kge <- function(R, AG, BETA) 1 - sqrt((R - 1) ^ 2 + (AG - 1) ^ 2 + (BETA - 1) ^ 2)


#------------------------------------------------------------------------------
#' compute the bias between simulated and observed data
#' 
#' This function compute the bias (unitless) between simulated and observed
#' streamflow data.
#' 
#' @param obs numeric vector. Observed streamflow vector.
#' @param sim numeric vector. Simulated streamflow vector.
#' @param na.rm logical. Whether missing values should be ignored.
#' @param sim_minus_obs logical.  should it be sim - obs? (or the other
#' way around)
#' @return Returns a float value corresponding to the bias.
#' @export
bias <- function(obs, sim, na.rm = TRUE, sim_minus_obs = TRUE) {
  if (length(obs) != length(sim)) stop("obs and sim must have the same length")
  if (na.rm) {
    isna <- is.na(obs) | is.na(sim)
    obs <- obs[!isna]
    sim <- sim[!isna]
  }
  if (sim_minus_obs) sum(sim - obs) / sum(obs) else sum(obs - sim) / sum(obs)
}

#------------------------------------------------------------------------------
#' compute the coefficient of determination between simulated and observed data
#' 
#' This function compute the compute the coefficient of determination using
#' the cor() function between simulated and observed data
#' 
#' @param obs numeric vector. Observed streamflow vector.
#' @param sim numeric vector. Simulated streamflow vector.
#' @param na.rm logical. Whether missing values should be ignored.
#' @param method character. see cor() function.
#' @return Returns a float value corresponding to the coefficient of 
#' determination.
#' @export
rsquare <- function(obs, sim, na.rm = TRUE,
 method = c("pearson", "kendall", "spearman")) {
  if (na.rm) use = "na.or.complete" else use = "everything"
  cor(obs, sim, use = use, method = method)
}

