
#------------------------------------------------------------------------------
#' Extract recession events from a streamflow time series
#' 
#' This function extracts the index of the start and end of each recession 
#' event which are returned in a two columns data.frame
#' 
#' @param Q numeric vector. Streamflow vector.
#' @param minduration integer. Minimum duration of a recession event in days
#' @param maxduration integer. Maximum duration of a recession event in days
#' @param minvalue float. Minimum streamflow value to consider a local maximum
#' as a potential start of a recession event
#' @param n_smooth integer. How many times the rolling mean is to be applied
#' on the streamflow timeseries before applying the recession events extraction
#' @param window_smooth integer. The size, in days, of the window over which
#' the rollogin mean is to be applied.
#' @param return_smoothed logical. Whether the smoothed streamflow timeseries 
#' should be returned as well.
#' @return Returns a two columns data.frame with the 'start' and 'end' index
#' of each recession event. If 'return_smoothed' is TRUE, a list with two
#' elements is returned: 'i_events' containing the the data.frame with the 
#' 'start' and 'end' index of the recession events and 'Q_smoothed' containing
#' the smoothed streamflow vector.
#' @export
recession_events_index <- function(Q, minduration, maxduration, minvalue,
  n_smooth, window_smooth, return_smoothed=FALSE) {
  # optional smooting of Q
  if (n_smooth > 0) {
    for (k in 1:n_smooth) {
      Q <- roll_mean(x = Q, n = window_smooth, fill = NA, align = "center")
    }
  }
  # find local maxima and minima
  minmax <- diff(sign(diff(Q)))
  ismax <- c(FALSE, minmax < 0, FALSE)
  # ismin <- c(FALSE, minmax > 0, FALSE)
  ismin <- c(FALSE, minmax > 0 | is.na(minmax), FALSE)
  # select those above a certain Q threshold
  ismaxenough <- Q > minvalue
  ismax <- ismax & ismaxenough
  # find end of recessions from minima and create event data.frame
  imax <- c(which(ismax), length(Q))
  imin <- which(ismin)
  # here some magic is done in C++ to find the minima following 
  # each local maxima (if any)
  ievents <- data.frame(Start = imax, End = rec_events(imax, imin))
  # remove events with no end
  ievents <- ievents[!is.na(ievents[, 2L]), ]
  # remove events that are too short 
  eventlength <- apply(ievents, 1, diff) + 1
  ievents <- ievents[eventlength >= minduration, ]
  # truncate events that are too long
  eventlength <- apply(ievents, 1, diff) + 1
  toolong <- eventlength >= maxduration
  ievents[toolong, 2] <- ievents[toolong, 1L] + maxduration - 1
  if (return_smoothed) return(list(Q_smoothed=Q, i_events=ievents))
  # return
  return(ievents)
}

#------------------------------------------------------------------------------
#' Computes the early and late recession times
#' 
#' This function fits exponential models to individual segments of the extracted
#' recessions and compute the correseponding early and late recession times.
#' 
#' By default the function considers two segments, the 5 first days and the 
#' segment from day 15 to day 30 from the start of each recession event. This is 
#' defined in the 'segments' argument as a list containing the indices (relative
#' to the start of each recession event) of each segments. The argument 'segment'
#' can be used to define 1, 2 or more segments.
#' 
#' @param Q numeric vector. Streamflow vector.
#' @param i_events numeric data.frame. A two-column data frame containing 
#' the index of the start and end of each recession events. See the function
#' recession_events_index().
#' @param segments list. A list defining the segments of recession events. 
#' See details.
#' @param n_min integer. The minimum number of valid (non-missing) streamflow
#' values to compute a valid recssion time. If, for a given segment, the number
#' of valid streamflow value is below 'n_min', the corresponding recession time
#' will be set to NA.
#' @param summary_fun function. The summary function to use aggregate the values
#' accross all recession events. Given a numeric vector as first argument, this
#' function must return a single value. It default to the 'median' function. If set
#' to NULL, a data.frame containing all the recession times for all recssion
#' events and all segments is returned.
#' @param ... Additional argument to pass to the summary function 'summary_fun'.
#' @return Returns a vector with as many elements as the defined segments containing
#' the aggregated recession times. If 'summary_fun' is not a function, a data.frame
#' with as many column as the defined segments and as many rows as the number of
#' recession events is returned.
#' @export
recession_times <- function(Q, i_events, segments = list(1L:5L, 15L:30L), n_min = 5L, summary_fun = median, ...) {
  n_seg <- length(segments)
  tau <- matrix(NA, nrow(i_events), n_seg, dimnames = list(NULL, if (is.null(names(segments))) paste0("Tau", 1:n_seg) else names(segments)))
  Q_events <- apply(i_events, 1, function(e) Q[e[1L]:e[2L]])
  for (k in 1:n_seg) {
    s <- segments[[k]]
    Q_e <- lapply(Q_events, function(e, s)  {
      e <- e[s]
      e <- e[e>0 & !is.na(e)]
      log(e)
    }, s = s)
    valid_Q_e <- lengths(Q_e) >= n_min
    Q_e <- Q_e[valid_Q_e]
    reg_res <- lapply(Q_e, function(e) .fast_lm(1:length(e), e))
    tau[valid_Q_e, k] <- -1 / unlist(lapply(reg_res, function(e) e$coefficients[2L]), use.names = FALSE)
  }
  if (is.function(summary_fun)) apply(tau, 2, summary_fun, ...) else as.data.frame(tau)
}


#------------------------------------------------------------------------------
#' Computes the recession time and recession concavity of a recession event.
#' 
#' This function fits a power law model ((1/tau)*Q^b) to a single recession provided
#' as numeric vector. It returns the computed recession time 'tau' and concavity 'b'.
#' 
#' @param x numeric vector. Vector containing the streamflow value for
#' a single recesion event.
#' @export
rec_power_law_model_fit <- function(x) {
  x <- log(x)
  invalid <- apply(is.na(x) | is.infinite(x), 1L, any)
  x <- x[!invalid, , drop = FALSE]
  y <- .fast_lm(x[, 1L], x[, 2L])$coefficients
  y[1] <- 1 / exp(y[1L])
  y
}

#------------------------------------------------------------------------------
#' Computes the recession time and recession concavity of multiple recessions.
#' 
#' This function fits power law models to all extracted recession events
#' and compute their correseponding recession times and concavity. 
#' One can choose to return all results (for all events) in a data.frame
#' or only a summary using the 'summary_fun' argument (typically 'median')
#' 
#' @param Q numeric vector. Streamflow vector.
#' @param i_events numeric data.frame. A two-column data frame containing 
#' the index of the start and end of each recession events. See the function
#' recession_events_index().
#' @param min_length integer. The minmum length (number of days) of a 
#' recession event to actually fit the power law model. If a recession is 
#' shorter, it is ignored (i.e. NA is returned for that particular event)
#' @param discard_length: integer. Number of days to ignore at the beginning
#' of each recession event.
#' @param summary_fun function. The summary function to use to aggregate
#' the values accross all recession events. Given a numeric vector as first
#' argument, this function must return a single value. It default to the 
#' 'median' function. If set to NULL, a data.frame containing all the 
#' recession times for all recssion events and all segments is returned.
#' @param ... Additional argument to pass to the summary function 'summary_fun'.
#' @return Returns a vector with as many elements as the defined segments 
#' containing the aggregated recession times. If 'summary_fun' is not a 
#' function, a data.frame with as many column as the defined segments and
#' as many rows as the number of recession events is returned.
#' @export
recession_concavity_time <- function(Q, i_events, 
    min_length = 0, discard_length = 0,
    summary_fun = median, ...) {
  dQ  <- diff(Q)
  l_I <- apply(i_events, 1L, function(e) e[2L] - e[1L])
  v_I <- l_I > min_length
  i_events <- i_events[v_I, ]
  n_I <- nrow(i_events)
  
  p <- matrix(NA_real_, n_I, 2L, dimnames = list(NULL, c("Tau", "b")))
  for (k in 1:n_I) {
    index <- (i_events[k, 1L] + discard_length):i_events[k, 2L]
    p[k, ] <- rec_power_law_model_fit(cbind(Q[index], -dQ[index]))
  }
  if (is.function(summary_fun)) apply(p, 2L, summary_fun, ...) else as.data.frame(p)
}
