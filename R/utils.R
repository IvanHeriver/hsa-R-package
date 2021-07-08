#------------------------------------------------------------------------------
#' Helper to synchronize timeseries
#' 
#' Find the overlapping period for time series with equal fixed time steps. 
#' Given two or more time vectors, the function returns a list containing
#' logical vectors which indicate, for each input time series and each time 
#' step whether it is included within the overlapping period.
#' 
#' If different classes types are found, the functions tries to convert the 
#' vectors using as.Date(), so make sure to have matching classes before using
#' this function to disable this behavior!
#' 
#' @param ... Time vectors (supporting the min/max methods)
#' @return A list of logical vectors with each component being the same length as
#' as the length of each input time vector. A TRUE value indicates that the
#' corresponding time step was found to be within the time range of all the other
#' time vectors.
#' @export
sync_timeseries <- function(...) {
  args <- list(...)
  n <- length(args)
  if (n == 1) return(rep(TRUE, length(args[[1]]))) # this doesn't make sens
  m <- vapply(args, length, numeric(1L))
  
  # deal with different time classes
  ts_class <- lapply(args, class)
  ts_first_class <- ts_class[[1]][1]
  different_classes <- FALSE
  for (k in 2:n) if (!identical(ts_class[[k]], ts_first_class)) different_classes <- TRUE
  if (different_classes) {
    warning("Vectors have different classes! They were all converted to 'Date' class")
    # try and convert everything to Date class
    args <- lapply(args, as.Date)
  }
  
  # get min/max of overlapping periods
  time_min <- max(vapply(args, min, numeric(1L)))
  time_max <- min(vapply(args, max, numeric(1L)))
  
  # check if time series do overlap
  if (time_max <=  time_min) warning("Vectors don't overlap!")
  
  # compute "overlap" logical vectors
  out <- list()
  for (k in 1:n) out[[k]] <- args[[k]] >= time_min & args[[k]] <= time_max
  
  return(out)
}

#------------------------------------------------------------------------------
#' Compute the hydrological year days and the seasons.
#' 
#' Compute a data.frame containing the hydrological year, and optionnally the
#' days of the hydrological year and the season which each day belongs to.
#' 
#' @param dates Date vector. A date/time vector (will be coerced to date using
#' as.Dates())
#' @param month numeric. Month number (1->12) defining the start of a 
#' hydrological year
#' @param day_of_the_year logical. should the days of the (hydrological) year be
#' computed?
#' @param seasons list. List defining the seasons in month number (1:12)
#' (or days if 'seasons_by_day' is TRUE).
#' @param seasons_by_day logical. are seasons defined by days (of the year)? in
#' the \code{seasons} argument.
#' @param minimal logical. Should only the minimal desired results be returned 
#' or intermediate results as well?
#' @return a data.frame containing the hydrological year, and optionnally the
#' days of the hydrological year and the season which each day belongs to.
#' @export
get_hydro_years_seasons <- function(dates, month = 9, day_of_the_year = TRUE, 
                                seasons = list("SummerFall" = 5:10, "WinterSpring" = c(11:12, 1:4)),
                                seasons_by_day = FALSE, minimal = FALSE) {
  dates <- as.Date(dates)
  
  # get hydrological year
  m <- as.numeric(format(dates, format = "%m"))
  hy <- y <- as.numeric(format(dates, format = "%Y"))
  m_prevy <- !m%in%c(month:12)
  hy[m_prevy] <- hy[m_prevy] - 1
  
  # get days of the year
  if (day_of_the_year) {
    j <- as.numeric(format(dates, format = "%j"))
    start_hy <- as.Date(paste0(y, "-", month, "-1"))
    start_y <- as.Date(paste0(y, "-1-1"))
    j_hy <- as.numeric(format(dates - start_hy + start_y, format = "%j"))
  } else {
    j_hy <- NA
  }
  
  # get seasons
  if (!is.null(seasons) && is.list(seasons)) {
    s_names <- names(seasons)
    s_names <- factor(s_names, levels = s_names)
    s_table <- data.frame(i = unname(unlist(seasons)), season = rep(s_names, unlist(lapply(seasons, length))))
    if (seasons_by_day && day_of_the_year){
      s <- s_table[match(j, s_table[, 1]), 2]
    } else {
      s <- s_table[match(m, s_table[, 1]), 2]
    }
  } else {
    s <- NA
  }

  # return resulting data.frame
  if (minimal) return(data.frame(hy = hy, j_hy = j_hy, s = s))
  data.frame(dates = dates, m = m, y = y, hy = hy, j = j, j_hy = j_hy, s = s)
}

#------------------------------------------------------------------------------
#' Compute the valid hydrological year
#' 
#' Given a hydrological year vector (a vector indicating the hydrological year
#' each time step belongs to) and a matrix with the same number
#' of rows (and any number of column), this function returns a logical
#' vector of the same length as the hydrological year vector that indicates
#' which are the time steps that are part of a valid hydrological year.
#' Validity of a hydrological year is assessed according  to the following
#' rules:
#' - is the year complete (i.e. is there at least \code{n} (e.g. 365) days)?
#' - is there less than \code{na.th} (proportion) missing values in the 
#' hydrological year?
#' 
#' @param hy numeric vector. a hydrological year vector (will be coerced to factor)
#' @param x numeric matrix. a matrix (or a vector) with as many row (elements) 
#' as the length of \code{hy} used to look for missing values
#' @param n integer. the minimum length of a complete hydrological year this is 
#' also use to compute the proportion of missing values.
#' @param na.th numeric. the minimal tolerated proportion of missing values
#' @return a logicalvector of the same length as the hydrological year vector
#' that indicates whether or not the time steps are part of a valid hydrological
#' year.
#' @export
get_valid_hydro_year <- function(hy, x = NULL, n = 365, na.th = 0.05) {
  # only full year
  rle_res <- rle(hy)
  unique_hy<- rle_res$values
  valid_hy_1 <- rep(rle_res$lengths >= n, rle_res$lengths)
 
  if (!is.null(x))  {
    # if x is provided, only year where the percentage of missing value is below na.th
    if (is.null(dim(x))) x <- matrix(x, length(x), 1)
    na <- apply(apply(x, 2, is.na), 1, any)
    valid_hy_2 <- rep(as.vector(tapply(na, hy, sum) / n <= na.th), rle_res$lengths)
    valid_hy_1 & valid_hy_2
  } else {
    valid_hy_1
  }
}

#------------------------------------------------------------------------------
#' What is within?
#' 
#' Given two vectors (typically date vectors), return a logical vector of the 
#' size of the first vector indicating whether or not the times/dates were 
#' found in the second vector
#' 
#' @param OriginalTime vector. first time/date vector
#' @param FilteredTime vector. second time/date vector
#' @return logical vector of the size of the
#' first vector indicating whether or not the times/dates were found in the
#' second vector
#' @export
get_what_is_within <- function(OriginalTime, FilteredTime) OriginalTime%in%FilteredTime

# FIXME: these functions could also be used for: PQ slopes, BFRmag!
#------------------------------------------------------------------------------
#' Compute a regime curve
#' 
#' Compute regime curves from the data provided in 'x'. It first apply 
#' smoothing (rolling mean) with a wondow 'n' of the data provided in 'x'. 
#' Then, it compute the inter-annual mean of each calendar day (provided in
#' vector 'hdays'). It recommended to ignore the last day using 
#' 'remove_last = TRUE'.
#' 
#' @param hdays numeric vector. vector of 'julian' date of length m which 
#' should be computed according to hydrological year (see function 
#' \code{get_hydro_years_seasons()})
#' @param x numeric vector. a vector (of length m) or a matrix (with m rows)
#' @param n integer. smoothing window
#' @param remove_last logical. should the last day of the year (366) be 
#' removed (recommended).
#' @return A 366 (remove_last = FALSE) or 365 (remove_last = TRUE) long 
#' numeric vector (or matrix with 366 or 365 rows) containing the regime curve(s) 
#' @export
regime_curve <- function(hdays, x, n = 1L, na.rm = TRUE, remove_last = TRUE) {
  if (is.null(dim(x))) x <- matrix(x, dimnames = list(NULL, as.character(substitute(x))))
  if (n > 1L) x <- apply(x, 2L, roll_mean, n = n, fill = NA)
  y <- data.frame(hdays, x)
  y <- dplyr::group_by(y, hdays)
  y <- dplyr::summarise_all(y, mean, na.rm = na.rm)
  y <- dplyr::data.frame(y)
  if (remove_last) return(y[-nrow(y), ]) else return(y)
}

#------------------------------------------------------------------------------
#' Compute cumulative curves on a data.frame of matrix
#' 
#' Given a vector or a matrix/data.frame 'x', compute the corresponding
#' cumulative curve. This function is only a wrapper around 
#' \code{cumsum_interpolate()} with a call to \code{apply()}.
#' 
#' @param x numeric vector, data.frame or matrix. a vector (of length m) or
#' a matrix/data.frame (with m rows)
#' @param na.action character. a character string indicating whether missing 
#' values should be interpolated or not (see \code{cumsum_interpolate()})
#' @return A data.frame with 'm' rows containing the cumulative version of 
#' 'x'.
#' @export
cumulative_curves <- function(x, na.action = "interpolate") {
  if (is.null(dim(x))) x <- matrix(x)
  y <- apply(x, 2L, cumsum_interpolate, na.action = na.action)
  as.data.frame(y)
}

# TODO: implement an option to verify the 'validity' of the interpolation
#       For example, it should not be longer than 'n' time step or I could
#       provide a value indicating how much the interpolated value weight
#       on the resulting cumulative curve
#------------------------------------------------------------------------------
#' Wrapper around the \code{cumsum} function where missing values can be 
#' interpolated.
#' 
#' @param x vector to use to compute cumulated values
#' @param na.action action to undertake when missing values are found (default:
#' 'interpolate', any other value leads to the original \code{cumsum} function)
#' @return the cumulative sum of x.
#' @export
cumsum_interpolate <- function(x, na.action = "interpolate") {
  isna <- is.na(x)
  if (all(isna)) stop("'x' cannot be only missing values!")
  if (any(isna)) {
    if (na.action == "interpolate") {
      y <- 1:length(x)
      x[isna] <- approx(x = y[!isna], y = x[!isna], xout = y[isna], method = "linear")$y
    } else if (is.numeric(na.action)) {
      x[isna] <- na.action
    }
    isna <- is.na(x)
    if (any(isna)) {
      warning("'x' has missing values at the start/end! They were treated as zeros.")
      x[isna] <- 0
      y <- cumsum(x)
      y[isna] <- NA
      return(y)
    }
  }
  cumsum(x)
}



# This function is here only to avoid code repetition (BFI, RC, etc.)
#------------------------------------------------------------------------------
#' Compute the ratio of two sums.
#' 
#' Compute the ratio between the sum of 'top' and the sum of 'bottom'.
#' 
#' @param top numeric vector. Its sum will be the nominator of the ratio
#' @param bottom numeric vector. Its sum will be the denominator of the ratio
#' @param na.rm logical. Should missing values be omited?
#' @return A single numeric value corresponding to the ratio of the sum of
#' 'top' and the sum of 'bottom'
#' @export
sum_ratio <- function(top, bottom, na.rm = TRUE) sum(top, na.rm = na.rm) / sum(bottom, na.rm = na.rm)


.fast_lm <- function(x, y) {
  x <- cbind(1, x)
  .lm.fit(x, y)
}
