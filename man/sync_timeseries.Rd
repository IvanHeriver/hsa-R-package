% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{sync_timeseries}
\alias{sync_timeseries}
\title{Helper to synchronize timeseries}
\usage{
sync_timeseries(...)
}
\arguments{
\item{...}{Time vectors (supporting the min/max methods)}
}
\value{
A list of logical vectors with each component being the same length as
as the length of each input time vector. A TRUE value indicates that the
corresponding time step was found to be within the time range of all the other
time vectors.
}
\description{
Find the overlapping period for time series with equal fixed time steps.
Given two or more time vectors, the function returns a list containing
logical vectors which indicate, for each input time series and each time
step whether it is included within the overlapping period.
}
\details{
If different classes types are found, the functions tries to convert the
vectors using as.Date(), so make sure to have matching classes before using
this function to disable this behavior!
}
