% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/baseflow.R
\name{baseflow_index}
\alias{baseflow_index}
\title{Compute the baseflow index}
\usage{
baseflow_index(Q, Qbf, na.rm = TRUE)
}
\arguments{
\item{Q}{numeric vector. Streamflow vector.}

\item{Qbf}{numeric vector. Baseflow vector.}

\item{na.rm}{logical. Should missing values be omited?}
}
\value{
A single value corresponding to the baseflow index
}
\description{
This function compute the baseflow index: the ratio between the total
baseflow volume and the total streamflow volume.
}
\details{
Warning: No checks are done on inputs except the time series length
of Q and Qbf that should match. You should have selected the proper
period before hand: both vector should be of same length and span over
the same period.
}
