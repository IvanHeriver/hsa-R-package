% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recessions_analysis.R
\name{recession_concavity_time}
\alias{recession_concavity_time}
\title{Computes the recession time and recession concavity of multiple recessions.}
\usage{
recession_concavity_time(
  Q,
  i_events,
  min_length = 0,
  discard_length = 0,
  summary_fun = median,
  ...
)
}
\arguments{
\item{Q}{numeric vector. Streamflow vector.}

\item{i_events}{numeric data.frame. A two-column data frame containing
the index of the start and end of each recession events. See the function
recession_events_index().}

\item{min_length}{integer. The minmum length (number of days) of a
recession event to actually fit the power law model. If a recession is
shorter, it is ignored (i.e. NA is returned for that particular event)}

\item{summary_fun}{function. The summary function to use to aggregate
the values accross all recession events. Given a numeric vector as first
argument, this function must return a single value. It default to the
'median' function. If set to NULL, a data.frame containing all the
recession times for all recssion events and all segments is returned.}

\item{...}{Additional argument to pass to the summary function 'summary_fun'.}

\item{discard_length:}{integer. Number of days to ignore at the beginning
of each recession event.}
}
\value{
Returns a vector with as many elements as the defined segments
containing the aggregated recession times. If 'summary_fun' is not a
function, a data.frame with as many column as the defined segments and
as many rows as the number of recession events is returned.
}
\description{
This function fits power law models to all extracted recession events
and compute their correseponding recession times and concavity.
One can choose to return all results (for all events) in a data.frame
or only a summary using the 'summary_fun' argument (typically 'median')
}
