% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pq_approach.R
\name{pq_curve_slopes}
\alias{pq_curve_slopes}
\title{P-Q signatures computations using P-Q cumulative curve}
\usage{
pq_curve_slopes(
  PQ,
  start = 15,
  end = 183,
  bp = mean(c(start, end)),
  intercept = TRUE
)
}
\arguments{
\item{PQ}{numeric vector. The cumulative precipitation minus
cumulative streamflow (P-Q) vector.}

\item{start}{numeric. The day number of the start of the period to
search for a change of trend in P-Q.}

\item{end}{numeric. The day number of the end of the period to
search for a change of trend in P-Q.}

\item{bp}{numeric value. The day number of the initial guess of threshold
date}

\item{intercept}{logical. Should the intercept be estimated (default: TRUE)
or fixed to c(0, 0) (FALSE)?}
}
\value{
A named vector of length 6 containing the breakpoint date ('bp'),
breakpoint strength signature ('bp_strength'), the slope of the first (dry)
period ('slp_dry') and its corresponding intercept ('b_dry'), the slope of
the second (wet) period ('slp_wet') and its corresponding intercept
('b_wet')
}
\description{
Compute the seasonal response change signatures: breackpoint date,
first period slope (dry), second period slope (wet) and the intercepts
associated with these two slopes (only there for plotting purposes).
}
\details{
This function uses a segmented regression to find the dry and wet slopes.
It also returns the intercepts of the two lines (which are not signatures)
that can be used to plot the lines.
}
