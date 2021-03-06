% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_hydro_years_seasons}
\alias{get_hydro_years_seasons}
\title{Compute the hydrological year days and the seasons.}
\usage{
get_hydro_years_seasons(
  dates,
  month = 9,
  day_of_the_year = TRUE,
  seasons = list(SummerFall = 5:10, WinterSpring = c(11:12, 1:4)),
  seasons_by_day = FALSE,
  minimal = FALSE
)
}
\arguments{
\item{dates}{Date vector. A date/time vector (will be coerced to date using
as.Dates())}

\item{month}{numeric. Month number (1->12) defining the start of a
hydrological year}

\item{day_of_the_year}{logical. should the days of the (hydrological) year be
computed?}

\item{seasons}{list. List defining the seasons in month number (1:12)
(or days if 'seasons_by_day' is TRUE).}

\item{seasons_by_day}{logical. are seasons defined by days (of the year)? in
the \code{seasons} argument.}

\item{minimal}{logical. Should only the minimal desired results be returned
or intermediate results as well?}
}
\value{
a data.frame containing the hydrological year, and optionnally the
days of the hydrological year and the season which each day belongs to.
}
\description{
Compute a data.frame containing the hydrological year, and optionnally the
days of the hydrological year and the season which each day belongs to.
}
