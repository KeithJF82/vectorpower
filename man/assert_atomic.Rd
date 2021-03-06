% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{assert_atomic}
\alias{assert_atomic}
\title{x is atomic}
\usage{
assert_atomic(
  x,
  message = "\%s must be atomic (see ?is.atomic)",
  name = deparse(substitute(x))
)
}
\arguments{
\item{x}{Value to check if atomic}

\item{message}{Message to output if x not atomic}

\item{name}{Name}
}
\description{
Check that value x is atomic
}
\details{
Check that x is atomic
}
