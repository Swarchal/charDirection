% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chdir.R
\name{check_datasets}
\alias{check_datasets}
\title{sanity check for datasets}
\usage{
check_datasets(ctrl, expm)
}
\arguments{
\item{ctrl}{matrix of control data}

\item{expm}{matrix of experiment data}
}
\value{
error message if error, otherwise silent
}
\description{
Will stop and return error message if rows have constant variance,
or if they have differing numbers of rows, or contain NA values.
}

