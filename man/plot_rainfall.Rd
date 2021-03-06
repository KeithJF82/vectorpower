% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_functions.R
\name{plot_rainfall}
\alias{plot_rainfall}
\title{Plot year-round rainfall based on parameter file}
\usage{
plot_rainfall(dataset_folder = "")
}
\arguments{
\item{dataset_folder}{Dataset folder}
}
\description{
Function for displaying year-round rainfall data based on parameter values in input file
}
\details{
Reads data from model_parameters.txt in dataset folder, calculates rainfall on each day of the year
         and outputs resulting values as a graph
}
