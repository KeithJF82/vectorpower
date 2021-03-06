% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_functions.R
\name{create_ento_start_data}
\alias{create_ento_start_data}
\title{Add entomological start data file to existing data folder}
\usage{
create_ento_start_data(dataset_folder = "")
}
\arguments{
\item{dataset_folder}{Dataset folder (must be an existing dataset folder containing model_parameters.txt
and start_data.txt files)}
}
\description{
Function which creates a starting data file for entomological modelling in existing data folder
}
\details{
Takes in the name of an existing data folder and creates start_data_ento.txt file from
             information in start_data.txt file.
}
