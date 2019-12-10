# This is an example of how to create a new dataset folder and populate it with equilibrium
# starting data for a range of baseline malaria level values. It uses data drawn from the online
# vectorpower repository files.

# Input name of existing folder where new dataset folder is to be created
output_location = getwd()  

# Input name of dataset folder to be created in output_location
dataset_location = "/Demo Folder 1/"

# Input vector of Values of steady-state annual EIR to generate data for (annual EIR under 
# variable rainfall will not be identical)
EIR_values=1.0*c(1:5)

# Input integer number of years (1-20) to run main-population model for to reach equilibrium
nyears=1

devtools::install_github("KeithJF82/vectorpower")
library(vectorpower)
param_file=url("https://raw.githubusercontent.com/KeithJF82/vectorpower/master/inst/extdata/param_examples/model_parameters_Mali.txt")
age_file=url("https://raw.githubusercontent.com/KeithJF82/vectorpower/master/inst/extdata/param_examples/age_data_145.txt")
het_file=url("https://raw.githubusercontent.com/KeithJF82/vectorpower/master/inst/extdata/param_examples/het_data_9.txt")
dataset_folder=paste(output_location,dataset_location,sep="")

data1 <- dataset_create(dataset_folder=dataset_folder,param_file=param_file,age_file=age_file,
                        het_file=het_file,EIR_values=EIR_values,nyears=nyears)

