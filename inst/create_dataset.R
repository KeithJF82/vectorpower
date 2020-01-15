# This is an example of how to create a new dataset folder and populate it with equilibrium
# starting data for a range of baseline malaria level values

# Input location of parameter files
input_location = "C:/Users/kjfras16/Documents/0 - Git repositories/vectorpower/inst/extdata/param_examples/"  

# Input name of existing folder where new dataset folder is to be created
output_location = "C:/Users/kjfras16/Documents/0 - Git repositories/vectorpower/inst/extdata/param_examples/"  

# Input name of dataset folder to be created in output_location
dataset_location = "C:/Users/kjfras16/Documents/0 - Git repositories/vectorpower/inst/extdata/DemoFolder1/"

# Input vector of Values of steady-state annual EIR to generate data for (annual EIR under 
# variable rainfall will not be identical)
EIR_values=0.1*c(1:5)

# Input integer number of years (1-20) to run main-population model for to reach equilibrium
nyears=1

library(devtools)
setwd(input_location)
load_all("vectorpower")
param_file=paste(input_location,"model_parameters_const.txt",sep="")
age_file=paste(input_location,"age_data_145.txt",sep="")
het_file=paste(input_location,"het_data_9.txt",sep="")

data1 <- dataset_create(dataset_folder=dataset_location,param_file=param_file,age_file=age_file,
                        het_file=het_file,EIR_values=EIR_values,nyears=nyears)

