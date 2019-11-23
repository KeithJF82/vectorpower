# This is an example of how to create a new dataset folder and populate it with equilibrium
# starting data for a range of baseline malaria level values

# Input location of 'vectorpower' folder, e.g. "C:/Users/Name/Documents/GitHub/"
vectorpower_location = "C:/Users/Keith/Documents/Work/"  

# Input name of existing folder where new dataset folder is to be created
output_location = "C:/Users/Keith/Documents/Work/vectorpower demo/"  

# Input name of dataset folder to be created in output_location
dataset_location = "Demo Folder 1/"

# Input vector of Values of steady-state annual EIR to generate data for (annual EIR under 
# variable rainfall will not be identical)
EIR_values=1.0*c(1:5)

# Input integer number of years (1-20) to run main-population model for to reach equilibrium
nyears=1

library(devtools)
setwd(vectorpower_location)
load_all("vectorpower")
param_file=paste(vectorpower_location,"vectorpower/inst/extdata/param_examples/model_parameters_Mali.txt",sep="")
age_file=paste(vectorpower_location,"vectorpower/inst/extdata/param_examples/age_data_145.txt",sep="")
het_file=paste(vectorpower_location,"vectorpower/inst/extdata/param_examples/het_data_9.txt",sep="")
dataset_folder=paste(output_location,dataset_location,sep="")

data1 <- dataset_create(dataset_folder=dataset_folder,param_file=param_file,age_file=age_file,
                        het_file=het_file,EIR_values=EIR_values,nyears=nyears)

