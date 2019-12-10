devtools::install_github("KeithJF82/vectorpower")
library(vectorpower)

## Worked example 1 - Creating a new data set

# This is an example of how to create a new dataset folder and populate it with equilibrium starting data for a range of baseline 
# malaria level values. The R code used here can also be found in the create_dataset.R file in the vectorpower/inst/ folder.
# 
# To create a new data set, it is necessary to specify the files to be used to load the parameter values governing the malaria 
# model. Three files are used: 
#   
# -param_file, which contains the main set of model parameters along with the parameters used to calculate
# rainfall as a function of time (which controls the seasonality of malaria)
# -age_file, which specifies the age intervals used for the human population
# -het_file, which specifies the parameter values controlling biting heterogeneity
# 
# In this example, the three parameter files are loaded from the online vectorpower repository. 

param_file=url("https://raw.githubusercontent.com/KeithJF82/vectorpower/master/inst/extdata/param_examples/model_parameters_Mali.txt")
age_file=url("https://raw.githubusercontent.com/KeithJF82/vectorpower/master/inst/extdata/param_examples/age_data_145.txt")
het_file=url("https://raw.githubusercontent.com/KeithJF82/vectorpower/master/inst/extdata/param_examples/het_data_9.txt")

# If you are running the example on your own computer, you may wish to change input_location to specify the location of the 
# example files in your offline copy of vectorpower, e.g.:
  
# input_folder="C:/Users/YourName/GitHub/vectorpower/inst/extdata/param_examples/"
# param_file=paste(input_folder,"model_parameters_Mali.txt",sep="")
# age_file=paste(input_folder,"age_data_145.txt",sep="")
# het_file=paste(input_folder,"het_data_9.txt",sep="")

# The location of the folder where the new data set is to be created must also be specified, along with the name of the new folder 
# to create inside that folder where the files making up the dataset will be placed. You may wish to change output_location if you 
# are running this example on your own computer.
# 
# If a folder with the same already exists and contains an existing dataset, the dataset files in that folder will be overwritten.

output_location = getwd() 
dataset_location = "/Demo Folder 1/"
dataset_folder=paste(output_location,dataset_location,sep="")

# The baseline malaria levels represented by the input data to go in the dataset are set by specifying a vector of values of 
# steady-state annual entomological inoculation rate (EIR). These values will not exactly match the year-round EIR values 
# represented by the data once the time-dependent model has been run; the final dataset will include a file specifying the 
# year-round EIR at equilibrium under the time-dependent model.
# 
# In this example, the EIR values are between 1 and 5 infectious bites per person per year; these represent low baseline malaria 
# levels. High malaria levels require longer computation times to avoid errors; the computation time increases with  approximately 
# the square root of EIR (so an EIR of 100 requires around 10 times the computation time of an EIR of 1).

EIR_values=1.0*c(1:5)

# Next, the integer number of years (from 1-20) for which the time-dependent model should be run to reach equilibrium must be 
# specified. The more years for which it is run, the more stable the equilibrium; if the equilibrium is less stable, the malaria 
# level will vary more significantly from year to year in calculations if the model is run for multiple years. For most purposes, 
# running for at least 5 years is sufficient. In this example, nyears is set to 1 for speed.

nyears=1

# With all the parameters specified, the dataset_create function is run. This will first set initial conditions using 
# steady-state calculations, then run the time-dependent model for the specified number of years to reach equilibrium. The results will be output as the dataset's starting data. The year-round data output from the starting data (assuming no new interventions) will then be calculated and output as an aid to users of the dataset.

data1 <- dataset_create(dataset_folder=dataset_folder,param_file=param_file,age_file=age_file,
                        het_file=het_file,EIR_values=EIR_values,nyears=nyears)

