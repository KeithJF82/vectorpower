library(vectorpower)

param_file=url("https://raw.githubusercontent.com/KeithJF82/vectorpower/master/inst/extdata/param_examples/model_parameters_const.txt")
age_file=url("https://raw.githubusercontent.com/KeithJF82/vectorpower/master/inst/extdata/param_examples/age_data_145.txt")
het_file=url("https://raw.githubusercontent.com/KeithJF82/vectorpower/master/inst/extdata/param_examples/het_data_9.txt")

output_location = system.file("extdata",package="vectorpower")
dataset_location = "DemoFolder1"
dataset_folder=paste(output_location,dataset_location,sep="/")

EIR_values=0.1*c(1:5)
nyears=1

data1 <- dataset_create(dataset_folder=dataset_folder,param_file=param_file,age_file=age_file,
                        het_file=het_file,EIR_values=EIR_values,nyears=nyears)

