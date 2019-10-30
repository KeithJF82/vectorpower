#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.
#' @useDynLib vectorpower
#' @import graphics
#' @import utils
NULL
#------------------------------------------------
# TODO - Create function for setting up or adjusting parameter and age data values in R?
#------------------------------------------------
#' @title Create new data set
#'
#' @description Function which sets up a folder containing new data in the correct formats which can be used for
#'              input into mainpop().
#'
#' @details     Takes in a folder name, designated parameter files and a vector of target values of annual entomological 
#'              inoculation rate (EIR) under constant rainfall conditions, creates a folder containing parameter, starting data
#'              and annual EIR data files.
#'
#' @param dataset_folder    Folder where new dataset to be located. Will be created if not already in existence
#' @param EIR_values        Vector of target values of annual entomological inoculation rate under constant rainfall conditions
#' @param param_file        File containing model parameters to use for new dataset
#' @param age_file          File containing age parameters to use for new dataset
#' @param nyears            Number of years to run model for to achieve equilibrium
#'
#' @export

dataset_create <- function (dataset_folder="",EIR_values=c(1.0),param_file="",age_file="",nyears=10)
{
  # Error checking (TODO - Finish)
  assert_string(dataset_folder)
  assert_numeric(EIR_values)
  assert_string(param_file) # TODO - Change to check that file exists
  assert_string(age_file)   # TODO - Change to check that file exists
  
  #
  if(dir.exists(dataset_folder)==FALSE){dir.create(dataset_folder)}
  
  # Copy parameter files from existing folder to new one
  cat("\nCreating new model_parameters.txt and age_data.txt files.")
  param_file_new=paste(dataset_folder,"model_parameters.txt",sep="")
  age_file_new=paste(dataset_folder,"age_data.txt",sep="")
  file.copy(from=param_file,to=param_file_new, overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)
  file.copy(from=age_file,to=age_file_new, overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)
  
  # Create new steady-state constant-rainfall data for desired EIR [TODO: or mv0] values using parameter [TODO: and age] files
  n_mv_values=length(EIR_values)
  file_endpoints = paste(dataset_folder,"start_data.txt",sep="")
  params <- read.table(paste(dataset_folder,"model_parameters.txt",sep=""), header=TRUE) # Read in model parameters
  cat("\nGenerating preliminary data (steady state at constant rainfall)")
  trial_params <- list(EIRy=EIR_values,file_endpoints=file_endpoints)
  raw_data <- rcpp_mainpop_ss(params,trial_params)
  mv_values <- raw_data$mv_values
  
  # Run model without interventions from steady-state data for nyears years in one go to achieve equilibrium
  cat("\nBeginning main population calculations to equilibrium.")
  n_mv_set=c(1:n_mv_values)
  output_folder=paste(dataset_folder,"Temp/",sep="")
  if(dir.exists(output_folder)==FALSE){dir.create(output_folder)}
  eq_data <- mainpop(input_folder = dataset_folder, output_folder = output_folder,n_mv_set = n_mv_set, int_v_varied = 0, 
                     int_values=c(0.0), start_interval = 0.0, time_values=365*c(0:nyears) )
  plot_mainpop_data(input_list=eq_data,set_n_int=1,benchmark = "EIR",age_start = 0,age_end = 65.0)
  
  # Transfer endpoints as new start data
  cat("\n\nCreating new starting_data.txt file from results.")
  file.copy(from=paste(output_folder,"endpoints.txt",sep=""),to=paste(dataset_folder,"start_data.txt",sep=""),
            overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)
  
  # Run for 1 year with daily data points to produce year-round data
  cat("\nRunning main population model for 1 year to establish year-round values\n(annual EIR and incidence, annual average prevalences).")
  annual_data <- mainpop(input_folder = dataset_folder,n_mv_set = n_mv_set, int_v_varied = 0, int_values=c(0.0),
                         start_interval = 0.0, time_values=1.0*c(0:364) )
  plot_mainpop_data(input_list=annual_data,set_n_int=1,benchmark = "EIR")
  
  # Output annual data for reference
  cat("\n\nOutputting annual data.\n")
  annual_EIR_file=paste(dataset_folder,"annual_data.txt",sep="")
  EIRy_values = rep(0,n_mv_values)
  slide_prevy_values=EIRy_values
  clin_incy_values=EIRy_values
  pcr_prevy_values=EIRy_values
  cat("n_mv\tmv0\tEIRy_ss\tEIRy\tslide_prev_y\tclin_inc_y\tpcr_prev_y\n",file=annual_EIR_file)
  for(i in n_mv_set){
    if(i==1){cat("n_mv\tmosq_density\tEIRy_ss\tEIRy\tslide_prev_y\tclin_inc_y\tpcr_prev_y\n")}
    EIRy_values[i]=sum(annual_data$EIR_benchmarks[,i])
    slide_prevy_values[i]=sum(annual_data$slide_prev_benchmarks[,,i])/365
    clin_incy_values[i]=sum(annual_data$clin_inc_benchmarks[,,i])
    pcr_prevy_values[i]=sum(annual_data$pcr_prev_benchmarks[,,i])/365
    cat(i,mv_values[i],EIR_values[i],EIRy_values[i],slide_prevy_values[i],clin_incy_values[i],pcr_prevy_values[i],"\n",sep="\t")
    cat(i,mv_values[i],EIR_values[i],EIRy_values[i],slide_prevy_values[i],clin_incy_values[i],pcr_prevy_values[i],"\n",
        file=annual_EIR_file,append=TRUE,sep="\t")
  }
  
  # TODO - Delete Temp folder
  
  return(annual_data)
}
