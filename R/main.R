#------------------------------------------------
#' @title Dummy function
#'
#' @description Simple test function that demonstrates some of the features of
#'   this package.
#'
#' @details Takes a vector of values, returns the square.
#'
#' @param x vector of values
#'
#' @export
#' @examples
#' # Find square of first 100 values
#' dummy1(1:100)

dummy1 <- function(x = 1:5) {
  
  # print message to console
  message("running R dummy1 function")
  
  # get arguments in list form
  args <- list(x = x)
  
  # run C++ function with these arguments
  output_raw <- dummy1_cpp(args)
  
  # some optional processing of output
  message("processing output")
  ret <- output_raw$x_squared
  
  # return
  return(ret)
}

#------------------------------------------------
#' @title Main population evaluation
#'
#' @description Function which takes in trial parameters, sends to C++ for computation, then collects and organizes results
#'
#' @details Takes a list of parameters, returns a list of raw data (data also saved to files as backup). 
#'
#' @param rain_file       File containing rainfall parameter values
#' @param parameter_file  File containing simulation parameter values
#' @param input_file      File containing mosquito and human input data
#' @param n_mv0_values    Number of lines of data in input file	
#' @param n_mv0_start     Number of first set of data to use in input file
#' @param n_mv0_end       Number of last set of data to use in input file (must be less than n_lines)
#' @param int_v_varied    Intervention parameter given variable value (0= None, 1=ATSB kill rate, 2=bednet coverage, 3=IRS coverage)
#' @param n_int_values    Number of values of varied intervention parameter to use (will be reset to 1 if int_v_varied=0)
#' @param int_v_min       Minimum intervention parameter value (sole value used if n_int_values=1; unused if int_v_varied=0)			
#' @param int_v_max       Maximum intervention parameter value (unused if n_int_values=1 or int_v_varied=0)			
#' @param date_start      Day of the year when simulation starts (should be set based on input_file)		
#' @param date_int        Day of the year when intervention trial starts (period from date_start to date_int used to equilibrate lagged FOI data)
#' @param time_interval   Length of time (in days) between measurements of benchmarks
#' @param n_pts           Number of data points to take (one taken at date_int, so end time = date_int + ((n_pts-1)*time_interval)) 
#' @param file_summary    Benchmark values at data capture points (sum over ages/heterogeneity categories where applicable)
#' @param file_benchmarks Benchmark values at data capture points (sum over heterogeneity categories where applicable)
#' @param file_endpoints  Detailed endpoint data (only used if int_v_varied=0)
#' @param file_EIRd       Daily EIR data during intervention period for each run
#' @param file_imm_start  Immunity data at start of intervention for each mv_num value
#'
#' @export
#' @examples

mainpop <- function (rain_file="Rainfall parameters.txt", parameter_file="Simulation parameters.txt",
                     input_file="Starting data.txt", n_mv0_values=600, n_mv0_start=25, n_mv0_end=50, int_v_varied=1, 
                     n_int_values=5, int_v_min=0.05, int_v_max=0.25, date_start=180.0, date_int=211.0, time_interval=7.0,
                     n_pts=13, file_summary="Benchmark_summary.txt", file_benchmarks="Benchmark_details.txt",
                     file_endpoints="endpoints.txt", file_EIRd="EIR.txt", file_imm_start="imm.txt") {
  
  trial_params <- list(rain_file=rain_file, parameter_file=parameter_file,input_file=input_file, n_mv0_values=n_mv0_values, 
                       n_mv0_start=n_mv0_start, n_mv0_end=n_mv0_end, int_v_varied=int_v_varied, n_int_values=n_int_values, 
                       int_v_min=int_v_min, int_v_max=int_v_max, date_start=date_start, date_int=date_int, 
                       time_interval=time_interval, n_pts=n_pts, file_summary=file_summary, file_benchmarks=file_benchmarks,
                       file_endpoints=file_endpoints, file_EIRd=file_EIRd, file_imm_start=file_imm_start)
  params <- utils::read.table(trial_params$parameter_file, header=TRUE)
  
  # Run simulation of main population
  raw_data <- rcpp_mainpop(params,trial_params)
  
  # placeholder
  results_data <- data.frame(raw_data)
  
  return(results_data)
}


