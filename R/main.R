#------------------------------------------------
#' @title Main population evaluation
#'
#' @description Function which takes in trial parameters for one or more populations, sends to C++ for computation 
#' of changes over time using deterministic continuum model, then collects and organizes results (EIR over time + 
#' starting immunity data for intervention computation or endpoint population data for control population, 
#' benchmark values for both)
#'
#' @details Takes a list of parameters, returns a list of raw data (data also saved to files as backup). 
#'
#' @param input_folder    Folder containing rainfall, model parameter and starting data
#' @param output_folder   Folder to send output files
#' @param n_lines         Number of lines of data in input file	
#' @param n_mv_start      Number of first set of data to use in input file
#' @param n_mv_end        Number of last set of data to use in input file (must be less than n_lines)
#' @param int_v_varied    Intervention parameter given variable value 
#'                        (0= None, 1=ATSB kill rate, 2=bednet coverage, 3=IRS coverage)
#' @param n_int_values    Number of values of varied intervention parameter to use (will be reset to 1 if int_v_varied=0)
#' @param int_v_min       Minimum intervention parameter value (sole value used if n_int_values=1; unused if int_v_varied=0)
#' @param int_v_max       Maximum intervention parameter value (unused if n_int_values=1 or int_v_varied=0)			
#' @param date_start      Day of the year when simulation starts (should be set based on input_file)		
#' @param date_int        Day of the year when intervention trial starts 
#'                        (period from date_start to date_int used to equilibrate lagged FOI data)
#' @param time_interval   Length of time (in days) between measurements of benchmarks
#' @param n_pts           Number of data points to take 
#'                        (one taken at date_int, so end time = date_int + ((n_pts-1)*time_interval)) 
#'
#' @export

mainpop <- function (input_folder = "inst/extdata/Constant/",output_folder = "inst/extdata/Constant/results_example/",
                     n_lines = 1, n_mv_set=c(0), int_v_varied = 1, n_int_values = 1, int_v_min = 0.0, 
                     int_v_max = 0.0, date_start = 0.0, date_int = 31.0, time_interval = 7.0, n_pts = 2)
  {

  na=145
  num_het=9
  n_cats=na*num_het
  n_days=((n_pts-1)*time_interval)+1
  n_mv_values=length(n_mv_set)
  n_mv_start=n_mv_set[1]
  n_mv_end=n_mv_set[n_mv_values]
  n_runs=n_mv_values*max(n_int_values,1)
  n_pts2=n_runs*n_pts
  
  parameter_file = paste(input_folder,"model_parameters.txt",sep="")
  params <- utils::read.table(parameter_file, header=TRUE)
  
  input_file = paste(input_folder,"start_data.txt",sep="")
  
  # TODO - Read in data from files. Automatically determine n_mv_values from no. lines in input_file, check if compatible
  #        with n_mv_start and n_mv_end
  # TODO - Take in list of mv values as vector so they don't need to be consecutive
  # TODO - Take in list of data collection dates as vector as well?
  
  file_benchmarks = paste(output_folder,"Benchmark_details.txt",sep="")
  file_EIRd = paste(output_folder,"EIR.txt",sep="")
  file_imm_start = paste(output_folder,"imm.txt",sep="")
  file_endpoints = paste(output_folder,"endpoints.txt",sep="")
  
  trial_params <- list(input_file=input_file, n_lines=n_lines, n_mv_set=n_mv_set,
                       n_mv_start=n_mv_start, n_mv_end=n_mv_end, int_v_varied=int_v_varied, n_int_values=n_int_values, 
                       int_v_min=int_v_min, int_v_max=int_v_max, date_start=date_start, date_int=date_int, 
                       time_interval=time_interval, n_pts=n_pts, file_benchmarks=file_benchmarks,
                       file_endpoints=file_endpoints, file_EIRd=file_EIRd, file_imm_start=file_imm_start)
  
  # Run simulation of main population
  raw_data <- rcpp_mainpop(params,trial_params)
  
  output_data <- list(EIR_output_values = raw_data$EIR_output_values,
                      slide_prev_values = matrix(data=raw_data$slide_prev_output_values,nrow=n_pts2,ncol=n_cats,byrow=TRUE), 
                      pcr_prev_values = matrix(data=raw_data$pcr_prev_output_values,nrow=n_pts2,ncol=na,byrow=TRUE),
                      clin_inc_values = matrix(data=raw_data$clin_inc_output_values,nrow=n_pts2,ncol=na,byrow=TRUE), 
                      EIR_daily_data = matrix(data=raw_data$EIR_daily_data,nrow=n_runs,ncol=n_days,byrow=TRUE),
                      IB_start_data = matrix(data=raw_data$IB_start_data,nrow=n_mv_values,ncol=n_cats,byrow=TRUE),
                      IC_start_data = matrix(data=raw_data$IC_start_data,nrow=n_mv_values,ncol=n_cats,byrow=TRUE),
                      ID_start_data = matrix(data=raw_data$ID_start_data,nrow=n_mv_values,ncol=n_cats,byrow=TRUE))
  
  return(output_data)
}

#------------------------------------------------
#' @title Cluster input setup
#'
#' @description Function for taking benchmark data output by mainpop() and selecting the desired benchmark
#'              (EIR, slide prevalence, PCR prevalence, clinical incidence) data to use for input in
#'              setting up clusters, plus creating list of intervention parameter values
#'
#' @details Takes in detailed benchmark data as a list (currently from a file) and outputs lists of chosen 
#'          benchmark values and intervention parameter values (currently as files)
#'
#' @param input_folder1       Folder containing age distribution data
#' @param input_folder2       Folder containing benchmark data
#' @param set_npt             Data point to use (1-max)
#' @param benchmark           Benchmark type to use in choosing clusters ("EIR", "slide_prev", "pcr_prev", or "clin_inc")
#' @param age_start           Starting age to use when calculating prevalence or incidence over age range (not used with EIR)
#' @param age_end             End age to use when calculating prevalence or incidence over age range (not used with EIR)
#'
#' @export

cluster_input_setup <- function(input_folder1 = "inst/extdata/Constant/",
                                input_folder2 = "inst/extdata/Constant/results_example/",
                                set_npt = 1,benchmark = "EIR",age_start = 0,age_end = 65.0){
  # TODO - Change input settings so that settings under different intervention conditions can be used
  
  age_file=paste(input_folder1,"age_data.txt",sep="")
  benchmark_file = paste(input_folder2,"Benchmark_details.txt",sep="")
  
  benchmark_data = utils::read.table(benchmark_file,header=TRUE,sep="\t")
  col_length = length(benchmark_data$run)
  n_mv_values = max(benchmark_data$n_mv)+1
  n_int_values = max(benchmark_data$int_num)+1
  n_pts = max(benchmark_data$npt)
    
  age_data = utils::read.table(age_file,header=TRUE,sep="\t")
  n_age_cats = length(age_data$age0)
  n_age_start = findInterval(age_start,age_data$age0)
  n_age_end = findInterval(age_end,age_data$age1)
  density_sum = 0
  benchmark_value = 0
  
  EIR = benchmark_data[[7]]
  slide_prev_by_age = list(rep(c(1:col_length),n_age_cats))
  pcr_prev_by_age = slide_prev_by_age
  clin_inc_by_age = slide_prev_by_age
  for(i in 1:n_age_cats){
    slide_prev_by_age[[i]] = benchmark_data[[7+i]]
    pcr_prev_by_age[[i]] = benchmark_data[[n_age_cats+7+i]]
    clin_inc_by_age[[i]] = benchmark_data[[(2*n_age_cats)+7+i]]
  }
  j = (c(0:(n_mv_values-1))*n_int_values*n_pts)+set_npt
  
  if(benchmark == "EIR"){
    benchmark_value = EIR[j]
  }else{
    if(benchmark == "slide_prev"){ benchmark_data = slide_prev_by_age}
    if(benchmark == "pcr_prev"){ benchmark_data = pcr_prev_by_age}
    if(benchmark == "clin_inc"){ benchmark_data = clin_inc_by_age }
    
    for(i in n_age_start:n_age_end){
      density_sum = density_sum + age_data$density[i]
      benchmark_value = benchmark_value + benchmark_data[[i]][j]
    }
    benchmark_value = benchmark_value/density_sum
  }
  
  utils::write.table(data.frame(n_mv=c(1:length(benchmark_value))-1,Benchmark=benchmark_value),
                     file=paste(input_folder2,"benchmark_list.txt",sep=""),sep="\t",
                     row.names=FALSE,col.names=TRUE)
  
  if(n_int_values==1){int_values=benchmark_data$param_int[1]}else{
    int_values=min(benchmark_data$param_int)+
      (c(1:(n_int_values)-1)*((max(benchmark_data$param_int)-min(benchmark_data$param_int))/(n_int_values-1)))
  }
  utils::write.table(data.frame(n_int=c(0:(n_int_values-1)),param_int=int_values),
                     file=paste(input_folder2,"int_list.txt",sep=""),sep="\t",
                     row.names=FALSE,col.names=TRUE)
}

#------------------------------------------------
#' @title Create clusters
#'
#' @description Function for creating cumulative probability distributions and generating clusters from them
#'
#' @details Takes in previously generated benchmark and intervention data (currently in files) and outputs cluster data
#'          (currently as files)
#'
#' @param input_folder        Folder containing input files
#' @param nclusters           Number of clusters to create
#' @param benchmark_mean      Mean of benchmark value distribution
#' @param benchmark_stdev     Standard deviation of benchmark value distribution
#' @param int_mean            Mean of intervention parameter value distribution
#' @param int_stdev           Standard deviation of intervention parameter value distribution
#'
#' @export

clusters_create <- function(input_folder = "inst/extdata/Constant/results_example/",nclusters=100,
                            benchmark_mean=0.25, benchmark_stdev=0.025, int_mean=0.15, int_stdev=0.05){
  
  benchmark_file=paste(input_folder,"benchmark_list.txt",sep="")
  int_file=paste(input_folder,"int_list.txt",sep="")
  
  data_bm <- utils::read.table(benchmark_file,header=TRUE,sep="\t")
  data_int <- utils::read.table(int_file,header=TRUE,sep="\t")
  nv_B=length(data_bm$n_mv)
  nv_I=length(data_int$n_int)
  graphics::par(mfrow=c(1,2))
  graphics::plot(data_bm$Benchmark,data_bm$mv_num)
  graphics::plot(data_int$param_int,data_int$int_num)
  
  nprobs=10000
  mid=nprobs/2
  sigma0=nprobs/10
  prob=rep(0,nprobs)
  cprob=prob
  v_bm1=prob
  v_bm2=prob
  v_i1=prob
  v_i2=prob
  mvn_index=prob
  int_index=prob
  for(i in 1:nprobs){
    v_bm1[i]=benchmark_mean+(((i-mid)*benchmark_stdev*10)/nprobs)
    v_i1[i]=int_mean+(((i-mid)*int_stdev*10)/nprobs)
    prob[i]=exp(-0.5*((i-mid)/sigma0)^2)
  }
  prob=prob/sum(prob)
  cprob[1]=prob[1]
  for(i in 2:nprobs){
    cprob[i]=cprob[i-1]+prob[i]
  }
  
  for(i in 1:nprobs){
    j=findInterval(v_bm1[i],data_bm$Benchmark)+1
    if(j>nv_B){j=nv_B}
    mvn_index[i]=j
    v_bm2[i]=data_bm$Benchmark[j]
    j=findInterval(v_i1[i],data_int$param_int)+1
    if(j>nv_I){j=nv_I}
    int_index[i]=j
    v_i2[i]=data_int$param_int[j]
  }
  
  clusters <- data.frame(Cluster=c(1:nclusters),CP_B=stats::runif(nclusters,0,1),n_B=rep(0,nclusters),B=rep(0,nclusters),
                         CP_I=stats::runif(nclusters,0,1),n_I=rep(0,nclusters),I=rep(0,nclusters),n_run=rep(0,nclusters))
  for(i in 1:nclusters){
    j=findInterval(clusters$CP_B[i],cprob)+1
    clusters$B[i]=v_bm2[j]
    clusters$n_B[i]=mvn_index[j]-1
    j=findInterval(clusters$CP_I[i],cprob)+1
    clusters$I[i]=v_i2[j]
    clusters$n_I[i]=int_index[j]-1
  }
  clusters$n_run=((clusters$n_B)*nv_I)+clusters$n_I
  
  graphics::par(mfrow=c(1,2))
  graphics::matplot(cprob,v_bm1,type="l",col=1)
  graphics::matplot(cprob,v_bm2,type="l",col=2,add=TRUE)
  graphics::matplot(clusters$CP_B,clusters$B,type="p",pch=1,col=3,add=TRUE)
  graphics::matplot(cprob,v_i1,type="l",col=1)
  graphics::matplot(cprob,v_i2,type="l",col=2,add=TRUE)
  graphics::matplot(clusters$CP_I,clusters$I,type="p",pch=1,col=3,add=TRUE)
  
  utils::write.table(clusters,file=paste(input_folder,"cluster_list.txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE)
  
}

#------------------------------------------------
#' @title Cohort evaluation
#'
#' @description Run stochastic individual model across trial cohort individuals in one or more clusters
#'
#' @details Takes a list of parameters, returns a list of raw data (data also saved to files as backup). 
#'
#' @param input_folder1   Folder containing model parameter data
#' @param input_folder2   Folder containing cluster data
#' @param n_mv_values     Number of mosquito density values in starting immunity data
#' @param n_EIR_values    Number of sets of EIR data
#' @param time_interval   Time interval between data points
#' @param n_divs          Number of data points (including 0)
#' @param n_clusters      Number of clusters in input data
#' @param prop_T_c        Treatment probability in trial cohort
#' @param n_patients      Number of trial cohort patients per cluster
#' @param age_c0          Minimum age of trial cohort patients
#' @param age_c1          Maximum age of trial cohort patients
#'
#' @export

cohort <- function(input_folder1 = "inst/extdata/Constant/",input_folder2 = "inst/extdata/Constant/results_example/",
                   n_mv_values = 1,n_EIR_values = 1,time_interval = 7.0,n_divs = 13,n_clusters = 1,prop_T_c = 0.9,
                   n_patients = 100,age_c0 = 0.5,age_c1 = 10.0){
  
  # TODO - 1) Reduce use of files
  #        2) Take some parameters directly from input data (n_mv_values, n_EIR_values, time_interval)
  
  parameter_file = paste(input_folder1,"model_parameters.txt",sep="")
  input_file1 = paste(input_folder2,"EIR.txt",sep="")
  input_file2 = paste(input_folder2,"imm.txt",sep="")
  input_file3 = paste(input_folder2,"cluster_list.txt",sep="")
  file_summary = paste(input_folder2,"summary.txt",sep="")
  file_frequency = paste(input_folder2,"frequency.txt",sep="")
  params <- utils::read.table(parameter_file, header=TRUE)
  
  cohort_params <- list(input_file1 = input_file1,input_file2 = input_file2,input_file3 = input_file3,
                        file_summary = file_summary,file_frequency = file_frequency,n_mv_values = n_mv_values,
                        n_EIR_values = n_EIR_values,time_interval = time_interval,n_divs = n_divs,n_clusters = n_clusters,
                        prop_T_c = prop_T_c,n_patients = n_patients,age_c0 = age_c0,age_c1 = age_c1)		
  
  raw_data <- rcpp_cohort(params,cohort_params)
  
  # placeholder; TODO - Once export of data in R list format added to cohort.cpp, add processing here
  results_data <- data.frame(raw_data)
  
  return(results_data)
}
