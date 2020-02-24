#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.
#' @useDynLib vectorpower
#' @import graphics
#' @import utils
NULL
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
#' @param input_data      List containing parameters, starting data etc. created by load_inputs function
#' @param output_folder   Folder to send output files (no files saved if set to NA)
#' @param n_mv_set        Vector of mosquito density number values to use (must be increasing order)
#' @param int_v_varied    Intervention parameter given variable value 
#'                        (0= None, 1=ATSB kill rate, 2=bednet coverage, 3=IRS coverage)
#' @param int_values      Vector of values of varied intervention parameter (will be unused if int_v_varied=0)
#' @param start_interval  Period from start date before intervention starts (used to equilibrate lagged FOI data)
#' @param time_values     Vector of time checkpoints
#'
#' @export

mainpop <- function (input_data = list(),output_folder = NA,n_mv_set=c(1), int_v_varied=0, int_values= c(0.0),
                     start_interval = 31.0, time_values=c(0.0,7.0))
{
  # Input error checking
  assert_list(input_data)
  if(is.na(output_folder)==0){assert_string(output_folder)}
  assert_int(int_v_varied)
  assert_in(int_v_varied,0:3)
  assert_numeric(int_values)
  assert_single_numeric(start_interval)
  assert_numeric(time_values)
  # TODO - Change these to something that works with URLs
  # assert_file_exists(input_files$age_file)
  # assert_file_exists(input_files$het_file)
  # assert_file_exists(input_files$param_file)
  # assert_file_exists(input_files$start_file)
  # assert_file_exists(input_files$annual_file)
  
  n_mv_set=input_data$n_mv_set
  n_pts=length(time_values)
  n_days=max(time_values)+1
  n_mv_values=length(n_mv_set)
  n_mv_end=n_mv_set[n_mv_values]
  if(int_v_varied==0) { int_values=c(0.0) }
  n_int_values=length(int_values)
  na=input_data$params$na
  num_het=input_data$params$num_het
  n_cats=na*num_het
  
  # TODO - Move creation of output file names to mainpop.cpp
  if(is.na(output_folder)==FALSE){
    file_benchmarks = paste(output_folder,"Benchmark_details.txt",sep="/")
    file_EIRd = paste(output_folder,"EIR.txt",sep="/")
    file_imm_start = paste(output_folder,"imm.txt",sep="/")
    file_endpoints = paste(output_folder,"endpoints.txt",sep="/")
    flag_file=1
  } else {
    file_benchmarks = NA
    file_EIRd = NA
    file_imm_start = NA
    file_endpoints = NA
    flag_file=0
  }
  
  # Organize trial parameters into list
  trial_params <- list(age_data=input_data$age_data, n_mv_values=n_mv_values, int_v_varied=int_v_varied, int_values=int_values,
                       start_interval=start_interval, time_values=time_values, n_pts=n_pts, flag_file=flag_file,
                       file_benchmarks=file_benchmarks,file_endpoints=file_endpoints, file_EIRd=file_EIRd, 
                       file_imm_start=file_imm_start)
  
  # Run simulation of main population
  raw_data <- rcpp_mainpop(params=input_data$params,inputs=input_data$start_data,trial_params)
  
  # process raw output data
  {
    n_pt_names=paste("n_pt",c(1:n_pts),sep="")
    na_names=paste("na",c(1:na),sep="")
    n_cat_names=paste("n_cat",c(1:n_cats),sep="")
    n_days_names=paste("n_days",c(1:n_days),sep="")
    n_mv_names=paste("n_mv",c(1:n_mv_values),sep="")
    n_int_names=paste("n_int",c(1:n_int_values),sep="")
    dimnames_list=list(na_names,n_pt_names,n_int_names,n_mv_names)
    dimnames_list2=list(n_pt_names,n_int_names,n_mv_names)
    EIR_benchmarks = array(data=raw_data$EIR_benchmarks,dim=c(n_pts,n_int_values,n_mv_values),dimnames=list(n_pt_names,n_int_names,n_mv_names))
    slide_prev_benchmarks = array(data=raw_data$slide_prev_benchmarks,dim=c(na,n_pts,n_int_values,n_mv_values),dimnames=dimnames_list)
    pcr_prev_benchmarks = array(data=raw_data$pcr_prev_benchmarks,dim=c(na,n_pts,n_int_values,n_mv_values),dimnames=dimnames_list)
    clin_inc_benchmarks = array(data=raw_data$clin_inc_benchmarks,dim=c(na,n_pts,n_int_values,n_mv_values),dimnames=dimnames_list)
    M_benchmarks = array(data=raw_data$M_benchmarks,dim=c(n_pts,n_int_values,n_mv_values),dimnames=dimnames_list2)
    M_spor_benchmarks = array(data=raw_data$M_spor_benchmarks,dim=c(n_pts,n_int_values,n_mv_values),dimnames=dimnames_list2)
    EIR_daily_data = array(data=raw_data$EIR_daily_data,dim=c(n_days,n_int_values,n_mv_values),dimnames=list(n_days_names,n_int_names,n_mv_names))
    IB_start_data = array(data=raw_data$IB_start_data,dim=c(n_cats,n_mv_values),dimnames=list(n_cat_names,n_mv_names))
    IC_start_data = array(data=raw_data$IC_start_data,dim=c(n_cats,n_mv_values),dimnames=list(n_cat_names,n_mv_names))
    ID_start_data = array(data=raw_data$ID_start_data,dim=c(n_cats,n_mv_values),dimnames=list(n_cat_names,n_mv_names))
  }
  
  output_data <- list(n_mv_values=n_mv_values,n_int_values=n_int_values,n_pts=n_pts,n_mv_set=n_mv_set,
                      time_values=time_values,int_values=int_values,
                      EIR_benchmarks=EIR_benchmarks,slide_prev_benchmarks=slide_prev_benchmarks,
                      pcr_prev_benchmarks=pcr_prev_benchmarks,clin_inc_benchmarks=clin_inc_benchmarks,
                      M_benchmarks=M_benchmarks,M_spor_benchmarks=M_spor_benchmarks,
                      EIR_daily_data=EIR_daily_data,annual_data=input_data$annual_data,
                      IB_start_data=IB_start_data,IC_start_data=IC_start_data,ID_start_data=ID_start_data,
                      params=input_data$params)
  
  return(output_data)
}

#------------------------------------------------
#' @title Cluster input setup
#'
#' @description Function for taking benchmark data output by mainpop() and selecting the desired benchmark
#'              (EIR, slide prevalence, PCR prevalence, clinical incidence) data to use for input in
#'              setting up clusters, plus creating list of intervention parameter values
#'
#' @details Takes in detailed benchmark data as a list and outputs lists of chosen 
#'          benchmark values and intervention parameter values (as both list and files)
#'
#' @param input_list          List containing mainpop output data
#' @param benchmark           Benchmark type to use in choosing clusters ("EIR", "slide_prev", "pcr_prev", or "clin_inc"
#'                            for new data, "EIR_annual", "slide_prev_annual", "pcr_prev_annual", "clin_inc_annual" for
#'                            pre-existing annual data in input folder)
#' @param set_n_pt            Data point to use (1-max)
#' @param set_n_int           Intervention number to use (1-max)
#' @param age_start           Starting age to use when calculating prevalence or incidence over age range (not used with EIR)
#' @param age_end             End age to use when calculating prevalence or incidence over age range (not used with EIR)
#' @param plot_flag           Logical operator indicating whether or not to plot graph of read values
#'
#' @export

cluster_input_setup <- function(input_list=list(), benchmark = "EIR",set_n_pt = 1,set_n_int=1,age_start = 0,age_end = 65.0,
                                plot_flag=TRUE){
  
  # Input error checking
  assert_list(input_list)
  assert_in(benchmark,c("EIR","slide_prev","pcr_prev","clin_inc","EIR_annual","slide_prev_annual","pcr_prev_annual","clin_inc_annual"))
  assert_int(set_n_pt)
  assert_in(set_n_pt,c(1:input_list$n_pts))
  assert_int(set_n_int)
  assert_in(set_n_int,c(1:input_list$n_int_values))
  assert_bounded(age_start,0.0,65.0)
  assert_bounded(age_end,age_start,65.0)
  assert_logical(plot_flag)
  
  n_age_start = findInterval(age_start,input_list$params$age_years)
  n_age_end = findInterval(age_end,input_list$params$age_years)
  
  n_EIR=switch(benchmark,"EIR"=1,"EIR_annual"=1,"slide_prev"=2,"pcr_prev"=2,"clin_inc"=2,
               "slide_prev_annual"=2,"pcr_prev_annual"=2,"clin_inc_annual"=2)
  n_annual=switch(benchmark,"EIR"=1,"EIR_annual"=2,"slide_prev"=1,"pcr_prev"=1,"clin_inc"=1,
                  "slide_prev_annual"=2,"pcr_prev_annual"=2,"clin_inc_annual"=2)
  n_AE=(10*n_EIR)+n_annual
  
  benchmark_values=0
  if(n_AE==11){benchmark_values = input_list$EIR_benchmarks[set_n_pt,set_n_int,]}
  if(n_AE==12){benchmark_values = input_list$annual_data$EIRy[input_list$n_mv_set]}
  if(n_AE==21){
    density_sum = 0
    if(benchmark == "slide_prev"){ benchmark_data = input_list$slide_prev_benchmarks[,set_n_pt,set_n_int,c(1:input_list$n_mv_values)]}
    if(benchmark == "pcr_prev"){ benchmark_data = input_list$pcr_prev_benchmarks[,set_n_pt,set_n_int,c(1:input_list$n_mv_values)]}
    if(benchmark == "clin_inc"){ benchmark_data = input_list$clin_inc_benchmarks[,set_n_pt,set_n_int,c(1:input_list$n_mv_values)] }
    for(i in n_age_start:n_age_end){
      density_sum = density_sum + input_list$params$den_norm[i]
      benchmark_values = benchmark_values + benchmark_data[i,]
    }
    benchmark_values = benchmark_values/density_sum
  }
  if(n_AE==22){
    if(benchmark == "slide_prev_annual"){ benchmark_values = input_list$annual_data$slide_prev_y[input_list$n_mv_set] }
    if(benchmark == "pcr_prev_annual"){ benchmark_values = input_list$annual_data$pcr_prev_y[input_list$n_mv_set] }
    if(benchmark == "clin_inc_annual"){ benchmark_values = input_list$annual_data$clin_inc_y[input_list$n_mv_set] }
  }
  
  if(plot_flag==TRUE){
    if(n_annual==1){ matplot(c(1:input_list$n_mv_values),benchmark_values,type="p",pch=2,col=2,xlab="N_M",ylab=benchmark)}
    else{ matplot(c(1:input_list$n_mv_values),benchmark_values,type="p",pch=2,col=2,xlab="N_M",ylab=benchmark) }
  }
  
  output <- list(benchmark=benchmark,set_n_pt = set_n_pt,set_n_int=set_n_int,age_start = age_start,age_end = age_end,
                 benchmark_values=benchmark_values,int_values=input_list$int_values,n_mv_set=c(1:input_list$n_mv_values))
  
  return(output)
}

#------------------------------------------------
#' @title Create clusters
#'
#' @description Function for creating cumulative probability distributions and generating clusters from them
#'
#' @details Takes in previously generated benchmark and intervention data (currently in files) and outputs cluster data
#'          (currently as files)
#'
#' @param input_list          List containing input_data, produced by cluster_input_setup()
#' @param n_clusters          Number of clusters to create
#' @param benchmark_mean      Mean of benchmark value distribution
#' @param benchmark_stdev     Standard deviation of benchmark value distribution
#' @param int_mean            Mean of intervention parameter value distribution
#' @param int_stdev           Standard deviation of intervention parameter value distribution
#'
#' @export

clusters_create <- function(input_list=list(),n_clusters=100,benchmark_mean=0.25, benchmark_stdev=0.025,
                            int_mean=0.15, int_stdev=0.05){
  
  # Input error checking (TODO - finish)
  assert_list(input_list)
  assert_single_int(n_clusters)
  assert_single_numeric(benchmark_mean)
  assert_single_numeric(benchmark_stdev)
  assert_single_numeric(int_mean)
  assert_single_numeric(int_stdev)
  
  nv_B=length(input_list$benchmark_values)
  nv_I=length(input_list$int_values)
  
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
    if(benchmark_stdev==0.0){benchmark_stdev=benchmark_mean/100.0}
    v_bm1[i]=benchmark_mean+(((i-mid)*benchmark_stdev*10)/nprobs)
    if(int_stdev==0.0){int_stdev=int_mean/100.0}
    v_i1[i]=int_mean+(((i-mid)*int_stdev*10)/nprobs)
    prob[i]=exp(-0.5*((i-mid)/sigma0)^2)
  }
  prob=prob/sum(prob)
  cprob[1]=prob[1]
  for(i in 2:nprobs){
    cprob[i]=cprob[i-1]+prob[i]
  }
  
  for(i in 1:nprobs){
    j=findPosition(v_bm1[i],input_list$benchmark_values)
    mvn_index[i]=j
    v_bm2[i]=input_list$benchmark_values[j]
    j=findPosition(v_i1[i],input_list$int_values)
    int_index[i]=j
    v_i2[i]=input_list$int_values[j]
  }
  
  clusters <- data.frame(CP_B=stats::runif(n_clusters,0,1),n_B=rep(0,n_clusters),B=rep(0,n_clusters),
                         CP_I=stats::runif(n_clusters,0,1),n_I=rep(0,n_clusters),I=rep(0,n_clusters),n_run=rep(0,n_clusters))
  for(i in 1:n_clusters){
    j=findPosition(clusters$CP_B[i],cprob)
    clusters$B[i]=v_bm2[j]
    clusters$n_B[i]=mvn_index[j]
    j=findPosition(clusters$CP_I[i],cprob)
    clusters$I[i]=v_i2[j]
    clusters$n_I[i]=int_index[j]
  }
  clusters$n_run=((clusters$n_B-1)*nv_I)+clusters$n_I-1
  
  par(mfrow=c(1,2))
  matplot(cprob,v_bm1,type="l",col=1,xlab="Cumulative probability",ylab=input_list$benchmark,ylim=c(min(v_bm2),max(v_bm2)))
  matplot(cprob,v_bm2,type="l",col=2,add=TRUE)
  matplot(clusters$CP_B,clusters$B,type="p",pch=1,col=3,add=TRUE)
  matplot(cprob,v_i1,type="l",col=1,xlab="Cumulative probability",ylab="Intervention parameter",ylim=c(min(v_i2),max(v_i2)))
  matplot(cprob,v_i2,type="l",col=2,add=TRUE)
  matplot(clusters$CP_I,clusters$I,type="p",pch=1,col=3,add=TRUE)
  par(mfrow=c(1,1))
  
  return(clusters)
}

#------------------------------------------------
#' @title Cohort evaluation
#'
#' @description Run stochastic individual model across trial cohort individuals in one or more clusters
#'
#' @details Takes a list of parameters, returns a list of raw data (data also saved to files as backup). 
#'
#' @param mainpop_data    List output by mainpop() containing main population data
#' @param cluster_data    List output by clusters_create() containing cluster data
#' @param n_patients      Number of trial cohort patients per cluster
#' @param prop_T_c        Treatment probability in cohort
#' @param age_start       Minimum age of trial cohort patients
#' @param age_end         Maximum age of trial cohort patients
#' @param flag_output     Integer indicating whether to show progress of simulation
#'
#' @export

cohort <- function(mainpop_data = list(), cluster_data=list(), n_patients = 1,
                   prop_T_c = 0.9,age_start = 0.5,age_end = 10.0, flag_output=0){
  
  # Input error checking (TODO - finish)
  assert_list(mainpop_data)
  assert_list(cluster_data)
  assert_single_int(n_patients)
  assert_single_bounded(prop_T_c,0.0,1.0)
  assert_single_bounded(age_start,0.0,65.0)
  assert_single_bounded(age_end,age_start,65.0)
  assert_single_int(flag_output)
  assert_in(flag_output,c(0,1))
  
  n_clusters=length(cluster_data$n_run)
  time_values=mainpop_data$time_values
  trial_params <- list(n_patients = n_patients,n_clusters=n_clusters, flag_output = flag_output,
                       time_values=time_values,tmax_i=length(mainpop_data$EIR_daily_data[,1,1]),
                       prop_T_c = prop_T_c,age_start = age_start,age_end = age_end,
                       EIR_daily_data = as.vector(mainpop_data$EIR_daily_data),
                       IB_start_data = as.vector(mainpop_data$IB_start_data),
                       IC_start_data = as.vector(mainpop_data$IC_start_data),
                       ID_start_data = as.vector(mainpop_data$ID_start_data))
  
  raw_data <- rcpp_cohort(mainpop_data$params,trial_params,cluster_data)
  results_data <- data.frame(raw_data)
  
  cluster_names=paste("cluster",c(1:n_clusters),sep="/")
  patient_names=paste("patient",c(1:n_patients),sep="/")
  n_pts=length(mainpop_data$time_values)
  n_pt_names=paste("n_pt",c(1:n_pts),sep="/")
  patients_age_outputs = array(data=raw_data$patients_age_outputs,dim=c(n_patients,n_clusters),
                               dimnames=list(patient_names,cluster_names))
  patients_het_outputs = array(data=raw_data$patients_het_outputs,dim=c(n_patients,n_clusters),
                               dimnames=list(patient_names,cluster_names))
  patients_status_outputs = array(data=raw_data$patients_status_outputs,dim=c(n_pts,n_patients,n_clusters),
                                  dimnames=list(n_pt_names,patient_names,cluster_names))
  p_det_outputs = array(data=raw_data$p_det_outputs,dim=c(n_pts,n_patients,n_clusters),
                                  dimnames=list(n_pt_names,patient_names,cluster_names))
  clin_inc_outputs = array(data=raw_data$clin_inc_outputs,dim=c(n_pts,n_clusters),
                        dimnames=list(n_pt_names,cluster_names))
  
  results=list(n_clusters=n_clusters,n_patients=n_patients,n_pts=n_pts,time_values=time_values,
               patients_age_outputs=patients_age_outputs,patients_het_outputs=patients_het_outputs,
               patients_status_outputs=patients_status_outputs,p_det_outputs=p_det_outputs,clin_inc_outputs=clin_inc_outputs)
  
  return(results)
}
