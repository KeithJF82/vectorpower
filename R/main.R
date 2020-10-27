#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.
#' @useDynLib vectorpower
#' @import graphics
#' @import utils
#' @import tidyverse
#' @import geepack
#' @import pbapply
#' @import stats
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
#' @param int_v_varied    Intervention parameter given variable value 
#'                        (0= None, 1=ATSB kill rate, 2=bednet coverage, 3=IRS coverage)
#' @param int_values      Vector of values of varied intervention parameter (will be unused if int_v_varied=0)
#' @param start_interval  Period from start date before intervention starts (used to equilibrate lagged FOI data)
#' @param time_values     Vector of time checkpoints
#' @param flag_dt_adjust  Integer indicating whether or not to adjust dt by mosquito density in rcpp_mainpop
#'                        (0 = No, 1 = Yes)
#'
#' @export

mainpop <- function (input_data = list(),output_folder = NA,int_v_varied=0, int_values= c(0.0),
                     start_interval = 31.0, time_values=c(0.0,7.0), flag_dt_adjust=1)
{
  # Input error checking
  assert_list(input_data)
  if(is.na(output_folder)==0){assert_string(output_folder)}
  assert_int(int_v_varied)
  assert_in(int_v_varied,0:4)
  assert_numeric(int_values)
  assert_single_numeric(start_interval)
  assert_numeric(time_values)
  assert_in(flag_dt_adjust,0:1)
  
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
  trial_params <- list(age_data=input_data$age_data, n_mv_values=n_mv_values, int_v_varied=int_v_varied, 
                       int_values=int_values,start_interval=start_interval, time_values=time_values, n_pts=n_pts, 
                       flag_file=flag_file,file_benchmarks=file_benchmarks,file_endpoints=file_endpoints, 
                       file_EIRd=file_EIRd, file_imm_start=file_imm_start,flag_dt_adjust=flag_dt_adjust)
  
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
    EIR_benchmarks = array(data=raw_data$EIR_benchmarks,dim=c(n_pts,n_int_values,n_mv_values),
                           dimnames=list(n_pt_names,n_int_names,n_mv_names))
    slide_prev_benchmarks = array(data=raw_data$slide_prev_benchmarks,dim=c(na,n_pts,n_int_values,n_mv_values),
                                  dimnames=dimnames_list)
    pcr_prev_benchmarks = array(data=raw_data$pcr_prev_benchmarks,dim=c(na,n_pts,n_int_values,n_mv_values),
                                dimnames=dimnames_list)
    clin_inc_benchmarks = array(data=raw_data$clin_inc_benchmarks,dim=c(na,n_pts,n_int_values,n_mv_values),
                                dimnames=dimnames_list)
    M_benchmarks = array(data=raw_data$M_benchmarks,dim=c(n_pts,n_int_values,n_mv_values),dimnames=dimnames_list2)
    M_spor_benchmarks = array(data=raw_data$M_spor_benchmarks,dim=c(n_pts,n_int_values,n_mv_values),
                              dimnames=dimnames_list2)
    EIR_daily_data = array(data=raw_data$EIR_daily_data,dim=c(n_days,n_int_values,n_mv_values),
                           dimnames=list(n_days_names,n_int_names,n_mv_names))
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
#' @param benchmark           Benchmark to use in choosing clusters ("EIR", "slide_prev", "pcr_prev", or "clin_inc"
#'                            for new data, "EIR_annual", "slide_prev_annual", "pcr_prev_annual", "clin_inc_annual" 
#'                            for pre-existing annual data in input folder)
#' @param set_n_pt            Data point to use (1-max)
#' @param set_n_int           Intervention number to use (1-max)
#' @param age_start           Starting age to use when calculating prevalence or incidence over age range 
#'                            (not used with EIR)
#' @param age_end             End age to use when calculating prevalence or incidence over age range 
#'                            (not used with EIR)
#' @param plot_flag           Logical operator indicating whether or not to plot graph of read values
#'
#' @export

cluster_input_setup <- function(input_list=list(), benchmark = "EIR",set_n_pt = 1,set_n_int=1,age_start = 0,
                                age_end = 65.0,plot_flag=TRUE){
  
  # Input error checking
  assert_list(input_list)
  assert_in(benchmark,c("EIR","slide_prev","pcr_prev","clin_inc","EIR_annual","slide_prev_annual","pcr_prev_annual",
                        "clin_inc_annual"))
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
    if(benchmark=="slide_prev"){ benchmark_data = input_list$slide_prev_benchmarks[,set_n_pt,set_n_int,
                                                                                     c(1:input_list$n_mv_values)]}
    if(benchmark=="pcr_prev"){ benchmark_data = input_list$pcr_prev_benchmarks[,set_n_pt,set_n_int,
                                                                                 c(1:input_list$n_mv_values)]}
    if(benchmark=="clin_inc"){ benchmark_data = input_list$clin_inc_benchmarks[,set_n_pt,set_n_int,
                                                                                 c(1:input_list$n_mv_values)] }
    for(i in n_age_start:n_age_end){
      density_sum = density_sum + input_list$params$den_norm[i]
      benchmark_values = benchmark_values + benchmark_data[i,]
    }
    benchmark_values = benchmark_values/density_sum
  }
  if(n_AE==22){
    if(benchmark=="slide_prev_annual"){ benchmark_values = input_list$annual_data$slide_prev_y[input_list$n_mv_set] }
    if(benchmark=="pcr_prev_annual"){ benchmark_values = input_list$annual_data$pcr_prev_y[input_list$n_mv_set] }
    if(benchmark=="clin_inc_annual"){ benchmark_values = input_list$annual_data$clin_inc_y[input_list$n_mv_set] }
  }
  
  if(plot_flag==TRUE){
    if(n_annual==1){ 
      matplot(c(1:input_list$n_mv_values),benchmark_values,type="p",pch=2,col=2,xlab="N_M",ylab=benchmark)}
    else{ matplot(c(1:input_list$n_mv_values),benchmark_values,type="p",pch=2,col=2,xlab="N_M",ylab=benchmark) }
  }
  
  output <- list(benchmark=benchmark,set_n_pt = set_n_pt,set_n_int=set_n_int,age_start = age_start,age_end = age_end,
                 benchmark_values=benchmark_values,int_values=input_list$int_values,
                 n_mv_set=c(1:input_list$n_mv_values))
  
  return(output)
}

#------------------------------------------------
#' @title Create clusters
#'
#' @description Function for creating cumulative probability distributions and generating clusters from them
#'
#' @details Takes in previously generated benchmark and intervention data (currently in files) 
#'          and outputs cluster data (currently as files)
#'
#' @param input_list          List containing input_data, produced by cluster_input_setup()
#' @param n_clusters          Number of clusters to create
#' @param benchmark_mean      Mean of benchmark value distribution
#' @param benchmark_stdev     Standard deviation of benchmark value distribution
#' @param int_mean            Mean of intervention parameter value distribution
#' @param int_stdev           Standard deviation of intervention parameter value distribution
#' @param plot_flag           True/false flag indicating whether to plot graphs of cluster data
#'
#' @export

clusters_create <- function(input_list=list(),n_clusters=100,benchmark_mean=0.25, benchmark_stdev=0.025,
                            int_mean=0.15, int_stdev=0.05, plot_flag=FALSE){
  
  # Input error checking (TODO - finish)
  assert_list(input_list)
  assert_single_int(n_clusters)
  assert_single_numeric(benchmark_mean)
  assert_single_numeric(benchmark_stdev)
  assert_single_numeric(int_mean)
  assert_single_numeric(int_stdev)
  assert_logical(plot_flag)
  
  nv_B=length(input_list$benchmark_values)
  nv_I=length(input_list$int_values)
  
  nprobs=10000
  mid=nprobs*0.5
  sigma0=nprobs*0.1
  prob=rep(0,nprobs)
  cprob=prob
  v_bm1=prob
  v_bm2=prob
  v_i1=prob
  v_i2=prob
  mvn_index=prob
  int_index=prob
  for(i in 1:nprobs){
    if(benchmark_stdev==0.0){benchmark_stdev=benchmark_mean*0.01}
    v_bm1[i]=benchmark_mean+(((i-mid)*benchmark_stdev*10)/nprobs)
    if(int_stdev==0.0){int_stdev=int_mean*0.01}
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
  
  clusters <- data.frame(CP_B=runif(n_clusters,0,1),n_B=rep(0,n_clusters),B=rep(0,n_clusters),
                         CP_I=runif(n_clusters,0,1),n_I=rep(0,n_clusters),I=rep(0,n_clusters),
                         n_run=rep(0,n_clusters))
  for(i in 1:n_clusters){
    j=findPosition(clusters$CP_B[i],cprob)
    clusters$B[i]=v_bm2[j]
    clusters$n_B[i]=mvn_index[j]
    j=findPosition(clusters$CP_I[i],cprob)
    clusters$I[i]=v_i2[j]
    clusters$n_I[i]=int_index[j]
  }
  clusters$n_run=((clusters$n_B-1)*nv_I)+clusters$n_I-1
  
  if(plot_flag==TRUE){
    par(mfrow=c(1,2))
    matplot(cprob,v_bm1,type="l",col=1,xlab="Cumulative probability",ylab=input_list$benchmark,
            ylim=c(min(v_bm2),max(v_bm2)))
    matplot(cprob,v_bm2,type="l",col=2,add=TRUE)
    matplot(clusters$CP_B,clusters$B,type="p",pch=1,col=3,add=TRUE)
    matplot(cprob,v_i1,type="l",col=1,xlab="Cumulative probability",ylab="Intervention parameter",
            ylim=c(min(v_i2),max(v_i2)))
    matplot(cprob,v_i2,type="l",col=2,add=TRUE)
    matplot(clusters$CP_I,clusters$I,type="p",pch=1,col=3,add=TRUE)
    par(mfrow=c(1,1))
  }
  
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
  
  cluster_names=paste("cluster",c(1:n_clusters),sep="")
  patient_names=paste("patient",c(1:n_patients),sep="")
  n_pts=length(mainpop_data$time_values)
  n_pt_names=paste("n_pt",c(1:n_pts),sep="")
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
               patients_status_outputs=patients_status_outputs,
               p_det_outputs=p_det_outputs,clin_inc_outputs=clin_inc_outputs)
  
  return(results)
}


#------------------------------------------------
#' @title Cohort evaluation (alternate)
#'
#' @description Run stochastic individual model across trial cohort individuals in one or more clusters
#'              (alternate version incorporating more details of testing and reactive treatment during trial,
#'              in progress)
#'
#' @details Takes a list of parameters, returns a list of raw data (data also saved to files as backup). 
#'
#' @param mainpop_data      List output by mainpop() containing main population data
#' @param cluster_data      List output by clusters_create() containing cluster data
#' @param n_patients        Number of trial cohort patients per cluster
#' @param prop_T_c          Treatment probability in cohort
#' @param age_start         Minimum age of trial cohort patients
#' @param age_end           Maximum age of trial cohort patients
#' @param flag_output       Integer indicating whether to show progress of simulation
#' @param test_time_values  Time points at which tests are administered
#' @param test_type         Type of test administered ("clin" = clinical, "RDT" = rapid diagnostic test)
#' @param flag_pre_clearing Integer indicating whether patients given pre-trial prophylaxis (if 1, place all patients
#'                          into prophylaxis category at start)
#' @param censor_period     Time period after a positive test during which a patient is not counted towards incidence
#' @param flag_reactive_treatment Integer indicating whether patients are automatically given prophylaxis after a
#'                                positive test (shifting them into treatment category if a clinical case,
#'                                prophylaxis category otherwise)
#'
#' @export

cohort2 <- function(mainpop_data = list(), cluster_data=list(), n_patients = 1, prop_T_c = 0.9,
                   age_start = 0.5,age_end = 10.0, flag_output=0, test_time_values = c(), test_type = "clin",
                   flag_pre_clearing = 1, censor_period = 0.0, flag_reactive_treatment = 1){
  
  # Input error checking (TODO - finish)
  assert_list(mainpop_data)
  assert_list(cluster_data)
  assert_single_int(n_patients)
  assert_single_bounded(prop_T_c,0.0,1.0)
  assert_single_bounded(age_start,0.0,65.0)
  assert_single_bounded(age_end,age_start,65.0)
  assert_single_int(flag_output)
  assert_in(flag_output,c(0,1))
  assert_in(test_type, c("clin","RDT"))
  assert_single_bounded(min(test_time_values),min(mainpop_data$time_values),max(mainpop_data$time_values))
  assert_single_bounded(max(test_time_values),min(mainpop_data$time_values),max(mainpop_data$time_values))
  assert_single_int(flag_pre_clearing)
  assert_in(flag_pre_clearing,c(0,1))
  assert_single_numeric(censor_period)
  assert_single_int(flag_reactive_treatment)
  assert_in(flag_reactive_treatment,c(0,1))
  
  if(test_type=="clin"){flag_test_type=1}
  if(test_type=="RDT"){flag_test_type=2}
  n_clusters=length(cluster_data$n_run)
  time_values=mainpop_data$time_values
  trial_params <- list(n_patients = n_patients,n_clusters = n_clusters, flag_output = flag_output,
                       flag_pre_clearing = flag_pre_clearing, flag_reactive_treatment = flag_reactive_treatment,
                       time_values = time_values, test_time_values = test_time_values, 
                       flag_test_type = flag_test_type, censor_period = censor_period, 
                       tmax_i=length(mainpop_data$EIR_daily_data[,1,1]), prop_T_c = prop_T_c,
                       age_start = age_start,age_end = age_end,
                       EIR_daily_data = as.vector(mainpop_data$EIR_daily_data),
                       IB_start_data = as.vector(mainpop_data$IB_start_data),
                       IC_start_data = as.vector(mainpop_data$IC_start_data),
                       ID_start_data = as.vector(mainpop_data$ID_start_data))
  
  raw_data <- rcpp_cohort2(mainpop_data$params,trial_params,cluster_data)
  results_data <- data.frame(raw_data)
  
  cluster_names=paste("cluster",c(1:n_clusters),sep="")
  patient_names=paste("patient",c(1:n_patients),sep="")
  n_pts=length(test_time_values)
  n_pt_names=paste("n_pt",c(1:n_pts),sep="")
  patients_age_outputs = array(data=raw_data$patients_age_outputs,dim=c(n_patients,n_clusters),
                               dimnames=list(patient_names,cluster_names))
  patients_het_outputs = array(data=raw_data$patients_het_outputs,dim=c(n_patients,n_clusters),
                               dimnames=list(patient_names,cluster_names))
  patients_status_outputs = array(data=raw_data$patients_status_outputs,dim=c(n_pts,n_patients,n_clusters),
                                  dimnames=list(n_pt_names,patient_names,cluster_names))
  patients_test_outputs = array(data=raw_data$patients_test_outputs,dim=c(n_pts,n_patients,n_clusters),
                                dimnames=list(n_pt_names,patient_names,cluster_names))
  p_det_outputs = array(data=raw_data$p_det_outputs,dim=c(n_pts,n_patients,n_clusters),
                        dimnames=list(n_pt_names,patient_names,cluster_names))
  test_incidence_outputs=array(data=results_data$test_incidence_outputs,dim=c(length(test_time_values),n_clusters))
  
  results <- list(n_clusters=n_clusters,n_patients=n_patients,n_pts=n_pts,time_values=time_values,
                  test_time_values=test_time_values,patients_age_outputs=patients_age_outputs,
                  patients_het_outputs=patients_het_outputs,test_incidence_outputs=test_incidence_outputs,
                  patients_status_outputs=patients_status_outputs,patients_test_outputs=patients_test_outputs)
  
  return(results)
}
#------------------------------------------------
#' @title Entomological-only calculations
#'
#' @description Function which takes in trial parameters for one or more mosquito populations, sends to C++ for 
#' computation of changes over time using deterministic continuum model, then collects and organizes results (TBA)
#'
#' @details Takes a list of parameters, returns a list of raw data (data also saved to files as backup). 
#'
#' @param input_folder    Folder containing parameter files
#' @param output_folder   Folder to send output files (no files saved if set to NA)
#' @param int_v_varied    Intervention parameter given variable value 
#'                        (0= None, 1=ATSB kill rate, 2=bednet coverage, 3=IRS coverage)
#' @param int_values      Vector of values of varied intervention parameter (will be unused if int_v_varied=0)
#' @param n_mv_set        Vector of mosquito density number values to use (must be increasing order)
#' @param start_interval  Period from start date before intervention starts (used to equilibrate lagged FOI data)
#' @param time_values     Vector of time checkpoints
#' @param flag_dt_adjust  Integer indicating whether or not to adjust dt by mosquito density in rcpp_mainpop
#'                        (0 = No, 1 = Yes)
#'
#' @export

ento <- function (input_folder = "",output_folder = NA,int_v_varied=0, int_values= c(0.0), n_mv_set=c(),
                     start_interval = 31.0, time_values=c(0.0,7.0), flag_dt_adjust=1)
{
  # Input error checking
  param_file=file=paste(input_folder,"model_parameters.txt",sep="/")
  assert_file_exists(param_file)
  start_data_file=file=paste(input_folder,"start_data_ento.txt",sep="/")
  assert_file_exists(start_data_file)
  if(is.na(output_folder)==0){assert_string(output_folder)}
  assert_int(int_v_varied)
  assert_in(int_v_varied,0:4)
  assert_numeric(int_values)
  assert_single_numeric(start_interval)
  assert_numeric(time_values)
  assert_in(flag_dt_adjust,0:1)
  
  n_pts=length(time_values)
  n_days=max(time_values)+1
  n_mv_values=length(n_mv_set)
  n_mv_end=n_mv_set[n_mv_values]
  if(int_v_varied==0) { int_values=c(0.0) }
  n_int_values=length(int_values)
  na=50 # TODO - Load from file
  
  params <- as.list(read.table(param_file, header=TRUE))
  start_data1 <- read.table(start_data_file,header=TRUE,nrows=n_mv_end)
  start_data2 <- list(mv_input=start_data1$mv0[n_mv_set], 
                      EL_input=start_data1$EL[n_mv_set], LL_input=start_data1$LL[n_mv_set], 
                      PL_input=start_data1$PL[n_mv_set], Sv_input=c(), Ev_input=c(), Iv_input=c())
  for(i in 1:na){
    start_data2$Sv_input <- append(start_data2$Sv_input,start_data1[[5+i]][n_mv_set])
    start_data2$Ev_input <- append(start_data2$Ev_input,start_data1[[5+i+na]][n_mv_set])
    start_data2$Iv_input <- append(start_data2$Iv_input,start_data1[[5+i+(2*na)]][n_mv_set])
  }
  
  if(is.na(output_folder)==FALSE){
    file_endpoints = paste(output_folder,"endpoints.txt",sep="/")
    flag_file=1
  } else {
    file_endpoints = NA
    flag_file=0
  }
  
  # Organize trial parameters into list
  trial_params <- list(n_mv_values=n_mv_values, int_v_varied=int_v_varied, int_values=int_values,
                       start_interval=start_interval, time_values=time_values, n_pts=n_pts, flag_file=flag_file,
                       file_endpoints=file_endpoints, flag_dt_adjust=flag_dt_adjust)
  
  # Run simulation of main population
  raw_data <- rcpp_ento(params=params,inputs=start_data2,trial_params)
  
  # process raw output data
  {
    n_pt_names=paste("n_pt",c(1:n_pts),sep="")
    na_names=paste("na",c(1:na),sep="")
    n_days_names=paste("n_days",c(1:n_days),sep="")
    n_mv_names=paste("n_mv",c(1:n_mv_values),sep="")
    n_int_names=paste("n_int",c(1:n_int_values),sep="")
    dimnames_list=list(na_names,n_pt_names,n_int_names,n_mv_names)
    dimnames_list2=list(n_pt_names,n_int_names,n_mv_names)
    M_benchmarks = array(data=raw_data$M_benchmarks,dim=c(n_pts,n_int_values,n_mv_values),
                         dimnames=dimnames_list2)
    M_spor_benchmarks = array(data=raw_data$M_spor_benchmarks,dim=c(n_pts,n_int_values,n_mv_values),
                              dimnames=dimnames_list2)
    M_3g_benchmarks = array(data=raw_data$M_3g_benchmarks,dim=c(n_pts,n_int_values,n_mv_values),
                            dimnames=dimnames_list2)
  }
  
  output_data <- list(n_mv_values=n_mv_values,n_int_values=n_int_values,n_pts=n_pts,n_mv_set=n_mv_set,
                      time_values=time_values,int_values=int_values,
                      M_benchmarks=M_benchmarks,M_spor_benchmarks=M_spor_benchmarks,M_3g_benchmarks=M_3g_benchmarks,
                      params=params)
  
  return(output_data)
}

#------------------------------------------------
#' @title Combined control/intervention cluster trial run
#'
#' @description Function which takes in trial parameters for main population modelling results and both control and 
#'              intervention clusters, and runs
#'
#' @details Takes a list of parameters, returns output data as list 
#'
# @param output_file_cluster Name of file to send cluster data; NA if data not to be saved
# @param output_file_indiv   Name of file to send individual data; NA if data not to be saved
#' @param mainpop_data       List containing main population data output by mainpop()
#' @param n_clusters         Number of clusters (control and intervention)
#' @param n_patients         Number of patients per cluster
#' @param benchmark          Benchmark type to use in choosing clusters ("EIR", "slide_prev", "pcr_prev", or 
#'                           "clin_inc" for new data, "EIR_annual", "slide_prev_annual", "pcr_prev_annual", 
#'                           "clin_inc_annual" for pre-existing annual data in input folder)
#' @param age_start          Starting age to use when calculating prevalence or incidence over age range 
#'                           (not used with EIR)
#' @param age_end            End age to use when calculating prevalence or incidence over age range 
#'                           (not used with EIR)     
#' @param test_time_values   Time points at which tests are administered
#' @param test_type          Type of test administered ("clin" = clinical, "RDT" = rapid diagnostic test)
#' @param flag_pre_clearing  Integer indicating whether patients given pre-trial prophylaxis (if 1, place all patients
#'                           into prophylaxis category at start)
#' @param censor_period      Time period after a positive test during which a patient is not counted towards incidence
#' @param flag_reactive_treatment Integer indicating whether patients are automatically given prophylaxis after a
#'                                positive test (shifting them into treatment category if a clinical case,
#'                                prophylaxis category otherwise)
#' @param prop_T_c            Treatment probability
#' @param benchmark_mean      Mean of benchmark value distribution
#' @param benchmark_stdev     Standard deviation of benchmark value distribution
#' @param int_mean            Mean of intervention parameter value distribution
#' @param int_stdev           Standard deviation of intervention parameter value distribution
#' @param data_level           Data type to output ("Cluster", "Individual" or "Both")
#'
#' @export

crt_combined <- function(#output_file_cluster=NA,output_file_indiv=NA,
                         mainpop_data=list(),n_clusters=1,n_patients=1,
                         benchmark="EIR_annual",age_start=0.0,age_end=65.0,test_time_values=c(0.0,1.0),
                         test_type="RDT",flag_pre_clearing=0,censor_period=0.0, flag_reactive_treatment=1,
                         prop_T_c=0.9,benchmark_mean=0.0,benchmark_stdev=0.0,int_mean=0.0,int_stdev=0.0,
                         data_level="Both"){
  
  assert_list(mainpop_data)
  assert_single_int(mainpop_data$n_mv_values)
  assert_single_bounded(max(test_time_values),0,max(mainpop_data$time_values))
  assert_in(data_level,c("Cluster","Individual","Both"))
  
  n_mv_values=mainpop_data$n_mv_values
  n_int_values=mainpop_data$n_int_values
  setup_list <- cluster_input_setup(input_list=mainpop_data, set_n_pt=1, set_n_int=1, benchmark = benchmark, 
                                    plot_flag=FALSE)
  n_pts=length(test_time_values)
  
  
  output <- list()
  
  # Create list of control clusters
  output$cluster_list_con <- clusters_create(input_list=setup_list,n_clusters=n_clusters,
                                             benchmark_mean=benchmark_mean,benchmark_stdev=benchmark_stdev, 
                                             int_mean=0.0, int_stdev=0.0)
  
  # Create list of intervention clusters
  output$cluster_list_int <- clusters_create(input_list=setup_list,n_clusters=n_clusters, 
                                             benchmark_mean=benchmark_mean, benchmark_stdev=benchmark_stdev, 
                                             int_mean=int_mean, int_stdev=int_stdev)
  
  # Simulate control clusters
  output$cohort_data_con <- cohort2(mainpop_data=mainpop_data,cluster_data=output$cluster_list_con, 
                                    n_patients = n_patients,prop_T_c = prop_T_c, 
                                    age_start = age_start, age_end = age_end, 
                                    test_time_values = test_time_values, test_type = test_type, 
                                    flag_pre_clearing = flag_pre_clearing, censor_period = censor_period, 
                                    flag_reactive_treatment = flag_reactive_treatment)
  
  # Simulate intervention clusters
  output$cohort_data_int <- cohort2(mainpop_data=mainpop_data,cluster_data=output$cluster_list_int, 
                                    n_patients = n_patients,prop_T_c = prop_T_c, 
                                    age_start = age_start, age_end = age_end, 
                                    test_time_values = test_time_values, test_type = test_type, 
                                    flag_pre_clearing = flag_pre_clearing, censor_period = censor_period, 
                                    flag_reactive_treatment = flag_reactive_treatment)
  
  if(data_level=="Cluster" || data_level=="Both"){
    fill_c=rep(0,2*n_clusters*n_pts)
    Date_fill_c=rep(test_time_values,2*n_clusters)
    Trt_fill_c=c(rep(0,n_clusters*n_pts),rep(1,n_clusters*n_pts)) 
    output$output_cluster <- data.frame(Date=Date_fill_c,Cluster=fill_c,Trt=Trt_fill_c,Numerator=fill_c,
                                             Denominator=fill_c,Incidence=fill_c)
  }
  if(data_level=="Individual" || data_level=="Both"){
    n_lines=2*n_clusters*n_patients*n_pts
    fill=rep(0,n_lines)
    output$output_indiv <- data.frame(Date=fill,Patient=fill,Cluster=fill,Trt=fill,Age=fill,
                                      Het_cat=fill,Status=fill,Test_status=fill)
    output$output_indiv$Date=rep(test_time_values,2*n_clusters*n_patients)
    output$output_indiv$Trt=c(rep(0,n_clusters*n_patients*n_pts),rep(1,n_clusters*n_patients*n_pts))
  }
  
  # Output cluster level data
  if(data_level=="Cluster" || data_level=="Both"){
    line=0
    for(i in 1:n_clusters){
      for(j in 1:n_pts){
        inc=output$cohort_data_con$test_incidence_outputs[j,i]
        if(j==1){ denominator=n_patients } else {denominator=n_patients-numerator}
        numerator=inc*denominator
        line=line+1
        output$output_cluster$Cluster[line]=i
        output$output_cluster$Numerator[line]=numerator
        output$output_cluster$Denominator[line]=denominator
        output$output_cluster$Incidence[line]=inc
      }
    }
    for(i in 1:n_clusters){
      i2=i+n_clusters
      for(j in 1:n_pts){
        inc=output$cohort_data_int$test_incidence_outputs[j,i]
        if(j==1){ denominator=n_patients } else {denominator=n_patients-numerator}
        numerator=inc*denominator
        line=line+1
        output$output_cluster$Cluster[line]=i
        output$output_cluster$Numerator[line]=numerator
        output$output_cluster$Denominator[line]=denominator
        output$output_cluster$Incidence[line]=inc
      }
    } 
  }
  
  if(data_level=="Individual" || data_level=="Both"){
    output$output_indiv$Status=c(as.vector(output$cohort_data_con$patients_status_outputs),
                                 as.vector(output$cohort_data_int$patients_status_outputs))
    output$output_indiv$Test_status=c(as.vector(output$cohort_data_con$patients_test_outputs),
                                      as.vector(output$cohort_data_int$patients_test_outputs))
    lines1=c(1:(n_patients*n_pts))
    lines2=c(1:n_pts)
    for(i in 1:n_clusters){
      output$output_indiv$Cluster[lines1]=rep(i,n_patients*n_pts)
      lines1=lines1+(n_patients*n_pts)
      for(j in 1:n_patients){
        age=output$cohort_data_con$patients_age_outputs[j,i]
        het_cat=output$cohort_data_con$patients_het_outputs[j,i]
        output$output_indiv$Patient[lines2]=rep(j,n_pts)
        output$output_indiv$Age[lines2]=rep(age,n_pts)
        output$output_indiv$Het_cat[lines2]=rep(het_cat,n_pts)
        lines2=lines2+n_pts
      }
    }
    for(i in 1:n_clusters){
      output$output_indiv$Cluster[lines1]=rep(i,n_patients*n_pts)
      lines1=lines1+(n_patients*n_pts)
      for(j in 1:n_patients){
        age=output$cohort_data_int$patients_age_outputs[j,i]
        het_cat=output$cohort_data_int$patients_het_outputs[j,i]
        output$output_indiv$Patient[lines2]=rep(j,n_pts)
        output$output_indiv$Age[lines2]=rep(age,n_pts)
        output$output_indiv$Het_cat[lines2]=rep(het_cat,n_pts)
        lines2=lines2+n_pts
      }
    } 
  }
  
  return(output)
}

#------------------------------------------------
#' @title Power calculation alt
#'
#' @description Function which takes in trial parameters for main population modelling results and both control and 
#'              intervention clusters, computes data using single instance of crt_combined and calculates power
#'
#' @details Takes a list of parameters, returns inf_boolen. 
#'
#' @param mainpop_data       List containing main population data output by mainpop()
#' @param n_clusters         Number of clusters (control and intervention)
#' @param n_patients         Number of patients per cluster
#' @param data_level         Data level to use ("Cluster" or "Individual")
#' @param alpha              Significance level
#' @param effect_size_sign   Sign of effect to check for ("Positive" or "Negative")
#' @param trial_type         Trial type (superiority, equivalence, non-inferiority - "sup", "eq", "ni")
#' @param delta_eq           Delta for equivalence trials
#' @param delta_ni           Delta for non-inferiority trials
#' @param n_sims             Number of iterations
#' @param benchmark          Benchmark type to use in choosing clusters ("EIR", "slide_prev", "pcr_prev", or 
#'                           "clin_inc" for new data, "EIR_annual", "slide_prev_annual", "pcr_prev_annual", 
#'                           "clin_inc_annual" for pre-existing annual data in input folder)
#' @param age_start          Starting age to use when calculating prevalence or incidence over age range 
#'                           (not used with EIR)
#' @param age_end            End age to use when calculating prevalence or incidence over age range 
#'                           (not used with EIR)     
#' @param test_time_values   Time points at which tests are administered
#' @param test_type          Type of test administered ("clin" = clinical, "RDT" = rapid diagnostic test)
#' @param flag_pre_clearing  Integer indicating whether patients given pre-trial prophylaxis (if 1, place all patients
#'                           into prophylaxis category at start)
#' @param censor_period      Time period after a positive test during which a patient is not counted towards incidence
#' @param flag_reactive_treatment Integer indicating whether patients are automatically given prophylaxis after a
#'                                positive test (shifting them into treatment category if a clinical case,
#'                                prophylaxis category otherwise)
#' @param prop_T_c            Treatment probability
#' @param benchmark_mean      Mean of benchmark value distribution
#' @param benchmark_stdev     Standard deviation of benchmark value distribution
#' @param int_mean            Mean of intervention parameter value distribution
#' @param int_stdev           Standard deviation of intervention parameter value distribution
#'
#' @export

power_compute <- 
  function(mainpop_data=list(),n_clusters=1,n_patients=1,data_level="Cluster",
           alpha=0.05,effect_size_sign="Positive",trial_type="sup",delta_eq=NULL,delta_ni=NULL,n_sims=1,
           benchmark="EIR_annual",age_start=0.0,age_end=65.0,test_time_values=c(0.0,1.0),
           test_type="RDT",flag_pre_clearing=0,censor_period=0.0, flag_reactive_treatment=1,prop_T_c=0.9,
           benchmark_mean=0.0,benchmark_stdev=0.0,int_mean=0.0,int_stdev=0.0)
  {
    
    assert_in(data_level,c("Cluster","Individual"))
    n_clusters_total=n_sims*n_clusters
    if(is.null(alpha)) {alpha <- 0.05}
    if(is.null(n_sims)){ n_sims <- 1000}
    if(is.null(trial_type)) {trial_type="sup"}
    if(trial_type=="ni")
    {
      if(is.null(delta_ni)) {stop("A value for delta for Non-inferiority trials MUST be specified")}
      if(delta_ni > 0 & effect_size_sign=="Positive") {stop("delta_ni and effect_size_sign must not have same sign")}
      if(delta_ni < 0 & effect_size_sign=="Negative") {stop("delta_ni and effect_size_sign must not have same sign")}
    }
    
    results <- list()
    
    ######### Generate data for use
    
    output <- 
      crt_combined(mainpop_data=mainpop_data,n_clusters=n_clusters_total,n_patients=n_patients,benchmark=benchmark,
                   age_start=age_start,age_end=age_end,test_time_values=test_time_values,
                   test_type=test_type,flag_pre_clearing=flag_pre_clearing,censor_period=censor_period,
                   flag_reactive_treatment=flag_reactive_treatment,prop_T_c=prop_T_c,benchmark_mean=benchmark_mean,
                   benchmark_stdev=benchmark_stdev,int_mean=int_mean,int_stdev=int_stdev,data_level=data_level)
    
    if(data_level=="Cluster"){ 
      data_full <- output$output_cluster 
    } else { 
      data_full <- output$output_indiv }
    
    
    ######### Ensure that treatment variable is specified first in formula
    
    data_full$Trt <- as.factor(data_full$Trt)
    data_full <- within(data_full, Trt<- relevel(Trt, ref = "1"))
    
    if(data_level=="Cluster"){
      Outcome <- "cbind(Numerator, Denominator)" 
      covariates=c("Trt","Date") 
      n_lines=n_clusters*length(test_time_values)
    } else {
      Outcome <- "Test_status"  
      covariates=c("Trt","Date","Age")
      n_lines=n_clusters*n_patients*length(test_time_values)
    }
    
    # This returns the formula:
    gee_formula <- as.formula(paste(Outcome,paste(covariates,collapse=" + "),sep=" ~ "))
    
    Cluster="Cluster"
    Patient="Patient"
    n_lines2=n_lines*n_sims
    inf_boolen=rep(NA,n_sims)
    
    for(sim in 1:n_sims){
      
      #Lines in data_full from which to take data
      values=c((((sim-1)*n_lines)+1):(sim*n_lines),
               (((sim-1)*n_lines)+1+n_lines2):((sim*n_lines)+n_lines2)) 
      if(data_level=="Cluster"){
        data<-data.frame(Date=data_full$Date[values],Cluster=data_full$Cluster[values],Trt=data_full$Trt[values],
                         Numerator=data_full$Numerator[values],Denominator=data_full$Denominator[values],
                         Incidence=data_full$Incidence[values])
        geepack_inf <- geeglm(formula=gee_formula,data = data,id = Cluster,
                              family = binomial(link="logit"),corstr = "exchangeable")
      } else {
        data<-data.frame(Date=data_full$Date[values],Patient=data_full$Patient[values],Cluster=data_full$Cluster[values],
                         Trt=data_full$Trt[values],Age=data_full$Age[values],Het_cat=data_full$Het_cat[values],
                         Status=data_full$Status[values],Test_status=data_full$Test_status[values])
        geepack_inf <- geeglm(formula=gee_formula,data = data,id = interaction(Cluster,Patient), 
                              family = binomial(link="logit"),corstr = "exchangeable")
      }
      
      ###### Estimate CI
      
      geepack_inf_solution <- coef(summary(geepack_inf))
      
      ci_table <- with(as.data.frame(geepack_inf_solution),
                       cbind(lwr= Estimate - qnorm(1-(0.5*alpha))*Std.err,
                             upr= Estimate + qnorm(1-(0.5*alpha))*Std.err))
      
      params_soln <- cbind.data.frame(geepack_inf_solution,ci_table)
      
      ################ SUPERIORITY TRIAL
      if(trial_type=="sup")
      {
        ###### Given direction of effect size, determine significance
        if (effect_size_sign=="Positive"){inf_boolen[sim] <- ifelse(as.vector(params_soln["Trt","lwr"]) > 0, 1, 0)}
        else{inf_boolen[sim] <- ifelse(as.vector(params_soln["Trt","upr"]) < 0, 1, 0)}
      }
      
      ################ EQUIVALENCE TRIAL
      if(trial_type=="eq")
      {
        ###### Regardless of  direction of effect size, equvalence will have CI containing zero and within abs(delta)
        inf_boolen[sim] <- ifelse(((as.vector(params_soln["Trt","lwr"])< 0 & 
                                      as.vector(params_soln["Trt","lwr"]) > - abs(delta_eq)) &
                                     (as.vector(params_soln["Trt","upr"]) > 0 & 
                                        as.vector(params_soln["Trt","upr"]) < + abs(delta_eq))),
                                  1,0)
      }
      
      ################ NON-INFERIORITY TRIAL
      if(trial_type=="ni")
      {
        ###### Determine the direction of effect size then define boundaries appropriately
        
        if(effect_size_sign=="Negative")
        {
          inf_boolen[sim] <- ifelse(((as.vector(params_soln["Trt","lwr"]) < delta_ni) & 
                                       (as.vector(params_soln["Trt","lwr"]) < 0 )) &
                                      ((as.vector(params_soln["Trt","upr"]) < delta_ni) & 
                                         (as.vector(params_soln["Trt","upr"]) > 0 )),
                                    1,0)
        } else {
          inf_boolen[sim] <- ifelse(((as.vector(params_soln["Trt","lwr"]) > delta_ni) & 
                                       (as.vector(params_soln["Trt","lwr"]) < 0 )) &
                                      ((as.vector(params_soln["Trt","upr"]) > delta_ni) & 
                                         (as.vector(params_soln["Trt","upr"]) > 0 )),
                                    1,0)
        }
      }
      
    }
    
    results$cluster_list_con = output$cluster_list_con
    results$cluster_list_int = output$cluster_list_int
    results$cohort_data_con = output$cohort_data_con
    results$cohort_data_int = output$cohort_data_int
    results$power_estimate=mean(inf_boolen)
    return(results)
  }
