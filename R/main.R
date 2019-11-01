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
#' @param input_folder    Folder containing rainfall, model parameter and starting data
#' @param output_folder   Folder to send output files (no files saved if set to NA)
#' @param n_mv_set        Vector of mosquito density number values to use (must be increasing order)
#' @param int_v_varied    Intervention parameter given variable value 
#'                        (0= None, 1=ATSB kill rate, 2=bednet coverage, 3=IRS coverage)
#' @param int_values      Vector of values of varied intervention parameter (will be unused if int_v_varied=0)
#' @param start_interval  Period from start date before intervention starts (used to equilibrate lagged FOI data)
#' @param time_values     Vector of time checkpoints
#'
#' @export

mainpop <- function (input_folder = "inst/extdata/Constant/",output_folder = NA,n_mv_set=c(1), 
                     int_v_varied=0, int_values= c(0.0),
                     start_interval = 31.0, time_values=c(0.0,7.0))
  {
  # Input error checking (TODO: finish)
  assert_string(input_folder)
  assert_in(int_v_varied,0:3)
  assert_numeric(n_mv_set)
  assert_numeric(int_values)
  assert_numeric(time_values)

  # Set up parameters 
  age_data=age_data_setup(read.table(paste(input_folder,"age_data.txt",sep=""),header=TRUE,sep="\t")[[1]])
  het_data = as.list(read.table(paste(input_folder,"het_data.txt",sep=""),header=TRUE,sep="\t"))  # Read in biting heterogeneity data
  params <- as.list(read.table(paste(input_folder,"model_parameters.txt",sep=""), header=TRUE))   # Read in model parameters
  na=length(age_data$age_width)
  num_het=length(het_data$het_x)
  params=c(na=na,num_het=num_het,params,age_data,het_data)
  n_pts=length(time_values)
  n_cats=na*num_het
  n_days=max(time_values)+1
  n_mv_values=length(n_mv_set)
  n_mv_end=n_mv_set[n_mv_values]
  if(int_v_varied==0) { int_values=c(0.0) }
  n_int_values=length(int_values)
  n_runs=n_mv_values*n_int_values
  
  # Read in data from input files
  input_file = paste(input_folder,"start_data.txt",sep="")
  input_data <- read.table(input_file,header=TRUE,nrows=n_mv_end)
  inputs <- list(mv_input=input_data$mv0[n_mv_set], 
                 EL_input=input_data$EL[n_mv_set], LL_input=input_data$LL[n_mv_set], PL_input=input_data$PL[n_mv_set],
                 Sv_input=input_data$Sv1[n_mv_set], Ev_input=input_data$Ev1[n_mv_set], Iv_input=input_data$Iv1[n_mv_set],
                 S_input=c(), T_input=c(), D_input=c(), A_input=c(), U_input=c(), P_input=c(),
                 ICA_input=c(), ICM_input=c(), IB_input=c(), ID_input=c())
  for(i in 1:n_cats){
    inputs$S_input <- append(inputs$S_input,input_data[[8+i]][n_mv_set])
    inputs$T_input <- append(inputs$T_input,input_data[[8+i+n_cats]][n_mv_set])
    inputs$D_input <- append(inputs$D_input,input_data[[8+i+(2*n_cats)]][n_mv_set])
    inputs$A_input <- append(inputs$A_input,input_data[[8+i+(3*n_cats)]][n_mv_set])
    inputs$U_input <- append(inputs$U_input,input_data[[8+i+(4*n_cats)]][n_mv_set])
    inputs$P_input <- append(inputs$P_input,input_data[[8+i+(5*n_cats)]][n_mv_set])
    inputs$ICA_input <- append(inputs$ICA_input,input_data[[8+i+(6*n_cats)]][n_mv_set])
    inputs$ICM_input <- append(inputs$ICM_input,input_data[[8+i+(7*n_cats)]][n_mv_set])
    inputs$IB_input <- append(inputs$IB_input,input_data[[8+i+(8*n_cats)]][n_mv_set])
    inputs$ID_input <- append(inputs$ID_input,input_data[[8+i+(9*n_cats)]][n_mv_set])
  }
  
  # TODO - Move the creation of the file names to mainpop.cpp
  if(is.na(output_folder)==FALSE){
    file_benchmarks = paste(output_folder,"Benchmark_details.txt",sep="")
    file_EIRd = paste(output_folder,"EIR.txt",sep="")
    file_imm_start = paste(output_folder,"imm.txt",sep="")
    file_endpoints = paste(output_folder,"endpoints.txt",sep="")
    flag_file=1
  } else {
    file_benchmarks = NA
    file_EIRd = NA
    file_imm_start = NA
    file_endpoints = NA
    flag_file=0
  }
  
  # Organize trial parameters into list
  trial_params <- list(age_data=age_data, n_mv_values=n_mv_values, int_v_varied=int_v_varied, int_values=int_values,
                       start_interval=start_interval, time_values=time_values, n_pts=n_pts, flag_file=flag_file,
                       file_benchmarks=file_benchmarks,file_endpoints=file_endpoints, file_EIRd=file_EIRd, 
                       file_imm_start=file_imm_start)
  
  # Run simulation of main population
  raw_data <- rcpp_mainpop(params,inputs,trial_params)
  
  # process raw output data
  {
  n_run_names=paste("n_run",c(1:n_runs),sep="")
  n_pt_names=paste("n_pt",c(1:n_pts),sep="")
  na_names=paste("na",c(1:na),sep="")
  n_cat_names=paste("n_cat",c(1:n_cats),sep="")
  n_days_names=paste("n_days",c(1:n_days),sep="")
  n_mv_names=paste("n_mv",c(1:n_mv_values),sep="")
  n_int_names=paste("n_int",c(1:n_int_values),sep="")
  slide_prev_names=paste("slide_prev",c(1:na),sep="")
  pcr_prev_names=paste("pcr_prev",c(1:na),sep="")
  clin_inc_names=paste("clin_inc",c(1:na),sep="")
  dimnames_list=list(na_names,n_pt_names,n_run_names)
  EIR_benchmarks = array(data=raw_data$EIR_benchmarks,dim=c(n_pts,n_runs),dimnames=list(n_pt_names,n_run_names))
  slide_prev_benchmarks = array(data=raw_data$slide_prev_benchmarks,dim=c(na,n_pts,n_runs),dimnames=dimnames_list)
  pcr_prev_benchmarks = array(data=raw_data$pcr_prev_benchmarks,dim=c(na,n_pts,n_runs),dimnames=dimnames_list)
  clin_inc_benchmarks = array(data=raw_data$clin_inc_benchmarks,dim=c(na,n_pts,n_runs),dimnames=dimnames_list)
  EIR_daily_data = array(data=raw_data$EIR_daily_data,dim=c(n_days,n_runs),dimnames=list(n_days_names,n_run_names))
  IB_start_data = array(data=raw_data$IB_start_data,dim=c(n_cats,n_mv_values),dimnames=list(n_cat_names,n_mv_names))
  IC_start_data = array(data=raw_data$IC_start_data,dim=c(n_cats,n_mv_values),dimnames=list(n_cat_names,n_mv_names))
  ID_start_data = array(data=raw_data$ID_start_data,dim=c(n_cats,n_mv_values),dimnames=list(n_cat_names,n_mv_names))
  }
  
  output_data <- list(params=params,n_mv_values=n_mv_values,n_int_values=n_int_values,n_pts=n_pts,
                      time_values=time_values,int_values=int_values,
                      EIR_benchmarks=EIR_benchmarks,slide_prev_benchmarks=slide_prev_benchmarks,
                      pcr_prev_benchmarks=pcr_prev_benchmarks,clin_inc_benchmarks=clin_inc_benchmarks,
                      EIR_daily_data=EIR_daily_data,
                      IB_start_data=IB_start_data,IC_start_data=IC_start_data,ID_start_data=ID_start_data)
  
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
#' @param benchmark           Benchmark type to use in choosing clusters ("EIR", "slide_prev", "pcr_prev", or "clin_inc")
#' @param set_n_pt            Data point to use (1-max)
#' @param set_n_int           Intervention number to use (1-max)
#' @param age_start           Starting age to use when calculating prevalence or incidence over age range (not used with EIR)
#' @param age_end             End age to use when calculating prevalence or incidence over age range (not used with EIR)
#'
#' @export

cluster_input_setup <- function(input_list=list(), benchmark = "EIR",set_n_pt = 1,set_n_int=1,age_start = 0,age_end = 65.0){
  
  # Input error checking (TODO - finish)
  assert_in(set_n_pt,1:input_list$n_pts)
  assert_in(benchmark,c("EIR","slide_prev","pcr_prev","clin_inc"))
  assert_list(input_list)
  assert_int(set_n_pt)
  assert_int(set_n_int)
  assert_bounded(age_start,0.0,65.0)
  assert_bounded(age_end,age_start,65.0)
  assert_in(set_n_pt,c(1:input_list$n_pts))
  assert_in(set_n_int,c(1:input_list$n_int_values))
  
  n_age_start = findInterval(age_start,input_list$params$age_years)
  n_age_end = findInterval(age_end,input_list$params$age_years)
  
  density_sum = 0
  benchmark_values = 0
  if(benchmark == "EIR"){
    benchmark_values = input_list$EIR_benchmarks[set_n_pt,(input_list$n_int_values*c(0:(input_list$n_mv_values-1)))+set_n_int]
  }else{
    j=input_list$n_int_values*c(1:input_list$n_mv_values)
    if(benchmark == "slide_prev"){ benchmark_data = input_list$slide_prev_benchmarks[,set_n_pt,j]}
    if(benchmark == "pcr_prev"){ benchmark_data = input_list$pcr_prev_benchmarks[,set_n_pt,j]}
    if(benchmark == "clin_inc"){ benchmark_data = input_list$clin_inc_benchmarks[,set_n_pt,j] }
    
    for(i in n_age_start:n_age_end){
      density_sum = density_sum + input_list$params$den_norm[i]
      benchmark_values = benchmark_values + benchmark_data[i,]
    }
    benchmark_values = benchmark_values/density_sum
  }
  int_values=input_list$int_values
  
  matplot(c(1:input_list$n_mv_values),benchmark_values,type="p",pch=2,col=2,xlab="N_M",ylab=benchmark)
  
  ret <- list(benchmark=benchmark,set_n_pt = set_n_pt,set_n_int=set_n_int,age_start = age_start,age_end = age_end,
              benchmark_values=benchmark_values,int_values=int_values)
  
  return(ret)
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
  assert_int(n_clusters)
  assert_numeric(benchmark_mean)
  assert_numeric(benchmark_stdev)
  assert_numeric(int_mean)
  assert_numeric(int_stdev)
  
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
    j=findInterval(v_bm1[i],input_list$benchmark_values)+1
    if(j>nv_B){j=nv_B}
    mvn_index[i]=j
    v_bm2[i]=input_list$benchmark_values[j]
    j=findInterval(v_i1[i],input_list$int_values)+1
    if(j>nv_I){j=nv_I}
    int_index[i]=j
    v_i2[i]=input_list$int_values[j]
  }
  
  clusters <- data.frame(CP_B=stats::runif(n_clusters,0,1),n_B=rep(0,n_clusters),B=rep(0,n_clusters),
                         CP_I=stats::runif(n_clusters,0,1),n_I=rep(0,n_clusters),I=rep(0,n_clusters),n_run=rep(0,n_clusters))
  for(i in 1:n_clusters){
    j=findInterval(clusters$CP_B[i],cprob)+1
    clusters$B[i]=v_bm2[j]
    clusters$n_B[i]=mvn_index[j]-1
    j=findInterval(clusters$CP_I[i],cprob)+1
    clusters$I[i]=v_i2[j]
    clusters$n_I[i]=int_index[j]-1
  }
  clusters$n_run=((clusters$n_B)*nv_I)+clusters$n_I
  
  par(mfrow=c(1,2))
  matplot(cprob,v_bm1,type="l",col=1,xlab="Cumulative probability",ylab=input_list$benchmark)
  matplot(cprob,v_bm2,type="l",col=2,add=TRUE)
  matplot(clusters$CP_B,clusters$B,type="p",pch=1,col=3,add=TRUE)
  matplot(cprob,v_i1,type="l",col=1,xlab="Cumulative probability",ylab="Intervention parameter")
  matplot(cprob,v_i2,type="l",col=2,add=TRUE)
  matplot(clusters$CP_I,clusters$I,type="p",pch=1,col=3,add=TRUE)
  
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
#' @param output_folder   Folder to send results data (set to NA to omit saving to file)
#' @param time_interval   Time interval between data points
#' @param n_pts           Number of data points (including 0)
#' @param prop_T_c        Treatment probability in trial cohort
#' @param age_start       Minimum age of trial cohort patients
#' @param age_end         Maximum age of trial cohort patients
#'
#' @export

cohort <- function(mainpop_data = list(), cluster_data=list(),n_patients = 100,output_folder = NA,
                   time_interval = 7.0,n_pts = 13,prop_T_c = 0.9,age_start = 0.5,age_end = 10.0){
  
  # Input error checking (TODO - finish)
  assert_list(mainpop_data)
  assert_list(cluster_data)
  assert_int(n_patients)
  assert_numeric(time_interval)
  assert_int(n_pts)
  assert_bounded(prop_T_c,0.0,1.0)
  assert_bounded(age_start,0.0,65.0)
  assert_bounded(age_end,age_start,65.0)
  
  n_clusters=length(cluster_data$n_B)
  if(is.na(output_folder)==FALSE){
    file_summary = paste(output_folder,"summary.txt",sep="")
    file_frequency = paste(output_folder,"frequency.txt",sep="")
    flag_file=1
  } else {
    file_summary = NA
    file_frequency = NA
    flag_file=0
  }

  trial_params <- list(file_summary = file_summary,file_frequency = file_frequency,n_patients = n_patients,n_clusters=n_clusters,
                       time_interval = time_interval,n_pts = n_pts,prop_T_c = prop_T_c,age_start = age_start,age_end = age_end,
                       flag_file = flag_file, EIR_daily_data = as.vector(mainpop_data$EIR_daily_data),
                       IB_start_data = as.vector(mainpop_data$IB_start_data),
                       IC_start_data = as.vector(mainpop_data$IC_start_data),
                       ID_start_data = as.vector(mainpop_data$ID_start_data))

  raw_data <- rcpp_cohort(mainpop_data$params,trial_params,cluster_data)
  results_data <- data.frame(raw_data)
  
  cluster_names=paste("cluster",c(1:n_clusters),sep="")
  patient_names=paste("patient",c(1:n_patients),sep="")
  n_pt_names=paste("n_pt",c(1:n_pts),sep="")
  patients_status_outputs = array(data=raw_data$patients_status_outputs,dim=c(n_pts,n_patients,n_clusters),
                                  dimnames=list(n_pt_names,patient_names,cluster_names))
  

  return(patients_status_outputs)
}
