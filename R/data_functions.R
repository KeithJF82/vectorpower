#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.
#' @useDynLib vectorpower
#' @import graphics
#' @import utils
NULL
#------------------------------------------------
# TODO - Create function for setting up or adjusting parameter and age data values in R?
#------------------------------------------------
#' @title Load data set from folder
#'
#' @description Function which takes the name of a dataset folder and sets the locations for the files 
#'
#' @details     Takes the location of a dataset folder, sets the locations for the files containing parameter values, 
#'              age and heterogeneity data, starting data and (if applicable) annual data and outputs them as a list for
#'              input to mainpop(). The filenames can also be set manually if they are not in a single dataset folder,
#'              e.g. if they are to be loaded from online
#'                            
#'
#' @param dataset_folder    Dataset folder location
#'
#' @export

load_dataset_folder <- function (dataset_folder="")
{
  
  assert_string(dataset_folder)
  
  file_list=list(age_file=NA,het_file=NA,param_file=NA,start_file=NA,annual_file=NA)
  if(dir.exists(dataset_folder)==1){
    file_list$age_file=paste(dataset_folder,"age_data.txt",sep="/")
    file_list$het_file=paste(dataset_folder,"het_data.txt",sep="/")
    file_list$param_file=paste(dataset_folder,"model_parameters.txt",sep="/")
    file_list$start_file=paste(dataset_folder,"start_data.txt",sep="/")
    file_list$annual_file=paste(dataset_folder,"annual_data.txt",sep="/")
  }
  else{cat("\nFolder does not exist.\n")}
  
  assert_file_exists(file_list$age_file)
  assert_file_exists(file_list$het_file)
  assert_file_exists(file_list$param_file)
  assert_file_exists(file_list$start_file)
  
  return(file_list)
}
#------------------------------------------------
#' @title Load input data from list of files
#'
#' @description Function which takes a list of files and creates list of input data
#'
#' @details     Takes a list of files (as file locations or URLS, created by load_dataset_folder or set manually)
#'              and creates list of input data, including starting model data and parameter values                            
#'
#' @param input_files     List of files containing parameter values, age and heterogeneity data, starting data, 
#'                        and (if applicable) annual data
#' @param n_mv_set        Vector of mosquito density number values to use (must be increasing order)
#'
#' @export

load_inputs <- function (input_files="",n_mv_set=c())
{
  assert_list(input_files)
  assert_numeric(n_mv_set)
  # TODO - Change these to something that works with URLs
  # assert_file_exists(input_files$age_file)
  # assert_file_exists(input_files$het_file)
  # assert_file_exists(input_files$param_file)
  # assert_file_exists(input_files$start_file)
  # assert_file_exists(input_files$annual_file)
  
  n_mv_values=length(n_mv_set)
  n_mv_end=n_mv_set[n_mv_values]
  
  # Read in parameter data from files
  age_data=age_data_setup(read.table(input_files$age_file,header=TRUE,sep="\t")[[1]]) # Read in age data
  het_data = as.list(read.table(input_files$het_file,header=TRUE,sep="\t"))   # Read in biting heterogeneity data
  params <- as.list(read.table(input_files$param_file, header=TRUE))          # Read in model parameters
  input_data <- read.table(input_files$start_file,header=TRUE,nrows=n_mv_end) # Read in starting data
  if(is.na(input_files$annual_file)==TRUE){ annual_data <- list()
  } else { annual_data <- as.list(read.table(input_files$annual_file, header=TRUE)) }
  
  na=length(age_data$age_width)
  num_het=length(het_data$het_x)
  n_cats=na*num_het
  params=c(na=na,num_het=num_het,params,age_data,het_data)
  
  start_data <- list(mv_input=input_data$mv0[n_mv_set], 
                 EL_input=input_data$EL[n_mv_set], LL_input=input_data$LL[n_mv_set], PL_input=input_data$PL[n_mv_set],
                 Sv_input=input_data$Sv1[n_mv_set], Ev_input=input_data$Ev1[n_mv_set], Iv_input=input_data$Iv1[n_mv_set],
                 S_input=c(), T_input=c(), D_input=c(), A_input=c(), U_input=c(), P_input=c(),
                 ICA_input=c(), ICM_input=c(), IB_input=c(), ID_input=c())
  for(i in 1:n_cats){
    start_data$S_input <- append(start_data$S_input,input_data[[8+i]][n_mv_set])
    start_data$T_input <- append(start_data$T_input,input_data[[8+i+n_cats]][n_mv_set])
    start_data$D_input <- append(start_data$D_input,input_data[[8+i+(2*n_cats)]][n_mv_set])
    start_data$A_input <- append(start_data$A_input,input_data[[8+i+(3*n_cats)]][n_mv_set])
    start_data$U_input <- append(start_data$U_input,input_data[[8+i+(4*n_cats)]][n_mv_set])
    start_data$P_input <- append(start_data$P_input,input_data[[8+i+(5*n_cats)]][n_mv_set])
    start_data$ICA_input <- append(start_data$ICA_input,input_data[[8+i+(6*n_cats)]][n_mv_set])
    start_data$ICM_input <- append(start_data$ICM_input,input_data[[8+i+(7*n_cats)]][n_mv_set])
    start_data$IB_input <- append(start_data$IB_input,input_data[[8+i+(8*n_cats)]][n_mv_set])
    start_data$ID_input <- append(start_data$ID_input,input_data[[8+i+(9*n_cats)]][n_mv_set])
  }
  
  inputs <- list(n_mv_set = n_mv_set, start_data = start_data,params = params, annual_data = annual_data)
  
  return(inputs)
  
}

#------------------------------------------------
#' @title Set up age data from age distribution data
#'
#' @description Function which uses the age distribution data (width of each age category in years) from an age_data.txt file
#'              to set up all static age-related parameter values
#'
#' @details     Takes in a vector of age width data and outputs a list of vectors of age parameter values
#'
#' @param age_width_years    Vector of age width data in years
#'
#' @export

age_data_setup <- function(age_width_years=c()){
  
  assert_numeric(age_width_years)
  assert_vector(age_width_years)
  
  dy = 365.0
  eta = 1.0/(21.0*dy)
  age_width=age_width_years*dy
  
  na=length(age_width)
  age=rep(0,na)
  age_years=age
  age[1] = 0.0
  age_years[1]=0.0
  for(i in 2:na){
    age[i] = age[i-1]+age_width[i-1]
    age_years[i] = age_years[i-1]+age_width_years[i-1]
  }
  
  age_rate = 1.0/age_width
  rem_rate = age_rate+eta
  den=rep(0,na)
  den[1] = 1.0 / (1.0 + (age_rate[1] / eta))
  for(i in 2:na){
    den[i] = (age_rate[i - 1]*den[i - 1]) / rem_rate[i]
  }
  den_norm=den/sum(den)
  
  age_data <- list(age_width=age_width,eta=eta,age=age,age_years=age_years,
                   age_rate=age_rate,rem_rate=rem_rate,den=den,den_norm=den_norm)
  
  return(age_data)
}

#------------------------------------------------
#' @title Plot year-round rainfall based on parameter file
#'
#' @description Function for displaying year-round rainfall data based on parameter values in input file
#'
#' @details Reads data from model_parameters.txt in dataset folder, calculates rainfall on each day of the year
#'          and outputs resulting values as a graph
#'
#' @param dataset_folder      Dataset folder
#' 
#' @export

plot_rainfall <- function(dataset_folder=""){
  
  assert_string(dataset_folder)
  
  if(dir.exists(dataset_folder)==FALSE){ cat("\nFolder does not exist.\n")} else {
    params <- read.table(paste(dataset_folder,"model_parameters.txt",sep="/"), header=TRUE) 
    
    Rnorm=params$Rnorm
    rconst=params$rconst
    Re0=params$Re0
    Re1=params$Re1
    Re2=params$Re2
    Re3=params$Re3
    Re4=params$Re4
    Im0=params$Im0
    Im1=params$Im1
    Im2=params$Im2
    Im3=params$Im3
    Im4=params$Im4 
    rain_coeff=((Rnorm) / (14.0*8.0*64.0*365.0));
    trig_coeff1 = (2.0*3.1415926536)/365; 
    trig_coeff2 = trig_coeff1*2.0;
    trig_coeff3 = trig_coeff1*3.0;
    trig_coeff4 = trig_coeff1*4.0;
    trig_coeff5 = trig_coeff1*5.0;
    
    days=c(1:365)
    rainfall=rep(0,365)
    
    for(i in days){
      rainfall[i]=rain_coeff*(rconst + (2.0*(
        (Re0*cos(trig_coeff1*i)) + (Im0*sin(trig_coeff1*i)) +
          (Re1*cos(trig_coeff2*i)) + (Im1*sin(trig_coeff2*i)) +
          (Re2*cos(trig_coeff3*i)) + (Im2*sin(trig_coeff3*i)) +
          (Re3*cos(trig_coeff4*i)) + (Im3*sin(trig_coeff4*i)) +
          (Re4*cos(trig_coeff5*i)) + (Im4*sin(trig_coeff5*i))
      )));
    }
    
    matplot(days,rainfall,type="l",col=1,lwd=1.5,xlab="Day",ylab="Rainfall (a.u.)")
    #matplot(params$date_start,rainfall[params$date_start],type="p",pch=2,col=2,add=TRUE)
    #legend("topright", inset=0.01, legend=c("Rainfall","Start date"),lwd=1.0,col=c(1:2),horiz=FALSE,bg='white',cex=1.0)
    return(rainfall)
  }
  
}

#------------------------------------------------
#' @title Plot pre-calculated data from dataset folder
#'
#' @description Function for displaying pre-calculated annual data or mosquito density parameter values from dataset folder
#'
#' @details Reads data from annual_data.txt in dataset folder and outputs chosen benchmark values as a graph if desired; 
#'          returns list of x and y values as a list
#'
#' @param input_folder      Dataset folder
#' @param xvalues           Data to display on x-axis ("N_M" for line number, "M" for mosquito density parameter)
#' @param yvalues           Data to display on y-axis ("M", "EIR", "slide_prev", "pcr_prev", or "clin_inc")
#'                          EIR represents annual total, prevalences and incidence represent year-round averages
#' @param plot_flag         Logical operator indicating whether or not to plot graph of read values
#' 
#' @export

get_folder_data <- function(input_folder="",xvalues="N_M",yvalues = "M", plot_flag = TRUE){
  
  # Input error checking (TODO - finish)
  assert_in(xvalues,c("N_M","M"))
  assert_in(yvalues,c("M","EIR","EIR_ss","slide_prev","pcr_prev","clin_inc"))
  assert_logical(plot_flag)
  
  annual_data <- as.list(read.table(paste(input_folder,"annual_data.txt",sep="/"), header=TRUE))   # Read in annual data
  n_mv_values=length(annual_data$n_mv)
  
  if(xvalues=="N_M"){ xdata = annual_data$n_mv } else { xdata = annual_data$mv0 }
  
  if(yvalues=="M"){ydata = annual_data$mv0}
  if(yvalues=="EIR"){ydata = annual_data$EIRy}
  if(yvalues=="EIR_ss"){ydata = annual_data$EIRy_ss}
  if(yvalues=="slide_prev"){ydata = annual_data$slide_prev_y}
  if(yvalues=="clin_inc"){ydata = annual_data$clin_inc_y}
  if(yvalues=="pcr_prev"){ydata = annual_data$pcr_prev_y}
  
  if(plot_flag==TRUE){
    matplot(xdata,ydata,type="p",pch=2,col=2,xlab=xvalues,ylab=yvalues,ylim=c(0,max(ydata)))
  }
  
  output <- list(xdata,ydata)
  names(output)[[1]]=xvalues
  names(output)[[2]]=yvalues
  
  return(output)
}

#------------------------------------------------
#' @title Get selected main population output data as a function of time
#'
#' @description Takes in detailed benchmark data as a list and  returns selected benchmark values as a list
#'
#' @details Function for taking output of mainpop(), selecting the desired benchmark
#'          (EIR, slide prevalence, PCR prevalence, clinical incidence) data, and outputting
#'          values of the benchmark as a function of time for the desired baseline malaria level(s)
#'          and intervention parameter value(s)
#'
#' @param input_list          List containing mainpop output data
#' @param benchmark           Benchmark type to use in choosing clusters ("EIR", "slide_prev", "pcr_prev", or "clin_inc")
#'                            Represents annual total in case of EIR, year-round average for others
#' @param set_n_mv            Baseline malaria level number(s) to output data for (from numbers in main population data, not
#'                            original data set)
#' @param set_n_int           Intervention number(s) to output data for (from numbers in main population data)
#' @param age_start           Starting age to use when calculating prevalence or incidence over age range (not used with EIR)
#' @param age_end             End age to use when calculating prevalence or incidence over age range (not used with EIR)
#' 
#' @export

get_mainpop_data <- function(input_list=list(),benchmark = "EIR", set_n_mv=1, set_n_int=1, age_start = 0, age_end = 65.0){
  
  # Input error checking
  assert_list(input_list)
  assert_in(benchmark,c("EIR","slide_prev","pcr_prev","clin_inc"))
  assert_int(set_n_mv)
  assert_int(set_n_int)
  assert_bounded(age_start,0.0,65.0)
  assert_bounded(age_end,age_start,65.0)
  assert_in(set_n_mv,c(1:length(input_list$n_mv_set)))
  assert_in(set_n_int,c(1:input_list$n_int_values))
  
  n_age_start = findInterval(age_start,input_list$params$age_years)
  n_age_end = findInterval(age_end,input_list$params$age_years)
  
  density_sum = 0
  benchmark_values = 0
  if(benchmark == "EIR"){
    benchmark_values = input_list$EIR_benchmarks[,set_n_int,set_n_mv]
  }else{
    if(benchmark == "slide_prev"){ benchmark_data = input_list$slide_prev_benchmarks}
    if(benchmark == "pcr_prev"){ benchmark_data = input_list$pcr_prev_benchmarks}
    if(benchmark == "clin_inc"){ benchmark_data = input_list$clin_inc_benchmarks }
    
    for(i in n_age_start:n_age_end){
      density_sum = density_sum + input_list$params$den_norm[i]
      benchmark_values = benchmark_values + benchmark_data[i,,set_n_int,set_n_mv]
    }
    benchmark_values = benchmark_values/density_sum
  }
  
  output <- list(data = list(), set_n_mv = set_n_mv, set_n_int = set_n_int)
  output$data <- list(input_list$time_values,benchmark_values)
  names(output$data)[[1]]="time (days)"
  names(output$data)[[2]]=benchmark
  
  return(output)
}

#------------------------------------------------
#' @title Plot graph of main population data
#'
#' @description Plot graph of main population data selected using get_mainpop_data function
#'
#' @details Plot graph of main population data selected using get_mainpop_data function
#'
#' @param input_list List containing main population benchmark data selected using get_mainpop_data
#' 
#' @export

plot_mainpop_data <- function(input_list=list()){
  
  # Input error checking
  assert_list(input_list)
  assert_list(input_list$data)
  benchmark=names(input_list$data[2])
  assert_in(benchmark,c("EIR","slide_prev","pcr_prev","clin_inc"))
  
  set_n_mv = input_list$set_n_mv
  set_n_int = input_list$set_n_int
  n_mv_values=length(set_n_mv)
  n_int_values=length(set_n_int)
  n_curves=n_mv_values*n_int_values
  time_values=input_list$data[[1]]
  benchmark_values=input_list$data[[2]]
  ylim=c(min(benchmark_values),max(benchmark_values))
  
  if(n_int_values==1){
    xlim=c(min(time_values),max(time_values)*1.25)
    titles=rep(" ",n_curves)
    colours=1+c(1:n_curves)
    for(i in 1:n_curves){
      titles[i]=paste("n_mv=",set_n_mv[i],sep="")
    }
    matplot(time_values,benchmark_values,type="p",pch=1,col=colours,xlab="time (days)",ylab=benchmark,
            xlim=xlim,ylim=ylim)
    legend("bottomright", inset=0.01, legend=titles, pch=1,col=colours,horiz=FALSE,bg='white',cex=1.0)
  } else {
    xlim=c(min(time_values),max(time_values)*1.4)
    titles=rep(" ",n_curves)
    colours=1+c(1:n_mv_values)
    points=rep(0,n_curves)
    for(i in 1:n_int_values){
      for(j in 1:n_mv_values){
        nt=((i-1)*n_mv_values)+j
        titles[nt]=paste("n_int=",set_n_int[i],", n_mv=",set_n_mv[j],sep="")
        points[nt]=i
      }
    }
    if(n_mv_values==1){
      matplot(time_values,benchmark_values[,1],type="p",pch=1,col=colours,xlab="time (days)",ylab=benchmark,xlim=xlim,ylim=ylim)
      for(i in 2:n_int_values){
        matplot(time_values,benchmark_values[,i],type="p",pch=i,col=1+c(1:n_mv_values),xlab="time (days)",ylab=benchmark,add=TRUE)
      }
      
    } else {
      matplot(time_values,benchmark_values[,1,],type="p",pch=1,col=colours,xlab="time (days)",ylab=benchmark,xlim=xlim,ylim=ylim)
      for(i in 2:n_int_values){
        matplot(time_values,benchmark_values[,i,],type="p",pch=i,col=1+c(1:n_mv_values),xlab="time (days)",ylab=benchmark,add=TRUE)
      }
      
    }
    legend("bottomright", inset=0.01, legend=titles, pch=points,col=rep(colours,n_int_values), horiz=FALSE,bg='white',cex=1.0)
  }
  
  return(NULL)
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
#' @param set_n_pt            Data point to use (1-max) if using new data
#' @param set_n_int           Intervention number to use (1-max) if using new data
#' @param age_start           Starting age to use when calculating prevalence or incidence over age range 
#'                            (not used with EIR)
#' @param age_end             End age to use when calculating prevalence or incidence over age range 
#'                            (not used with EIR)
#' @param plot_flag           Logical operator indicating whether or not to plot graph of read values
#'
#' @export

cluster_input_setup <- function(input_list=list(), benchmark = "EIR",set_n_pt = 1,set_n_int=1,age_start = 0,
                                age_end = 65.0,plot_flag=FALSE){
  
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
#' @title Display results of cohort simulations
#'
#' @description Function which takes the output of the cohort() function and displays selected results
#'
#' @details Takes in a data frame of the form output by the cohort() function, calculates selected 
#' benchmark values and plots them on a graph
#'
#' @param cohort_data   List of form output by cohort() function
#' @param benchmark     Benchmark type to output ("slide_prev", "pcr_prev", "clin_inc" or "test_inc")
#' @param flag_output   Integer indicating how to present data (1: all clusters, 2: average over all clusters)
#'
#' @export

plot_cohort_data <- function(cohort_data = list(), benchmark="slide_prev",flag_output=1){

  assert_list(cohort_data)
  assert_in(benchmark,c("slide_prev","pcr_prev","clin_inc","test_inc"))
  assert_single_int(flag_output)
  assert_in(flag_output,c(1,2))
  
  n_pts=cohort_data$n_pts
  n_patients=cohort_data$n_patients
  n_clusters=cohort_data$n_clusters
  nvalues=n_pts*n_clusters
  dv=1.0/n_patients
  
  # Extract/calculate values of chosen benchmark
  benchmark_values=array(data=rep(0.0,nvalues),dim=c(n_pts,n_clusters))
  for(i in 1:n_pts){
    for(j in 1:n_clusters){
      #patient_status_outputs[npt,patient,cluster]
      if(benchmark=="clin_inc"){ benchmark_values[i,j]=cohort_data$clin_inc_outputs[i,j] } else {
        
        if(benchmark=="test_inc"){
          benchmark_values[i,j]=cohort_data$test_incidence_outputs[i,j]
        } else {
          frequencies=rep(0,6)
          for(k in 1:n_patients){
            status=cohort_data$patients_status_outputs[i,k,j]+1
            frequencies[status]=frequencies[status]+1
          }
          
          if(benchmark=="pcr_prev"){benchmark_values[i,j]=sum(frequencies[2:5])*dv}
          
          if(benchmark=="slide_prev"){benchmark_values[i,j]=(sum(frequencies[2:3])+sum(cohort_data$p_det_outputs[i,,j]))*dv
          } 
        }
      }
    }
  }
  
  # Plot graph
  if(benchmark=="test_inc"){time_values=cohort_data$test_time_values} else {time_values=cohort_data$time_values}
  if(flag_output==1){
    if(n_clusters>1){
      matplot(time_values,benchmark_values[,1],type="p",pch=2,col=2,xlab="time (days)",ylab=benchmark,
              ylim=c(0,max(benchmark_values)))
      for(i in 2:n_clusters){
        matplot(time_values,benchmark_values[,i],type="p",pch=1+i,col=1+i, xaxt="n",xlab="",ylab="",add=TRUE)
      }
    }else{
      matplot(time_values,benchmark_values,type="p",pch=2,col=2,xlab="time (days)",ylab=benchmark,
              ylim=c(0,max(benchmark_values)))
    }
    legend("bottomleft", inset=0.01, legend=c(1:n_clusters),pch=1+c(1:n_clusters), col=1+c(1:n_clusters),
           horiz=FALSE,bg='white',cex=1.0)
    output=benchmark_values
  } else {
    if(n_clusters>1){
      benchmark_values_mean=rep(0,n_pts)
      for(i in 1:n_pts){
        benchmark_values_mean[i]=mean(benchmark_values[i,])
      }
    }else{
      benchmark_values_mean=benchmark_values
    }
    
    matplot(time_values,benchmark_values_mean,type="p",pch=2,col=2,xlab="time (days)",ylab=benchmark,
            ylim=c(0,max(benchmark_values_mean)))
    output=benchmark_values_mean
  }
  
  return(output)
}

#------------------------------------------------
#' @title Add entomological start data file to existing data folder
#'
#' @description Function which creates a starting data file for entomological modelling in existing data folder
#'
#' @details     Takes in the name of an existing data folder and creates start_data_ento.txt file from
#'              information in start_data.txt file.
#'
#' @param dataset_folder    Dataset folder (must be an existing dataset folder containing model_parameters.txt
#'                          and start_data.txt files)
#' 
#' @export

create_ento_start_data <- function (dataset_folder=""){
  assert_file_exists(paste(dataset_folder,"model_parameters.txt",sep="/"))
  start_data_file=paste(dataset_folder,"start_data.txt",sep="/")
  assert_file_exists(start_data_file)

  start_data=read.table(start_data_file,header=TRUE)
  start_data_file_new=paste(dataset_folder,"start_data_ento.txt",sep="/")
  
  cat(file=start_data_file_new,"n_run\tmv0\tEL\tLL\tPL")
  na_m=50
  for(i in 1:na_m){cat(file=start_data_file_new,"\tSv",i,sep="",append=TRUE)}
  for(i in 1:na_m){cat(file=start_data_file_new,"\tEv",i,sep="",append=TRUE)}
  for(i in 1:na_m){cat(file=start_data_file_new,"\tIv",i,sep="",append=TRUE)}
  for(j in 1:length(start_data$n_run)){
    cat(file=start_data_file_new,"\n",append=TRUE)
    cat(file=start_data_file_new,j,start_data$mv0[j],start_data$EL[j],start_data$LL[j],start_data$PL[j],sep="\t",append=TRUE) 
    cat(file=start_data_file_new,"\t",start_data$Sv1[j],sep="",append=TRUE)
    for(i in 2:na_m){cat(file=start_data_file_new,"\t0",sep="",append=TRUE)}
    cat(file=start_data_file_new,"\t",start_data$Ev1[j],sep="",append=TRUE)
    for(i in 2:na_m){cat(file=start_data_file_new,"\t0",sep="",append=TRUE)}
    cat(file=start_data_file_new,"\t",start_data$Iv1[j],sep="",append=TRUE)
    for(i in 2:na_m){cat(file=start_data_file_new,"\t0",sep="",append=TRUE)}
  }
  
  
  return(NULL)
}
