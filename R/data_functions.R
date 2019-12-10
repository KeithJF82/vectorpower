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
#' @param age_file          File containing age distribution parameters to use for new dataset
#' @param het_file          File containing biting heterogeneity parameters to use for new dataset
#' @param nyears            Number of years to run model for to achieve equilibrium
#'
#' @export

dataset_create <- function (dataset_folder="",EIR_values=c(1.0),param_file="",age_file="",het_file="",nyears=10)
{
  # Error checking (TODO - Finish)
  assert_string(dataset_folder)
  assert_numeric(EIR_values)
  # assert_string(param_file) # TODO - Change to check that file exists
  # assert_string(age_file)   # TODO - Change to check that file exists
  # assert_string(het_file)   # TODO - Change to check that file exists
  assert_single_int(nyears)
  assert_single_bounded(nyears,1,20)
  
  #
  if(dir.exists(dataset_folder)==FALSE){dir.create(dataset_folder)}
  
  # Copy parameter files from existing folder to new one
  cat("\nCreating new model_parameters.txt and age_data.txt files.")
  param_file_new=paste(dataset_folder,"model_parameters.txt",sep="")
  age_file_new=paste(dataset_folder,"age_data.txt",sep="")
  het_file_new=paste(dataset_folder,"het_data.txt",sep="")
  file.copy(from=param_file,to=param_file_new, overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)
  file.copy(from=age_file,to=age_file_new, overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)
  file.copy(from=het_file,to=het_file_new, overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)
  
  # Create new steady-state constant-rainfall data for desired EIR [TODO: or mv0] values using parameter [TODO: and age] files
  n_mv_values=length(EIR_values)
  file_endpoints = paste(dataset_folder,"start_data.txt",sep="")
  age_data=age_data_setup(read.table(age_file_new,header=TRUE,sep="\t")[[1]])
  het_data = as.list(read.table(paste(dataset_folder,"het_data.txt",sep=""),header=TRUE,sep="\t"))
  params <- read.table(paste(dataset_folder,"model_parameters.txt",sep=""), header=TRUE) 
  na=length(age_data$age_width)
  num_het=length(het_data$het_x)
  params=c(na=na,num_het=num_het,params,age_data,het_data)
  
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
    EIRy_values[i]=sum(annual_data$EIR_benchmarks[,1,i])
    slide_prevy_values[i]=sum(annual_data$slide_prev_benchmarks[,,1,i])/365
    clin_incy_values[i]=sum(annual_data$clin_inc_benchmarks[,,1,i])
    pcr_prevy_values[i]=sum(annual_data$pcr_prev_benchmarks[,,1,i])/365
    cat(i,mv_values[i],EIR_values[i],EIRy_values[i],slide_prevy_values[i],clin_incy_values[i],pcr_prevy_values[i],"\n",sep="\t")
    cat(i,mv_values[i],EIR_values[i],EIRy_values[i],slide_prevy_values[i],clin_incy_values[i],pcr_prevy_values[i],"\n",
        file=annual_EIR_file,append=TRUE,sep="\t")
  }
  
  # Delete Temp folder
  setwd(dataset_folder)
  unlink("Temp", recursive=TRUE,force=TRUE)
  
  return(annual_data)
}
#------------------------------------------------
#' @title Plot selected main population output data as a function of time
#'
#' @description Function for taking output of mainpop(), selecting the desired benchmark
#'              (EIR, slide prevalence, PCR prevalence, clinical incidence) data, and 
#'              plotting it on a graph as an aid to setting up clusters
#'
#' @details Takes in detailed benchmark data as a list and plots a graph based on the
#'          input parameters ; returns graph x and y values as a list
#'
#' @param input_list          List containing mainpop output data
#' @param benchmark           Benchmark type to use in choosing clusters ("EIR", "slide_prev", "pcr_prev", or "clin_inc")
#' @param set_n_int           Intervention number to use (1-max)
#' @param age_start           Starting age to use when calculating prevalence or incidence over age range (not used with EIR)
#' @param age_end             End age to use when calculating prevalence or incidence over age range (not used with EIR)

plot_mainpop_data <- function(input_list=list(),benchmark = "EIR", set_n_int=1, age_start = 0, age_end = 65.0){
  
  # Input error checking (TODO - finish)
  assert_list(input_list)
  assert_in(benchmark,c("EIR","slide_prev","pcr_prev","clin_inc"))
  assert_int(set_n_int)
  assert_bounded(age_start,0.0,65.0)
  assert_bounded(age_end,age_start,65.0)
  assert_in(set_n_int,c(1:input_list$n_int_values))
  
  n_age_start = findInterval(age_start,input_list$params$age_years)
  n_age_end = findInterval(age_end,input_list$params$age_years)
  
  density_sum = 0
  benchmark_values = 0
  if(benchmark == "EIR"){
    benchmark_values = input_list$EIR_benchmarks[,set_n_int,input_list$n_mv_set]
  }else{
    if(benchmark == "slide_prev"){ benchmark_data = input_list$slide_prev_benchmarks[,,set_n_int,input_list$n_mv_set]}
    if(benchmark == "pcr_prev"){ benchmark_data = input_list$pcr_prev_benchmarks[,,set_n_int,input_list$n_mv_set]}
    if(benchmark == "clin_inc"){ benchmark_data = input_list$clin_inc_benchmarks[,,set_n_int,input_list$n_mv_set] }
    
    for(i in n_age_start:n_age_end){
      density_sum = density_sum + input_list$params$den_norm[i]
      benchmark_values = benchmark_values + benchmark_data[i,,]
    }
    benchmark_values = benchmark_values/density_sum
  }
  int_values=input_list$int_values
  
  if(input_list$n_mv_values>1){
    matplot(input_list$time_values,benchmark_values[,1],type="p",pch=2,col=2,xlab="time (days)",ylab=benchmark,
            ylim=c(0,max(benchmark_values)))
    for(i in 2:input_list$n_mv_values){
      matplot(input_list$time_values,benchmark_values[,i],type="p",pch=2,col=1+i, xaxt="n",xlab="",ylab="",add=TRUE)
    }
  }else{
    matplot(input_list$time_values,benchmark_values,type="p",pch=2,col=2,xlab="time (days)",ylab=benchmark,
            ylim=c(0,max(benchmark_values)))
  }
  legend("bottomleft", inset=0.01, legend=c(1:input_list$n_mv_values), lwd=1.0,col=1+c(1:input_list$n_mv_values),
         horiz=FALSE,bg='white',cex=1.0)
  
  output <- list(input_list$time_values,benchmark_values)
  names(output)[[1]]="time (days)"
  names(output)[[2]]=benchmark
  
  return(output)
}

#------------------------------------------------
#' @title Plot year-round rainfall based on parameter file
#'
#' @description Function for displaying year-round rainfall data based on parameter values in input file
#'
#' @details Reads data from model_parameters.txt in dataset folder, calculates rainfall on each day of the year
#'          and outputs resulting values as a graph
#'
#' @param input_folder      Dataset folder

plot_rainfall <- function(input_folder=""){
  
  assert_string(input_folder)
  
  if(dir.exists(input_folder)==FALSE){ cat("\nFolder does not exist.\n")} else {
    params <- read.table(paste(input_folder,"model_parameters.txt",sep=""), header=TRUE) 
    
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
    
    days=c(0:364)
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
    matplot(params$date_start,rainfall[params$date_start],type="p",pch=2,col=2,add=TRUE)
    legend("topright", inset=0.01, legend=c("Rainfall","Start date"),lwd=1.0,col=c(1:2),horiz=FALSE,bg='white',cex=1.0)
    
  }
  
}
#------------------------------------------------
#' @title Plot pre-calculated data from dataset folder
#'
#' @description Function for displaying pre-calculated annual data or mosquito density parameter values from dataset folder
#'
#' @details Reads data from annual_data.txt in dataset folder and outputs chosen benchmark values as a graph; returns graph 
#'          x and y values as a list
#'
#' @param input_folder      Dataset folder
#' @param xvalues           Data to display on x-axis ("N_M" for line number, "M" for mosquito density parameter)
#' @param yvalues           Data to display on y-axis ("M", "EIR", "slide_prev", "pcr_prev", or "clin_inc")

plot_folder_data <- function(input_folder="",xvalues="N_M",yvalues = "M"){
  
  # Input error checking (TODO - finish)
  assert_in(xvalues,c("N_M","M"))
  assert_in(yvalues,c("M","EIR","EIR_ss","slide_prev","pcr_prev","clin_inc"))
  
  annual_data <- as.list(read.table(paste(input_folder,"annual_data.txt",sep=""), header=TRUE))   # Read in annual data
  n_mv_values=length(annual_data$n_mv)
  
  if(xvalues=="N_M"){ xdata = annual_data$n_mv } else { xdata = annual_data$mv0 }
  
  if(yvalues=="M"){ydata = annual_data$mv0}
  if(yvalues=="EIR"){ydata = annual_data$EIRy}
  if(yvalues=="EIR_ss"){ydata = annual_data$EIRy_ss}
  if(yvalues=="slide_prev"){ydata = annual_data$slide_prev_y}
  if(yvalues=="clin_inc"){ydata = annual_data$clin_inc_y}
  if(yvalues=="pcr_prev"){ydata = annual_data$pcr_prev_y}
  
  matplot(xdata,ydata,type="p",pch=2,col=2,xlab=xvalues,ylab=yvalues,ylim=c(0,max(ydata)))
  
  output <- list(xdata,ydata)
  names(output)[[1]]=xvalues
  names(output)[[2]]=yvalues
  
  return(output)
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
