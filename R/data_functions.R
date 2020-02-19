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

load_dataset <- function (dataset_folder="")
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
  # Error checking
  assert_string(dataset_folder)
  assert_numeric(EIR_values)
  assert_single_int(nyears)
  assert_single_bounded(nyears,1,20)
  # TODO - Find way to check input data exists when files are identified using URLs
  # assert_file_exists(age_file)
  # assert_file_exists(het_file)
  # assert_file_exists(param_file)
  
  if(dir.exists(dataset_folder)==FALSE){dir.create(dataset_folder)}
  
  # Copy parameter data to new files
  cat("\nCreating new model_parameters.txt, het_data.txt and age_data.txt files.\n")
  param_data=read.table(param_file,header=TRUE)
  age_data=read.table(age_file,header=TRUE)
  het_data=read.table(het_file,header=TRUE)
  age_file_new=paste(dataset_folder,"age_data.txt",sep="/")
  het_file_new=paste(dataset_folder,"het_data.txt",sep="/")
  param_file_new=paste(dataset_folder,"model_parameters.txt",sep="/")
  start_file_new=paste(dataset_folder,"start_data.txt",sep="/")
  write.table(param_data,file=param_file_new,col.names=TRUE,sep="\t",quote=FALSE)
  write.table(age_data,file=age_file_new,col.names=TRUE,sep="\t",quote=FALSE)
  write.table(het_data,file=het_file_new,col.names=TRUE,sep="\t",quote=FALSE)
  
  # Create new steady-state constant-rainfall data for desired EIR [TODO: or mv0] values using parameter [TODO: and age] files
  n_mv_values=length(EIR_values)
  age_data=age_data_setup(read.table(age_file_new,header=TRUE,sep="\t")[[1]])
  het_data = read.table(het_file_new,header=TRUE,sep="\t")
  params <- read.table(param_file_new, header=TRUE) 
  na=length(age_data$age_width)
  num_het=length(het_data$het_x)
  params=c(na=na,num_het=num_het,params,age_data,het_data)
  
  cat("\nGenerating preliminary data (steady state at constant rainfall)")
  trial_params <- list(EIRy=EIR_values,file_endpoints=start_file_new)
  raw_data <- rcpp_mainpop_ss(params,trial_params)
  mv_values <- raw_data$mv_values
  
  # Run model without interventions from steady-state data for nyears years in one go to achieve equilibrium
  cat("\nBeginning main population calculations to equilibrium.")
  n_mv_set=c(1:n_mv_values)
  output_folder=paste(dataset_folder,"Temp",sep="/")
  if(dir.exists(output_folder)==FALSE){dir.create(output_folder)}
  input_files=list(age_file=age_file_new,het_file=het_file_new,param_file=param_file_new,
                   start_file=start_file_new,annual_file=NA)
  eq_data <- mainpop(input_files = input_files, output_folder = output_folder,n_mv_set = n_mv_set, int_v_varied = 0, 
                     int_values=c(0.0), start_interval = 0.0, time_values=365*c(0:nyears) )
  plot_mainpop_data(input_list=eq_data,set_n_int=1,benchmark = "EIR",age_start = 0,age_end = 65.0)
  
  # Transfer endpoints as new start data
  cat("\n\nCreating new starting_data.txt file from results.")
  file.copy(from=paste(output_folder,"endpoints.txt",sep="/"),to=start_file_new,
            overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)
  
  # Run for 1 year with daily data points to produce year-round data
  cat("\nRunning main population model for 1 year to establish year-round values")
  cat("\n(annual EIR and incidence, annual average prevalences).")
  input_files=list(age_file=age_file_new,het_file=het_file_new,param_file=param_file_new,
                   start_file=start_file_new,annual_file=NA)
  annual_data <- mainpop(input_files = input_files,n_mv_set = n_mv_set, int_v_varied = 0, int_values=c(0.0),
                         start_interval = 0.0, time_values=1.0*c(0:364) )
  plot_mainpop_data(input_list=annual_data,set_n_int=1,benchmark = "EIR")
  
  # Output annual data for reference
  cat("\n\nOutputting annual data.\n")
  annual_EIR_file=paste(dataset_folder,"annual_data.txt",sep="/")
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
#'                            Represents annual total in case of EIR, year-round average for others
#' @param set_n_int           Intervention number to use (1-max)
#' @param age_start           Starting age to use when calculating prevalence or incidence over age range (not used with EIR)
#' @param age_end             End age to use when calculating prevalence or incidence over age range (not used with EIR)
#' 
#' @export

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
    benchmark_values = input_list$EIR_benchmarks[,set_n_int,]
  }else{
    if(benchmark == "slide_prev"){ benchmark_data = input_list$slide_prev_benchmarks}
    if(benchmark == "pcr_prev"){ benchmark_data = input_list$pcr_prev_benchmarks}
    if(benchmark == "clin_inc"){ benchmark_data = input_list$clin_inc_benchmarks }
    
    for(i in n_age_start:n_age_end){
      density_sum = density_sum + input_list$params$den_norm[i]
      benchmark_values = benchmark_values + benchmark_data[i,,set_n_int,]
    }
    benchmark_values = benchmark_values/density_sum
  }
  
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
#' @title Display results of cohort simulations
#'
#' @description Function which takes the output of the cohort() function and displays selected results
#'
#' @details Takes in a data frame of the form output by the cohort() function, calculates selected 
#' benchmark values and plots them on a graph
#'
#' @param cohort_data List of form output by cohort() function
#' @param benchmark Benchmark type to output ("slide_prev", "pcr_prev", or "clin_inc")
#' @param flag_output Integer indicating how to present data (1: all clusters, 
#'                            2: average over all clusters)
#'
#' @export

plot_cohort_data <- function(cohort_data = list(), benchmark="slide_prev",flag_output=1){

  assert_list(cohort_data)
  assert_in(benchmark,c("slide_prev","pcr_prev","clin_inc"))
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
  
  # Plot graph
  if(flag_output==1){
    if(n_clusters>1){
      matplot(cohort_data$time_values,benchmark_values[,1],type="p",pch=2,col=2,xlab="time (days)",ylab=benchmark,
              ylim=c(0,max(benchmark_values)))
      for(i in 2:n_clusters){
        matplot(cohort_data$time_values,benchmark_values[,i],type="p",pch=2,col=1+i, xaxt="n",xlab="",ylab="",add=TRUE)
      }
    }else{
      matplot(cohort_data$time_values,benchmark_values,type="p",pch=2,col=2,xlab="time (days)",ylab=benchmark,
              ylim=c(0,max(benchmark_values)))
    }
    legend("bottomleft", inset=0.01, legend=c(1:n_clusters), lwd=1.0,col=1+c(1:n_clusters),
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
    matplot(cohort_data$time_values,benchmark_values_mean,type="p",pch=2,col=2,xlab="time (days)",ylab=benchmark,
            ylim=c(0,max(benchmark_values_mean)))
    output=benchmark_values_mean
  }
  
  return(output)
}

#------------------------------------------------
#' @title Tabulate results of cohort simulations
#'
#' @description Function which takes the output of the cohort() function and tabulates results
#'
#' @details Takes the output of the cohort() function and tabulates results by cluster, patient and time point for 
#'          processing using statistical models. Each patient across all clusters will have a unique ID number
#'          (patient_ID = (patient's number in their cluster) + (patient's cluster number - 1)*(number of patients per cluster))
#'
#' @param cohort_data List of form output by cohort() function
#' @param file        Location/name of file to send results
#'
#' @export

tabulate_cohort_data <- function(cohort_data = list(), file="dummy.txt"){
  assert_list(cohort_data)
  
  # Set up data frame
  npts_time=length(cohort_data$time_values)
  n_patients=cohort_data$n_patients
  n_clusters=cohort_data$n_clusters
  n_lines=n_clusters*npts_time*n_patients
  cluster_output=data.frame(n_cluster=rep(0,n_lines),n_patient=rep(0,n_lines),patient_ID=rep(0,n_lines),time=rep(0,n_lines),
                            age=rep(0,n_lines),het_category=rep(0,n_lines),S=rep(0,n_lines),T=rep(0,n_lines),D=rep(0,n_lines),
                            A=rep(0,n_lines),U=rep(0,n_lines),P=rep(0,n_lines),p_det=rep(0,n_lines))
  
  # Transfer cohort results data to data frame
  line=0
  for(n_c in 1:n_clusters){
    for(i in 1:n_patients){
      patient_ID=((n_c-1)*n_patients)+i
      for(j in 1:npts_time){
        line=line+1
        cluster_output$patient_ID[line]=patient_ID
        cluster_output$time[line]=cohort_data$time_values[j]
        state=cohort_data$patients_status_outputs[j,i,n_c]
        if(state==0){cluster_output$S[line]=1}
        else{
          if(state==1){cluster_output$T[line]=1}
          else{
            if(state==2){cluster_output$D[line]=1}
            else{
              if(state==3){
                cluster_output$A[line]=1
                cluster_output$p_det[line]=cohort_data$p_det_outputs[j,i,n_c]
              }
              else{
                if(state==4){cluster_output$U[line]=1}
                else{cluster_output$P[line]=1}
              }
            }
          }
        }
        cluster_output$n_cluster[line]=n_c
        cluster_output$n_patient[line]=i
        cluster_output$age[line]=cohort_data$patients_age_outputs[i,n_c]
        cluster_output$het_category[line]=cohort_data$patients_het_outputs[i,n_c]
      }
    }
  }
  
  # Write data frame to file
  write.table(cluster_output,row.names=FALSE,file=file)
  
  return(NULL)
}
