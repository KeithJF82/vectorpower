library(vectorpower)
data_location=system.file("extdata/",package="vectorpower")

# Input location of dataset folder to use and view annual EIR data
dataset_folder = paste(data_location,"/DemoFolder1/",sep="")
basedata = plot_folder_data(input_folder=dataset_folder,xvalues="N_M",yvalues="EIR")

# Set inputs
{
  start_interval = 0.0  # Length of time to run main population model from starting date
                        # before intervention begins
                        # Default start date for Mali parameter data is day 180 - Jul 1st
  time_values = c(0.0,31.0,62.0,92.0,123.0,153.0) # Time points at which benchmark data output
  n_mv_set = c(1:5) # List of sets of values to load from starting data
  int_values = 0.1*c(0:1) # List of intervention parameter values
  control_folder=paste(dataset_folder,"Control/",sep="")   # Folder (inside dataset folder) to send control cluster output
  int_folder=paste(dataset_folder,"Intervention/",sep="") # Folder to send intervention cluster output
  benchmark1="EIR" # Benchmark to display from main population data after simulations complete
  benchmark2="EIR_annual" # Malaria benchmark to use to select data for clusters
  age_start_coh=0.5 # Minimum age of cohort patients
  age_end_coh=10.0  # Maximum age of cohort patients
  n_clusters=10    # Number of clusters to generate
  n_patients=100    # Number of patients per cluster
  benchmark_mean=0.36  # Mean value of benchmark2 to use to create distribution
  benchmark_stdev=0.085  # Standard deviation of benchmark2 to use to create distribution
  data_file=paste(data_location,"mainpop_data.Rds",sep="")
  cluster_file1=paste(data_location,"cluster_list_con.Rds",sep="")
  cluster_file2=paste(data_location,"cluster_list_int.Rds",sep="")
}

# Simulate main population
mainpop_data <- mainpop(input_files = load_dataset(dataset_folder), n_mv_set = n_mv_set, int_v_varied = 1, int_values=int_values,
                         start_interval = start_interval, time_values=time_values)

# Save main population data to .rds file
saveRDS(mainpop_data,file=data_file)

# Load previously generated data from .rds file
mainpop_data <- readRDS(data_file)

# Plot control results - EIR over time for different baseline malaria levels
plot_con <- plot_mainpop_data(input_list=mainpop_data,set_n_int=1,benchmark=benchmark1)

# Plot intervention results - EIR over time for different baseline malaria levels
plot_int <- plot_mainpop_data(input_list=mainpop_data,set_n_int=2,benchmark=benchmark1)

# Create list of benchmark and intervention parameter values from which to draw clusters
setup_list <- cluster_input_setup(input_list=mainpop_data, set_n_pt=1, set_n_int=1, benchmark = benchmark2)

# ---------------------------------------------------------------------------------------------------------------------------

# Create list of control clusters (note setting of intervention parameter mean 
# and standard deviation to 0)
cluster_list_con <- clusters_create(input_list=setup_list,n_clusters=n_clusters,
                                    benchmark_mean=benchmark_mean, 
                                    benchmark_stdev=benchmark_stdev, 
                                    int_mean=0.0, int_stdev=0.0)

# Create list of intervention clusters
cluster_list_int <- clusters_create(input_list=setup_list,n_clusters=n_clusters,
                                    benchmark_mean=benchmark_mean, 
                                    benchmark_stdev=benchmark_stdev, 
                                    int_mean=0.1, int_stdev=0.0)

# Save cluster sets to files
saveRDS(cluster_list_con,file=cluster_file1) 
saveRDS(cluster_list_int,file=cluster_file2)

# Load previously generated cluster data from .rds files
cluster_list_con <- readRDS(cluster_file1)
cluster_list_int <- readRDS(cluster_file2)

# ---------------------------------------------------------------------------------------------------------------------------

# Simulate control clusters
cohort_data_con <- cohort(mainpop_data=mainpop_data,cluster_data=cluster_list_con,prop_T_c = 0.9,n_patients = n_patients,
                          age_start = age_start_coh,age_end = age_end_coh,flag_output=0)

# Simulate intervention clusters
cohort_data_int <- cohort(mainpop_data=mainpop_data,cluster_data=cluster_list_int,prop_T_c = 0.9,n_patients = n_patients,
                          age_start = age_start_coh,age_end = age_end_coh,flag_output=0)
