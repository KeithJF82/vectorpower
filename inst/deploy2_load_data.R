library(vectorpower)
data_location=system.file("extdata/",package="vectorpower")

# Set inputs
{
  data_file=paste(data_location,"/mainpop_data.Rds",sep="")
  cluster_file1=paste(data_location,"/cluster_list_con.Rds",sep="")
  cluster_file2=paste(data_location,"/cluster_list_int.Rds",sep="")
}

# Load previously generated data from .rds file
mainpop_data <- readRDS(data_file)
cluster_list_con <- readRDS(cluster_file1)
cluster_list_int <- readRDS(cluster_file2)

# Simulate control clusters
cohort_data_con <- cohort(mainpop_data=mainpop_data,cluster_data=cluster_list_con,
                          flag_output=0,prop_T_c = 0.9,n_patients = 100,
                          age_start = 0.5,age_end = 10.0)
plot_cohort_data(cohort_data_con,"slide_prev")

# Simulate intervention clusters
cohort_data_int <- cohort(mainpop_data=mainpop_data,cluster_data=cluster_list_int,
                          flag_output=0,prop_T_c = 0.9,n_patients = 100,
                          age_start = 0.5,age_end = 10.0)
