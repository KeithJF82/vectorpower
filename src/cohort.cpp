#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

struct patients //Structure containing individual patient data for cohort
{
	int status; //Infection group - 0=S, 1=T, 2=D, 3=A, 4=U, 5=P
	int na, num_het; //Age and heterogeneity groups
	double IB, IC, ID; //Immunities
	int infected; //Flag indicating patient in early stage of infection
	double delay; //Infection delay time
}
patient[1000];

//Global variables-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double dy = 365.0;			//Days in a year
double tinterval1 = 1.0;	//Interval between calculations of current prevalence/incidence values
int na = 145;				//Number of age categories in main population
int num_het = 9;			//Number of heterogeneity categories
double het_x[] = { -4.5127459, -3.2054291, -2.076848, -1.0232557, 0.0, 1.0232557, 2.076848, 3.2054291, 4.5127459 }; // Biting heterogeneity
double het_wt[] = { 0.00002235, 0.00278914, 0.04991641, 0.2440975, 0.40634921, 0.2440975, 0.04991641, 0.00278914, 0.00002235 }; // Fractions in each heterogeneity category

//Additional functions----------------------------------------------------------------------------------------------------------------------------------------------------------------------

int intdiv(double x, double y);
double randgen(int n);	//Generate random number between 0 and 1 to n decimal places
int rcpp_to_int(SEXP x);
double rcpp_to_double(SEXP x);
string rcpp_to_string(SEXP x);

// [[Rcpp::export]]
int rcpp_cohort(List params, List cohort_params)
{ 
	int num = 0; 
	int na_c0, na_c1, na_c2, n_c, mvn, int_num, data_num, i, j, ij, n, nt, div, pos, ntmax, tflag, interval, infected; 
	double v, dt, t, EIRd, EIRt, t_mark1, t_mark2, p_multiplier, prob, cprobmax, p_inf_bite, p_inf_from_bite, p_clin_inf, IC_cur, IB_cur, ID_cur, rate_db_t, rate_dc_t, rate_dd_t;
	FILE *output1;
	FILE *output2;
	FILE *input_EIR; 
	FILE *input_immunity; 
	FILE *input_clusters; 
	srand(time(NULL));
 
	output1 = NULL; 
	output2 = NULL; 
	input_EIR = NULL; 
	input_immunity = NULL; 
	input_clusters = NULL;
 
	//Load input parameter data from R------------------------------------------------------------------------------------------------------------------------------------------
 	
	Rcout << "\nLoading input parameter data from R\n";
	R_FlushConsole();
	string input_filename1 = rcpp_to_string(cohort_params["input_file1"]);
 	string input_filename2 = rcpp_to_string(cohort_params["input_file2"]);
 	string input_filename3 = rcpp_to_string(cohort_params["input_file3"]);
 	string output_filename1 = rcpp_to_string(cohort_params["file_summary"]);	//Individual data (fraction in each category) every tinterval2 days over intervention period for each iteration (single run) or average cohort prevalence/indidence/PCR test positivity rate every tinterval2 days after intervention start (multi-run)
 	string output_filename2 = rcpp_to_string(cohort_params["file_endpoints"]);	//Positive PCR test frequency data for each iteration (single run) or average across all iterations (multi-run)
 	//string output_filename3 = rcpp_to_string(cohort_params["output_file3"]);	//PCR test results (0/1) for each individual every tinterval2 days over intervention period for each iteration (single run)
 	int n_mv_values = rcpp_to_int(cohort_params["n_mv_values"]);				//Number of lines in immunity data input file (
 	int n_EIR_values = rcpp_to_int(cohort_params["n_EIR_values"]);				//Number of lines in EIR data input file (number of unique runs)
 	double tinterval2 = rcpp_to_double(cohort_params["time_interval"]);			//Interval between calculations of regular test probabilities and averaged prevalence/incidence values
 	int n_divs = rcpp_to_int(cohort_params["n_divs"]);							//Number of tinterval2-length divisions to run after intervention starts
	double tmax = (n_divs - 1) * tinterval2;									//Duration of intervention trial 
	int tmax_i = intdiv(tmax, 1.0);
	int tmax_i2 = tmax_i + 1;
 	int n_clusters= rcpp_to_int(cohort_params["n_clusters"]);					//Number of lines in cluster data input file (number of clusters)
	double prop_T_c = rcpp_to_int(cohort_params["prop_T_c"]);					//Proportion of clinical cases successfully treated in cohort
	int n_patients = rcpp_to_int(cohort_params["n_patients"]);					//Number of patients in cohort
	double age_c0 = rcpp_to_int(cohort_params["age_c0"]);						//Minimum age in cohort
	double age_c1 = rcpp_to_int(cohort_params["age_c1"]);						//Maximum age in cohort
	double age_c2 = age_c1 + (tmax / dy);										//Maximum age in cohort taking into account ageing during trial
 
	Rcout << "\nEIR data input file: " << input_filename1 << "\nImmunity data input file: " << input_filename2 << "\nCluster data input file: " << input_filename3;
	Rcout << "\nn_mv_values: " << n_mv_values << "\nn_EIR_values: " << n_EIR_values << "\nn_divs: " << n_divs << "\ntinterval2: " << tinterval2 << "\nn_clusters: " << n_clusters << "\ntmax: " << tmax << "\ttmax_i: " << tmax_i;
	R_FlushConsole();
	 
	//Set up constant parameters------------------------------------------------------------------------------------------------------------------------------------------------
	 
	Rcout <<"\nInitializing constants.\n";
	R_FlushConsole();
	int n_cats = na * num_het; //Total number of age/heterogeneity categories in main population 
	double inv_dy = 1.0 / dy;
	double dv_p1 = 1.0 / n_patients;
	double dv_p2 = dv_p1*(dy / tinterval2);
	double rho = rcpp_to_double(params["rho"]);										// Age-dependent biting parameter
	double a0 = rcpp_to_double(params["a0"]);										// Age-dependent biting parameter
	double sigma2 = rcpp_to_double(params["sigma2"]);								// Variance of log heterogeneity in biting
	double dur_E = rcpp_to_double(params["dur_E"]);									// Latent period
	double dur_T = rcpp_to_double(params["dur_T"]);									// Average time to recover from disease and parasitaemia when treated
	double dur_D = rcpp_to_double(params["dur_D"]);									// Average time to move from clinical disease to asymptomatic when not successfully treated
	double dur_A = rcpp_to_double(params["dur_A"]);									// Average time to move from asymptomatic patent infection to sub-patent infection
	double dur_U = rcpp_to_double(params["dur_U"]);									// Average time to recover from sub-patent infection
	double dur_P = rcpp_to_double(params["dur_P"]);									// Average time to leave protected state
 
	// Immunity parameters - infection
	double bh = rcpp_to_double(params["bh"]);										// Maximum probability due to no immunity
	double bmin = rcpp_to_double(params["bmin"]);									// Maximum relative reduction due to immunity
	double db = rcpp_to_double(params["db"]);										// Inverse of decay rate
	double kb = rcpp_to_double(params["kb"]);										// Shape parameter
	double ub = rcpp_to_double(params["ub"]);										// Duration in which immunity is not boosted
	double inv_IB0 = 1.0 / rcpp_to_double(params["IB0"]);							// Inverse of scale parameter
	double bmin_rev = 1.0 - bmin;
	double IB_boost = 1.0 / (ub + 1.0);

	// Immunity parameters - detection
	double dmin = rcpp_to_double(params["dmin"]);									// Minimum probability due to maximum immunity
	double dd = rcpp_to_double(params["dd"]);										// Inverse of decay rate
	double kd = rcpp_to_double(params["kd"]);										// Shape parameter
	double ud = rcpp_to_double(params["ud"]);										// Duration in which immunity is not boosted
	double ad0 = rcpp_to_double(params["ad0"]);										// Scale parameter relating age to immunity
	double fd0 = rcpp_to_double(params["fd0"]);										// Time scale at which immunity changes with age
	double gammad = rcpp_to_double(params["gammad"]);								// Shape parameter relating age to immunity
	double inv_ID0 = 1.0 / rcpp_to_double(params["ID0"]);							// Scale parameter
	double dmin_rev = 1.0 - dmin;
	double ID_boost = 1.0 / (ud + 1.0);

	// Immunity parameters - clinical disease
	double phi0 = rcpp_to_double(params["phi0"]);									// Maximum probability due to no immunity
	double phi1 = rcpp_to_double(params["phi1"]);									// Maximum relative reduction due to no immunity
	double dc = rcpp_to_double(params["dc"]);										// Inverse of decay rate
	double kc = rcpp_to_double(params["kc"]);										// Shape parameter
	double uc = rcpp_to_double(params["uc"]);										// Duration in which immunity is not boosted
	double inv_IC0 = 1.0 / rcpp_to_double(params["IC0"]);							// Scale parameter
	double phi1_rev = 1.0 - phi1;
	double IC_boost = 1.0 / (uc + 1.0);
	
	//Constant derived values--------------------------------------------------------------------------------------------------------------------------

	Rcout << "\nInitializing constant derived values.\n";
	R_FlushConsole();

	double rT = 1.0 / dur_T;			//Rate of recovery from disease and parasitaemia when treated
	double rD = 1.0 / dur_D;			//Rate of moving from clinical disease to asymptomatic when not successfully treated
	double rA0 = 1.0 / dur_A;			//Rate of moving from asymptomatic patent infection to sub-patent infection
	double rU = 1.0 / dur_U;			//Rate of recovery from sub - patent infection
	double rP = 1.0 / (dur_P - dur_T);	//Rate of leaving protected state
	double rate_dc = 1.0 / dc;
	double rate_db = 1.0 / db;
	double rate_dd = 1.0 / dd;
	double *age_width = (double *)malloc(na * sizeof(double));
	double *age_rate = (double *)malloc(na * sizeof(double));
	double *rem_rate = (double *)malloc(na * sizeof(double));
	double *age = (double *)malloc(na * sizeof(double));
	double *agey0 = (double *)malloc(na * sizeof(double));
	double *agey1 = (double *)malloc(na * sizeof(double));
	double eta = 1.0 / (21.0*dy);		//death rate, for exponential population distribution
	int agefinding0 = 1;				//Flags indicating which cohort age parameter is being sought while running through ages
	int agefinding1 = 0;
	int agefinding2 = 0;

	na_c0 = 0;							//Minimum age category number in main population corresponding to first age category in cohort; determined from age_c0
	na_c1 = na - 1;						//Maximum age category number in main population corresponding to last age category in cohort at start; determined from age_c1
	na_c2 = na - 1;						//Maximum age category number in main population corresponding to last age category in cohort at end; determined from age_c2
	for (i = 0; i < na; i++)
	{
		if (i < 25) { age_width[i] = dy * 0.04; }
		if (i > 24 && i < 45) { age_width[i] = dy * 0.05; }
		if (i > 44 && i < 60) { age_width[i] = dy / 15.0; }
		if (i > 59 && i < 70) { age_width[i] = dy * 0.1; }
		if (i > 69 && i < 78) { age_width[i] = dy * 0.125; }
		if (i > 77 && i < 83) { age_width[i] = dy * 0.2; }
		if (i > 82 && i < 99) { age_width[i] = dy * 0.25; }
		if (i > 98 && i < 119) { age_width[i] = dy * 0.5; }
		if (i > 118 && i < 129) { age_width[i] = dy; }
		if (i > 128) { age_width[i] = dy * 2.0; }
		age_rate[i] = 1.0 / age_width[i];
		rem_rate[i] = age_rate[i] + eta;
		age[i] = i == 0 ? 0.0 : (age[i - 1] + age_width[i - 1]);
		if (agefinding0 == 1 && age[i] >= age_c0 * dy)
		{
			na_c0 = i;
			agefinding0 = 0;
			agefinding1 = 1;
		}
		if (agefinding1 == 1 && age[i] >= age_c1 * dy)
		{
			na_c1 = i;
			agefinding1 = 0;
			agefinding2 = 1;
		}
		if (agefinding2 == 1 && age[i] >= age_c2 * dy)
		{
			na_c2 = i;
			agefinding2 = 0;
		}
		agey0[i] = age[i] * inv_dy;
	}

	int na_c = na_c2 - na_c0 + 1;			//Total number of age categories within cohort
	int n_cats_c = na_c * num_het;			//Total number of age and heterogeneity categories in cohort
	for (i = 0; i < na - 1; i++) { agey1[i] = age[i + 1] * inv_dy; }

	//Calculate the proportion in each age group using discrete time, to match the equilibrium from the differential eqns
	double *den = (double *)malloc(na * sizeof(double));
	double *foi_age = (double *)malloc(na * sizeof(double));
	double densum = 0.0;
	for (i = 0; i < na; i++)
	{
		den[i] = i == 0 ? (1.0 / (1.0 + (age_rate[0] / eta))) : (age_rate[i - 1] * den[i - 1] / (rem_rate[i]));
		densum += den[i];
		foi_age[i] = 1.0 - (rho*exp(-age[i] / a0));
	}
	double inv_densum = 1.0 / densum;

	//Heterogeneity
	double *rel_foi = (double *)malloc(num_het * sizeof(double));
	for (j = 0; j < num_het; j++) { rel_foi[j] = num_het == 0 ? 1.0 : exp(-(sigma2*0.5) + (pow(sigma2, 0.5)*het_x[j])); }
	double *fd = (double *)malloc(na * sizeof(double));
	for (i = 0; i < na; i++) { fd[i] = 1.0 - ((1.0 - fd0) / (1.0 + (pow(age[i] / ad0, gammad)))); }

	//Set up additional arrays-------------------------------------------------------------------------------------------------------------------------------------

	double* EIR_data = (double*)malloc(tmax_i2 * sizeof(double));			//Daily entomological inoculation rate values saved for cohort calculations
	double *IB_data = (double *)malloc(n_cats * sizeof(double));
	double *IC_data = (double *)malloc(n_cats * sizeof(double));
	double *ID_data = (double *)malloc(n_cats * sizeof(double));
	double *cprob = (double *)malloc(n_cats_c * sizeof(double));				//Cumulative probability distribution used to determine age/heterogenity category to place randomly generated patients
	int *patients_status = (int *)malloc(6 * sizeof(int));						//Tally of number of patients in each status category (S,T,D,A,U,P)
	int *pcr_test_results = (int *)malloc(n_divs*n_patients * sizeof(int));		//PCR test results for individual patients at tinterval2 checkpoints
	double *pcr_distribution = (double *)malloc((n_divs + 1) * sizeof(double)); //Frequency distribution of number of positive PCR test results over entire trial
	double *slide_prev_coh_values = (double *)malloc(n_divs * sizeof(double));	//Values of slide prevalence (cohort) at tinterval2 checkpoints
	double *clin_inc_coh_values = (double *)malloc(n_divs * sizeof(double));	//Values of clinical incidence (cohort) at tinterval2 checkpoints
	double *pcr_fraction_values = (double *)malloc(n_divs * sizeof(double));	//Values of proportion of patients in cohort giving positive PCR test results at tinterval2 checkpoints

	//Load data from file------------------------------------------------------------------------------------------------------------------------------------------

	double *IB_input = (double *)malloc(n_cats* n_mv_values * sizeof(double));
	double *IC_input = (double *)malloc(n_cats* n_mv_values * sizeof(double));
	double *ID_input = (double *)malloc(n_cats* n_mv_values * sizeof(double));
	double* EIR_input = (double*)malloc(n_EIR_values * tmax_i2 * sizeof(double));
	int *run_mvn = (int *)malloc(n_clusters * sizeof(int));
	int *run_int = (int *)malloc(n_clusters * sizeof(int));
	int *run_data_num = (int *)malloc(n_clusters * sizeof(int));

	input_EIR = fopen(input_filename1.c_str(), "r");//fopen_s(&input_EIR, input_filename1.c_str(), "r"); //
	input_immunity = fopen(input_filename2.c_str(), "r");//fopen_s(&input_immunity, input_filename2.c_str(), "r"); //
	input_clusters = fopen(input_filename3.c_str(), "r");//fopen_s(&input_clusters, input_filename3.c_str(), "r"); //
	if (input_EIR == NULL || input_immunity == NULL || input_clusters == NULL)
	{
		Rcout << "Missing input file.";
		goto finish;
	}
	
	Rcout << "\nLoading EIR data.\n";
	R_FlushConsole();
	fseek(input_EIR, 8 + (6 * tmax_i2) - ((tmax_i2 < 10 ? tmax_i2 : 10)) + (tmax_i2 < 100 ? 0 : (tmax_i2 - 100)) + (tmax_i2 < 1000 ? 0 : (tmax_i2 - 1000)), SEEK_SET); //Skip over header of EIR input file
	pos = 0;
	for (i = 0; i < n_EIR_values; i++)
	{
		fscanf(input_EIR, "%d", &num);
		for (j = 0; j < tmax_i2; j++)
		{
			fscanf(input_EIR, "%lf", &EIR_input[pos]);
			pos++;
		}
	}
	fclose(input_EIR);

	Rcout << "\nLoading immunity data.\n";
	R_FlushConsole();
	fseek(input_immunity, 8 + (15 * n_cats) - ((n_cats < 10 ? n_cats : 10) * 3) + (n_cats < 100 ? 0 : (n_cats - 100) * 3) + (n_cats < 1000 ? 0 : (n_cats - 1000) * 3), SEEK_CUR); //Skip over header of immunity input file
	for (i = 0; i < n_mv_values; i++)
	{
		fscanf(input_immunity, "%d", &num);
		ij = i*n_cats;
		for (j = 0; j < n_cats; j++) { fscanf(input_immunity, "%lf", &IB_input[ij + j]); }
		for (j = 0; j < n_cats; j++) { fscanf(input_immunity, "%lf", &IC_input[ij + j]); }
		for (j = 0; j < n_cats; j++) { fscanf(input_immunity, "%lf", &ID_input[ij + j]); }
	}
	fclose(input_immunity);

	Rcout << "\nLoading cluster data.\n";
	R_FlushConsole();
	fseek(input_clusters, 51, SEEK_CUR); //Skip over header of cluster input file
	for (n_c = 0; n_c < n_clusters; n_c++)
	{
		fscanf(input_clusters, "%d", &num);
		fscanf(input_clusters, "%lf", &v);
		fscanf(input_clusters, "%d", &run_mvn[n_c]);
		fscanf(input_clusters, "%lf", &v);
		fscanf(input_clusters, "%lf", &v);
		fscanf(input_clusters, "%d", &run_int[n_c]);
		fscanf(input_clusters, "%lf", &v);
		fscanf(input_clusters, "%d", &run_data_num[n_c]);
		Rcout << "\nCluster " << n_c << ":\tmvn=" << run_mvn[n_c] << "\tint_num=" << run_int[n_c] << "\trun_data_num=" << run_data_num[n_c];
	}
	fclose(input_clusters);

	//Set up conditions for one or more runs-------------------------------------------------------------------------------------------------------------------------------------

	Rcout << "\nSetting up conditions.\n";
	R_FlushConsole();
	if (n_clusters > 1)
	{
		output1 = fopen(output_filename1.c_str(), "w");//fopen_s(&output1, output_filename1.c_str(), "w"); //
		fprintf(output1, "n_patients\t%i\nCluster\tmvn\tint_num", n_patients);
		for (div = 1; div <= n_divs; div++) { fprintf(output1, "\tprev_c%i", div); }
		for (div = 1; div <= n_divs; div++) { fprintf(output1, "\tinc_c%i", div); }
		for (div = 1; div <= n_divs; div++) { fprintf(output1, "\tpcr_p%i", div); }
		fclose(output1);
		output2 = fopen(output_filename2.c_str(), "w");//fopen_s(&output2, output_filename2.c_str(), "w"); //
		fprintf(output2, "n_patients\t%i\nCluster\tmvn\tint_num", n_patients);
		for (div = 0; div <= n_divs; div++) { fprintf(output2, "\tF%i", div); }
		fclose(output2);
	}

	for (n_c = 0; n_c < n_clusters; n_c++) //This for() loop runs all the simulations, jumping to start_run for each one and jumping back to run_complete when it is finished
	{
		mvn = run_mvn[n_c];
		int_num = run_int[n_c];
		data_num = run_data_num[n_c];
		Rcout << "\nProcessing cluster " << n_c << ":\tmvn=" << mvn << "\tint_num=" << int_num << "\tdata_num=" << data_num;
		goto start_run;
	run_complete:
		Rcout << "\nCluster " << n_c << " complete.\n";
	}

	goto finish;

	//-------------------------------------------------------------------------------------Compute over time--------------------------------------------------------------------------------------------------------

start_run:

	ij = data_num*tmax_i2;
	for (i = 0; i < tmax_i2; i++) { EIR_data[i] = EIR_input[ij + i]; }
	ij = mvn*n_cats;
	for (i = 0; i < n_cats; i++)
	{
		pos = ij + i;
		IB_data[i] = IB_input[pos];
		IC_data[i] = IC_input[pos];
		ID_data[i] = ID_input[pos];
	}

	//Set up cumulative probability distribution from which to draw randomly generated patients

	for (pos = 0; pos < n_cats_c; pos++) { cprob[pos] = 0.0; }
	for (i = na_c0; i <= na_c1; i++)
	{
		for (j = 0; j < num_het; j++)
		{
			pos = ((i - na_c0)*num_het) + j;
			cprob[pos] = (den[i] * het_wt[j])*inv_densum;
			if (pos > 0) { cprob[pos] += cprob[pos - 1]; }
		}
	}

	cprobmax = cprob[((na_c1 - na_c0)*num_het) + num_het - 1];
	for (pos = 0; pos < n_cats_c; pos++) { if (cprob[pos] > 0.0) { cprob[pos] /= cprobmax; } }
	for (i = 0; i < 6; i++) { patients_status[i] = 0; }

	dt = 0.1; //Time increment for cohort data. Like that used for the main population, this may be adjusted automatically to prevent errors
	tflag = 0;

	for (n = 0; n < n_patients; n++) //Randomly determine patient's age and heterogeneity category using cprob
	{
		prob = randgen(4);		
		for (i = na_c0; i <= na_c1; i++)
		{
			for (j = 0; j < num_het; j++)
			{
				if (prob <= cprob[((i - na_c0)*num_het) + j]) { goto found; }
			}
		}
		found:
		patient[n].na = i;
		patient[n].num_het = j;
	}

	restart_coh:
	for (div = 0; div < n_divs; div++)
	{
		clin_inc_coh_values[div] = 0.0;
		slide_prev_coh_values[div] = 0.0;
		pcr_fraction_values[div] = 0.0;
		pcr_distribution[div] = 0.0;
	}
	pcr_distribution[n_divs] = 0.0;

	if (n_clusters == 1)
	{
		output1 = fopen(output_filename1.c_str(), "w");//fopen_s(&output1, output_filename1.c_str(), "w"); //
		output2 = fopen(output_filename2.c_str(), "w");//fopen_s(&output2, output_filename2.c_str(), "w"); //
		fprintf(output1, "\n\nCluster:\t%i\tmvn:\t%i\tint_num\t%i\nt\tS\tT\tD\tA\tU\tP\tEIR\tprev\tclin_inc\tpcr_fraction", n_c, mvn, int_num);
		fprintf(output2, "\n\nCluster:\t%i\tmvn:\t%i\tint_num\t%i\npatient", n_c, mvn, int_num);
		for (i = 0; i < n_divs; i++) { fprintf(output2, "\tPCR_test_%i", i); }
		fclose(output1);
		fclose(output2);
	}
	p_multiplier = 1.0 / dt;
	ntmax = intdiv(tmax, dt);
	Rcout << "\nRunning individual model; dt=" << dt << "; ntmax=" << ntmax;
	R_FlushConsole();

	for (n = 0; n < n_patients; n++)
	{
		patient[n].status = 0; //All patients start in group S, representing pre-intervention treatment
		pos = (mvn*n_cats) + (patient[n].na*num_het) + patient[n].num_het;
		patient[n].IB = IB_input[pos];
		patient[n].IC = IC_input[pos];
		patient[n].ID = ID_input[pos];
		patient[n].infected = 0;
		patient[n].delay = 0.0;
	}

	t_mark1 = tinterval1;
	t_mark2 = 0.0;// tinterval2;
	div = 0;
	interval = 0;
	rate_db_t = rate_db*dt;
	rate_dc_t = rate_dc*dt;
	rate_dd_t = rate_dd*dt;
	for (nt = 0; nt <= ntmax; nt++)
	{
		t = nt*dt;
		EIRd = EIR_data[interval];
		EIRt = EIRd*dt;
		for (n = 0; n < n_patients; n++)
		{
			infected = 0;
			i = patient[n].na;
			j = patient[n].num_het;
			IB_cur = patient[n].IB;
			IC_cur = patient[n].IC;
			ID_cur = patient[n].ID;
			p_inf_bite = EIRt * rel_foi[j] * foi_age[i]; //Probability of an infectious bite

			//Check that dt has not been set too high
			if (p_inf_bite > 0.9) { tflag = 1; }
			if (tflag > 0)
			{
				tflag = 0;
				if (dt <= 0.0025) { goto end_coh; }
				dt = max(0.0025, dt*0.5);
				goto restart_coh;
			}
			if (randgen(4) <= p_inf_bite)
			{
				patient[n].IB += IB_boost;
				p_inf_from_bite = bh*(IB_cur > 0.0 ? ((bmin_rev / (1.0 + pow(IB_cur * inv_IB0, kb))) + bmin) : 1.0); //Probability of an infection from an infectious bite
				if (randgen(4) <= p_inf_from_bite)
				{
					infected = 1;
					patient[n].IC += IC_boost;
					patient[n].ID += ID_boost;
				}
			}

			if (patient[n].infected == 1)
			{
				patient[n].delay -= dt;
				if (patient[n].delay <= 0.0)
				{
					patient[n].infected = 0;
					patient[n].delay = 0.0;
					p_clin_inf = phi0 * ((phi1_rev / (1.0 + pow(IC_cur * inv_IC0, kc))) + phi1); //Probability of a clinical infection
					if (randgen(4) <= p_clin_inf)
					{/*clinical infection*/
						clin_inc_coh_values[div] += dv_p2;
						if (randgen(2) <= prop_T_c) {/*to T*/ patient[n].status = 1; }
						else {/*to D*/ patient[n].status = 2; }
					}
					else {/*to A*/ patient[n].status = 3; }
				}
			}
			else
			{
				switch (patient[n].status)
				{
				case 0: //S
				{
					if (infected == 1)
					{
						patient[n].infected = 1;
						patient[n].delay = dur_E;
					}
				}
					break;
				case 1: //T
				{
					if (randgen(4)*p_multiplier <= rT) {/*to P*/ patient[n].status = 5; }
				}
					break;
				case 2: //D
				{
					if (randgen(4)*p_multiplier <= rD) {/*to A*/ patient[n].status = 3; }
				}
					break;
				case 3: //A
				{
					if (infected == 1)
					{
						patient[n].infected = 1;
						patient[n].delay = dur_E;
					}
					else
					{
						if (randgen(4)*p_multiplier <= rA0) {/*to U*/ patient[n].status = 4; }
					}
				}
					break;
				case 4: //U
				{
					if (infected == 1)
					{
						patient[n].infected = 1;
						patient[n].delay = dur_E;
					}
					else
					{
						if (randgen(4)*p_multiplier <= rU) {/*to S*/ patient[n].status = 0; }
					}
				}
					break;
				case 5: //P
				{
					if (randgen(4)*p_multiplier <= rP) {/*to S*/ patient[n].status = 0; }
				}
					break;
				default: { Rcout << "\nError in switch 1!\n"; }
				}
			}
			//Natural decrease in immunity
			patient[n].IB -= IB_cur*rate_db_t;
			patient[n].IC -= IC_cur*rate_dc_t;
			patient[n].ID -= ID_cur*rate_dd_t;
		}

		if (t >= t_mark1)
		{
			interval++;
			t_mark1 += tinterval1;
		}

		//Output patient data and benchmarks
		if (t >= t_mark2)
		{
			for (i = 0; i < 6; i++) { patients_status[i] = 0; }
			pos = div*n_patients;
			for (n = 0; n < n_patients; n++)
			{
				patients_status[patient[n].status]++;
				switch (patient[n].status)
				{
				case 1: //T
				{
					pcr_test_results[pos] = 1;
					pcr_fraction_values[div] += dv_p1;
					slide_prev_coh_values[div] += dv_p1;
				}
					break;
				case 2: //D
				{
					pcr_test_results[pos] = 1;
					pcr_fraction_values[div] += dv_p1;
					slide_prev_coh_values[div] += dv_p1;
				}
					break;
				case 3: //A
				{
					pcr_test_results[pos] = 1;
					pcr_fraction_values[div] += dv_p1;
					if (randgen(4) <= dmin + (dmin_rev / (1.0 + (fd[patient[n].na] * pow(patient[n].ID * inv_ID0, kd))))) { slide_prev_coh_values[div] += dv_p1; }
				}
					break;
				case 4: //U
				{
					pcr_test_results[pos] = 1;
					pcr_fraction_values[div] += dv_p1;
				}
					break;
				default: {pcr_test_results[pos] = 0; }
				}
				pos++;
			}

			if (n_clusters == 1)
			{
				output1 = fopen(output_filename1.c_str(), "a");//fopen_s(&output1, output_filename1.c_str(), "a"); //
				fprintf(output1, "\n%.0f", t);
				for (i = 0; i < 6; i++) { fprintf(output1, "\t%i", patients_status[i]); }
				fprintf(output1, "\t%.3e\t%.3e\t%.3e\t%.3e", EIRd, slide_prev_coh_values[div], clin_inc_coh_values[div], pcr_fraction_values[div]);
				Rcout << "\n\tt = " << t << " EIR = " << EIRd << " prev = " << slide_prev_coh_values[div] << " inc = " << clin_inc_coh_values[div] << " pcr = " << pcr_fraction_values[div];	
				fclose(output1);
			}
			//Rcout << "\nt = " << t << "\tEIR = " << EIRd << "\tprev = " << slide_prev_coh_values[div];
			t_mark2 += tinterval2;
			div++;
			for (n = 0; n < n_patients; n++)
			{
				j = 0;
				for (i = 0; i < n_divs; i++) { if (pcr_test_results[(i*n_patients) + n] == 1) { j++; } }
				pcr_distribution[j] += dv_p1;
			}
		}

	}
	
	output2 = fopen(output_filename2.c_str(), "a");//fopen_s(&output2, output_filename2.c_str(), "a"); //
	if (n_clusters > 1)
	{
		output1 = fopen(output_filename1.c_str(), "a");//fopen_s(&output1, output_filename1.c_str(), "a"); //
		fprintf(output1, "\n%i\t%i\t%i", n_c, mvn, int_num);
		for (div = 0; div < n_divs; div++) { fprintf(output1, "\t%.3e", slide_prev_coh_values[div]); }
		for (div = 0; div < n_divs; div++) { fprintf(output1, "\t%.3e", clin_inc_coh_values[div]); }
		for (div = 0; div < n_divs; div++) { fprintf(output1, "\t%.3e", pcr_fraction_values[div]); }
		fclose(output1);	
		fprintf(output2, "\n%i\t%i\t%i", n_c, mvn, int_num);
		for (div = 0; div <= n_divs; div++) { fprintf(output2, "\t%.3e", pcr_distribution[div] / n_divs); }
	}
	else
	{
		for (n = 0; n < n_patients; n++)
		{
			fprintf(output2, "\n%i", n);
			for (div = 0; div < n_divs; div++) { fprintf(output2, "\t%i", pcr_test_results[(div*n_patients) + n]); }
		}
	}
	fclose(output2);

	end_coh:

	goto run_complete;

	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	finish:
	if (input_EIR != NULL) { fclose(input_EIR); }
	if (input_immunity != NULL) { fclose(input_immunity); }
	if (input_clusters != NULL) { fclose(input_clusters); }
	if (output1 != NULL) { fclose(output1); }
	if (output2 != NULL) { fclose(output2); }
	Rcout << "\nProgram complete\n";
	return 0;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int intdiv(double x,double y) //Outputs product of two doubles, rounding down, as an integer
{
	int i;
	double z = x / y;
	i = z - fmod(z, 1.0);
	return i;
}

double randgen(int n)
{
	int m;
	double value = 0.0;
	for (m = 1; m <= n; m++) { value += (rand() % 10)*pow(10.0, -m); }
	return value;
}

// converts input from List format to int format.
int rcpp_to_int(SEXP x) { return as<int>(x); }

// converts input from List format to double format.
double rcpp_to_double(SEXP x) { return as<double>(x); }

// converts input from List format to string format.
string rcpp_to_string(SEXP x) { return as<string>(x); }