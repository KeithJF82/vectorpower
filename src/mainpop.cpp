// CPP file for simulation of malaria intervention trial - main population portion
#include <Rcpp.h>
#include "main.h"
using namespace Rcpp;
using namespace std;

// Main program------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// [[Rcpp::export]]
int rcpp_mainpop(List params, List trial_params)
{
	int n_run, mv_num, int_num, i, j, ij, n, nt, nt_E, div, nt_latgam, nt_latmosq, pos, pos2, ntmax, dur_Ei, latgami, latmosqi, tflag, interval, increment, data_saved, flag_int, np;
	double mv0, t, t_int, mu_net_irs, prevent_net_irs, prop_T_rev, av0, muv1, K0, year_day, rain, KL, beta, Surv1, EIRd, mv, incv, incv0, FOIv, FOIv0, age0, age1, t_mark1, t_mark2, EIR_sum, dt2;
	double mu_atsb, cov_nets, cov_irs, prop_T, rN, rNW, dNW, rI, rIW, dIW, dIF, EL, LL, PL, Sv1, Ev1, Iv1, Rnorm, rconst, Re0, Re1, Re2, Re3, Re4, Im0, Im1, Im2, Im3, Im4;
	double S_sum, T_sum, D_sum, A_sum, U_sum, P_sum, H, H_inv, delS, delT, delD, delA1, delA2, delP, delU1, delU2, inv_x, inv_KL;
	double S_cur, T_cur, D_cur, A_cur, U_cur, P_cur, FOI_cur, foi_age_cur, clin_inc_cur, ICA_cur, ICM_cur, IB_cur, ID_cur, EIR_cur, delSv1, EL_cur, LL_cur, PL_cur, dEL, dLL, Sv_cur, Ev_cur, Iv_cur;
	FILE* benchmark_summary = NULL;
	FILE* benchmarks_by_age = NULL;
	FILE* endpoint_data = NULL;
	FILE* data_EIR = NULL;
	FILE* data_immunity = NULL;
	FILE* input = NULL;
	FILE* frain = NULL;
	int flag_error = 0;

	//Constants (TODO: Make global)----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	double dy = 365.0;			// Days in a year
	double tinterval1 = 1.0;	// Interval between calculations of current prevalence / incidence values
	int na = 145;				// Number of age categories in main population
	int num_het = 9;			// Number of heterogeneity categories
	double het_x[] = { -4.5127459, -3.2054291, -2.076848, -1.0232557, 0.0, 1.0232557, 2.076848, 3.2054291, 4.5127459 }; // Biting heterogeneity
	double het_wt[] = { 0.00002235, 0.00278914, 0.04991641, 0.2440975, 0.40634921, 0.2440975, 0.04991641, 0.00278914, 0.00002235 }; // Fractions in each heterogeneity category
	double KMIN = 0.0005;		// Minimum value of larval birth rate coefficient K0 to prevent zero or negative mosquito numbers

	// Load environment/vector/parasite/etc. parameters from R------------------------------------------------------------------------------------------------------------------------------------------------

	string input_filename = rcpp_to_string(trial_params["input_file"]);				// Name of file containing mosquito and human input data
	string rain_filename = rcpp_to_string(trial_params["rain_file"]);				// Name of file containing rainfall parameter values
	int n_mv0_values = rcpp_to_int(trial_params["n_mv0_values"]);								// Number of lines of data in input file

	double mu_atsb_def = rcpp_to_double(params["mu_atsb_def"]);						// Attractive targeted sugar bait (ATSB) killing rate associated with data
	double cov_nets_def = rcpp_to_double(params["cov_nets_def"]);					// Proportion of people protected by bednets associated with data
	double cov_irs_def = rcpp_to_double(params["cov_irs_def"]);						// Proportion of people protected by interior residual spraying (IRS) associated with data
	double prop_T_def = rcpp_to_double(params["prop_T_def"]);						// Proportion of clinical cases successfully treated in main population associated with data

	double Q0 = rcpp_to_double(params["Q0"]);										// Default anthropophagy
	double muv0 = rcpp_to_double(params["muv0"]);									// Default mosquito birth / death rate
	double phi_bite_bed = rcpp_to_double(params["phi_bite_bed"]);					// Default proportion of bites delivered to people in bed
	double p_repel_net_entry = rcpp_to_double(params["p_repel_net_entry"]);			// Default probability of mosquito being repelled by bednet before entering house
	double p_repel_net_bed = rcpp_to_double(params["p_repel_net_bed"]);				// Default probability of mosquito being repelled by bednet after entering house but before biting
	double p_kill_net = rcpp_to_double(params["p_kill_net"]);						// Default probability of mosquito being killed by bednet while trying to bite
	double phi_bite_house = rcpp_to_double(params["phi_bite_house"]);				// Default proportion of bites delivered to people in houses
	double p_repel_irs_entry = rcpp_to_double(params["p_repel_irs_entry"]);			// Default probability of mosquito being repelled by IRS before entering house
	double p_repel_irs_bed = rcpp_to_double(params["p_repel_irs_bed"]);				// Default probability of mosquito being repelled by IRS after entering house but before biting
	double p_kill_irs1 = rcpp_to_double(params["p_kill_irs1"]);						// Default probability of mosquito being killed by IRS while trying to bite
	double p_kill_irs2 = rcpp_to_double(params["p_kill_irs2"]);						// Default probability of mosquito being killed by IRS after biting

	double rho = rcpp_to_double(params["rho"]);										// Age-dependent biting parameter
	double a0 = rcpp_to_double(params["a0"]);										// Age-dependent biting parameter
	double sigma2 = rcpp_to_double(params["sigma2"]);								// Variance of log heterogeneity in biting
	double dur_E = rcpp_to_double(params["dur_E"]);									// Latent period
	double dur_T = rcpp_to_double(params["dur_T"]);									// Average time to recover from disease and parasitaemia when treated
	double dur_D = rcpp_to_double(params["dur_D"]);									// Average time to move from clinical disease to asymptomatic when not successfully treated
	double dur_A = rcpp_to_double(params["dur_A"]);									// Average time to move from asymptomatic patent infection to sub-patent infection
	double dur_U = rcpp_to_double(params["dur_U"]);									// Average time to recover from sub-patent infection
	double dur_P = rcpp_to_double(params["dur_P"]);									// Average time to leave protected state
	double latmosq = rcpp_to_double(params["latmosq"]);								// Time lag for infection in mosquitoes
	double latgam = rcpp_to_double(params["latgam"]);								// Time lag mosquitoes to become infectious

	// Infectivity to mosquitoes
	double cD = rcpp_to_double(params["cD"]);										// Untreated disease
	double cU = rcpp_to_double(params["cU"]);										// Sub-patent infection
	double cT = rcpp_to_double(params["cT"]);										// Treated disease
	double gamma_inf = rcpp_to_double(params["gamma_inf"]);							// Parameter for infectiousness of state A

	// Immunity parameters - infection
	double bh = rcpp_to_double(params["bh"]);										// Maximum probability due to no immunity
	double bmin = rcpp_to_double(params["bmin"]);									// Maximum relative reduction due to immunity
	double db = rcpp_to_double(params["db"]);										// Inverse of decay rate
	double kb = rcpp_to_double(params["kb"]);										// Shape parameter
	double ub = rcpp_to_double(params["ub"]);										// Duration in which immunity is not boosted
	double inv_IB0 = 1.0 / rcpp_to_double(params["IB0"]);							// Inverse of scale parameter

	// Immunity parameters - detection
	double dmin = rcpp_to_double(params["dmin"]);									// Minimum probability due to maximum immunity
	double dd = rcpp_to_double(params["dd"]);										// Inverse of decay rate
	double kd = rcpp_to_double(params["kd"]);										// Shape parameter
	double ud = rcpp_to_double(params["ud"]);										// Duration in which immunity is not boosted
	double ad0 = rcpp_to_double(params["ad0"]);										// Scale parameter relating age to immunity
	double fd0 = rcpp_to_double(params["fd0"]);										// Time scale at which immunity changes with age
	double gammad = rcpp_to_double(params["gammad"]);								// Shape parameter relating age to immunity
	double inv_ID0 = 1.0 / rcpp_to_double(params["ID0"]);							// Scale parameter

	// Immunity parameters - clinical disease
	double phi0 = rcpp_to_double(params["phi0"]);									// Maximum probability due to no immunity
	double phi1 = rcpp_to_double(params["phi1"]);									// Maximum relative reduction due to no immunity
	double dc = rcpp_to_double(params["dc"]);										// Inverse of decay rate
	double kc = rcpp_to_double(params["kc"]);										// Shape parameter
	double uc = rcpp_to_double(params["uc"]);										// Duration in which immunity is not boosted
	double P_IC_M = rcpp_to_double(params["P_IC_M"]);								// 
	double dm = rcpp_to_double(params["dm"]);										// Inverse of decay rate of maternal immunity
	double inv_IC0 = 1.0 / rcpp_to_double(params["IC0"]);							// Scale parameter

	// Larval parameters	
	double mue = rcpp_to_double(params["mue"]);										// Death rate of early instar
	double mul = rcpp_to_double(params["mul"]);										// Death rate of late instar
	double mup = rcpp_to_double(params["mup"]);										// Death rate of pupae
	double de = rcpp_to_double(params["de"]);										// Duration of early instar
	double dl = rcpp_to_double(params["dl"]);										// Duration of late instar
	double dp = rcpp_to_double(params["dp"]);										// Duration of pupae
	double eov = rcpp_to_double(params["eov"]);										// Eggs per ovulation
	double gammal = rcpp_to_double(params["gammal"]);								// Density dependence term
	double dgon = rcpp_to_double(params["dgon"]);									// Gonotrophic cycle length

	// Load trial parameter data from R--------------------------------------------------------------------------------------------------------------------------------------------

	int n_mv0_start = rcpp_to_int(trial_params["n_mv0_start"]);								// Number of first set of data to use in input file (must be less than or equal to n_mv0_end)
	int n_mv0_end = rcpp_to_int(trial_params["n_mv0_end"]);									// Number of last set of data to use in input file (must be less than n_mv0_values)
	if (n_mv0_start > n_mv0_end || n_mv0_end >= n_mv0_values)
	{
		Rcout << "\nData set number values error: n_mv0_values=" << n_mv0_values << " n_mv0_start=" << n_mv0_start << " n_mv0_end=" << n_mv0_end << "\n";
		flag_error = 1;
	}
	int int_v_varied = rcpp_to_int(trial_params["int_v_varied"]);					// Intervention parameter given variable value (0=model parameter set, 1=ATSB kill rate, 2=bednet coverage, 3=IRS coverage)
	if (int_v_varied < 0 || int_v_varied > 3)
	{
		Rcout << "\nError in intervention parameter choice: int_v_varied=" << int_v_varied << "\n";
		flag_error = 1;
	}
	int n_int_values = rcpp_to_int(trial_params["n_int_values"]);					// Number of intervention parameter values to use
	double int_v_min = rcpp_to_double(trial_params["int_v_min"]);					// Minimum intervention parameter value to use (sole value used if n_int_values=1)
	double int_v_max = rcpp_to_double(trial_params["int_v_max"]);					// Maximum intervention parameter value to use (unused if n_int_values=1)
	if (int_v_max < int_v_min || int_v_min < 0.0)
	{
		Rcout << "\nError in intervention parameter values: int_v_min=" << int_v_min << " int_v_max=" << int_v_max << "\n";
		flag_error = 1;
	}
	double date_start = rcpp_to_double(trial_params["date_start"]);					// Day of the year when simulation starts (should be set based on input_file)						 
	double date_int = rcpp_to_double(trial_params["date_int"]);						// Day of the year when intervention trial starts (period from date_start to date_atsb used to equilibrate lagged FOI data)
	if (date_int < date_start)
	{
		Rcout << "\nError in start dates: date_start=" << date_start << " date_int=" << date_int << "\n";
		flag_error = 1;
	}
	double tinterval2 = rcpp_to_double(trial_params["time_interval"]);				// Interval between calculations of regular test probabilities and averaged prevalence/incidence values
	int n_pts = rcpp_to_int(trial_params["n_pts"]);									// Number of data points to be taken (starting at time 0 after start of intervention and thereafter every tinterval2 days) 
	string file_summary = rcpp_to_string(trial_params["file_summary"]);
	string file_benchmarks = rcpp_to_string(trial_params["file_benchmarks"]);
	string file_endpoints = rcpp_to_string(trial_params["file_endpoints"]);
	string file_EIRd = rcpp_to_string(trial_params["file_EIRd"]);
	string file_imm_start = rcpp_to_string(trial_params["file_imm_start"]);

	// Constant derived values--------------------------------------------------------------------------------------------------------------------------

	int n_mv_values = n_mv0_end - n_mv0_start + 1;									// Number of mv0 values to be used
	int n_cats = na * num_het;												// Total number of age/heterogeneity categories in main population
	int size1 = na * sizeof(double);										// Array dimensions
	int size2 = num_het * sizeof(double);
	int size3 = n_cats * sizeof(double);
	int size4 = n_mv_values * sizeof(double);
	int size5 = n_cats * size4;
	int n_runmax = n_int_values * n_mv_values;								// Total number of simulations to run
	double int_start_time = date_int - date_start;							// Time after start of run at which ATSB starts being used
	double tmax_c = (n_pts - 1) * tinterval2;								// Duration of intervention use in days
	int tmax_c_i = intdiv(tmax_c, 1.0);										// Duration of intervention use as integer
	double tmax = int_start_time + tmax_c;									// Total duration of run in days
	double inv_dy = 1.0 / dy;												// Inverse of days in a year (fraction of a year in a day)
	double trig_coeff1 = 2.0 * 3.1415926536 * inv_dy;						// Coefficient used in calculation of rainfall
	double trig_coeff2 = trig_coeff1 * 2.0;
	double trig_coeff3 = trig_coeff1 * 3.0;
	double trig_coeff4 = trig_coeff1 * 4.0;
	double trig_coeff5 = trig_coeff1 * 5.0;
	double rT = 1.0 / dur_T;												// Rate of recovery from disease and parasitaemia when treated
	double rD = 1.0 / dur_D;												// Rate of moving from clinical disease to asymptomatic when not successfully treated
	double rA0 = 1.0 / dur_A;												// Rate of moving from asymptomatic patent infection to sub-patent infection
	double rU = 1.0 / dur_U;												// Rate of recovery from sub - patent infection
	double rP = 1.0 / (dur_P - dur_T);										// Rate of leaving protected state
	double cDU = cD - cU;
	double bmin_rev = 1.0 - bmin;
	double dmin_rev = 1.0 - dmin;
	double phi1_rev = 1.0 - phi1;
	double rate_dc = 1.0 / dc;												// Decay rate of clinical disease immunity
	double rate_dcm = 1.0 / dm;												// Decay rate of maternal immunity
	double rate_db = 1.0 / db;												// Decay rate of infection immunity
	double rate_dd = 1.0 / dd;												// Decay rate of detection immunity
	double rate_de = 1.0 / de;
	double rate_dl = 1.0 / dl;
	double* slide_prev_age= (double*)malloc(size1);
	double* pcr_prev_age = (double*)malloc(size1);
	double* clin_inc_age = (double*)malloc(size1);
	double* age_width = (double*)malloc(size1);								// Widths of age categories in days
	double* age_rate = (double*)malloc(size1);								// Rates of transition between age categories
	double* rem_rate = (double*)malloc(size1);								// Total loss rate from each age category (ageing + death)
	double* age = (double*)malloc(size1);									// Age values corresponding to starts of age categories (days)
	double* agey0 = (double*)malloc(size1);									// Age values corresponding to starts of age categories (years)
	double* agey1 = (double*)malloc(size1);									// Age values corresponding to ends of age categories (years)
	double eta = 1.0 / (21.0 * dy);											// Death rate, for exponential population distribution
	double* den = (double*)malloc(size1);
	double* foi_age = (double*)malloc(size1);
	int* age20i = (int*)malloc(na * sizeof(int));

	// Set age category parameters
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
		if (i > 118 && i < 129) { age_width[i] = dy * 1.0; }
		if (i > 128) { age_width[i] = dy * 2.0; }
		age_rate[i] = 1.0 / age_width[i];
		rem_rate[i] = age_rate[i] + eta;
		age[i] = i == 0 ? 0.0 : (age[i - 1] + age_width[i - 1]);
		agey0[i] = age[i] * inv_dy;
	}
	for (i = 0; i < na - 1; i++) { agey1[i] = age[i + 1] * inv_dy; }
	for (i = 0; i < na; i++)
	{
		den[i] = i == 0 ? (1.0 / (1.0 + (age_rate[0] / eta))) : (age_rate[i - 1] * den[i - 1] / (rem_rate[i]));
		foi_age[i] = 1.0 - (rho * exp(-age[i] / a0));
		if (i == 0) { age20i[i] = 0; }
		else
		{
			if (age[i] >= 20.0 * dy && age[i - 1] < 20.0 * dy) { age20i[i] = i; }
			else { age20i[i] = age20i[i - 1]; }
		}
	}
	int age20u = age20i[na - 1];
	int age20l = age20u - 1;
	double age_20_factor = (((20 * dy) - age[age20l] - (0.5 * age_width[age20l])) * 2.0) / (age_width[age20l] + age_width[age20u]);

	double* rel_foi = (double*)malloc(size2);
	double* b = (double*)malloc(size3);					// Probability of infection from an infectious bite
	double* FOI = (double*)malloc(size3);				// Force of infection (mosquito -> human)
	double* FOI0 = (double*)malloc(size3);				// Force of infection at start of trial
	double* p_det = (double*)malloc(size3);				// Probability of asymptomatic infection being detected
	double* fd = (double*)malloc(size1);				//
	double* phi = (double*)malloc(size3);				// Probability of infection becoming clinical
	double* ICM_init = (double*)malloc(size2);			// Value of maternally inherited immunity in newborns
	double* rate_ibaq = (double*)malloc(size3);	//
	double* rate_clinaq = (double*)malloc(size3);//
	double* rate_detaq = (double*)malloc(size3);//
	double* cA = (double*)malloc(size3);				//
	double* FOIvij = (double*)malloc(size3);			//Force of infection (human -> mosquito)

	for (j = 0; j < num_het; j++) { rel_foi[j] = num_het == 0 ? 1.0 : exp(-(sigma2 / 2.0) + (pow(sigma2, 0.5) * het_x[j])); }	// Relative biting rate in heterogeneity category i
	double omega = 1.0 - ((rho * eta) / (eta + (1.0 / a0)));																	// Normalising constant for biting rate by age
	double* x_I = (double*)malloc(size1);																						// Intermediate variables for calculating immunity functions
	double* inv_x_I = (double*)malloc(size1);
	double inv_omega = 1.0 / omega;
	x_I[0] = den[0] / eta;
	inv_x_I[0] = 1.0 / x_I[0];
	fd[0] = 1.0 - ((1.0 - fd0) / (1.0 + pow(age[0] / ad0, gammad)));
	for (i = 1; i < na; i++)
	{
		x_I[i] = den[i] / (den[i - 1] * age_rate[i - 1]);
		inv_x_I[i] = 1.0 / x_I[i];
		fd[i] = 1.0 - ((1.0 - fd0) / (1.0 + pow(age[i] / ad0, gammad)));
	}

	double rnyears = 8.0;
	double ppyear = 64.0;
	double beta_larval = (eov * muv0 * exp(-muv0 * dgon)) / (1.0 - exp(-muv0 * dgon));
	double b_lambda = ((gammal * mul) / mue) - (de / dl) + ((gammal - 1.0) * mul * de);
	double lambda = (-0.5 * b_lambda) + pow((0.25 * pow(b_lambda, 2.0)) + ((gammal * beta_larval * mul * de) / (2.0 * mue * muv0 * dl * (1.0 + (dp * mup)))), 0.5);
	double term = (lambda / (mue * de)) - (1.0 / (mul * dl)) - 1.0;
	double rain_coeff = 0.0;
	double dt = 0.1;																			// Time increment (default initial value)
	double dp_int = n_int_values == 1 ? 0.0 : (int_v_max - int_v_min) / (n_int_values - 1);		// Increment between variable parameter values
	
	// Set up additional parameters-------------------------------------------------------------------------------------------------------------------------------------

	double* S = (double*)malloc(size3);													// Susceptible humans
	double* T = (double*)malloc(size3);													// Humans in treatment for clinical symptoms
	double* D = (double*)malloc(size3);													// Humans with clinical symptoms not receiving treatment
	double* A = (double*)malloc(size3);													// Humans with asymptomatic patent infections
	double* U = (double*)malloc(size3);													// Humans with asymptomatic sub-patent infections
	double* P = (double*)malloc(size3);													// Humans in post-treatment stage
	double* slide_prev = (double*)malloc(size3);										// Prevalence
	double* pcr_prev = (double*)malloc(size3);											// Prevalence (PCR test)
	double* clin_inc = (double*)malloc(size3);											// Clinical incidence
	double* ICA = (double*)malloc(size3);												// Immunity level against clinical infection (adult)
	double* ICM = (double*)malloc(size3);												// Immunity level against clinical infection (maternally inherited)
	double* IB = (double*)malloc(size3);												// Immunity level against infection from infectious bite
	double* ID = (double*)malloc(size3);												// Immunity level determining p_det

	double* FOI_lag = (double*)malloc(n_cats * 5000 * sizeof(double));					// Time-lagged force of mosquito -> human infection (to account for incubation period)
	double* FOIv_lag = (double*)malloc(5000 * sizeof(double));							// Time-lagged force of human -> mosquito infection (to account for incubation period)
	double* incv_lag = (double*)malloc(5000 * sizeof(double));							// Time-lagged incubation of malaria in infected mosquito

	double* EIR_output_values = (double*)malloc(n_pts * sizeof(double));				// Values of entomological inoculation rate at tinterval2 checkpoints
	double* slide_prev_output_values = (double*)malloc(n_pts * sizeof(double));			// Values of slide prevalence at tinterval2 checkpoints
	double* pcr_prev_output_values = (double*)malloc(n_pts * sizeof(double));			// Values of PCR prevalence at tinterval2 checkpoints
	double* clin_inc_output_values = (double*)malloc(n_pts * sizeof(double));			// Values of clinical incidence at tinterval2 checkpoints
	double* EIR_data = (double*)malloc((tmax_c_i + 1) * sizeof(double));				// Daily entomological inoculation rate values saved for cohort calculations

	double* mv_input = (double*)malloc(size4);
	double* EL_input = (double*)malloc(size4);
	double* LL_input = (double*)malloc(size4);
	double* PL_input = (double*)malloc(size4);
	double* Sv_input = (double*)malloc(size4);
	double* Ev_input = (double*)malloc(size4);
	double* Iv_input = (double*)malloc(size4);
	double* S_input = (double*)malloc(size5);
	double* T_input = (double*)malloc(size5);
	double* D_input = (double*)malloc(size5);
	double* A_input = (double*)malloc(size5);
	double* U_input = (double*)malloc(size5);
	double* P_input = (double*)malloc(size5);
	double* ICA_input = (double*)malloc(size5);
	double* ICM_input = (double*)malloc(size5);
	double* IB_input = (double*)malloc(size5);
	double* ID_input = (double*)malloc(size5);
	double* IB_output = (double*)malloc(size3);
	double* IC_output = (double*)malloc(size3);
	double* ID_output = (double*)malloc(size3);

	//Load data from files------------------------------------------------------------------------------------------------------------------------------------------

	fopen_s(&frain, rain_filename.c_str(), "r");
	if (frain == NULL)
	{
		Rcout << "Rain file " << frain << " not found.\n";
		flag_error = 1;
	}
	else
	{
		fscanf(frain, "%lf", &Rnorm);
		fscanf(frain, "%lf", &rconst);
		fscanf(frain, "%lf", &Re0);
		fscanf(frain, "%lf", &Re1);
		fscanf(frain, "%lf", &Re2);
		fscanf(frain, "%lf", &Re3);
		fscanf(frain, "%lf", &Re4);
		fscanf(frain, "%lf", &Im0);
		fscanf(frain, "%lf", &Im1);
		fscanf(frain, "%lf", &Im2);
		fscanf(frain, "%lf", &Im3);
		fscanf(frain, "%lf", &Im4);
		rain_coeff = ((dy * Rnorm) / (14.0 * rnyears * ppyear));
		fclose(frain);
	}

	fopen_s(&input, input_filename.c_str(), "r");
	if (input == NULL)
	{
		Rcout << "Input file " << input << " not found.\n";
		flag_error = 1;
	}
	else
	{
		fseek(input, 31 + (46 * n_cats) - (10 * (n_cats < 10 ? n_cats : 10)) + (n_cats < 100 ? 0 : (n_cats - 100) * 10) + (n_cats < 1000 ? 0 : (n_cats - 1000) * 10), SEEK_CUR);// Skip over header of input file
		ij = 7 + (10 * n_cats);
		if (n_mv0_start > 0) //Skip over unused lines
		{
			for (i = 0; i < n_mv0_start; i++)
			{
				fscanf(input, "%d", &i);
				for (j = 0; j < ij; j++) { fscanf(input, "%lf", &t); }
			}
		}
		for (i = 0; i < n_mv_values; i++)
		{
			fscanf(input, "%d", &n);
			fscanf(input, "%lf", &mv_input[i]);
			fscanf(input, "%lf", &EL_input[i]);
			fscanf(input, "%lf", &LL_input[i]);
			fscanf(input, "%lf", &PL_input[i]);
			fscanf(input, "%lf", &Sv_input[i]);
			fscanf(input, "%lf", &Ev_input[i]);
			fscanf(input, "%lf", &Iv_input[i]);
			ij = i * n_cats;
			for (j = 0; j < n_cats; j++) { fscanf(input, "%lf", &S_input[ij + j]); }
			for (j = 0; j < n_cats; j++) { fscanf(input, "%lf", &T_input[ij + j]); }
			for (j = 0; j < n_cats; j++) { fscanf(input, "%lf", &D_input[ij + j]); }
			for (j = 0; j < n_cats; j++) { fscanf(input, "%lf", &A_input[ij + j]); }
			for (j = 0; j < n_cats; j++) { fscanf(input, "%lf", &U_input[ij + j]); }
			for (j = 0; j < n_cats; j++) { fscanf(input, "%lf", &P_input[ij + j]); }
			for (j = 0; j < n_cats; j++) { fscanf(input, "%lf", &ICA_input[ij + j]); }
			for (j = 0; j < n_cats; j++) { fscanf(input, "%lf", &ICM_input[ij + j]); }
			for (j = 0; j < n_cats; j++) { fscanf(input, "%lf", &IB_input[ij + j]); }
			for (j = 0; j < n_cats; j++) { fscanf(input, "%lf", &ID_input[ij + j]); }
		}
		fclose(input);
	}

	// Set up conditions for one or more runs-------------------------------------------------------------------------------------------------------------------------------------

		if (flag_error == 1)//If one or more errors were found during setup, skip simulations and go to end.
	{
		Rcout << "\nOne or more errors found. Ending program.\n";
		goto finish;
	}
	else
	{
		Rcout << "\ninput file name: " << input_filename << "\nn_mv0_start/max: " << n_mv0_start << "/" << n_mv0_end << "\nint_v_varied: " << int_v_varied << "\nn_int_values: " << n_int_values << "\nint_v_min/max: " << int_v_min << "/" << int_v_max;
		Rcout << "\nSimulation start date: " << date_start << "\nIntervention start date: " << date_int;
		Rcout << "\nBenchmarks summary file: " << file_summary << "\nBenchmark details file: " << file_benchmarks << "\nEndpoints file: " << file_endpoints << "\nDaily EIR file: " << file_EIRd << "\nStarting immunity file: " << file_imm_start;
		R_FlushConsole();
	}

	fopen_s(&benchmark_summary, file_summary.c_str(), "w");
	fprintf(benchmark_summary, "run\tmv_num\tmv0\tint_num\tmu_atsb\tcov_nets\tcov_irs\tprop_T");
	for (i = 0; i < 4; i++)
	{
		fprintf(benchmark_summary, "\t-");
		for (j = 0; j < n_pts; j++) { fprintf(benchmark_summary, "\t%.0f", j * tinterval2); }
	}
	fclose(benchmark_summary);

	fopen_s(&benchmarks_by_age, file_benchmarks.c_str(), "w");
	fprintf(benchmarks_by_age, "run\tmv_num\tint_num\tnpt\tEIR");
	for (i = 0; i < na; i++) { fprintf(benchmarks_by_age, "\tslide_prev%i",i); }
	for (i = 0; i < na; i++) { fprintf(benchmarks_by_age, "\tpcr_prev%i", i); }
	for (i = 0; i < na; i++) { fprintf(benchmarks_by_age, "\tclin_inc%i", i); }
	fclose(benchmarks_by_age);

	fopen_s(&data_EIR, file_EIRd.c_str(), "w");
	fprintf(data_EIR, "n_run");
	for (t = 0.0; t <= tmax_c; t += tinterval1) { fprintf(data_EIR, "\tEIR%.0f", t); }
	fclose(data_EIR);
	fopen_s(&data_immunity, file_imm_start.c_str(), "w");
	fprintf(data_immunity, "mv_num");
	for (i = 0; i < n_cats; i++) { fprintf(data_immunity, "\tIB%i", i); }
	for (i = 0; i < n_cats; i++) { fprintf(data_immunity, "\tIC%i", i); }
	for (i = 0; i < n_cats; i++) { fprintf(data_immunity, "\tID%i", i); }
	fclose(data_immunity);
	if (int_v_varied == 0) // Generate header for detailed endpoint data file
	{
		fopen_s(&endpoint_data, file_endpoints.c_str(), "w");
		fprintf(endpoint_data, "n_run\tmv0\tEL\tLL\tPL\tSv1\tEv1\tIv1");
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\tS%i", i); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\tT%i", i); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\tD%i", i); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\tA%i", i); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\tU%i", i); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\tP%i", i); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\tICA%i", i); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\tICM%i", i); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\tIB%i", i); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\tID%i", i); }
		fclose(endpoint_data);
	}

	Rcout << "\n\nDefault trial parameter values:\nCase treatment rate=" << prop_T_def << "\nATSB killing rate=" << mu_atsb_def << "\nBednet coverage=" << cov_nets_def << "\nIRS coverage=" << cov_irs_def;
	for (n_run = 0; n_run < n_runmax; n_run++)// This for() loop runs all the simulations, jumping to start_run for each one and jumping back to run_complete when it is finished
	{
		int_num = n_run % n_int_values;
		mv_num = (n_run - int_num) / n_int_values;
		mv0 = mv_input[mv_num];
		if (int_num == 0) { dt = 0.05 / sqrt(mv0); }

		mu_atsb = mu_atsb_def;
		cov_nets = cov_nets_def;
		cov_irs = cov_irs_def;
		prop_T = prop_T_def;

		goto start_run;
	run_complete:
		Rcout << "\nRun " << n_run << " complete.\n";
		R_FlushConsole();
	}

	goto finish;

	//-----------------------------------------------------------------------------Set up initial data for run--------------------------------------------------------------------

start_run:

	prop_T_rev = 1.0 - prop_T;											// Probability of treatment not being received by clinical cases		
	rN = cov_nets * p_repel_net_entry;
	rNW = cov_nets * p_repel_net_bed;
	dNW = cov_nets * p_kill_net;
	rI = cov_irs * p_repel_irs_entry;
	rIW = cov_irs * p_repel_irs_bed;
	dIW = cov_irs * p_kill_irs1;
	dIF = cov_irs * p_kill_irs2;
	K0 = (mv0 * 2.0 * dl * muv0 * (1.0 + (dp * mup)) * (gammal * (lambda + 1.0))) / term;																				//Larval birth rate coefficient, set based on mv0
	prevent_net_irs = phi_bite_house - ((phi_bite_house - phi_bite_bed) * (1.0 - rI) * (1.0 - rIW - dIW) * (1.0 - dIF)) -// Proportion of bites prevented by bednets and/or interior residual spraying
		(phi_bite_bed * (1.0 - rI) * (1.0 - rIW - dIW) * (1.0 - dIF) * (1.0 - rN) * (1.0 - rNW - dNW));
	mu_net_irs = 0.3333 * Q0 * (((phi_bite_house - phi_bite_bed) * (1.0 - rI) * (dIW + ((1.0 - rIW - dIW) * dIF))) +				// Additional death rate due to bednets and/or interior residual spraying
		(phi_bite_bed * (1.0 - rI) * (1.0 - rN) * (dIW + ((1.0 - rIW - dIW) * (dNW + ((1.0 - rNW - dNW) * dIF))))));
	av0 = 0.3333 * Q0 * (1.0 - prevent_net_irs);																																				// Biting rate on humans / mosquito
	muv1 = muv0 + mu_atsb + mu_net_irs;																																							// Total mosquito death rate
	Rcout << "\nRun " << n_run << ": Data line=" << (n_mv0_start + mv_num) << " Mosquito density=" << mv0 << " Intervention number=" << int_num;

restart_run:
	dt2 = dt * 0.1;
	
	ntmax = intdiv(tmax, dt) + 1;
	dur_Ei = intdiv(dur_E, dt) + 1;
	latgami = intdiv(latgam, dt) + 1;
	latmosqi = intdiv(latmosq, dt) + 1;

	year_day = fmod(date_start, dy); // Day of the year at the start of the run
	rain = rain_coeff * (rconst + (2.0 * (
		(Re0 * cos(trig_coeff1 * year_day)) + (Im0 * sin(trig_coeff1 * year_day)) +
		(Re1 * cos(trig_coeff2 * year_day)) + (Im1 * sin(trig_coeff2 * year_day)) +
		(Re2 * cos(trig_coeff3 * year_day)) + (Im2 * sin(trig_coeff3 * year_day)) +
		(Re3 * cos(trig_coeff4 * year_day)) + (Im3 * sin(trig_coeff4 * year_day)) +
		(Re4 * cos(trig_coeff5 * year_day)) + (Im4 * sin(trig_coeff5 * year_day))
		)));
	KL = max(KMIN, K0 * rain);

	EL = EL_input[mv_num];
	LL = LL_input[mv_num];
	PL = PL_input[mv_num];
	Sv1 = Sv_input[mv_num];
	Ev1 = Ev_input[mv_num];
	Iv1 = Iv_input[mv_num];

	EIRd = av0 * Iv1 * inv_omega;
	pos = mv_num * n_cats;
	for (i = 0; i < n_cats; i++)
	{
		S[i] = S_input[pos];
		T[i] = T_input[pos];
		D[i] = D_input[pos];
		A[i] = A_input[pos];
		U[i] = U_input[pos];
		P[i] = P_input[pos];
		ICA[i] = ICA_input[pos];
		ICM[i] = ICM_input[pos];
		IB[i] = IB_input[pos];
		ID[i] = ID_input[pos];
		pos++;
	}

	H = arraysum(S, n_cats) + arraysum(T, n_cats) + arraysum(D, n_cats) + arraysum(A, n_cats) + arraysum(U, n_cats) + arraysum(P, n_cats);
	H_inv = 1.0 / H;
	pos = 0;
	FOIv0 = 0.0;
	for (i = 0; i < na; i++) //Run through age categories
	{
		foi_age_cur = foi_age[i];
		for (j = 0; j < num_het; j++) //Run through heterogeneity categories
		{
			// Normalize S,T,D,A,U,P to sum to 1.0
			S[pos] *= H_inv;
			T[pos] *= H_inv;
			D[pos] *= H_inv;
			A[pos] *= H_inv;
			U[pos] *= H_inv;
			P[pos] *= H_inv;
			EIR_cur = EIRd * foi_age_cur * rel_foi[j];
			b[pos] = bh * ((bmin_rev / (1.0 + pow(IB[pos] * inv_IB0, kb))) + bmin);
			FOI0[pos] = EIR_cur * (IB[pos] > 0.0 ? b[pos] : bh);
			p_det[pos] = dmin + (dmin_rev / (1.0 + (fd[i] * pow(ID[pos] * inv_ID0, kd))));
			cA[pos] = cU + (cDU * pow(p_det[pos], gamma_inf));
			FOIv0 += foi_age_cur * av0 * ((cT * T[pos]) + (cD * D[pos]) + (cA[pos] * A[pos]) + (cU * U[pos])) * rel_foi[j] * inv_omega;
			pos++;
		}
	}

	Surv1 = exp(-muv1 * latmosq);
	incv0 = FOIv0 * Sv1 * Surv1;
	for (i = 0; i < 5000; i++)
	{
		FOIv_lag[i] = FOIv0;
		incv_lag[i] = incv0;
		ij = i * n_cats;
		for (pos = 0; pos < n_cats; pos++) { FOI_lag[ij + pos] = FOI0[pos]; }
	}

	//-------------------------------------------------------------------------------------Compute over time--------------------------------------------------------------------------------------------------------

	div = 0;
	interval = 0;
	increment = 0;
	t_mark1 = tinterval1;
	t_mark2 = 0.0;
	EIR_sum = 0.0;
	for (n = 0; n <= tmax_c_i; n++) { EIR_data[n] = 0.0; }
	data_saved = 0;
	flag_int = 0;
	for (nt = 0; nt <= ntmax; nt++)
	{
		tflag = 0;											//Flag indicates whether dt needs to be altered
		t = nt * dt;										// Time since start of simulation
		t_int = t - int_start_time;							// Time since start of intervention use
		nt_E = nt % dur_Ei;
		nt_latgam = nt % latgami;
		nt_latmosq = nt % latmosqi;
		year_day = fmod(t + date_start, dy);				// Day of the year
		mv = Sv1 + Ev1 + Iv1;								// Total number of mosquitoes

		if (flag_int == 0 && t_int >= 0.0)					// Value of varied trial parameters set after start of intervention
		{
			flag_int = 1;
			switch (int_v_varied)
			{
			case 0: {}
					break;
			case 1:
			{
				mu_atsb = int_v_min + (int_num * dp_int);
				Rcout << "\nIntervention begun. ATSB kill rate=" << mu_atsb;
			}
			break;
			case 2:
			{
				cov_nets = int_v_min + (int_num * dp_int);
				rN = cov_nets * p_repel_net_entry;
				rNW = cov_nets * p_repel_net_bed;
				dNW = cov_nets * p_kill_net;
				Rcout << "\nIntervention begun. Bednet coverage=" << cov_nets;
			}
			break;
			case 3:
			{
				cov_irs = int_v_min + (int_num * dp_int);
				rI = cov_irs * p_repel_irs_entry;
				rIW = cov_irs * p_repel_irs_bed;
				dIW = cov_irs * p_kill_irs1;
				dIF = cov_irs * p_kill_irs2;
				Rcout << "\nIntervention begun. IRS coverage=" << cov_irs;
			}
			break;
			default: {Rcout << "\n\tint_v_varied error\n"; }
			}
			prevent_net_irs = phi_bite_house - ((phi_bite_house - phi_bite_bed) * (1.0 - rI) * (1.0 - rIW - dIW) * (1.0 - dIF)) -
				(phi_bite_bed * (1.0 - rI) * (1.0 - rIW - dIW) * (1.0 - dIF) * (1.0 - rN) * (1.0 - rNW - dNW));
			mu_net_irs = 0.3333 * Q0 * (((phi_bite_house - phi_bite_bed) * (1.0 - rI) * (dIW + ((1.0 - rIW - dIW) * dIF))) +
				(phi_bite_bed * (1.0 - rI) * (1.0 - rN) * (dIW + ((1.0 - rIW - dIW) * (dNW + ((1.0 - rNW - dNW) * dIF))))));
			av0 = 0.3333 * Q0 * (1.0 - prevent_net_irs);
			muv1 = muv0 + mu_atsb + mu_net_irs;
		}

		rain = rain_coeff * (rconst + (2.0 * (
			(Re0 * cos(trig_coeff1 * year_day)) + (Im0 * sin(trig_coeff1 * year_day)) +
			(Re1 * cos(trig_coeff2 * year_day)) + (Im1 * sin(trig_coeff2 * year_day)) +
			(Re2 * cos(trig_coeff3 * year_day)) + (Im2 * sin(trig_coeff3 * year_day)) +
			(Re3 * cos(trig_coeff4 * year_day)) + (Im3 * sin(trig_coeff4 * year_day)) +
			(Re4 * cos(trig_coeff5 * year_day)) + (Im4 * sin(trig_coeff5 * year_day))
			)));
		KL = max(K0 * rain, KMIN);
		beta = ((0.5 * PL) / dp);								// Rate of spawning of adult mosquitoes
		EIRd = av0 * Iv1 * inv_omega;							// Daily entomological inoculation rate

		for (j = 0; j < num_het; j++) { ICM_init[j] = P_IC_M * (ICA[(age20l * num_het) + j] + age_20_factor * (ICA[(age20u * num_het) + j] - ICA[(age20l * num_het) + j])); }
		pos = 0;
		ij = nt_E * n_cats;
		for (i = 0; i < na; i++) //Run through age categories
		{
			foi_age_cur = foi_age[i];
			for (j = 0; j < num_het; j++) //Run through heterogeneity categories
			{
				EIR_cur = EIRd * foi_age_cur * rel_foi[j];																							// Entomological inoculation rate
				b[pos] = bh * (bmin_rev / (1.0 + pow(IB[pos] * inv_IB0, kb)) + bmin);
				FOI[pos] = t > dur_E ? FOI_lag[ij + pos] : FOI0[pos];																				// Force of infection (mosquito -> human)
				FOI_lag[ij + pos] = EIR_cur * (IB[pos] > 0.0 ? b[pos] : bh);
				FOI_cur = FOI[pos];
				p_det[pos] = dmin + (dmin_rev / (1.0 + (fd[i] * pow(ID[pos] * inv_ID0, kd))));														// Probability of case detection
				rate_ibaq[pos] = EIR_cur / ((EIR_cur * ub) + 1.0);																					// Rate of change of
				rate_detaq[pos] = FOI_cur / ((FOI_cur * ud) + 1.0);																					// Rate of change of
				rate_clinaq[pos] = FOI_cur / ((FOI_cur * uc) + 1.0);																				// Rate of change of
				phi[pos] = phi0 * ((phi1_rev / (1.0 + pow((ICM[pos] + ICA[pos]) * inv_IC0, kc))) + phi1);											// Probability of clinical disease
				clin_inc[pos] = phi[pos] * FOI_cur * (S[pos] + A[pos] + U[pos]);																	// Clinical incidence rate
				cA[pos] = cU + ((cD - cU) * pow(p_det[pos], gamma_inf));
				FOIvij[pos] = foi_age_cur * av0 * ((cT * T[pos]) + (cD * D[pos]) + (cA[pos] * A[pos]) + (cU * U[pos])) * rel_foi[j] * inv_omega;	// Force of infection (human -> mosquito)
				pos++;
			}
		}

		//Key human data - main population
		for (i = na - 1; i >= 0; i--) //Run through age categories (backwards, so that forward ageing does not affect the current time point)
		{
			age1 = (i < na - 1 ? age_rate[i] : 0.0) + eta;
			inv_x = inv_x_I[i];
			pos = i * num_het;
			if (i > 0)
			{
				age0 = age_rate[i - 1];
				pos2 = pos - num_het;
			}
			else
			{
				age0 = 0.0;
				pos2 = 0;
			}
			for (j = 0; j < num_het; j++) //Run through heterogeneity categories
			{
				// Current parameter values
				S_cur = S[pos];
				T_cur = T[pos];
				D_cur = D[pos];
				A_cur = A[pos];
				U_cur = U[pos];
				P_cur = P[pos];
				FOI_cur = FOI[pos];
				clin_inc_cur = clin_inc[pos];
				ICA_cur = ICA[pos];
				ICM_cur = ICM[pos];
				IB_cur = IB[pos];
				ID_cur = ID[pos];

				// Changes within time increment 
				delS = FOI_cur * S_cur;
				delT = rT * T_cur;
				delD = rD * D_cur;
				delA1 = rA0 * A_cur;
				delA2 = FOI_cur * A_cur;
				delP = rP * P_cur;
				delU1 = rU * U_cur;
				delU2 = FOI_cur * U_cur;

				//Apply changes
				S[pos] += dt * (-delS + delP + delU1 - (age1 * S_cur) + (age0 * S[pos2]) + (i == 0 ? eta * het_wt[j] : 0.0));
				T[pos] += dt * ((prop_T * clin_inc_cur) - delT - (age1 * T_cur) + (age0 * T[pos2]));
				D[pos] += dt * ((prop_T_rev * clin_inc_cur) - delD - (age1 * D_cur) + (age0 * D[pos2]));
				A[pos] += dt * (((1.0 - phi[pos]) * (delS + delA2 + delU2)) + delD - delA1 - delA2 - (age1 * A_cur) + (age0 * A[pos2]));
				U[pos] += dt * (delA1 - delU1 - delU2 - (age1 * U_cur) + (age0 * U[pos2]));
				P[pos] += dt * (delT - delP - (age1 * P_cur) + (age0 * P[pos2]));
				ICA[pos] += dt * (rate_clinaq[pos] - (ICA_cur * rate_dc) - (inv_x * (i == 0 ? ICA_cur : (ICA_cur - ICA[pos2]))));
				ICM[pos] += dt * (-(ICM_cur * rate_dcm) - (inv_x * (i == 0 ? (ICM_cur - ICM_init[j]) : (ICM_cur - ICM[pos2]))));
				IB[pos] += dt * (rate_ibaq[pos] - (IB_cur * rate_db) - (inv_x * (i == 0 ? IB_cur : (IB_cur - IB[pos2]))));
				ID[pos] += dt * (rate_detaq[pos] - (ID_cur * rate_dd) - (inv_x * (i == 0 ? ID[j] : (ID_cur - ID[pos2]))));

				pos++;
				pos2++;
			}
		}

		S_sum = arraysum(S, n_cats);
		T_sum = arraysum(T, n_cats);
		D_sum = arraysum(D, n_cats);
		A_sum = arraysum(A, n_cats);
		U_sum = arraysum(U, n_cats);
		P_sum = arraysum(P, n_cats);
		H = S_sum + T_sum + D_sum + A_sum + U_sum + P_sum;// Total normalized number of humans; this should always sum to 1.0

		// Save daily immunity and EIR data for cohort calculations
		if (t_int >= 0.0)
		{
			EIR_sum += EIRd;
			increment++;
			if (data_saved == 0)
			{
				for (pos = 0; pos < n_cats; pos++)
				{
					IB_output[pos] = IB[pos];
					IC_output[pos] = ICA[pos] + ICM[pos];
					ID_output[pos] = ID[pos];
				}
				data_saved = 1;
			}
		}
		
		FOIv = t > latgam ? FOIv_lag[nt_latgam] : FOIv0;
		FOIv_lag[nt_latgam] = arraysum(FOIvij, n_cats);
		incv = t > latmosq ? incv_lag[nt_latmosq] : incv0;
		Surv1 = exp(-muv1 * latmosq);
		incv_lag[nt_latmosq] = FOIv * Sv1 * Surv1;
		inv_KL = 1.0 / KL;

		//Larval data calculation; calculated in smaller time increments
		for (np = 0; np < 10; np++)
		{
			EL_cur = EL;
			LL_cur = LL;
			PL_cur = PL;
			dEL = EL * rate_de;
			dLL = LL * rate_dl;
			EL += dt2 * ((beta_larval * mv) - (mue * (1.0 + ((EL_cur + LL_cur) * inv_KL)) * EL_cur) - dEL);
			LL += dt2 * (dEL - (mul * (1.0 + ((gammal * (EL_cur + LL_cur)) * inv_KL)) * LL_cur) - dLL);
			PL += dt2 * (dLL - (mup * PL_cur) - (PL_cur / dp));
		}

		//Adult mosquito data calculation
		Sv_cur = Sv1;
		Ev_cur = Ev1;
		Iv_cur = Iv1;
		delSv1 = FOIv * Sv_cur;
		Sv1 += dt * (-(muv1 * Sv_cur) - delSv1 + beta);
		Ev1 += dt * (-(muv1 * Ev_cur) + delSv1 - incv);
		Iv1 += dt * (-(muv1 * Iv_cur) + incv);

		// Check that time increment has not been set too large
		if ((H - 1.0) * (H - 1.0) > 0.01) { tflag = 1; }
		if (EL < 0.0) { tflag = 2; }
		if (Sv1 < 0.0) { tflag = 3; }
		if (tflag > 0)
		{
			if (dt <= 0.0025) { goto end_run; }
			dt = max(0.0025, dt * 0.9);
			Rcout << "\ntflag = " << tflag << " Adjusting time increment. New dt = " << dt << "\n";
			goto restart_run;
		}

		// Outputs
		if (t >= t_mark1)
		{

			if (t_int > 0.0)
			{
				if (interval <= tmax_c_i) { EIR_data[interval] = EIR_sum / (increment * 1.0); }
				EIR_sum = 0.0;
				increment = 0;
				interval++;
			}

			t_mark1 += tinterval1;
		}
		if (t_int >= t_mark2)
		{
			pos = 0;
			for (i = 0; i < na; i++) //Run through age categories
			{
				slide_prev_age[i] = 0.0;
				pcr_prev_age[i] = 0.0;
				clin_inc_age[i] = 0.0;
				for (j = 0; j < num_het; j++)
				{
					slide_prev[pos] = T[pos] + D[pos] + (p_det[pos] * A[pos]);
					pcr_prev[pos] = T[pos] + D[pos] + A[pos] + U[pos];
					slide_prev_age[i] += slide_prev[pos];
					pcr_prev_age[i] += pcr_prev[pos];
					clin_inc_age[i] += clin_inc[pos];
					pos++;
				}				
			}

			EIR_output_values[div] = EIRd;
			slide_prev_output_values[div] = arraysum(slide_prev_age, na);
			pcr_prev_output_values[div] = arraysum(pcr_prev_age, na);
			clin_inc_output_values[div] = dy * arraysum(clin_inc_age, na);
			div++;
			t_mark2 += tinterval2;
			fopen_s(&benchmarks_by_age, file_benchmarks.c_str(), "a");
			fprintf(benchmarks_by_age, "\n%i\t%i\t%i\t%i\t%.3e", n_run, mv_num, int_num, div, EIRd);
			for (i = 0; i < na; i++) { fprintf(benchmarks_by_age, "\t%.3e", slide_prev_age[i]); }
			for (i = 0; i < na; i++) { fprintf(benchmarks_by_age, "\t%.3e", pcr_prev_age[i]); }
			for (i = 0; i < na; i++) { fprintf(benchmarks_by_age, "\t%.3e", clin_inc_age[i]); }
			fclose(benchmarks_by_age);
		}

	}

end_run:

	fopen_s(&benchmark_summary, file_summary.c_str(), "a");
	fprintf(benchmark_summary, "\n%i\t%i\t%.3e\t%i\t%.3f\t%.3f\t%.3f\t%.3f\tEIR%i", n_run, mv_num, mv0, int_num, mu_atsb, cov_nets, cov_irs, prop_T, n_run);
	for (j = 0; j < n_pts; j++) { fprintf(benchmark_summary, "\t%.3e", EIR_output_values[j]); }
	fprintf(benchmark_summary, "\tslide_prev_all%i", n_run);
	for (j = 0; j < n_pts; j++) { fprintf(benchmark_summary, "\t%.3e", slide_prev_output_values[j]); }
	fprintf(benchmark_summary, "\tpcr_prev_all%i", n_run);
	for (j = 0; j < n_pts; j++) { fprintf(benchmark_summary, "\t%.3e", pcr_prev_output_values[j]); }
	fprintf(benchmark_summary, "\tclin_inc_all%i", n_run);
	for (j = 0; j < n_pts; j++) { fprintf(benchmark_summary, "\t%.3e", clin_inc_output_values[j]); }
	fclose(benchmark_summary);

	// Output detailed data to use as input for future computation
	fopen_s(&data_EIR, file_EIRd.c_str(), "a");
	fprintf(data_EIR, "\n%i", n_run);
	for (i = 0; i <= tmax_c_i; i++) { fprintf(data_EIR, "\t%.3e", EIR_data[i]); }
	fclose(data_EIR);
	if (int_num == 0)
	{
		fopen_s(&data_immunity, file_imm_start.c_str(), "a");
		fprintf(data_immunity, "\n%i", mv_num);
		for (i = 0; i < n_cats; i++) { fprintf(data_immunity, "\t%.3e", IB_output[i]); }
		for (i = 0; i < n_cats; i++) { fprintf(data_immunity, "\t%.3e", IC_output[i]); }
		for (i = 0; i < n_cats; i++) { fprintf(data_immunity, "\t%.3e", ID_output[i]); }
		fclose(data_immunity);
	}
	if (int_v_varied == 0)
	{
		fopen_s(&endpoint_data, file_endpoints.c_str(), "a");
		fprintf(endpoint_data, "\n%i\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e", n_run, mv0, EL, LL, PL, Sv1, Ev1, Iv1);
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\t%.3e", S[i]); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\t%.3e", T[i]); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\t%.3e", D[i]); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\t%.3e", A[i]); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\t%.3e", U[i]); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\t%.3e", P[i]); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\t%.3e", ICA[i]); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\t%.3e", ICM[i]); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\t%.3e", IB[i]); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\t%.3e", ID[i]); }
		fclose(endpoint_data);
	}

	goto run_complete;

	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

finish:
	if (input != NULL) { fclose(input); }
	if (frain != NULL) { fclose(frain); }
	if (benchmark_summary != NULL) { fclose(benchmark_summary); }
	if (benchmarks_by_age != NULL) { fclose(benchmarks_by_age); }
	if (endpoint_data != NULL) { fclose(endpoint_data); }
	if (data_EIR != NULL) { fclose(data_EIR); }
	if (data_immunity != NULL) { fclose(data_immunity); }
	Rcout << "\nMain population computations complete.\n";
	return 0;
}
