// CPP file for simulation of malaria intervention trial - main population portion
#include <Rcpp.h>
#include "main.h"
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
Rcpp::List rcpp_mainpop(List params, List inputs, List trial_params)
{
	int n_run, n_mv, n_int, i, j, ij, nt, nt_E, div, nt_latgam, nt_latmosq, pos, pos2, ntmax, dur_Ei, latgami, latmosqi, tflag, interval, increment, data_saved, flag_int, np;
	double mv0, param_int, t, t_int, mu_protections, r_bite_prevent, prop_T_rev, av0, K0, year_day, rain, KL, beta, EIRd, mv, incv, incv0, FOIv, FOIv0, age0, age1, t_mark1, t_mark2, EIR_sum, dt2;
	double muv1, Surv1, beta_larval1;
	double mu_atsb, cov_nets, cov_irs, cov_set, prop_T, rN, rNW, dNW, rI, rIW, dIW, dIF, EL, LL, PL, Sv1, Ev1, Iv1, H, /*H_inv, */delS, delT, delD, delA1, delA2, delP, delU1, delU2, inv_x, inv_KL;
	double S_cur, T_cur, D_cur, A_cur, U_cur, P_cur, FOI_cur, foi_age_cur, clin_inc_cur, ICA_cur, ICM_cur, IB_cur, ID_cur, EIR_cur, delSv1, EL_cur, LL_cur, PL_cur, dEL, dLL, Sv_cur, Ev_cur, Iv_cur;
	FILE* benchmarks_by_age = NULL;
	FILE* endpoint_data = NULL;
	FILE* data_EIR = NULL;
	FILE* data_immunity = NULL;
	double dSET = 0.0;
	double rSET = 0.0;

	// Constants ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	double dy = 365.0;			// Days in a year
	//double tinterval1 = 1.0;	// Interval between calculations of current prevalence / incidence values
	double KMIN = 0.0005;		// Minimum value of larval birth rate coefficient K0 to prevent zero or negative mosquito numbers

	// Load environment/vector/parasite/etc. parameters from R------------------------------------------------------------------------------------------------------------------------------------------------
	
	// Rainfall parameters
	double date_start = rcpp_to_double(params["date_start"]);
	double Rnorm = rcpp_to_double(params["Rnorm"]);
	double rconst = rcpp_to_double(params["rconst"]);
	double Re0 = rcpp_to_double(params["Re0"]);
	double Re1 = rcpp_to_double(params["Re1"]);
	double Re2 = rcpp_to_double(params["Re2"]);
	double Re3 = rcpp_to_double(params["Re3"]);
	double Re4 = rcpp_to_double(params["Re4"]);
	double Im0 = rcpp_to_double(params["Im0"]);
	double Im1 = rcpp_to_double(params["Im1"]);
	double Im2 = rcpp_to_double(params["Im2"]);
	double Im3 = rcpp_to_double(params["Im3"]);
	double Im4 = rcpp_to_double(params["Im4"]);

	// Heterogeneity parameters
	int num_het = rcpp_to_int(params["num_het"]);								// Number of age categories in main population
	vector<double> het_x = rcpp_to_vector_double(params["het_x"]);				// Biting heterogeneity
	vector<double> het_wt = rcpp_to_vector_double(params["het_wt"]);			// Fractions in each heterogeneity category

	// Age parameters
	int na = rcpp_to_int(params["na"]);											// Number of age categories in main population
	vector<double> age_width = rcpp_to_vector_double(params["age_width"]);		// Widths of age categories in days
	vector<double> age_rate = rcpp_to_vector_double(params["age_rate"]);		// Rates of transition between age categories
	vector<double> rem_rate = rcpp_to_vector_double(params["rem_rate"]);		// Total loss rate from each age category (ageing + death)
	vector<double> age = rcpp_to_vector_double(params["age"]);					// Age values corresponding to starts of age categories (days)
	vector<double> den = rcpp_to_vector_double(params["den"]);
	double eta = rcpp_to_double(params["eta"]);									// Death rate, for exponential population distribution

	double mu_atsb_def = rcpp_to_double(params["mu_atsb_def"]);					// Attractive targeted sugar bait (ATSB) killing rate associated with data
	double cov_nets_def = rcpp_to_double(params["cov_nets_def"]);				// Proportion of people protected by bednets associated with data
	double cov_irs_def = rcpp_to_double(params["cov_irs_def"]);					// Proportion of people protected by interior residual spraying (IRS) associated with data
	double cov_set_def = 0.0;													// Proportion of people protected by SET
	double prop_T_def = rcpp_to_double(params["prop_T_def"]);					// Proportion of clinical cases successfully treated in main population associated with data

	double Q0 = rcpp_to_double(params["Q0"]);									// Default anthropophagy
	double muv0 = rcpp_to_double(params["muv0"]);								// Default mosquito birth / death rate
	double phi_bite_bed = rcpp_to_double(params["phi_bite_bed"]);				// Default proportion of bites delivered to people in bed
	double p_repel_net_entry = rcpp_to_double(params["p_repel_net_entry"]);		// Default probability of mosquito being repelled by bednet before entering house
	double p_repel_net_bed = rcpp_to_double(params["p_repel_net_bed"]);			// Default probability of mosquito being repelled by bednet after entering house but before biting
	double p_kill_net = rcpp_to_double(params["p_kill_net"]);					// Default probability of mosquito being killed by bednet while trying to bite
	double phi_bite_house = rcpp_to_double(params["phi_bite_house"]);			// Default proportion of bites delivered to people in houses
	double p_repel_irs_entry = rcpp_to_double(params["p_repel_irs_entry"]);		// Default probability of mosquito being repelled by IRS before entering house
	double p_repel_irs_bed = rcpp_to_double(params["p_repel_irs_bed"]);			// Default probability of mosquito being repelled by IRS after entering house but before biting
	double p_kill_irs1 = rcpp_to_double(params["p_kill_irs1"]);					// Default probability of mosquito being killed by IRS while trying to bite
	double p_kill_irs2 = rcpp_to_double(params["p_kill_irs2"]);					// Default probability of mosquito being killed by IRS after biting
	double p_repel_set = 0.229;
	double p_kill_set = 0.118;

	double rho = rcpp_to_double(params["rho"]);									// Age-dependent biting parameter
	double a0 = rcpp_to_double(params["a0"]);									// Age-dependent biting parameter
	double sigma2 = rcpp_to_double(params["sigma2"]);							// Variance of log heterogeneity in biting
	double dur_E = rcpp_to_double(params["dur_E"]);								// Latent period
	double dur_T = rcpp_to_double(params["dur_T"]);								// Average time to recover from disease and parasitaemia when treated
	double dur_D = rcpp_to_double(params["dur_D"]);								// Average time to move from clinical disease to asymptomatic when not successfully treated
	double dur_A = rcpp_to_double(params["dur_A"]);								// Average time to move from asymptomatic patent infection to sub-patent infection
	double dur_U = rcpp_to_double(params["dur_U"]);								// Average time to recover from sub-patent infection
	double dur_P = rcpp_to_double(params["dur_P"]);								// Average time to leave protected state
	double latmosq = rcpp_to_double(params["latmosq"]);							// Time lag for infection in mosquitoes
	double latgam = rcpp_to_double(params["latgam"]);							// Time lag mosquitoes to become infectious

	// Infectivity to mosquitoes
	double cD = rcpp_to_double(params["cD"]);									// Untreated disease
	double cU = rcpp_to_double(params["cU"]);									// Sub-patent infection
	double cT = rcpp_to_double(params["cT"]);									// Treated disease
	double gamma_inf = rcpp_to_double(params["gamma_inf"]);						// Parameter for infectiousness of state A

	// Immunity parameters - infection
	double bh = rcpp_to_double(params["bh"]);									// Maximum probability due to no immunity
	double bmin = rcpp_to_double(params["bmin"]);								// Maximum relative reduction due to immunity
	double db = rcpp_to_double(params["db"]);									// Inverse of decay rate
	double kb = rcpp_to_double(params["kb"]);									// Shape parameter
	double ub = rcpp_to_double(params["ub"]);									// Duration in which immunity is not boosted
	double inv_IB0 = 1.0 / rcpp_to_double(params["IB0"]);						// Inverse of scale parameter

	// Immunity parameters - detection
	double dmin = rcpp_to_double(params["dmin"]);								// Minimum probability due to maximum immunity
	double dd = rcpp_to_double(params["dd"]);									// Inverse of decay rate
	double kd = rcpp_to_double(params["kd"]);									// Shape parameter
	double ud = rcpp_to_double(params["ud"]);									// Duration in which immunity is not boosted
	double ad0 = rcpp_to_double(params["ad0"]);									// Scale parameter relating age to immunity
	double fd0 = rcpp_to_double(params["fd0"]);									// Time scale at which immunity changes with age
	double gammad = rcpp_to_double(params["gammad"]);							// Shape parameter relating age to immunity
	double inv_ID0 = 1.0 / rcpp_to_double(params["ID0"]);						// Scale parameter

	// Immunity parameters - clinical disease
	double phi0 = rcpp_to_double(params["phi0"]);								// Maximum probability due to no immunity
	double phi1 = rcpp_to_double(params["phi1"]);								// Maximum relative reduction due to no immunity
	double dc = rcpp_to_double(params["dc"]);									// Inverse of decay rate
	double kc = rcpp_to_double(params["kc"]);									// Shape parameter
	double uc = rcpp_to_double(params["uc"]);									// Duration in which immunity is not boosted
	double P_IC_M = rcpp_to_double(params["P_IC_M"]);							// 
	double dm = rcpp_to_double(params["dm"]);									// Inverse of decay rate of maternal immunity
	double inv_IC0 = 1.0 / rcpp_to_double(params["IC0"]);						// Scale parameter

	// Larval parameters	
	double mue = rcpp_to_double(params["mue"]);									// Death rate of early instar
	double mul = rcpp_to_double(params["mul"]);									// Death rate of late instar
	double mup = rcpp_to_double(params["mup"]);									// Death rate of pupae
	double de = rcpp_to_double(params["de"]);									// Duration of early instar
	double dl = rcpp_to_double(params["dl"]);									// Duration of late instar
	double dp = rcpp_to_double(params["dp"]);									// Duration of pupae
	double eov = rcpp_to_double(params["eov"]);									// Eggs per ovulation
	double gammal = rcpp_to_double(params["gammal"]);							// Density dependence term
	double dgon = rcpp_to_double(params["dgon"]);								// Gonotrophic cycle length

	// Load trial parameter data from R--------------------------------------------------------------------------------------------------------------------------------------------		
	
	// Variable mu_atsb for special runs - TODO: Integrate into package
	/*vector<double> mu_atsb_var = {};*/
	int flag_dt_adjust = rcpp_to_int(trial_params["flag_dt_adjust"]);					// Flag indicating whether or not to adjust dt by mosquito density (0=no,1=yes)
	int flag_file = rcpp_to_int(trial_params["flag_file"]);								// Flag indicating whether results data to be saved to files
	int int_v_varied = rcpp_to_int(trial_params["int_v_varied"]);						// Intervention parameter given variable value (0=model parameter set, 1=ATSB kill rate, 2=bednet coverage, 3=IRS coverage, 4=SET coverage)
	int n_mv_values = rcpp_to_int(trial_params["n_mv_values"]);							// Number of mosquito density values to be used		
	vector<double> int_values = rcpp_to_vector_double(trial_params["int_values"]);		// Vector of intervention parameter values
	int n_int_values = int_values.size();												// Number of intervention parameter values to use
	double date_int = date_start + rcpp_to_double(trial_params["start_interval"]);		// Day of the year when intervention trial starts (period from date_start to date_int used to equilibrate lagged FOI data)
	vector<double> time_values = rcpp_to_vector_double(trial_params["time_values"]);	// Vector of time benchmark points
	int n_pts = rcpp_to_int(trial_params["n_pts"]);										// Number of data points to be taken
	string file_benchmarks = flag_file == 1 ? rcpp_to_string(trial_params["file_benchmarks"]) : "na";
	string file_endpoints = flag_file == 1 ? rcpp_to_string(trial_params["file_endpoints"]) : "na";
	string file_EIRd = flag_file == 1 ? rcpp_to_string(trial_params["file_EIRd"]) : "na";
	string file_imm_start = flag_file == 1 ? rcpp_to_string(trial_params["file_imm_start"]) : "na";

	// Constant derived values--------------------------------------------------------------------------------------------------------------------------

	int n_cats = na * num_het;													// Total number of age/heterogeneity categories in main population
	int size1 = na * sizeof(double);											// Array dimensions
	int size2 = num_het * sizeof(double);
	int size3 = n_cats * sizeof(double);
	int n_runmax = n_int_values * n_mv_values;									// Total number of simulations to run
	double int_start_time = date_int - date_start;								// Time after start of run at which ATSB starts being used
	double tmax_c = time_values[n_pts-1];										// Duration of intervention use in days
	int tmax_c_i = intdiv(tmax_c, 1.0);											// Duration of intervention use as integer
	double tmax = int_start_time + tmax_c;										// Total duration of run in days
	double inv_dy = 1.0 / dy;													// Inverse of days in a year (fraction of a year in a day)
	double trig_coeff1 = 2.0 * 3.1415926536 * inv_dy;							// Coefficient used in calculation of rainfall
	double trig_coeff2 = trig_coeff1 * 2.0;
	double trig_coeff3 = trig_coeff1 * 3.0;
	double trig_coeff4 = trig_coeff1 * 4.0;
	double trig_coeff5 = trig_coeff1 * 5.0;
	double rT = 1.0 / dur_T;													// Rate of recovery from disease and parasitaemia when treated
	double rD = 1.0 / dur_D;													// Rate of moving from clinical disease to asymptomatic when not successfully treated
	double rA0 = 1.0 / dur_A;													// Rate of moving from asymptomatic patent infection to sub-patent infection
	double rU = 1.0 / dur_U;													// Rate of recovery from sub - patent infection
	double rP = 1.0 / (dur_P - dur_T);											// Rate of leaving protected state
	double cDU = cD - cU;
	double bmin_rev = 1.0 - bmin;
	double dmin_rev = 1.0 - dmin;
	double phi1_rev = 1.0 - phi1;
	double rate_dc = 1.0 / dc;													// Decay rate of clinical disease immunity
	double rate_dcm = 1.0 / dm;													// Decay rate of maternal immunity
	double rate_db = 1.0 / db;													// Decay rate of infection immunity
	double rate_dd = 1.0 / dd;													// Decay rate of detection immunity
	double rate_de = 1.0 / de;
	double rate_dl = 1.0 / dl;
	double* slide_prev_age = (double*)malloc(size1);
	double* pcr_prev_age = (double*)malloc(size1);
	double* clin_inc_age = (double*)malloc(size1);
	double* foi_age = (double*)malloc(size1);
	int* age20i = (int*)malloc(na * sizeof(int));

	for (i = 0; i < na; i++)
	{
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
	double* b = (double*)malloc(size3);							// Probability of infection from an infectious bite
	double* FOI = (double*)malloc(size3);						// Force of infection (mosquito -> human)
	double* FOI0 = (double*)malloc(size3);						// Force of infection at start of trial
	double* p_det = (double*)malloc(size3);						// Probability of asymptomatic infection being detected
	double* fd = (double*)malloc(size1);						//
	double* phi = (double*)malloc(size3);						// Probability of infection becoming clinical
	double* ICM_init = (double*)malloc(size2);					// Value of maternally inherited immunity in newborns
	double* rate_ibaq = (double*)malloc(size3);					//
	double* rate_clinaq = (double*)malloc(size3);				//
	double* rate_detaq = (double*)malloc(size3);				//
	double* cA = (double*)malloc(size3);						//
	double* FOIvij = (double*)malloc(size3);					// Force of infection (human -> mosquito)

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
	double beta_larval0 = (eov * muv0 * exp(-muv0 * dgon)) / (1.0 - exp(-muv0 * dgon));
	double b_lambda = ((gammal * mul) / mue) - (de / dl) + ((gammal - 1.0) * mul * de);
	double lambda = (-0.5 * b_lambda) + pow((0.25 * pow(b_lambda, 2.0)) + ((gammal * beta_larval0 * mul * de) / (2.0 * mue * muv0 * dl * (1.0 + (dp * mup)))), 0.5);
	double term = (lambda / (mue * de)) - (1.0 / (mul * dl)) - 1.0;
	double dt = 0.1;												// Time increment (default initial value)
	double rain_coeff = ((dy * Rnorm) / (14.0 * rnyears * ppyear));	// Rainfall coefficient
		
	// Set up additional parameters-------------------------------------------------------------------------------------------------------------------------------------

	double* S = (double*)malloc(size3);									// Susceptible humans
	double* T = (double*)malloc(size3);									// Humans in treatment for clinical symptoms
	double* D = (double*)malloc(size3);									// Humans with clinical symptoms not receiving treatment
	double* A = (double*)malloc(size3);									// Humans with asymptomatic patent infections
	double* U = (double*)malloc(size3);									// Humans with asymptomatic sub-patent infections
	double* P = (double*)malloc(size3);									// Humans in post-treatment stage
	double* slide_prev = (double*)malloc(size3);						// Prevalence
	double* pcr_prev = (double*)malloc(size3);							// Prevalence (PCR test)
	double* clin_inc = (double*)malloc(size3);							// Clinical incidence
	double* ICA = (double*)malloc(size3);								// Immunity level against clinical infection (adult)
	double* ICM = (double*)malloc(size3);								// Immunity level against clinical infection (maternally inherited)
	double* IB = (double*)malloc(size3);								// Immunity level against infection from infectious bite
	double* ID = (double*)malloc(size3);								// Immunity level determining p_det

	double* FOI_lag = (double*)malloc(n_cats * 12500 * sizeof(double));	// Time-lagged force of mosquito -> human infection (to account for incubation period)
	double* FOIv_lag = (double*)malloc(12500 * sizeof(double));			// Time-lagged force of human -> mosquito infection (to account for incubation period)
	double* incv_lag = (double*)malloc(12500 * sizeof(double));			// Time-lagged incubation of malaria in infected mosquito

	vector<double> EIR_benchmarks(n_pts* n_runmax, 0.0);				// Values of entomological inoculation rate at checkpoints
	vector<double> slide_prev_benchmarks(n_pts* n_runmax* na, 0.0);		// Values of slide prevalence at checkpoints
	vector<double> pcr_prev_benchmarks(n_pts* n_runmax* na, 0.0);		// Values of PCR prevalence at checkpoints
	vector<double> clin_inc_benchmarks(n_pts* n_runmax* na, 0.0);		// Values of clinical incidence at checkpoints
	vector<double> M_benchmarks(n_pts * n_runmax, 0.0);					// Values of total relative mosquito population (Sv+Ev+Iv) at checkpoints
	vector<double> M_spor_benchmarks(n_pts * n_runmax, 0.0);			// Values of sporozoite rate in mosquitoes (Iv/(Sv+Ev+Iv)) at checkpoints
	vector<double> EIR_daily_data(tmax_c_i*n_runmax, 0.0);		// Daily entomological inoculation rate values saved for cohort calculations
	vector<double> IB_start_data(n_mv_values* n_cats, 0.0);				// Starting IB values saved for cohort calculations
	vector<double> IC_start_data(n_mv_values* n_cats, 0.0);				// Starting IC values saved for cohort calculations
	vector<double> ID_start_data(n_mv_values* n_cats, 0.0);				// Starting ID values saved for cohort calculations

	vector<double> mv_input = rcpp_to_vector_double(inputs["mv_input"]);
	vector<double> EL_input = rcpp_to_vector_double(inputs["EL_input"]);
	vector<double> LL_input = rcpp_to_vector_double(inputs["LL_input"]);
	vector<double> PL_input = rcpp_to_vector_double(inputs["PL_input"]);
	vector<double> Sv_input = rcpp_to_vector_double(inputs["Sv_input"]);
	vector<double> Ev_input = rcpp_to_vector_double(inputs["Ev_input"]);
	vector<double> Iv_input = rcpp_to_vector_double(inputs["Iv_input"]);
	vector<double> S_input = rcpp_to_vector_double(inputs["S_input"]);
	vector<double> T_input = rcpp_to_vector_double(inputs["T_input"]);
	vector<double> D_input = rcpp_to_vector_double(inputs["D_input"]);
	vector<double> A_input = rcpp_to_vector_double(inputs["A_input"]);
	vector<double> U_input = rcpp_to_vector_double(inputs["U_input"]);
	vector<double> P_input = rcpp_to_vector_double(inputs["P_input"]);
	vector<double> ICA_input = rcpp_to_vector_double(inputs["ICA_input"]);
	vector<double> ICM_input = rcpp_to_vector_double(inputs["ICM_input"]);
	vector<double> IB_input = rcpp_to_vector_double(inputs["IB_input"]);
	vector<double> ID_input = rcpp_to_vector_double(inputs["ID_input"]);

	// Set up conditions for one or more runs-------------------------------------------------------------------------------------------------------------------------------------
			
	if (flag_file == 1) 
	{
		benchmarks_by_age = fopen(file_benchmarks.c_str(), "w");
		fprintf(benchmarks_by_age, "run\tn_mv\tmv0\tn_int\tparam_int\tnpt\tEIR");
		for (i = 0; i < na; i++) { fprintf(benchmarks_by_age, "\tslide_prev%i", i); }
		for (i = 0; i < na; i++) { fprintf(benchmarks_by_age, "\tpcr_prev%i", i); }
		for (i = 0; i < na; i++) { fprintf(benchmarks_by_age, "\tclin_inc%i", i); }
		fclose(benchmarks_by_age);

		data_EIR = fopen(file_EIRd.c_str(), "w");
		fprintf(data_EIR, "n_run");
		for (i = 1; i <= tmax_c_i; i++) { fprintf(data_EIR, "\tEIR%i", i); }
		fclose(data_EIR);
		data_immunity = fopen(file_imm_start.c_str(), "w");
		fprintf(data_immunity, "n_mv");
		for (i = 0; i < n_cats; i++) { fprintf(data_immunity, "\tIB%i", i); }
		for (i = 0; i < n_cats; i++) { fprintf(data_immunity, "\tIC%i", i); }
		for (i = 0; i < n_cats; i++) { fprintf(data_immunity, "\tID%i", i); }
		fclose(data_immunity);
		if (int_v_varied == 0) // Generate header for detailed endpoint data file
		{
			endpoint_data = fopen(file_endpoints.c_str(), "w");
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
			fprintf(endpoint_data, "\n");
			fclose(endpoint_data);
		}
	}

	// Rcout << "\n\nDefault trial parameter values:\nCase treatment rate=" << prop_T_def << "\nATSB killing rate=" << mu_atsb_def << "\nBednet coverage=" << cov_nets_def << "\nIRS coverage=" << cov_irs_def;
	for (n_run = 0; n_run < n_runmax; n_run++) // This for() loop runs all the simulations, jumping to start_run for each one and jumping back to run_complete when it is finished
	{
		n_int = n_run % n_int_values;
		n_mv = (n_run - n_int) / n_int_values;
		mv0 = mv_input[n_mv];
		if (flag_dt_adjust==1) 
		{ 
			if(n_int == 0) {dt = min(0.1,0.1 / sqrt(mv0)); }
		}
		else { dt=0.05; }
		mu_atsb = mu_atsb_def;
		cov_nets = cov_nets_def;
		cov_irs = cov_irs_def;
		cov_set = cov_set_def;
		prop_T = prop_T_def;

		goto start_run;
		run_complete:
		Rcout << "\nRun " << n_run + 1 << " complete. dt=" << dt << "\n";
		R_FlushConsole();
	}

	goto finish;

	//-----------------------------------------------------------------------------Set up initial data for run--------------------------------------------------------------------

	start_run:

	param_int = int_v_varied == 0 ? 0.0 : int_values[n_int];
	prop_T_rev = 1.0 - prop_T;																								// Probability of treatment not being received by clinical cases				
	restart_run:
	Rcout << "\nRun " << n_run + 1 << " Mosquito density=" << mv0 << " Intervention number=" << n_int << " dt=" << dt;
	rN = cov_nets * p_repel_net_entry;
	rNW = cov_nets * p_repel_net_bed;
	dNW = cov_nets * p_kill_net;
	rI = cov_irs * p_repel_irs_entry;
	rIW = cov_irs * p_repel_irs_bed;
	dIW = cov_irs * p_kill_irs1;
	dIF = cov_irs * p_kill_irs2;
	K0 = (mv0 * 2.0 * dl * muv0 * (1.0 + (dp * mup)) * (gammal * (lambda + 1.0))) / term;									// Larval birth rate coefficient, set based on mv0
	r_bite_prevent = phi_bite_house - ((phi_bite_house - phi_bite_bed) * (1.0 - rI) * (1.0 - rIW - dIW) * (1.0 - dIF)) -	// Proportion of bites prevented by bednets and/or interior residual spraying
			(phi_bite_bed * (1.0 - rI) * (1.0 - rIW - dIW) * (1.0 - dIF) * (1.0 - rN) * (1.0 - rNW - dNW));
	mu_protections = 0.3333 * Q0 * (((phi_bite_house - phi_bite_bed) * (1.0 - rI) * (dIW + ((1.0 - rIW - dIW) * dIF))) +		// Additional death rate due to bednets and/or interior residual spraying
			(phi_bite_bed * (1.0 - rI) * (1.0 - rN) * (dIW + ((1.0 - rIW - dIW) * (dNW + ((1.0 - rNW - dNW) * dIF))))));

	av0 = 0.3333 * Q0 * (1.0 - r_bite_prevent);																			// Biting rate on humans / mosquito
	muv1 = muv0 + mu_atsb + mu_protections;																						// Total mosquito death rate

	dt2 = dt * 0.1;		
	ntmax = intdiv(tmax, dt);
	dur_Ei = intdiv(dur_E, dt);
	latgami = intdiv(latgam, dt);
	latmosqi = intdiv(latmosq, dt);
	year_day = fmod(date_start, dy); // Day of the year at the start of the run
	rain = rain_coeff * (rconst + (2.0 * (
			(Re0 * cos(trig_coeff1 * year_day)) + (Im0 * sin(trig_coeff1 * year_day)) +
			(Re1 * cos(trig_coeff2 * year_day)) + (Im1 * sin(trig_coeff2 * year_day)) +
			(Re2 * cos(trig_coeff3 * year_day)) + (Im2 * sin(trig_coeff3 * year_day)) +
			(Re3 * cos(trig_coeff4 * year_day)) + (Im3 * sin(trig_coeff4 * year_day)) +
			(Re4 * cos(trig_coeff5 * year_day)) + (Im4 * sin(trig_coeff5 * year_day))
			)));
	KL = max(KMIN, K0 * rain);

	EL = EL_input[n_mv];
	LL = LL_input[n_mv];
	PL = PL_input[n_mv];
	Sv1 = Sv_input[n_mv];
	Ev1 = Ev_input[n_mv];
	Iv1 = Iv_input[n_mv];
	EIRd = av0 * Iv1 * inv_omega;
	for (i = 0; i < n_cats; i++)
	{
		pos = (i * n_mv_values) + n_mv;
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
	}

	pos = 0;
	FOIv0 = 0.0;
	for (i = 0; i < na; i++) // Run through age categories
	{
		foi_age_cur = foi_age[i];
		for (j = 0; j < num_het; j++) // Run through heterogeneity categories
		{
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
	for (i = 0; i < 12500; i++)
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
	t_mark1 = 1.0;
	t_mark2 = time_values[0];
	EIR_sum = 0.0;
	data_saved = 0;
	flag_int = 0;
	for (nt = 0; nt <= ntmax; nt++)
	{
		tflag = 0;								// Flag indicates whether dt needs to be altered
		t = nt * dt;							// Time since start of simulation
		t_int = t - int_start_time;				// Time since start of intervention use
		nt_E = nt % dur_Ei;
		nt_latgam = nt % latgami;
		nt_latmosq = nt % latmosqi;
		year_day = fmod(t + date_start, dy);	// Day of the year
		mv = Sv1 + Ev1 + Iv1;					// Total number of mosquitoes
		if (flag_int == 0 && t_int >= 0.0)		// Value of varied trial parameters set after start of intervention
		{
			flag_int = 1;
			switch (int_v_varied)
			{
				case 0: {}
					break;
				case 1:
				{
					mu_atsb = param_int;
					Rcout << "\nIntervention begun. ATSB kill rate=" << mu_atsb;
				}
					break;
				case 2:
				{
					cov_nets = param_int;
					rN = cov_nets * p_repel_net_entry;
					rNW = cov_nets * p_repel_net_bed;
					dNW = cov_nets * p_kill_net;
					Rcout << "\nIntervention begun. Bednet coverage=" << cov_nets;
				}
					break;
				case 3:
				{
					cov_irs = param_int;
					rI = cov_irs * p_repel_irs_entry;
					rIW = cov_irs * p_repel_irs_bed;
					dIW = cov_irs * p_kill_irs1;
					dIF = cov_irs * p_kill_irs2;
					Rcout << "\nIntervention begun. IRS coverage=" << cov_irs;
				}
					break;
				case 4:
				{
					cov_set = param_int;
					rSET = cov_set * p_repel_set;
					dSET = cov_set * p_kill_set;
					Rcout << "\nIntervention begun. SET coverage=" << cov_set;
				}
					break;
				default: { Rcout << "\n\tint_v_varied error\n"; }
			}
			R_FlushConsole();
			if(int_v_varied==4)
			{				
				r_bite_prevent = phi_bite_house-((phi_bite_house-phi_bite_bed)*(1.0-rSET-dSET))-(phi_bite_bed*(1.0-rSET-dSET)*(1.0-rNW-dNW));
				mu_protections = 0.3333 * Q0 * (((phi_bite_house-phi_bite_bed)*dSET)+(phi_bite_bed*(dSET+(dNW*(1.0-rSET-dSET)))));
			}
			else
			{
				r_bite_prevent = phi_bite_house - ((phi_bite_house - phi_bite_bed) * (1.0 - rI) * (1.0 - rIW - dIW) * (1.0 - dIF)) -
								(phi_bite_bed * (1.0 - rI) * (1.0 - rIW - dIW) * (1.0 - dIF) * (1.0 - rN) * (1.0 - rNW - dNW));
				mu_protections = 0.3333 * Q0 * (((phi_bite_house - phi_bite_bed) * (1.0 - rI) * (dIW + ((1.0 - rIW - dIW) * dIF))) +
								(phi_bite_bed * (1.0 - rI) * (1.0 - rN) * (dIW + ((1.0 - rIW - dIW) * (dNW + ((1.0 - rNW - dNW) * dIF))))));
			}
			av0 = 0.3333 * Q0 * (1.0 - r_bite_prevent);
			muv1 = muv0 + mu_atsb + mu_protections;
		}

		/*if (t_int >= 0.0 && n_int > 0) // Variable muv1 for special runs - TODO: Integrate into package
		{
			int day = t_int;
			muv1 = mu_atsb_var[day] + muv0 + mu_protections;
		}*/

		rain = rain_coeff * (rconst + (2.0 * (
				(Re0 * cos(trig_coeff1 * year_day)) + (Im0 * sin(trig_coeff1 * year_day)) +
				(Re1 * cos(trig_coeff2 * year_day)) + (Im1 * sin(trig_coeff2 * year_day)) +
				(Re2 * cos(trig_coeff3 * year_day)) + (Im2 * sin(trig_coeff3 * year_day)) +
				(Re3 * cos(trig_coeff4 * year_day)) + (Im3 * sin(trig_coeff4 * year_day)) +
				(Re4 * cos(trig_coeff5 * year_day)) + (Im4 * sin(trig_coeff5 * year_day))
				)));
		KL = max(K0 * rain, KMIN);
		beta = ((0.5 * PL) / dp);		// Rate of spawning of adult mosquitoes
		EIRd = av0 * Iv1 * inv_omega;	// Daily entomological inoculation rate

		for (j = 0; j < num_het; j++) { ICM_init[j] = P_IC_M * (ICA[(age20l * num_het) + j] + age_20_factor * (ICA[(age20u * num_het) + j] - ICA[(age20l * num_het) + j])); }
		pos = 0;
		ij = nt_E * n_cats;
		for (i = 0; i < na; i++) // Run through age categories
		{
			foi_age_cur = foi_age[i];
			for (j = 0; j < num_het; j++) // Run through heterogeneity categories
			{
				EIR_cur = EIRd * foi_age_cur * rel_foi[j];													// Entomological inoculation rate
				b[pos] = bh * (bmin_rev / (1.0 + pow(IB[pos] * inv_IB0, kb)) + bmin);
				FOI[pos] = t > dur_E ? FOI_lag[ij + pos] : FOI0[pos];										// Force of infection (mosquito -> human)
				FOI_lag[ij + pos] = EIR_cur * (IB[pos] > 0.0 ? b[pos] : bh);
				FOI_cur = FOI[pos];
				p_det[pos] = dmin + (dmin_rev / (1.0 + (fd[i] * pow(ID[pos] * inv_ID0, kd))));				// Probability of case detection
				rate_ibaq[pos] = EIR_cur / ((EIR_cur * ub) + 1.0);											// Rate of change of
				rate_detaq[pos] = FOI_cur / ((FOI_cur * ud) + 1.0);											// Rate of change of
				rate_clinaq[pos] = FOI_cur / ((FOI_cur * uc) + 1.0);										// Rate of change of
				phi[pos] = phi0 * ((phi1_rev / (1.0 + pow((ICM[pos] + ICA[pos]) * inv_IC0, kc))) + phi1);	// Probability of clinical disease
				clin_inc[pos] = phi[pos] * FOI_cur * (S[pos] + A[pos] + U[pos]);							// Clinical incidence rate
				cA[pos] = cU + ((cD - cU) * pow(p_det[pos], gamma_inf));
				FOIvij[pos] = foi_age_cur * av0 * ((cT * T[pos]) + (cD * D[pos]) + (cA[pos] * A[pos]) + (cU * U[pos])) * rel_foi[j] * inv_omega;		// Force of infection (human -> mosquito)
				pos++;
			}
		}

		// Key human data - main population
		for (i = na - 1; i >= 0; i--) // Run through age categories (backwards, so that forward ageing does not affect the current time point)
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
			for (j = 0; j < num_het; j++) // Run through heterogeneity categories
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

				// Apply changes
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
				
		// Save daily immunity and EIR data for cohort calculations
		if (t_int >= 0.0)
		{
			EIR_sum += EIRd;
			increment++;
			if (data_saved == 0)
			{
				data_saved = 1;
				ij = n_mv * n_cats;
				for (pos = 0; pos < n_cats; pos++)
				{
					pos2 = ij + pos;
					IB_start_data[pos2] = IB[pos];
					IC_start_data[pos2] = ICA[pos] + ICM[pos];
					ID_start_data[pos2] = ID[pos];
				}
			}
		}

		FOIv = t > latgam ? FOIv_lag[nt_latgam] : FOIv0;
		FOIv_lag[nt_latgam] = arraysum(FOIvij, n_cats);
		incv = t > latmosq ? incv_lag[nt_latmosq] : incv0;
		Surv1 = exp(-muv1 * latmosq);
		beta_larval1 = (eov * muv1 * exp(-muv1 * dgon)) / (1.0 - exp(-muv1 * dgon));
		incv_lag[nt_latmosq] = FOIv * Sv1 * Surv1;
		inv_KL = 1.0 / KL;

		// Larval data calculation; calculated in smaller time increments
		for (np = 0; np < 10; np++)
		{
			EL_cur = EL;
			LL_cur = LL;
			PL_cur = PL;
			dEL = EL * rate_de;
			dLL = LL * rate_dl;
			EL += dt2 * ((beta_larval1 * mv) - (mue * (1.0 + ((EL_cur + LL_cur) * inv_KL)) * EL_cur) - dEL);
			LL += dt2 * (dEL - (mul * (1.0 + ((gammal * (EL_cur + LL_cur)) * inv_KL)) * LL_cur) - dLL);
			PL += dt2 * (dLL - (mup * PL_cur) - (PL_cur / dp));
		}

		// Adult mosquito data calculation
		Sv_cur = Sv1;
		Ev_cur = Ev1;
		Iv_cur = Iv1;
		delSv1 = FOIv * Sv_cur;
		Sv1 += dt * (-(muv1 * Sv_cur) - delSv1 + beta);
		Ev1 += dt * (-(muv1 * Ev_cur) + delSv1 - incv);
		Iv1 += dt * (-(muv1 * Iv_cur) + incv);
		
		// Outputs
		if (t >= t_mark1)
		{
			t_mark1 += 1.0;
			H = arraysum(S, n_cats) + arraysum(T, n_cats) + arraysum(D, n_cats) + arraysum(A, n_cats) + arraysum(U, n_cats) + arraysum(P, n_cats); // Total normalized number of humans; this should always sum to 1.0
			/*S_cur=arraysum(S, n_cats);
			T_cur=arraysum(T, n_cats);
			D_cur=arraysum(D, n_cats);
			A_cur=arraysum(A, n_cats);
			U_cur=arraysum(U, n_cats);
			P_cur=arraysum(P, n_cats);
			H=S_cur+T_cur+D_cur+A_cur+U_cur+P_cur;
			Rcout << "\nt = " << t << "\nS = " << S_cur << " T = " << T_cur << " D = " << D_cur << " A = " << A_cur << " U = " << U_cur << " P = " << P_cur;
			Rcout << "\nH =" << H << " mosq_sum = " << Sv1+Ev1+Iv1 << " larv_sum = " << EL+LL+PL;*/

			// Check that time increment has not been set too large
			if ((H - 1.0) * (H - 1.0) > 0.01) { tflag = 1; }
			if (min(EL,min(LL,PL))<0.0) { tflag = 2; }
			if (min(Sv1,min(Ev1,Iv1))<0.0) { tflag = 3; }
			if (isnan(EIRd)==1){ tflag = 4;}
			if (tflag > 0)
			{
				if (dt <= 0.001) { goto end_run; }
				dt = max(0.001, dt * 0.5);
				Rcout << "\ntflag = " << tflag << " t = " << t << " Adjusting time increment. New dt = " << dt << "\n";
				R_FlushConsole();
				goto restart_run;
				goto end_run;
			}
			if (t_int > 0.0)
			{
				if (interval < tmax_c_i) { EIR_daily_data[(tmax_c_i* n_run) + interval] = EIR_sum / (increment * 1.0); }
				EIR_sum = 0.0;
				increment = 0;
				interval++;
			}
		}

		if (t_int >= t_mark2)
		{
			pos = 0;
			pos2 = (n_pts * n_run) + div;
			M_benchmarks[pos2] = Sv1 + Ev1 + Iv1;
			M_spor_benchmarks[pos2] = Iv1 / M_benchmarks[pos2];
			ij = (n_pts * n_run * na) + (div * na);
			for (i = 0; i < na; i++) // Run through age categories
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
				pos2 = ij + i;
				slide_prev_benchmarks[pos2] = slide_prev_age[i];
				pcr_prev_benchmarks[pos2] = pcr_prev_age[i];
				clin_inc_benchmarks[pos2] = clin_inc_age[i];
			}
			EIR_benchmarks[(n_run * n_pts) + div] = EIRd;
			div++;
			if (div < n_pts) { t_mark2 = time_values[div]; }
			if (flag_file == 1)
			{
				benchmarks_by_age = fopen(file_benchmarks.c_str(), "a");
				fprintf(benchmarks_by_age, "\n%i\t%i\t%.3e\t%i\t%.3e\t%i\t%.3e", n_run, n_mv, mv0, n_int, param_int, div, EIRd);
				for (i = 0; i < na; i++) { fprintf(benchmarks_by_age, "\t%.3e", slide_prev_benchmarks[ij + i]); }
				for (i = 0; i < na; i++) { fprintf(benchmarks_by_age, "\t%.3e", pcr_prev_benchmarks[ij + i]); }
				for (i = 0; i < na; i++) { fprintf(benchmarks_by_age, "\t%.3e", clin_inc_benchmarks[ij + i]); }
				fclose(benchmarks_by_age);
			}
			//Rcout << "\nt = " << t << " H = " << H << " EIR = " << EIRd;
			R_FlushConsole();
		}
	}

	end_run:
		
	if (flag_file == 1)
	{
		data_EIR = fopen(file_EIRd.c_str(), "a");
		fprintf(data_EIR, "\n%i", n_run);
		for (i = 0; i < tmax_c_i; i++) { fprintf(data_EIR, "\t%.3e", EIR_daily_data[(tmax_c_i * n_run) + i]); }
		fclose(data_EIR);
		if (n_int == 0)
		{
			data_immunity = fopen(file_imm_start.c_str(), "a");
			fprintf(data_immunity, "\n%i", n_mv);
			ij = n_mv * n_cats;
			for (i = 0; i < n_cats; i++) { fprintf(data_immunity, "\t%.3e", IB_start_data[ij + i]); }
			for (i = 0; i < n_cats; i++) { fprintf(data_immunity, "\t%.3e", IC_start_data[ij + i]); }
			for (i = 0; i < n_cats; i++) { fprintf(data_immunity, "\t%.3e", ID_start_data[ij + i]); }
			fclose(data_immunity);
		}
		if (int_v_varied == 0)
		{
			// Output detailed data to use as input for future computation
			endpoint_data = fopen(file_endpoints.c_str(), "a");
			fprintf(endpoint_data, "%i\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e", n_run + 1, mv0, EL, LL, PL, Sv1, Ev1, Iv1);
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
			fprintf(endpoint_data, "\n");
			fclose(endpoint_data);
		}
	}

	goto run_complete;

	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	finish:
	//Rcout << "\nComputations complete 1.\n";
	//R_FlushConsole();

	if (benchmarks_by_age != NULL) { fclose(benchmarks_by_age); }
	if (endpoint_data != NULL) { fclose(endpoint_data); }
	if (data_EIR != NULL) { fclose(data_EIR); }
	if (data_immunity != NULL) { fclose(data_immunity); }
	Rcout << "\nComputations complete.\n";
	R_FlushConsole();

	// Return list
	//return List::create(Named("dummy") = 1.0);
	return List::create(Named("EIR_benchmarks") = EIR_benchmarks, Named("slide_prev_benchmarks")= slide_prev_benchmarks, 
						Named("pcr_prev_benchmarks")= pcr_prev_benchmarks, Named("clin_inc_benchmarks")= clin_inc_benchmarks, 
						Named("M_benchmarks") = M_benchmarks, Named("M_spor_benchmarks") = M_spor_benchmarks,
						Named("EIR_daily_data")= EIR_daily_data, Named("IB_start_data") = IB_start_data, 
						Named("IC_start_data") = IC_start_data, Named("ID_start_data") = ID_start_data);
}