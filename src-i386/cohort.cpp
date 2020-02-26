#include <Rcpp.h>
#include "main.h"
using namespace Rcpp;
using namespace std;

struct patients//Structure containing individual patient data for cohort
{
	int status;//Infection group - 0=S, 1=T, 2=D, 3=A, 4=U, 5=P
	int na, num_het;//Age and heterogeneity groups
	double IB, IC, ID;//Immunities
	int infected;//Flag indicating patient in early stage of infection
	double delay;//Infection delay time
}
patient[1000];

// [[Rcpp::export]]
Rcpp::List rcpp_cohort(List params, List trial_params, List cluster_data)
{
	int na_c0, na_c1, na_c2, n_c, mvn, int_num, data_num, i, j, n, nt, div, pos, ntmax, tflag, interval, infected;
	double dt, t, EIRd, EIRt, t_mark1, t_mark2, p_multiplier, prob, cprobmax, p_inf_bite, p_inf_from_bite, p_clin_inf, IC_cur, IB_cur, ID_cur, rate_db_t, rate_dc_t, rate_dd_t, random_value;

	//Constants (TODO: Make global)----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	double dy = 365.0;											// Days in a year
	double inv_dy = 1.0 / dy;
	double tinterval1 = 1.0;									// Interval between calculations of current prevalence / incidence values

	//Load input parameter data from R------------------------------------------------------------------------------------------------------------------------------------------

	int flag_output = rcpp_to_int(trial_params["flag_output"]);
	Rcout << "\nBeginning cluster calculations\n";
	R_FlushConsole();
	vector<double> time_values = rcpp_to_vector_double(trial_params["time_values"]);	// Vector of time benchmark points
	int n_divs = time_values.size();
	vector<double> tinterval2(n_divs, 0.0);
	tinterval2[0] = time_values[0] = 0.0 ? 1.0 : time_values[0];
	for (div = 1; div < n_divs; div++) { tinterval2[div] = time_values[div] - time_values[div-1]; }
	int tmax_i = rcpp_to_int(trial_params["tmax_i"]);
	double tmax = tmax_i * 1.0;
	int n_clusters = rcpp_to_int(trial_params["n_clusters"]);					// Number of lines in cluster data input file (number of clusters)
	double prop_T_c = rcpp_to_double(trial_params["prop_T_c"]);					// Proportion of clinical cases successfully treated in cohort
	int n_patients = rcpp_to_int(trial_params["n_patients"]);					// Number of patients in cohort
	double age_start = rcpp_to_double(trial_params["age_start"]);					// Minimum age in cohort
	double age_end = rcpp_to_double(trial_params["age_end"]);						// Maximum age in cohort
	double age_c2 = age_end + (tmax * inv_dy);									// Maximum age in cohort taking into account ageing during trial

	//Set up constant parameters------------------------------------------------------------------------------------------------------------------------------------------------

	// Heterogeneity parameters
	int num_het = rcpp_to_int(params["num_het"]);						// Number of age categories in main population
	vector<double> het_x = rcpp_to_vector_double(params["het_x"]);		// Biting heterogeneity
	vector<double> het_wt = rcpp_to_vector_double(params["het_wt"]);	// Fractions in each heterogeneity category

	// Age parameters
	int na = rcpp_to_int(params["na"]);										// Number of age categories in main population
	vector<double> age_width = rcpp_to_vector_double(params["age_width"]);	// Widths of age categories in days
	vector<double> age_rate = rcpp_to_vector_double(params["age_rate"]);	// Rates of transition between age categories
	vector<double> rem_rate = rcpp_to_vector_double(params["rem_rate"]);	// Total loss rate from each age category (ageing + death)
	vector<double> age = rcpp_to_vector_double(params["age"]);				// Age values corresponding to starts of age categories (days)
	vector<double> den = rcpp_to_vector_double(params["den"]);

	int n_cats = na * num_het;									//Total number of age/heterogeneity categories in main population 
	double dv_p1 = 1.0 / n_patients;
	double rho = rcpp_to_double(params["rho"]);					// Age-dependent biting parameter
	double a0 = rcpp_to_double(params["a0"]);					// Age-dependent biting parameter
	double sigma2 = rcpp_to_double(params["sigma2"]);			// Variance of log heterogeneity in biting
	double dur_E = rcpp_to_double(params["dur_E"]);				// Latent period
	double dur_T = rcpp_to_double(params["dur_T"]);				// Average time to recover from disease and parasitaemia when treated
	double dur_D = rcpp_to_double(params["dur_D"]);				// Average time to move from clinical disease to asymptomatic when not successfully treated
	double dur_A = rcpp_to_double(params["dur_A"]);				// Average time to move from asymptomatic patent infection to sub-patent infection
	double dur_U = rcpp_to_double(params["dur_U"]);				// Average time to recover from sub-patent infection
	double dur_P = rcpp_to_double(params["dur_P"]);				// Average time to leave protected state

	// Immunity parameters - infection
	double bh = rcpp_to_double(params["bh"]);					// Maximum probability due to no immunity
	double bmin = rcpp_to_double(params["bmin"]);				// Maximum relative reduction due to immunity
	double db = rcpp_to_double(params["db"]);					// Inverse of decay rate
	double kb = rcpp_to_double(params["kb"]);					// Shape parameter
	double ub = rcpp_to_double(params["ub"]);					// Duration in which immunity is not boosted
	double inv_IB0 = 1.0 / rcpp_to_double(params["IB0"]);		// Inverse of scale parameter
	double bmin_rev = 1.0 - bmin;
	double IB_boost = 1.0 / (ub + 1.0);

	// Immunity parameters - detection
	double dmin = rcpp_to_double(params["dmin"]);						// Minimum probability due to maximum immunity
	double dd = rcpp_to_double(params["dd"]);							// Inverse of decay rate
	double kd = rcpp_to_double(params["kd"]);							// Shape parameter
	double ud = rcpp_to_double(params["ud"]);							// Duration in which immunity is not boosted
	double ad0 = rcpp_to_double(params["ad0"]);							// Scale parameter relating age to immunity
	double fd0 = rcpp_to_double(params["fd0"]);							// Time scale at which immunity changes with age
	double gammad = rcpp_to_double(params["gammad"]);					// Shape parameter relating age to immunity
	double inv_ID0 = 1.0 / rcpp_to_double(params["ID0"]);				// Scale parameter
	double dmin_rev = 1.0 - dmin;
	double ID_boost = 1.0 / (ud + 1.0);

	// Immunity parameters - clinical disease
	double phi0 = rcpp_to_double(params["phi0"]);						// Maximum probability due to no immunity
	double phi1 = rcpp_to_double(params["phi1"]);						// Maximum relative reduction due to no immunity
	double dc = rcpp_to_double(params["dc"]);							// Inverse of decay rate
	double kc = rcpp_to_double(params["kc"]);							// Shape parameter
	double uc = rcpp_to_double(params["uc"]);							// Duration in which immunity is not boosted
	double inv_IC0 = 1.0 / rcpp_to_double(params["IC0"]);				// Scale parameter
	double phi1_rev = 1.0 - phi1;
	double IC_boost = 1.0 / (uc + 1.0);

	//Constant derived values--------------------------------------------------------------------------------------------------------------------------

	R_FlushConsole();

	double rT = 1.0 / dur_T;			//Rate of recovery from disease and parasitaemia when treated
	double rD = 1.0 / dur_D;			//Rate of moving from clinical disease to asymptomatic when not successfully treated
	double rA0 = 1.0 / dur_A;			//Rate of moving from asymptomatic patent infection to sub-patent infection
	double rU = 1.0 / dur_U;			//Rate of recovery from sub - patent infection
	double rP = 1.0 / (dur_P - dur_T);	//Rate of leaving protected state
	double rate_dc = 1.0 / dc;
	double rate_db = 1.0 / db;
	double rate_dd = 1.0 / dd;
	int agefinding0 = 1;				//Flags indicating which cohort age parameter is being sought while running through ages
	int agefinding1 = 0;
	int agefinding2 = 0;

	na_c0 = 0;					//Minimum age category number in main population corresponding to first age category in cohort; determined from age_start
	na_c1 = na - 1;				//Maximum age category number in main population corresponding to last age category in cohort at start; determined from age_end
	na_c2 = na - 1;				//Maximum age category number in main population corresponding to last age category in cohort at end; determined from age_c2
	for (i = 0; i < na; i++)
	{
		if (agefinding0 == 1 && age[i] >= age_start * dy)
		{
			na_c0 = i;
			agefinding0 = 0;
			agefinding1 = 1;
		}
		if (agefinding1 == 1 && age[i] >= age_end * dy)
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
	}

	int na_c = na_c2 - na_c0 + 1;			//Total number of age categories within cohort
	int n_cats_c = na_c * num_het;			//Total number of age and heterogeneity categories in cohort

	//Calculate the proportion in each age group using discrete time, to match the equilibrium from the differential eqns
	double* foi_age = (double*)malloc(na * sizeof(double));
	double densum = 0.0;
	for (i = 0; i < na; i++)
	{
		densum += den[i];
		foi_age[i] = 1.0 - (rho * exp(-age[i] / a0));
	}
	double inv_densum = 1.0 / densum;

	//Heterogeneity
	double* rel_foi = (double*)malloc(num_het * sizeof(double));
	for (j = 0; j < num_het; j++) { rel_foi[j] = num_het == 0 ? 1.0 : exp(-(sigma2 * 0.5) + (pow(sigma2, 0.5) * het_x[j])); }
	double* fd = (double*)malloc(na * sizeof(double));
	for (i = 0; i < na; i++) { fd[i] = 1.0 - ((1.0 - fd0) / (1.0 + (pow(age[i] / ad0, gammad)))); }

	//Set up additional arrays-------------------------------------------------------------------------------------------------------------------------------------

	double* cprob = (double*)malloc(n_cats_c * sizeof(double));						//Cumulative probability distribution used to determine age/heterogenity category to place randomly generated patients
	double* clin_inc_values = (double*)malloc(n_divs * sizeof(double));             //Values of clinical incidence (cohort) between checkpoints
	vector<int> patients_status_outputs(n_clusters * n_patients * n_divs, 0);		//All patient statuses across all clusters at each time point, for output to R
	vector<double> p_det_outputs(n_clusters * n_patients * n_divs, 0.0);			//Probability of asymptomatic infection being detected, where relevant, for all patients across all clusters at each time point, for output to R
	vector<double> clin_inc_outputs(n_clusters * n_divs, 0.0);						//Clinical incidence per person per day in each cluster at each time point (averaged over time up to that point from previous point), for output to R
	vector<double> patients_age_outputs(n_clusters * n_patients, 0.0);				//Patient ages (at start of trial period) across all clusters
	vector<double> patients_het_outputs(n_clusters * n_patients, 0.0);				//Patient heterogeneity group numbers across all clusters

	//Load input data------------------------------------------------------------------------------------------------------------------------------------------

	vector<double> EIR_input = rcpp_to_vector_double(trial_params["EIR_daily_data"]);
	vector<double> IB_input = rcpp_to_vector_double(trial_params["IB_start_data"]);
	vector<double> IC_input = rcpp_to_vector_double(trial_params["IC_start_data"]);
	vector<double> ID_input = rcpp_to_vector_double(trial_params["ID_start_data"]);
	vector<int> run_mvn = rcpp_to_vector_int(cluster_data["n_B"]);
	vector<int> run_int = rcpp_to_vector_int(cluster_data["n_I"]);
	vector<int> run_data_num = rcpp_to_vector_int(cluster_data["n_run"]);

	//Set up conditions for one or more runs-------------------------------------------------------------------------------------------------------------------------------------
	
	for (n_c = 0; n_c < n_clusters; n_c++) //This for() loop runs all the simulations, jumping to start_run for each one and jumping back to run_complete when it is finished
	{
		mvn = run_mvn[n_c];
		int_num = run_int[n_c];
		data_num = run_data_num[n_c];
		if (flag_output == 1) { Rcout << "Processing cluster " << n_c << ":\tmvn=" << mvn << "\tint_num=" << int_num << "\tdata_num=" << data_num; }
		goto start_run;
	run_complete:
		if (flag_output == 1) { Rcout << "\n"; }
	}

	goto finish;

	//-------------------------------------------------------------------------------------Compute over time--------------------------------------------------------------------------------------------------------

start_run:

	//Set up cumulative probability distribution from which to draw randomly generated patients

	for (pos = 0; pos < n_cats_c; pos++) { cprob[pos] = 0.0; }
	for (i = na_c0; i <= na_c1; i++)
	{
		for (j = 0; j < num_het; j++)
		{
			pos = ((i - na_c0) * num_het) + j;
			cprob[pos] = (den[i] * het_wt[j]) * inv_densum;
			if (pos > 0) { cprob[pos] += cprob[pos - 1]; }
		}
	}

	cprobmax = cprob[((na_c1 - na_c0) * num_het) + num_het - 1];
	for (pos = 0; pos < n_cats_c; pos++) { if (cprob[pos] > 0.0) { cprob[pos] /= cprobmax; } }

	dt = 0.1;//Time increment for cohort data. Like that used for the main population, this may be adjusted automatically to prevent errors
	tflag = 0;

	for (n = 0; n < n_patients; n++)//Randomly determine patient's age and heterogeneity category using cprob
	{
		prob = runif1();
		for (i = na_c0; i <= na_c1; i++)
		{
			for (j = 0; j < num_het; j++)
			{
				if (prob <= cprob[((i - na_c0) * num_het) + j]) { goto found; }
			}
		}
	found:
		patient[n].na = i;
		patient[n].num_het = j;
		pos = ((n_c)*n_patients) + n;
		patients_age_outputs[pos] = age[i] * inv_dy;
		patients_het_outputs[pos] = j;
	}

restart:
	for (div = 0; div < n_divs; div++)
	{
		clin_inc_values[div] = 0.0;
	}

	p_multiplier = 1.0 / dt;
	ntmax = intdiv(tmax - 1.0, dt);
	R_FlushConsole();

	for (n = 0; n < n_patients; n++)
	{
		patient[n].status = 0;	//All patients start in group S, representing pre-intervention treatment
		pos = ((mvn - 1) * n_cats) + (patient[n].na * num_het) + patient[n].num_het;
		patient[n].IB = IB_input[pos];
		patient[n].IC = IC_input[pos];
		patient[n].ID = ID_input[pos];
		patient[n].infected = 0;
		patient[n].delay = 0.0;
	}

	t_mark1 = tinterval1;
	t_mark2 = time_values[0];
	div = 0;
	interval = 0;
	rate_db_t = rate_db * dt;
	rate_dc_t = rate_dc * dt;
	rate_dd_t = rate_dd * dt;
	for (nt = 0; nt <= ntmax; nt++)
	{
		t = nt * dt;

		//Output patient data and benchmarks
		if (t >= t_mark2)
		{
			pos = (n_c * n_patients * n_divs) + div;
			clin_inc_outputs[(n_c * n_divs) + div] = clin_inc_values[div];
			for (n = 0; n < n_patients; n++)
			{
				patients_status_outputs[pos] = patient[n].status;
				if (patient[n].status == 3) { p_det_outputs[pos] = dmin + (dmin_rev / (1.0 + (fd[patient[n].na] * pow(patient[n].ID * inv_ID0, kd)))); }
				pos += n_divs;
			}

			div++;
			if (div < n_divs) { t_mark2 = time_values[div]; }
		}

		EIRd = EIR_input[(data_num * tmax_i) + interval];
		EIRt = EIRd * dt;
		for (n = 0; n < n_patients; n++)
		{
			infected = 0;
			i = patient[n].na;
			j = patient[n].num_het;
			IB_cur = patient[n].IB;
			IC_cur = patient[n].IC;
			ID_cur = patient[n].ID;
			p_inf_bite = EIRt * rel_foi[j] * foi_age[i];//Probability of an infectious bite

			//Check that dt has not been set too high
			if (p_inf_bite > 0.9) { tflag = 1; }
			if (tflag > 0)
			{
				tflag = 0;
				if (dt <= 0.0025) { goto end; }
				dt = max(0.0025, dt * 0.5);
				goto restart;
			}
			if (runif1() <= p_inf_bite)
			{
				patient[n].IB += IB_boost;
				p_inf_from_bite = bh * (IB_cur > 0.0 ? ((bmin_rev / (1.0 + pow(IB_cur * inv_IB0, kb))) + bmin) : 1.0);	//Probability of an infection from an infectious bite
				if (runif1() <= p_inf_from_bite)
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
					p_clin_inf = phi0 * ((phi1_rev / (1.0 + pow(IC_cur * inv_IC0, kc))) + phi1);//Probability of a clinical infection
					if (runif1() <= p_clin_inf)
					{/*clinical infection*/
						clin_inc_values[div] += dv_p1 / tinterval2[div];
						random_value = runif1();
						if (random_value <= prop_T_c) {/*to T*/ patient[n].status = 1; }
						else {/*to D*/  patient[n].status = 2; }
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
					if (runif1() * p_multiplier <= rT) {/*to P*/ patient[n].status = 5; }
				}
				break;
				case 2: //D
				{
					if (runif1() * p_multiplier <= rD) {/*to A*/ patient[n].status = 3; }
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
						if (runif1() * p_multiplier <= rA0) {/*to U*/ patient[n].status = 4; }
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
						if (runif1() * p_multiplier <= rU) {/*to S*/ patient[n].status = 0; }
					}
				}
				break;
				case 5: //P
				{
					if (runif1() * p_multiplier <= rP) {/*to S*/ patient[n].status = 0; }
				}
				break;
				default: { Rcout << "\nError in switch 1!\n"; }
				}
			}
			//Natural decrease in immunity
			patient[n].IB -= IB_cur * rate_db_t;
			patient[n].IC -= IC_cur * rate_dc_t;
			patient[n].ID -= ID_cur * rate_dd_t;
		}

		if (t >= t_mark1)
		{
			interval++;
			t_mark1 += tinterval1;
		}
	}
	
end:

	goto run_complete;

	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

finish:

	Rcout << "\nCluster calculations complete\n";

	// Return list
	return List::create(Named("patients_age_outputs") = patients_age_outputs, Named("patients_het_outputs") = patients_het_outputs, Named("patients_status_outputs") = patients_status_outputs, 
						Named("p_det_outputs") = p_det_outputs, Named("clin_inc_outputs") = clin_inc_outputs);
}
