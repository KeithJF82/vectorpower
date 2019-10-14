#include <Rcpp.h>
#include "main.h"
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

// [[Rcpp::export]]
Rcpp::List rcpp_cohort(List params, List cohort_params, List mainpop_data, List cluster_data)
{
	int na_c0, na_c1, na_c2, n_c, mvn, int_num, data_num, i, j, n, nt, div, pos, ntmax, tflag, interval, infected;
	double dt, t, EIRd, EIRt, t_mark1, t_mark2, p_multiplier, prob, cprobmax, p_inf_bite, p_inf_from_bite, p_clin_inf, IC_cur, IB_cur, ID_cur, rate_db_t, rate_dc_t, rate_dd_t;
	FILE* output1;
	FILE* output2;

	output1 = NULL;
	output2 = NULL;

	//Constants (TODO: Make global)----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	double dy = 365.0;			// Days in a year
	double tinterval1 = 1.0;	// Interval between calculations of current prevalence / incidence values
	int na = 145;				// Number of age categories in main population
	int num_het = 9;			// Number of heterogeneity categories
	double het_x[] = { -4.5127459, -3.2054291, -2.076848, -1.0232557, 0.0, 1.0232557, 2.076848, 3.2054291, 4.5127459 }; // Biting heterogeneity
	double het_wt[] = { 0.00002235, 0.00278914, 0.04991641, 0.2440975, 0.40634921, 0.2440975, 0.04991641, 0.00278914, 0.00002235 }; // Fractions in each heterogeneity category

	//Load input parameter data from R------------------------------------------------------------------------------------------------------------------------------------------

	Rcout << "\nLoading input parameter data from R\n";
	R_FlushConsole();
	string output_filename1 = rcpp_to_string(cohort_params["file_summary"]);	//Individual data (fraction in each category) every tinterval2 days over intervention period for each iteration (single run) or average cohort prevalence/indidence/PCR test positivity rate every tinterval2 days after intervention start (multi-run)
	string output_filename2 = rcpp_to_string(cohort_params["file_frequency"]);	//Positive PCR test frequency data for each iteration (single run) or average across all iterations (multi-run)
	double tinterval2 = rcpp_to_double(cohort_params["time_interval"]);			//Interval between calculations of regular test probabilities and averaged prevalence/incidence values
	int n_divs = rcpp_to_int(cohort_params["n_divs"]);							//Number of tinterval2-length divisions to run after intervention starts
	double tmax = (n_divs - 1) * tinterval2;									//Duration of intervention trial 
	int n_clusters = rcpp_to_int(cohort_params["n_clusters"]);					//Number of lines in cluster data input file (number of clusters)
	double prop_T_c = rcpp_to_int(cohort_params["prop_T_c"]);					//Proportion of clinical cases successfully treated in cohort
	int n_patients = rcpp_to_int(cohort_params["n_patients"]);					//Number of patients in cohort
	double age_c0 = rcpp_to_int(cohort_params["age_c0"]);						//Minimum age in cohort
	double age_c1 = rcpp_to_int(cohort_params["age_c1"]);						//Maximum age in cohort
	double age_c2 = age_c1 + (tmax / dy);										//Maximum age in cohort taking into account ageing during trial

	//Set up constant parameters------------------------------------------------------------------------------------------------------------------------------------------------

	Rcout << "\nInitializing constants.\n";
	R_FlushConsole();
	double inv_dy = 1.0 / dy;
	double dv_p1 = 1.0 / n_patients;
	double dv_p2 = dv_p1 * (dy / tinterval2);
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
	double* age_width = (double*)malloc(na * sizeof(double));
	double* age_rate = (double*)malloc(na * sizeof(double));
	double* rem_rate = (double*)malloc(na * sizeof(double));
	double* age = (double*)malloc(na * sizeof(double));
	double* agey0 = (double*)malloc(na * sizeof(double));
	double* agey1 = (double*)malloc(na * sizeof(double));
	double eta = 1.0 / (21.0 * dy);		//death rate, for exponential population distribution
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
	double* den = (double*)malloc(na * sizeof(double));
	double* foi_age = (double*)malloc(na * sizeof(double));
	double densum = 0.0;
	for (i = 0; i < na; i++)
	{
		den[i] = i == 0 ? (1.0 / (1.0 + (age_rate[0] / eta))) : (age_rate[i - 1] * den[i - 1] / (rem_rate[i]));
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

	double* cprob = (double*)malloc(n_cats_c * sizeof(double));				//Cumulative probability distribution used to determine age/heterogenity category to place randomly generated patients
	int* patients_status = (int*)malloc(6 * sizeof(int));						//Tally of number of patients in each status category (S,T,D,A,U,P)
	int* pcr_test_results = (int*)malloc(n_divs * n_patients * sizeof(int));		//PCR test results for individual patients at tinterval2 checkpoints
	double* pcr_distribution = (double*)malloc((n_divs + 1) * sizeof(double)); //Frequency distribution of number of positive PCR test results over entire trial
	double* slide_prev_values = (double*)malloc(n_divs * sizeof(double));	//Values of slide prevalence (cohort) at tinterval2 checkpoints
	double* clin_inc_values = (double*)malloc(n_divs * sizeof(double));	//Values of clinical incidence (cohort) at tinterval2 checkpoints
	double* pcr_prev_values = (double*)malloc(n_divs * sizeof(double));	//Values of proportion of patients in cohort giving positive PCR test results at tinterval2 checkpoints

	//Load input data------------------------------------------------------------------------------------------------------------------------------------------

	vector<vector<double>> EIR_input = rcpp_to_matrix_double(mainpop_data["EIR_daily_data"]);
	vector<vector<double>> IB_input = rcpp_to_matrix_double(mainpop_data["IB_start_data"]);
	vector<vector<double>> IC_input = rcpp_to_matrix_double(mainpop_data["IC_start_data"]);
	vector<vector<double>> ID_input = rcpp_to_matrix_double(mainpop_data["ID_start_data"]);
	vector<int> run_mvn = rcpp_to_vector_int(cluster_data["n_B"]);
	vector<int> run_int = rcpp_to_vector_int(cluster_data["n_I"]);
	vector<int> run_data_num = rcpp_to_vector_int(cluster_data["n_run"]);

	//Set up conditions for one or more runs-------------------------------------------------------------------------------------------------------------------------------------

	Rcout << "\nSetting up conditions.\n";
	R_FlushConsole();
	if (n_clusters > 1)
	{
		output1 = fopen(output_filename1.c_str(), "w");//fopen_s(&output1, output_filename1.c_str(), "w"); //
		fprintf(output1, "n_patients\t%i\nCluster\tmvn\tint_num", n_patients);
		for (div = 1; div <= n_divs; div++) { fprintf(output1, "\tslide_prev%i", div); }
		for (div = 1; div <= n_divs; div++) { fprintf(output1, "\tpcr_p%i", div); }
		for (div = 1; div <= n_divs; div++) { fprintf(output1, "\tclin_incinc_%i", div); }
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
		Rcout << "\nProcessing cluster " << n_c << ":\tmvn=" << mvn << "\tint_num=" << int_num << "\tdata_num=" << data_num << "\tEIR[0]=" << EIR_input[0][data_num];
		goto start_run;
	run_complete:
		Rcout << "\nCluster " << n_c << " complete.\n";
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
	for (i = 0; i < 6; i++) { patients_status[i] = 0; }

	dt = 0.1; //Time increment for cohort data. Like that used for the main population, this may be adjusted automatically to prevent errors
	tflag = 0;

	for (n = 0; n < n_patients; n++) //Randomly determine patient's age and heterogeneity category using cprob
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
	}

restart:
	for (div = 0; div < n_divs; div++)
	{
		clin_inc_values[div] = 0.0;
		slide_prev_values[div] = 0.0;
		pcr_prev_values[div] = 0.0;
		pcr_distribution[div] = 0.0;
	}
	pcr_distribution[n_divs] = 0.0;

	if (n_clusters == 1)
	{
		output1 = fopen(output_filename1.c_str(), "w");//fopen_s(&output1, output_filename1.c_str(), "w"); //
		output2 = fopen(output_filename2.c_str(), "w");//fopen_s(&output2, output_filename2.c_str(), "w"); //
		fprintf(output1, "\n\nCluster:\t%i\tmvn:\t%i\tint_num\t%i\nt\tS\tT\tD\tA\tU\tP\tEIR\tprev\tclin_inc\tpcr_prev", n_c, mvn, int_num);
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
		pos = (patient[n].na * num_het) + patient[n].num_het;
		patient[n].IB = IB_input[pos][mvn];
		patient[n].IC = IC_input[pos][mvn];
		patient[n].ID = ID_input[pos][mvn];
		patient[n].infected = 0;
		patient[n].delay = 0.0;
	}

	t_mark1 = tinterval1;
	t_mark2 = 0.0;// tinterval2;
	div = 0;
	interval = 0;
	rate_db_t = rate_db * dt;
	rate_dc_t = rate_dc * dt;
	rate_dd_t = rate_dd * dt;
	for (nt = 0; nt <= ntmax; nt++)
	{
		t = nt * dt;
		EIRd = EIR_input[interval][data_num];
		EIRt = EIRd * dt;
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
				if (dt <= 0.0025) { goto end; }
				dt = max(0.0025, dt * 0.5);
				goto restart;
			}
			if (runif1() <= p_inf_bite)
			{
				patient[n].IB += IB_boost;
				p_inf_from_bite = bh * (IB_cur > 0.0 ? ((bmin_rev / (1.0 + pow(IB_cur * inv_IB0, kb))) + bmin) : 1.0); //Probability of an infection from an infectious bite
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
					p_clin_inf = phi0 * ((phi1_rev / (1.0 + pow(IC_cur * inv_IC0, kc))) + phi1); //Probability of a clinical infection
					if (runif1() <= p_clin_inf)
					{/*clinical infection*/
						clin_inc_values[div] += dv_p2;
						if (runif1() <= prop_T_c) {/*to T*/ patient[n].status = 1; }
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

		//Output patient data and benchmarks
		if (t >= t_mark2)
		{
			for (i = 0; i < 6; i++) { patients_status[i] = 0; }
			pos = div * n_patients;
			for (n = 0; n < n_patients; n++)
			{
				patients_status[patient[n].status]++;
				switch (patient[n].status)
				{
				case 1: //T
				{
					pcr_test_results[pos] = 1;
					pcr_prev_values[div] += dv_p1;
					slide_prev_values[div] += dv_p1;
				}
				break;
				case 2: //D
				{
					pcr_test_results[pos] = 1;
					pcr_prev_values[div] += dv_p1;
					slide_prev_values[div] += dv_p1;
				}
				break;
				case 3: //A
				{
					pcr_test_results[pos] = 1;
					pcr_prev_values[div] += dv_p1;
					if (runif1() <= dmin + (dmin_rev / (1.0 + (fd[patient[n].na] * pow(patient[n].ID * inv_ID0, kd))))) { slide_prev_values[div] += dv_p1; }
				}
				break;
				case 4: //U
				{
					pcr_test_results[pos] = 1;
					pcr_prev_values[div] += dv_p1;
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
				fprintf(output1, "\t%.3e\t%.3e\t%.3e\t%.3e", EIRd, slide_prev_values[div], clin_inc_values[div], pcr_prev_values[div]);
				Rcout << "\n\tt = " << t << " EIR = " << EIRd << " prev = " << slide_prev_values[div] << " inc = " << clin_inc_values[div] << " pcr = " << pcr_prev_values[div];
				fclose(output1);
			}
			//Rcout << "\nt = " << t << "\tEIR = " << EIRd << "\tprev = " << slide_prev_values[div];
			t_mark2 += tinterval2;
			div++;
			for (n = 0; n < n_patients; n++)
			{
				j = 0;
				for (i = 0; i < n_divs; i++) { if (pcr_test_results[(i * n_patients) + n] == 1) { j++; } }
				pcr_distribution[j] += dv_p1;
			}
		}

	}

	output2 = fopen(output_filename2.c_str(), "a");//fopen_s(&output2, output_filename2.c_str(), "a"); //
	if (n_clusters > 1)
	{
		output1 = fopen(output_filename1.c_str(), "a");//fopen_s(&output1, output_filename1.c_str(), "a"); //
		fprintf(output1, "\n%i\t%i\t%i", n_c, mvn, int_num);
		for (div = 0; div < n_divs; div++) { fprintf(output1, "\t%.3e", slide_prev_values[div]); }
		for (div = 0; div < n_divs; div++) { fprintf(output1, "\t%.3e", pcr_prev_values[div]); }
		for (div = 0; div < n_divs; div++) { fprintf(output1, "\t%.3e", clin_inc_values[div]); }
		fclose(output1);
		fprintf(output2, "\n%i\t%i\t%i", n_c, mvn, int_num);
		for (div = 0; div <= n_divs; div++) { fprintf(output2, "\t%.3e", pcr_distribution[div] / n_divs); }
	}
	else
	{
		for (n = 0; n < n_patients; n++)
		{
			fprintf(output2, "\n%i", n);
			for (div = 0; div < n_divs; div++) { fprintf(output2, "\t%i", pcr_test_results[(div * n_patients) + n]); }
		}
	}
	fclose(output2);

end:

	goto run_complete;

	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

finish:

	vector<double> dummy(1, 0.0);

	if (output1 != NULL) { fclose(output1); }
	if (output2 != NULL) { fclose(output2); }
	Rcout << "\nProgram complete\n";

	// Return list
	return List::create(Named("dummy") = dummy);
}
