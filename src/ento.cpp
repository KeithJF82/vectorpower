// CPP file for simulation of malaria intervention trial - entomological run only
#include <Rcpp.h>
#include "main.h"
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
Rcpp::List rcpp_ento(List params, List inputs, List trial_params)
{
	int i, n_run, n_mv, n_int, nt, ntmax, div, /*latmosqi, */tflag, flag_int, np, pos2;
	double mv0, param_int, t, t_int, mu_protections, /*r_bite_prevent, */muv1, K0, year_day, rain, KL, beta, mv, t_mark1, t_mark2, dt2;
	double mu_atsb, cov_nets, cov_irs, cov_set, rN, rNW, dNW, rI, rIW, dIW, dIF, EL, LL, PL, inv_KL;
	double EL_cur, LL_cur, PL_cur, Sv_cur, Ev_cur, Iv_cur, dEL, dLL, delSv, incv;
	FILE* endpoint_data = NULL;
	double dSET=0.0;
	double rSET=0.0;

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

	// Entomological parameters
	double mu_atsb_def = rcpp_to_double(params["mu_atsb_def"]);					// Attractive targeted sugar bait (ATSB) killing rate associated with data
	double cov_nets_def = rcpp_to_double(params["cov_nets_def"]);				// Proportion of people protected by bednets associated with data
	double cov_irs_def = rcpp_to_double(params["cov_irs_def"]);					// Proportion of people protected by interior residual spraying (IRS) associated with data
	double cov_set_def = 0.0;													// Proportion of people protected by SET
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
	int na_m = 50;																// Number of age levels for adult mosquitoes	
	int days_3g = 12;															// Number of days for 3 gonotrophic cycles
	double* Sv= (double*)malloc(na_m*sizeof(double));
	double* Ev= (double*)malloc(na_m*sizeof(double));
	double* Iv= (double*)malloc(na_m*sizeof(double));

	// Adult mosquito parameters
	double Q0 = rcpp_to_double(params["Q0"]);									// Default anthropophagy
	double muv0 = rcpp_to_double(params["muv0"]);								// Default mosquito birth / death rate
	//double rho = rcpp_to_double(params["rho"]);								// Age-dependent biting parameter
	//double a0 = rcpp_to_double(params["a0"]);									// Age-dependent biting parameter
	//double sigma2 = rcpp_to_double(params["sigma2"]);							// Variance of log heterogeneity in biting
	//double latmosq = rcpp_to_double(params["latmosq"]);							// Time lag for infection in mosquitoes
	//double latgam = rcpp_to_double(params["latgam"]);							// Time lag mosquitoes to become infectious

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
	
	//vector<double> mu_atsb_var = rcpp_to_vector_double(trial_params["mu_atsb_var"]);	// Variable mu_atsb for special runs - TODO: Integrate into package
	int flag_dt_adjust = rcpp_to_int(trial_params["flag_dt_adjust"]);					// Flag indicating whether or not to adjust dt by mosquito density (0=no,1=yes)
	int flag_file = rcpp_to_int(trial_params["flag_file"]);								// Flag indicating whether results data to be saved to files
	int int_v_varied = rcpp_to_int(trial_params["int_v_varied"]);						// Intervention parameter given variable value (0=model parameter set, 1=ATSB kill rate, 2=bednet coverage, 3=IRS coverage, 4=SET coverage)
	int n_mv_values = rcpp_to_int(trial_params["n_mv_values"]);							// Number of mosquito density values to be used		
	vector<double> int_values = rcpp_to_vector_double(trial_params["int_values"]);		// Vector of intervention parameter values
	int n_int_values = int_values.size();												// Number of intervention parameter values to use
	double date_int = date_start + rcpp_to_double(trial_params["start_interval"]);		// Day of the year when intervention trial starts (period from date_start to date_int used to equilibrate lagged FOI data)
	vector<double> time_values = rcpp_to_vector_double(trial_params["time_values"]);	// Vector of time benchmark points
	int n_pts = rcpp_to_int(trial_params["n_pts"]);										// Number of data points to be taken
	//string file_benchmarks = flag_file == 1 ? rcpp_to_string(trial_params["file_benchmarks"]) : "na";
	string file_endpoints = flag_file == 1 ? rcpp_to_string(trial_params["file_endpoints"]) : "na";

	// Constant derived values--------------------------------------------------------------------------------------------------------------------------

	int n_runmax = n_int_values * n_mv_values;									// Total number of simulations to run
	double int_start_time = date_int - date_start;								// Time after start of run at which ATSB starts being used
	double tmax_c = time_values[n_pts-1];										// Duration of intervention use in days
	//int tmax_c_i = intdiv(tmax_c, 1.0);											// Duration of intervention use as integer
	double tmax = int_start_time + tmax_c;										// Total duration of run in days
	double inv_dy = 1.0 / dy;													// Inverse of days in a year (fraction of a year in a day)
	double trig_coeff1 = 2.0 * 3.1415926536 * inv_dy;							// Coefficient used in calculation of rainfall
	double trig_coeff2 = trig_coeff1 * 2.0;
	double trig_coeff3 = trig_coeff1 * 3.0;
	double trig_coeff4 = trig_coeff1 * 4.0;
	double trig_coeff5 = trig_coeff1 * 5.0;
	double rate_de = 1.0 / de;
	double rate_dl = 1.0 / dl;
	//double omega = 1.0 - ((rho * eta) / (eta + (1.0 / a0)));	
	//double inv_omega = 1.0 / omega;

	double rnyears = 8.0;
	double ppyear = 64.0;
	double beta_larval = (eov * muv0 * exp(-muv0 * dgon)) / (1.0 - exp(-muv0 * dgon));
	double b_lambda = ((gammal * mul) / mue) - (de / dl) + ((gammal - 1.0) * mul * de);
	double lambda = (-0.5 * b_lambda) + pow((0.25 * pow(b_lambda, 2.0)) + ((gammal * beta_larval * mul * de) / (2.0 * mue * muv0 * dl * (1.0 + (dp * mup)))), 0.5);
	double term = (lambda / (mue * de)) - (1.0 / (mul * dl)) - 1.0;
	double dt = 0.1;												// Time increment (default initial value)
	double rain_coeff = ((dy * Rnorm) / (14.0 * rnyears * ppyear));	// Rainfall coefficient
		
	// Set up additional parameters-------------------------------------------------------------------------------------------------------------------------------------

	//double* FOIv_lag = (double*)malloc(12500 * sizeof(double));			// Time-lagged force of human -> mosquito infection (to account for incubation period)
	//double* incv_lag = (double*)malloc(12500 * sizeof(double));			// Time-lagged incubation of malaria in infected mosquito
	vector<double> M_benchmarks(n_pts * n_runmax, 0.0);					// Values of total relative mosquito population (Sv+Ev+Iv) at checkpoints
	vector<double> M_spor_benchmarks(n_pts * n_runmax, 0.0);			// Sporozoite-positive mosquito population (Iv) at checkpoints
	vector<double> M_3g_benchmarks(n_pts * n_runmax, 0.0);				// Population of mosquitoes over 3 gonotrophic cycles at checkpoints

	vector<double> mv_input = rcpp_to_vector_double(inputs["mv_input"]);
	vector<double> EL_input = rcpp_to_vector_double(inputs["EL_input"]);
	vector<double> LL_input = rcpp_to_vector_double(inputs["LL_input"]);
	vector<double> PL_input = rcpp_to_vector_double(inputs["PL_input"]);
	vector<double> Sv_input = rcpp_to_vector_double(inputs["Sv_input"]);
	vector<double> Ev_input = rcpp_to_vector_double(inputs["Ev_input"]);
	vector<double> Iv_input = rcpp_to_vector_double(inputs["Iv_input"]);
	
	// Set up conditions for one or more runs-------------------------------------------------------------------------------------------------------------------------------------
			
	if (flag_file == 1) 
	{
		if (int_v_varied == 0) // Generate header for detailed endpoint data file
		{
			endpoint_data = fopen(file_endpoints.c_str(), "w");
			fprintf(endpoint_data, "n_run\tmv0\tEL\tLL\tPL");
			for (i = 0; i < na_m; i++) { fprintf(endpoint_data, "\tSv%i", i); }
			for (i = 0; i < na_m; i++) { fprintf(endpoint_data, "\tEv%i", i); }
			for (i = 0; i < na_m; i++) { fprintf(endpoint_data, "\tIv%i", i); }
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
		
		goto start_run;
		run_complete:
		Rcout << "\nRun " << n_run + 1 << " complete. dt=" << dt << "\n";
		R_FlushConsole();
	}

	goto finish;

	//-----------------------------------------------------------------------------Set up initial data for run--------------------------------------------------------------------

	start_run:

	param_int = int_v_varied == 0 ? 0.0 : int_values[n_int];															// Probability of treatment not being received by clinical cases				
	rN = cov_nets * p_repel_net_entry;
	rNW = cov_nets * p_repel_net_bed;
	dNW = cov_nets * p_kill_net;
	rI = cov_irs * p_repel_irs_entry;
	rIW = cov_irs * p_repel_irs_bed;
	dIW = cov_irs * p_kill_irs1;
	dIF = cov_irs * p_kill_irs2;
	K0 = (mv0 * 2.0 * dl * muv0 * (1.0 + (dp * mup)) * (gammal * (lambda + 1.0))) / term;									// Larval birth rate coefficient, set based on mv0
	//r_bite_prevent = phi_bite_house - ((phi_bite_house - phi_bite_bed) * (1.0 - rI) * (1.0 - rIW - dIW) * (1.0 - dIF)) -	// Proportion of bites prevented by bednets and/or interior residual spraying
	//		(phi_bite_bed * (1.0 - rI) * (1.0 - rIW - dIW) * (1.0 - dIF) * (1.0 - rN) * (1.0 - rNW - dNW));
	mu_protections = 0.3333 * Q0 * (((phi_bite_house - phi_bite_bed) * (1.0 - rI) * (dIW + ((1.0 - rIW - dIW) * dIF))) +		// Additional death rate due to bednets and/or interior residual spraying
			(phi_bite_bed * (1.0 - rI) * (1.0 - rN) * (dIW + ((1.0 - rIW - dIW) * (dNW + ((1.0 - rNW - dNW) * dIF))))));
	//Values set manually when situation not accommodated by model
	//r_bite_prevent=0.707;
	//mu_protections=0.0865;
	//av0 = 0.3333 * Q0 * (1.0 - r_bite_prevent);																			// Biting rate on humans / mosquito
	muv1 = muv0 + mu_atsb + mu_protections;																						// Total mosquito death rate
	// muv1 = mu_atsb_var[0] + muv0 + mu_protections;																			// Variable muv1 for special runs - TODO: Integrate with package
	Rcout << "\nRun " << n_run + 1 << " Mosquito density=" << mv0 << " Intervention number=" << n_int << " dt=" << dt;

	restart_run:
	dt2 = dt * 0.1;		
	ntmax = intdiv(tmax, dt);
	//dur_Ei = intdiv(dur_E, dt);
	//latgami = intdiv(latgam, dt);
	//latmosqi = intdiv(latmosq, dt);
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
	for(i=0; i<na_m; i++)
	{
		Sv[i]=Sv_input[i];
		Ev[i]=Ev_input[i];
		Iv[i]=Iv_input[i];
	}
	//EIRd = av0 * Iv1 * inv_omega;

	/*Surv1 = exp(-muv1 * latmosq);
	incv0 = FOIv0 * Sv1 * Surv1;
	for (i = 0; i < 12500; i++)
	{
		FOIv_lag[i] = FOIv0;
		incv_lag[i] = incv0;
		ij = i * n_cats;
		for (pos = 0; pos < n_cats; pos++) { FOI_lag[ij + pos] = FOI0[pos]; }
	}*/

	//-------------------------------------------------------------------------------------Compute over time--------------------------------------------------------------------------------------------------------

	div = 0;
	t_mark1 = 1.0;
	t_mark2 = time_values[0];
	flag_int = 0;
	for (nt = 0; nt <= ntmax; nt++)
	{
		tflag = 0;								// Flag indicates whether dt needs to be altered
		t = nt * dt;							// Time since start of simulation
		t_int = t - int_start_time;				// Time since start of intervention use
		//nt_E = nt % dur_Ei;
		//nt_latgam = nt % latgami;
		//nt_latmosq = nt % latmosqi;
		year_day = fmod(t + date_start, dy);	// Day of the year
		mv = arraysum(Sv, na_m) + arraysum(Ev, na_m) + arraysum(Iv, na_m);					// Total number of mosquitoes
		if (flag_int == 0 && t_int >= 0.0)		// Value of varied trial parameters set after start of intervention
		{
			flag_int = 1;
			switch (int_v_varied)
			{
				case 0: {}
					break;
				case 1:
				{
					/*cov_nets=0.6;
					rN = cov_nets * p_repel_net_entry;
					rNW = cov_nets * p_repel_net_bed;
					dNW = cov_nets * p_kill_net;*/
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
				//r_bite_prevent = phi_bite_house-((phi_bite_house-phi_bite_bed)*(1.0-rSET-dSET))-(phi_bite_bed*(1.0-rSET-dSET)*(1.0-rNW-dNW));
				mu_protections = 0.3333 * Q0 * (((phi_bite_house-phi_bite_bed)*dSET)+(phi_bite_bed*(dSET+(dNW*(1.0-rSET-dSET)))));
			}
			else
			{
				//r_bite_prevent = phi_bite_house - ((phi_bite_house - phi_bite_bed) * (1.0 - rI) * (1.0 - rIW - dIW) * (1.0 - dIF)) -
				//				(phi_bite_bed * (1.0 - rI) * (1.0 - rIW - dIW) * (1.0 - dIF) * (1.0 - rN) * (1.0 - rNW - dNW));
				mu_protections = 0.3333 * Q0 * (((phi_bite_house - phi_bite_bed) * (1.0 - rI) * (dIW + ((1.0 - rIW - dIW) * dIF))) +
								(phi_bite_bed * (1.0 - rI) * (1.0 - rN) * (dIW + ((1.0 - rIW - dIW) * (dNW + ((1.0 - rNW - dNW) * dIF))))));
				//Values set manually when situation not accommodated by model
				//r_bite_prevent=0.707;
				//mu_protections=0.0865;
			}
			//av0 = 0.3333 * Q0 * (1.0 - r_bite_prevent);
			/*if (t_int >= 0.0 && n_int > 0) // Variable muv1 for special runs - TODO: Integrate into package
			{
				int day = t_int;
				muv1 = mu_atsb_var[day] + muv0 + mu_protections;
			}*/
			muv1 = muv0 + mu_atsb + mu_protections;

		}

		rain = rain_coeff * (rconst + (2.0 * (
				(Re0 * cos(trig_coeff1 * year_day)) + (Im0 * sin(trig_coeff1 * year_day)) +
				(Re1 * cos(trig_coeff2 * year_day)) + (Im1 * sin(trig_coeff2 * year_day)) +
				(Re2 * cos(trig_coeff3 * year_day)) + (Im2 * sin(trig_coeff3 * year_day)) +
				(Re3 * cos(trig_coeff4 * year_day)) + (Im3 * sin(trig_coeff4 * year_day)) +
				(Re4 * cos(trig_coeff5 * year_day)) + (Im4 * sin(trig_coeff5 * year_day))
				)));
		KL = max(K0 * rain, KMIN);
		beta = ((0.5 * PL) / dp);		// Rate of spawning of adult mosquitoes
		//EIRd = av0 * Iv1 * inv_omega;	// Daily entomological inoculation rate
		
		// TODO - Re-introduce infection of mosquitoes
		/*FOIv = t > latgam ? FOIv_lag[nt_latgam] : FOIv0;
		FOIv_lag[nt_latgam] = arraysum(FOIvij, n_cats);
		incv = t > latmosq ? incv_lag[nt_latmosq] : incv0;
		Surv1 = exp(-muv1 * latmosq);
		incv_lag[nt_latmosq] = FOIv * Sv1 * Surv1;*/
		inv_KL = 1.0 / KL;

		// Larval data calculation; calculated in smaller time increments
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

		// Adult mosquito data calculation
		i=na_m-1;
		Sv_cur = Sv[i];
		Ev_cur = Ev[i];
		Iv_cur = Iv[i];
		delSv = 0.0;//TODO - FOIv * Sv_cur;
		incv=0.0; //TODO - incv_lag[]
		Sv[i] += dt * (-(muv1 * Sv_cur) - delSv + Sv[i-1]);
		Ev[i] += dt * (-(muv1 * Ev_cur) + delSv - incv + Ev[i-1]); //TODO - incv
		Iv[i] += dt * (-(muv1 * Iv_cur) + incv + Iv[i-1]); //TODO - incv	
		for(i=na_m-1; i>0; i--)
		{
			Sv_cur = Sv[i];
			Ev_cur = Ev[i];
			Iv_cur = Iv[i];
			delSv = 0.0;//TODO - FOIv * Sv_cur;
			incv=0.0; //TODO - incv_lag[]
			Sv[i] += dt * (-(muv1 * Sv_cur) - delSv + Sv[i-1] - Sv_cur);
			Ev[i] += dt * (-(muv1 * Ev_cur) + delSv - incv + Ev[i-1] - Ev_cur); //TODO - incv
			Iv[i] += dt * (-(muv1 * Iv_cur) + incv + Iv[i-1] - Iv_cur); //TODO - incv		
		}		
		Sv_cur=Sv[0];
		Ev_cur = Ev[0];
		Iv_cur = Iv[0];
		delSv = 0.0;//TODO - FOIv * Sv_cur;
		incv=0.0; //TODO - incv_lag[]
		Sv[0] += dt*(beta - (muv1*Sv_cur) - Sv_cur);
		Ev[0] += dt*( -(muv1*Ev_cur) + delSv - incv - Ev_cur); //TODO - incv
		Iv[0] += dt*( -(muv1*Iv_cur) +incv - Iv_cur); //TODO - incv		
		
		// Outputs
		if (t >= t_mark1)
		{
			t_mark1 += 1.0;
			// Check that time increment has not been set too large
			if(arraysum(Sv, na_m) + arraysum(Ev, na_m) + arraysum(Iv, na_m)<=0.0){ tflag = 1; }
			if (min(EL,min(LL,PL))<0.0) { tflag = 2; }
			if (min(arraymin(Sv,na_m),min(arraymin(Ev,na_m),arraymin(Iv,na_m)))<0.0) { tflag = 3; }
			if (tflag > 0)
			{
				if (dt <= 0.001) { goto end_run; }
				dt = max(0.001, dt * 0.9);
				Rcout << "\ntflag = " << tflag << " t = " << t << " Adjusting time increment. New dt = " << dt << "\n";
				R_FlushConsole();
				goto restart_run;
				goto end_run;
			}
		}

		if (t_int >= t_mark2)
		{
			pos2 = (n_pts * n_run) + div;
			M_spor_benchmarks[pos2] = arraysum(Iv, na_m);
			M_benchmarks[pos2] = arraysum(Sv, na_m) + arraysum(Ev, na_m) + M_spor_benchmarks[pos2];
			M_3g_benchmarks[pos2] = M_benchmarks[pos2] - arraysum(Sv, days_3g) - arraysum(Ev, days_3g) - arraysum(Iv, days_3g);
			div++;
			if (div < n_pts) { t_mark2 = time_values[div]; }
			R_FlushConsole();
		}
	}

	end_run:
		
	if (flag_file == 1)
	{
		if (int_v_varied == 0)
		{
			// Output detailed data to use as input for future computation
			endpoint_data = fopen(file_endpoints.c_str(), "a");
			fprintf(endpoint_data, "%i\t%.3e\t%.3e\t%.3e\t%.3e", n_run + 1, mv0, EL, LL, PL);
			for (i = 0; i < na_m; i++) { fprintf(endpoint_data, "\t%.3e", Sv[i]); }
			for (i = 0; i < na_m; i++) { fprintf(endpoint_data, "\t%.3e", Ev[i]); }
			for (i = 0; i < na_m; i++) { fprintf(endpoint_data, "\t%.3e", Iv[i]); }
			fprintf(endpoint_data,"\n");
			fclose(endpoint_data);
		}
	}

	goto run_complete;

	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	finish:
	//Rcout << "\nComputations complete 1.\n";
	//R_FlushConsole();

	if (endpoint_data != NULL) { fclose(endpoint_data); }
	Rcout << "\nComputations complete.\n";
	R_FlushConsole();

	// Return list
	//return List::create(Named("dummy") = 1.0);
	return List::create(Named("M_benchmarks") = M_benchmarks, Named("M_spor_benchmarks") = M_spor_benchmarks,
						Named("M_3g_benchmarks") = M_3g_benchmarks);
}