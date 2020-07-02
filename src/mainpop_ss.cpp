// CPP file for simulation of malaria intervention trial - main population portion
#include <Rcpp.h>
#include "main.h"
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
Rcpp::List rcpp_mainpop_ss(List params, List trial_params)
{
	int n_mv, i, j, pos, pos2, n_repeats;
	double mv0, mu_net_irs, prevent_net_irs, av0, muv1, KL, Surv1, EIRd, FOIv0, FOIv1;
	double rN, rNW, dNW, rI, rIW, dIW, dIF, EL, LL, PL, Sv1, Ev1, Iv1, H, H_inv, dev, EIR_cur, FOI_cur, phi_cur, slide_prev,pcr_prev;
	FILE* endpoint_data = NULL;

	//Constants ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	double dy = 365.0;		// Days in a year

	// Load environment/vector/parasite/etc. parameters from R------------------------------------------------------------------------------------------------------------------------------------------------

	string file_endpoints = rcpp_to_string(trial_params["file_endpoints"]);
	vector<double> EIRy = rcpp_to_vector_double(trial_params["EIRy"]);
	int n_mv_values = EIRy.size();
	vector<double> mv_values(n_mv_values, 0.0);
	
	// Heterogeneity parameters
	int num_het = rcpp_to_int(params["num_het"]);							// Number of age categories in main population
	vector<double> het_x = rcpp_to_vector_double(params["het_x"]);			// Biting heterogeneity
	vector<double> het_wt = rcpp_to_vector_double(params["het_wt"]);		// Fractions in each heterogeneity category

	// Age parameters
	int na = rcpp_to_int(params["na"]);										// Number of age categories in main population
	vector<double> age_width = rcpp_to_vector_double(params["age_width"]);	// Widths of age categories in days
	vector<double> age_rate = rcpp_to_vector_double(params["age_rate"]);	// Rates of transition between age categories
	vector<double> rem_rate = rcpp_to_vector_double(params["rem_rate"]);	// Total loss rate from each age category (ageing + death)
	vector<double> age = rcpp_to_vector_double(params["age"]);				// Age values corresponding to starts of age categories (days)
	vector<double> den = rcpp_to_vector_double(params["den"]);
	double eta = rcpp_to_double(params["eta"]);								// Death rate, for exponential population distribution

	double mu_atsb = rcpp_to_double(params["mu_atsb_def"]);					// Attractive targeted sugar bait (ATSB) killing rate associated with data
	double cov_nets = rcpp_to_double(params["cov_nets_def"]);				// Proportion of people protected by bednets associated with data
	double cov_irs = rcpp_to_double(params["cov_irs_def"]);					 // Proportion of people protected by interior residual spraying (IRS) associated with data
	double prop_T = rcpp_to_double(params["prop_T_def"]);					// Proportion of clinical cases successfully treated in main population associated with data

	double Q0 = rcpp_to_double(params["Q0"]);								// Default anthropophagy
	double muv0 = rcpp_to_double(params["muv0"]);							// Default mosquito birth / death rate
	double phi_bite_bed = rcpp_to_double(params["phi_bite_bed"]);			// Default proportion of bites delivered to people in bed
	double p_repel_net_entry = rcpp_to_double(params["p_repel_net_entry"]);	// Default probability of mosquito being repelled by bednet before entering house
	double p_repel_net_bed = rcpp_to_double(params["p_repel_net_bed"]);		// Default probability of mosquito being repelled by bednet after entering house but before biting
	double p_kill_net = rcpp_to_double(params["p_kill_net"]);				// Default probability of mosquito being killed by bednet while trying to bite
	double phi_bite_house = rcpp_to_double(params["phi_bite_house"]);		// Default proportion of bites delivered to people in houses
	double p_repel_irs_entry = rcpp_to_double(params["p_repel_irs_entry"]);	// Default probability of mosquito being repelled by IRS before entering house
	double p_repel_irs_bed = rcpp_to_double(params["p_repel_irs_bed"]);		// Default probability of mosquito being repelled by IRS after entering house but before biting
	double p_kill_irs1 = rcpp_to_double(params["p_kill_irs1"]);				// Default probability of mosquito being killed by IRS while trying to bite
	double p_kill_irs2 = rcpp_to_double(params["p_kill_irs2"]);				// Default probability of mosquito being killed by IRS after biting

	double rho = rcpp_to_double(params["rho"]);								// Age-dependent biting parameter
	double a0 = rcpp_to_double(params["a0"]);								// Age-dependent biting parameter
	double sigma2 = rcpp_to_double(params["sigma2"]);						// Variance of log heterogeneity in biting
	double dur_T = rcpp_to_double(params["dur_T"]);							// Average time to recover from disease and parasitaemia when treated
	double dur_D = rcpp_to_double(params["dur_D"]);							// Average time to move from clinical disease to asymptomatic when not successfully treated
	double dur_A = rcpp_to_double(params["dur_A"]);							// Average time to move from asymptomatic patent infection to sub-patent infection
	double dur_U = rcpp_to_double(params["dur_U"]);							// Average time to recover from sub-patent infection
	double dur_P = rcpp_to_double(params["dur_P"]);							// Average time to leave protected state
	double latmosq = rcpp_to_double(params["latmosq"]);						// Time lag for infection in mosquitoes

	// Infectivity to mosquitoes
	double cD = rcpp_to_double(params["cD"]);								// Untreated disease
	double cU = rcpp_to_double(params["cU"]);								// Sub-patent infection
	double cT = rcpp_to_double(params["cT"]);								// Treated disease
	double gamma_inf = rcpp_to_double(params["gamma_inf"]);					// Parameter for infectiousness of state A

	// Immunity parameters - infection
	double bh = rcpp_to_double(params["bh"]); 
	double bmin = rcpp_to_double(params["bmin"]);							// Maximum relative reduction due to immunity
	double db = rcpp_to_double(params["db"]);								// Inverse of decay rate
	double kb = rcpp_to_double(params["kb"]);								// Shape parameter
	double ub = rcpp_to_double(params["ub"]);								// Duration in which immunity is not boosted
	double inv_IB0 = 1.0 / rcpp_to_double(params["IB0"]);					// Inverse of scale parameter

	// Immunity parameters - detection
	double dmin = rcpp_to_double(params["dmin"]);							// Minimum probability due to maximum immunity
	double dd = rcpp_to_double(params["dd"]);								// Inverse of decay rate
	double kd = rcpp_to_double(params["kd"]);								// Shape parameter
	double ud = rcpp_to_double(params["ud"]);								// Duration in which immunity is not boosted
	double ad0 = rcpp_to_double(params["ad0"]);								// Scale parameter relating age to immunity
	double fd0 = rcpp_to_double(params["fd0"]);								// Time scale at which immunity changes with age
	double gammad = rcpp_to_double(params["gammad"]);						// Shape parameter relating age to immunity
	double inv_ID0 = 1.0 / rcpp_to_double(params["ID0"]);					// Scale parameter

	// Immunity parameters - clinical disease
	double phi0 = rcpp_to_double(params["phi0"]);							// Maximum probability due to no immunity
	double phi1 = rcpp_to_double(params["phi1"]);							// Maximum relative reduction due to no immunity
	double dc = rcpp_to_double(params["dc"]);								// Inverse of decay rate
	double kc = rcpp_to_double(params["kc"]);								// Shape parameter
	double uc = rcpp_to_double(params["uc"]);								// Duration in which immunity is not boosted
	double P_IC_M = rcpp_to_double(params["P_IC_M"]);						// 
	double dm = rcpp_to_double(params["dm"]);								// Inverse of decay rate of maternal immunity
	double inv_IC0 = 1.0 / rcpp_to_double(params["IC0"]);					// Scale parameter

	// Larval parameters	
	double mue = rcpp_to_double(params["mue"]);								// Death rate of early instar
	double mul = rcpp_to_double(params["mul"]);								// Death rate of late instar
	double mup = rcpp_to_double(params["mup"]);								// Death rate of pupae
	double de = rcpp_to_double(params["de"]);								// Duration of early instar
	double dl = rcpp_to_double(params["dl"]);								// Duration of late instar
	double dp = rcpp_to_double(params["dp"]);								// Duration of pupae
	double eov = rcpp_to_double(params["eov"]);								// Eggs per ovulation
	double gammal = rcpp_to_double(params["gammal"]);						// Density dependence term
	double dgon = rcpp_to_double(params["dgon"]);							// Gonotrophic cycle length

	// Set up conditions for one or more runs-------------------------------------------------------------------------------------------------------------------------------------
	
	int n_cats = na * num_het;												// Total number of age/heterogeneity categories in main population
	endpoint_data = fopen(file_endpoints.c_str(), "w");
	fprintf(endpoint_data, "n_mv\tmv0\tEL\tLL\tPL\tSv1\tEv1\tIv1");
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

	// Set up arrays-------------------------------------------------------------------------------------------------------------------------------------

	int size1 = na * sizeof(double);										// Array dimensions
	int size2 = num_het * sizeof(double);
	int size3 = n_cats * sizeof(double);

	//Output data
	double* S = (double*)malloc(size3);										// Susceptible humans
	double* T = (double*)malloc(size3);										// Humans in treatment for clinical symptoms
	double* D = (double*)malloc(size3);										// Humans with clinical symptoms not receiving treatment
	double* A = (double*)malloc(size3);										// Humans with asymptomatic patent infections
	double* U = (double*)malloc(size3);										// Humans with asymptomatic sub-patent infections
	double* P = (double*)malloc(size3);										// Humans in post-treatment stage
	double* ICA = (double*)malloc(size3);									// Immunity level against clinical infection (adult)
	double* ICM = (double*)malloc(size3);									// Immunity level against clinical infection (maternally inherited)
	double* IB = (double*)malloc(size3);									// Immunity level against infection from infectious bite
	double* ID = (double*)malloc(size3);									// Immunity level determining p_det

	//Other data
	double* Y = (double*)malloc(size3);
	double* betaS = (double*)malloc(size3);
	double* betaT = (double*)malloc(size3);
	double* betaD = (double*)malloc(size3);
	double* betaA = (double*)malloc(size3);
	double* betaU = (double*)malloc(size3);
	double* betaP = (double*)malloc(size3);
	double* aT = (double*)malloc(size3);
	double* aD = (double*)malloc(size3);
	double* aP = (double*)malloc(size3);
	double* den_het = (double*)malloc(size3);

	// Constant derived values--------------------------------------------------------------------------------------------------------------------------

	double prop_T_rev = 1.0 - prop_T;										// Probability of treatment not being received by clinical cases	
	double inv_dy = 1.0 / dy;												// Inverse of days in a year (fraction of a year in a day)
	double cDU = cD - cU;
	double rT = 1.0 / dur_T;												// Rate of recovery from disease and parasitaemia when treated
	double rD = 1.0 / dur_D;												// Rate of moving from clinical disease to asymptomatic when not successfully treated
	double rA0 = 1.0 / dur_A;												// Rate of moving from asymptomatic patent infection to sub-patent infection
	double rU = 1.0 / dur_U;												// Rate of recovery from sub - patent infection
	double rP = 1.0 / (dur_P - dur_T);										// Rate of leaving protected state
	double bmin_rev = 1.0 - bmin;
	double dmin_rev = 1.0 - dmin;
	double phi1_rev = 1.0 - phi1;
	double* delta = (double*)malloc(size1);
	double* gamma = (double*)malloc(size1);
	double* foi_age = (double*)malloc(size1);
	int* age20i = (int*)malloc(na * sizeof(int));

	// Set age category parameters
	pos = 0;
	for (i = 0; i < na; i++)
	{
		gamma[i] = i == na - 1 ? 0.0 : rem_rate[i];
		foi_age[i] = 1.0 - (rho * exp(-age[i] / a0));
		if (i == 0) { age20i[i] = 0; }
		else
		{
			if (age[i] >= 20.0 * dy && age[i - 1] < 20.0 * dy) { age20i[i] = i; }
			else { age20i[i] = age20i[i - 1]; }
		}
		for (j = 0; j < num_het; j++) 
		{ 
			den_het[pos] = den[i] * het_wt[j];
			betaT[pos] = rT + gamma[i];
			betaD[pos] = rD + gamma[i];
			betaP[pos] = rP + gamma[i];
			pos++;
		}
	}

	int age20u = age20i[na - 1];
	int age20l = age20u - 1;
	double age_20_factor = (((20 * dy) - age[age20l] - (0.5 * age_width[age20l])) * 2.0) / (age_width[age20l] + age_width[age20u]);

	double* rel_foi = (double*)malloc(size2);
	double* b = (double*)malloc(size3);					// Probability of infection from an infectious bite
	double* p_det = (double*)malloc(size3);				// Probability of asymptomatic infection being detected
	double* fd = (double*)malloc(size1);				//
	double* ICM_init = (double*)malloc(size2);			// Value of maternally inherited immunity in newborns
	double* rate_ibaq = (double*)malloc(size3);			//
	double* rate_clinaq = (double*)malloc(size3);		//
	double* rate_detaq = (double*)malloc(size3);		//
	double* cA = (double*)malloc(size3);				//
	double* FOIvij = (double*)malloc(size3);			//Force of infection (human -> mosquito)

	for (j = 0; j < num_het; j++) { rel_foi[j] = num_het == 0 ? 1.0 : (exp((-sigma2*0.5) + (pow(sigma2, 0.5) * het_x[j]))); }	// Relative biting rate in heterogeneity category i
	double omega = 1.0 - ((rho * eta) / (eta + (1.0 / a0)));																	// Normalising constant for biting rate by age
	double* x_I = (double*)malloc(size1);
	double inv_omega = 1.0 / omega;
	x_I[0] = den[0] / eta;
	fd[0] = 1.0 - ((1.0 - fd0) / (1.0 + pow(age[0] / ad0, gammad)));
	delta[0] = 0.0;
	for (i = 1; i < na; i++)
	{
		delta[i] = age_rate[i - 1];
		x_I[i] = den[i] / (den[i - 1] * age_rate[i - 1]);
		fd[i] = 1.0 - ((1.0 - fd0) / (1.0 + pow(age[i] / ad0, gammad)));
	}

	double beta_larval = (eov * muv0 * exp(-muv0 * dgon)) / (1.0 - exp(-muv0 * dgon));
	double b_lambda = ((gammal * mul) / mue) - (de / dl) + ((gammal - 1.0) * mul * de);
	double lambda = (-0.5 * b_lambda) + pow((0.25 * pow(b_lambda, 2.0)) + ((gammal * beta_larval * mul * de) / (2.0 * mue * muv0 * dl * (1.0 + (dp * mup)))), 0.5);
	double term = (lambda / (mue * de)) - (1.0 / (mul * dl)) - 1.0;

	//------------------------------------------------------------------------------------------------------------------------------------------------
		
	rN = cov_nets * p_repel_net_entry;
	rNW = cov_nets * p_repel_net_bed;
	dNW = cov_nets * p_kill_net;
	rI = cov_irs * p_repel_irs_entry;
	rIW = cov_irs * p_repel_irs_bed;
	dIW = cov_irs * p_kill_irs1;
	dIF = cov_irs * p_kill_irs2;
	prevent_net_irs = phi_bite_house - ((phi_bite_house - phi_bite_bed) * (1.0 - rI) * (1.0 - rIW - dIW) * (1.0 - dIF)) -	// Proportion of bites prevented by bednets and/or interior residual spraying
		(phi_bite_bed * (1.0 - rI) * (1.0 - rIW - dIW) * (1.0 - dIF) * (1.0 - rN) * (1.0 - rNW - dNW));
	mu_net_irs = 0.3333 * Q0 * (((phi_bite_house - phi_bite_bed) * (1.0 - rI) * (dIW + ((1.0 - rIW - dIW) * dIF))) +		// Additional death rate due to bednets and/or interior residual spraying
		(phi_bite_bed * (1.0 - rI) * (1.0 - rN) * (dIW + ((1.0 - rIW - dIW) * (dNW + ((1.0 - rNW - dNW) * dIF))))));
	//Values set manually when situation not accommodated by model
	//prevent_net_irs=0.707;
	//mu_net_irs=0.0865;
	av0 = 0.3333 * Q0 * (1.0 - prevent_net_irs);																			// Biting rate on humans / mosquito
	muv1 = muv0 + mu_atsb + mu_net_irs;																						// Total mosquito death rate
	Surv1 = exp(-muv1 * latmosq); 
	
	for (n_mv = 0; n_mv < n_mv_values; n_mv++)
	{
		EIRd = EIRy[n_mv] * inv_dy;

		// Initial human characteristics
		pos = 0;
		for (i = 0; i < na; i++) //Run through age categories
		{
			for (j = 0; j < num_het; j++) //Run through heterogeneity categories
			{
				pos2=pos-num_het;
				EIR_cur = EIRd * foi_age[i] * rel_foi[j];
				FOI_cur = EIR_cur * bh;
				rate_ibaq[pos] = EIR_cur / ((EIR_cur * ub) + 1.0);
				rate_detaq[pos] = FOI_cur / ((FOI_cur * ud) + 1.0);
				rate_clinaq[pos] = FOI_cur / ((FOI_cur * uc) + 1.0);
				S[pos] = den_het[pos];
				T[pos] = 0.0;
				D[pos] = 0.0;
				A[pos] = 0.0;
				U[pos] = 0.0;
				P[pos] = 0.0;
				ICA[pos] = ((i == 0 ? 0.0 : ICA[pos2]) + (rate_clinaq[pos] * x_I[i])) / (1.0 + (x_I[i] / dc));
				ICM[pos] = 1.0;
				IB[pos] = ((i == 0 ? 0.0 : IB[pos2]) + (rate_ibaq[pos] * x_I[i])) / (1.0 + (x_I[i] / db));
				ID[pos] = ((i == 0 ? 0.0 : ID[pos2]) + (rate_detaq[pos] * x_I[i])) / (1.0 + (x_I[i] / dd));				 
				pos++;
			}
		}
		//Initial mosquito density
		mv0 = (0.2345 * EIRy[n_mv]) + 0.1216;
		
		//Calculations repeated until steady state reached
		FOIv1 = 0.0;
		n_repeats=0;
		repeat:
		n_repeats++;
		pos = 0;
		for (i = 0; i < na; i++) //Run through age categories
		{
			for (j = 0; j < num_het; j++) //Run through heterogeneity categories
			{
				b[pos] = bh * ((bmin_rev / (1.0 + pow(IB[pos] * inv_IB0, kb))) + bmin);
				p_det[pos] = dmin + (dmin_rev / (1.0 + (fd[i] * pow(ID[pos] * inv_ID0, kd))));
				cA[pos] = cU + (cDU * pow(p_det[pos], gamma_inf));
				FOIvij[pos] = foi_age[i] * av0 * ((cT * T[pos]) + (cD * D[pos])+(cA[pos] * A[pos]) + (cU * U[pos])) * rel_foi[j] * inv_omega;
				pos++;
			}
		}
		FOIv0 = FOIv1;
		FOIv1 = arraysum(FOIvij,n_cats);

		// Mosquito characteristics
		Iv1 = (mv0 * FOIv1 * Surv1) / (FOIv1 + muv1);
		Sv1 = (muv1 * Iv1) / (FOIv1 * Surv1);
		Ev1 = mv0 - Sv1 - Iv1;
		mv0 = (omega * EIRd * (FOIv1 + muv1)) / ((FOIv1 * Surv1) * av0);
		KL = (mv0 * (2.0 * dl * muv1 * (1.0 + (dp * mup))) * (gammal * (lambda + 1.0))) / term;
		
		// Larval characteristics
		EL = KL * (lambda / (gammal * (lambda + 1.0))) * term;
		LL = (KL / (gammal * (lambda + 1.0))) * term;
		PL = KL * (dp / (dl * (1.0 + (dp * mup)))) * (1.0 / (gammal * (lambda + 1.0))) * term;
	
		// Human characteristics
		for (j = 0; j < num_het; j++) { ICM_init[j] = P_IC_M * (ICA[(age20l * num_het) + j] + age_20_factor * (ICA[(age20u * num_het) + j] - ICA[(age20l * num_het) + j])); }
		pos = 0;
		for (i = 0; i < na; i++) //Run through age categories
		{
			for (j = 0; j < num_het; j++) //Run through heterogeneity categories
			{
				pos2 = pos - num_het; 
				EIR_cur = EIRd * foi_age[i] * rel_foi[j];
				FOI_cur = EIR_cur * (IB[pos] > 0.0 ? b[pos] : bh);
				phi_cur = (phi0 *phi1_rev) / (1.0 + pow((ICM[pos] + ICA[pos]) * inv_IC0, kc) + phi1);
				rate_ibaq[pos] = EIR_cur / ((EIR_cur * ub) + 1.0);
				rate_detaq[pos] = FOI_cur / ((FOI_cur * ud) + 1.0);
				rate_clinaq[pos] = FOI_cur / ((FOI_cur * uc) + 1.0);

				betaS[pos] = FOI_cur + gamma[i];
				betaA[pos] = (FOI_cur * phi_cur) + rA0 + gamma[i];
				betaU[pos] = FOI_cur + rU + gamma[i];
				aT[pos] = (FOI_cur * phi_cur * prop_T) / betaT[pos];
				aD[pos] = (FOI_cur * phi_cur * prop_T_rev) / betaD[pos];
				aP[pos] = (rT * aT[pos]) / betaP[pos];

				if (i == 0)
				{
					Y[pos] = den_het[pos] / (1.0 + aT[pos] + aD[pos] + aP[pos]);
					T[pos] = aT[pos] * Y[pos];
					D[pos] = aD[pos] * Y[pos];
					P[pos] = aP[pos] * Y[pos];
				}
				else
				{
					Y[pos] = (den_het[pos] - (delta[i] * ((T[pos2] / betaT[pos]) + (D[pos2] / betaD[pos]) + ((((rT * T[pos2]) / betaT[pos]) + (P[pos2])/ betaP[pos]))))) / (1.0 + aT[pos] + aD[pos] + aP[pos]);
					T[pos] = (aT[pos] * Y[pos]) + ((delta[i] * T[pos2]) / betaT[pos]);
					D[pos] = (aD[pos] * Y[pos]) + ((delta[i] * D[pos2]) / betaD[pos]);
					P[pos] = (aP[pos] * Y[pos]) + (delta[i] * (((rT * T[pos2]) / betaT[pos]) + P[pos2]) / betaP[pos]);
				}

				A[pos] = ((i == 0 ? 0.0 : (delta[i] * A[pos2])) + (FOI_cur * (1.0 - phi_cur) * Y[pos]) + (rD * D[pos])) / (betaA[pos] + (FOI_cur * (1 - phi_cur)));
				U[pos] = ((rA0 * A[pos]) + (i == 0 ? 0.0 : delta[i] * U[pos2])) / betaU[pos];
				S[pos] = max(0.0, Y[pos] - A[pos] - U[pos]);
				ICA[pos] = ((i == 0 ? 0.0 : ICA[pos2]) + (rate_clinaq[pos] * x_I[i])) / (1.0 + (x_I[i] / dc));
				ICM[pos] = i == na ? 0.0 : ICM_init[j]*((exp(-age[i]/dm)) - (exp(-age[i + 1]/dm)))*(dm / (age[i + 1] - age[i]));
				IB[pos] = ((i == 0 ? 0.0 : IB[pos2]) + (rate_ibaq[pos] * x_I[i])) / (1.0 + (x_I[i] / db));
				ID[pos] = ((i == 0 ? 0.0 : ID[pos2]) + (rate_detaq[pos] * x_I[i])) / (1.0 + (x_I[i] / dd));
				pos++;
			}
		}

		pos=0;
		slide_prev = 0.0;
		pcr_prev = 0.0;
		H = arraysum(S, n_cats) + arraysum(T, n_cats) + arraysum(D, n_cats) + arraysum(A, n_cats) + arraysum(U, n_cats) + arraysum(P, n_cats);
		H_inv = 1.0 / H;
		for (i = 0; i < na; i++) // Run through age categories
		{
			for (j = 0; j < num_het; j++)
			{
				slide_prev+= T[pos] + D[pos] + (p_det[pos] * A[pos]);
				pcr_prev+= T[pos] + D[pos] + A[pos] + U[pos];
				pos++;
			}
		}
		slide_prev*=H_inv;
		pcr_prev*=H_inv;

		if (FOIv1 <= 0.0) { dev = 1.0; } 
		else { dev = abs(FOIv1 - FOIv0) / (FOIv1 + FOIv0); }
		if (dev > 1.0e-8) { goto repeat; }

		//finish:
		mv_values[n_mv] = mv0;
		Rcout << "\nn_mv = " << n_mv + 1 << "\tMosquito density = " << mv0 << "\tslide_prev = " << slide_prev << "\tpcr_prev = " << pcr_prev;
		R_FlushConsole();
		H = arraysum(S, n_cats) + arraysum(T, n_cats) + arraysum(D, n_cats) + arraysum(A, n_cats) + arraysum(U, n_cats) + arraysum(P, n_cats);
		H_inv = 1.0 / H;
		pos = 0;
		for (i = 0; i < na; i++) //Run through age categories
		{
			for (j = 0; j < num_het; j++) //Run through heterogeneity categories
			{
				// Normalize S,T,D,A,U,P to sum to 1.0
				S[pos] *= H_inv;
				T[pos] *= H_inv;
				D[pos] *= H_inv;
				A[pos] *= H_inv;
				U[pos] *= H_inv;
				P[pos] *= H_inv;
				pos++;
			}
		}

		// Output data to use as input for future computation
		endpoint_data = fopen(file_endpoints.c_str(), "a");
		fprintf(endpoint_data, "%i\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e", n_mv + 1, mv0, EL, LL, PL, Sv1, Ev1, Iv1);
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\t%.4e", S[i]); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\t%.4e", T[i]); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\t%.4e", D[i]); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\t%.4e", A[i]); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\t%.4e", U[i]); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\t%.4e", P[i]); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\t%.4e", ICA[i]); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\t%.4e", ICM[i]); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\t%.4e", IB[i]); }
		for (i = 0; i < n_cats; i++) { fprintf(endpoint_data, "\t%.4e", ID[i]); }
		fprintf(endpoint_data, "\n");
		fclose(endpoint_data);
	}

	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	//finish2:
	if (endpoint_data != NULL) { fclose(endpoint_data); }
	Rcout << "\nSteady-state calculations complete.\n";

	// Return list
	return List::create(Named("mv_values") = mv_values);
}