#include "main.h"
#include <string>

using namespace std;

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double arraysum(double* V, int dim1)// Outputs sum of first dim1 entries in array
{
	int i;
	double sum = 0.0;
	for (i = 0; i < dim1; i++) { sum += V[i]; }
	return sum;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int intdiv(double x, double y)// Outputs product of two doubles, rounding down, as an integer
{
	int i;
	double z = x / y;
	i = z - fmod(z, 1.0);
	return i;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// converts input from List format to int format.
int rcpp_to_int(SEXP x) { return Rcpp::as<int>(x); }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// converts input from List format to double format.
double rcpp_to_double(SEXP x) { return Rcpp::as<double>(x); }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// converts input from List format to string format.
string rcpp_to_string(SEXP x) { return Rcpp::as<string>(x); }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// converts input from Rcpp::List format to vector<int> format.
vector<int> rcpp_to_vector_int(SEXP x) { return Rcpp::as<vector<int> >(x); }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// converts input from Rcpp::List format to vector<double> format.
vector<double> rcpp_to_vector_double(SEXP x) { return Rcpp::as<vector<double> >(x); }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// converts input from Rcpp::List format to vector<vector<double>> format.
vector<vector<double>> rcpp_to_matrix_double(Rcpp::List x) {
	int nrow = int(x.size());
	vector< vector<double> > x_mat(nrow);
	for (int i = 0; i < nrow; i++) {
		x_mat[i] = Rcpp::as<vector<double> >(x[i]);
	}
	return x_mat;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double runif1(double a, double b) {
  return R::runif(a,b);
}

