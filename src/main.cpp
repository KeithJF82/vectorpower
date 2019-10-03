#include "main.h"
#include <string>

using namespace std;

//------------------------------------------------
// Dummy function to test Rcpp working as expected
// [[Rcpp::export]]
Rcpp::List dummy1_cpp(Rcpp::List args) {
  
  // print message to console
  Rcpp::Rcout << "running C++ dummy1_cpp function\n";
  
  // get inputs from Rcpp format to base C++ format
  vector<double> x = Rcpp::as<vector<double>>(args("x"));
  
  // square values
  for (int i=0; i<int(x.size()); i++) {
    x[i] *= x[i];
  }
  
  // return as Rcpp list
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("x_squared") = x);
  return ret;
}

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

double randgen(int n)
{
	int m;
	double value = 0.0;
	for (m = 1; m <= n; m++) { value += (rand() % 10) * pow(10.0, -m); }
	return value;
}
