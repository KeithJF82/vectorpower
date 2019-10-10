
#include <Rcpp.h>
#include <string>

using namespace std;

//------------------------------------------------
Rcpp::List dummy1_cpp(Rcpp::List args);
double arraysum(double* V, int dim1);
int intdiv(double x, double y);
int rcpp_to_int(SEXP x);
double rcpp_to_double(SEXP x);
std::string rcpp_to_string(SEXP x);
double runif1(double a = 0.0, double b = 1.0);
