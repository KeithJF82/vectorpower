
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
vector<int> rcpp_to_vector_int(SEXP x);
vector<double> rcpp_to_vector_double(SEXP x);
vector<vector<double>> rcpp_to_matrix_double(Rcpp::List x);
double runif1(double a = 0.0, double b = 1.0);
