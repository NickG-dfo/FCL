//Function Call Test//
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericVector cppF_using_RF(Rcpp::NumericVector X){
	Rcpp::Function tF("temp_func");
	return tF(X);
}