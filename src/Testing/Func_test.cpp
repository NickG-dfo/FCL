//Function Call Test//
#include <Rcpp.h>

Rcpp::NumericVector cppF_using_RF(Rcpp::NumericVector X){
	Rcpp::Function tF("temp_func");
	return tF(X);
}

template<typename Type> Type MPtest(Type x)
{
	Rcpp::Function MPFunc("MP_f");
	auto value = MPFunc( Rcpp::wrap(x) );
	return Rcpp::as<Type>( value );
}

// [[Rcpp::export]]
double retFunc(double y)
{
	return MPtest<double>( y );
}