#include <math.h> //this is for std::pow, std::exp, and std::log
#include <stdio.h> //this is for std::printf
#include <string> //these are for std::string in 'str_is'
//#include <cstring>

// #include <Rcpp.h> //not needed if including Armadillo
#include <RcppArmadillo.h> //this is for use of Rcpp with mvnorm
#include <mvnorm.h> // this is a sub-header within Arma
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
//This ^ is required for Arma, and is essential for rmvnorm
using namespace Rcpp;

#include <FCLutils.hpp> //simple functions for conversions and containers
using namespace FCLutils;

#include <Mortality.hpp> //class object for natural mortality
#include <Ncalc.hpp> //objects and functions for abundances, indices, lengths, and data simulation
#include <Recruitment.hpp> //objects and functions for recruitment
#include <HCR.hpp> //objects and classes for MPs
#include <Fishing.hpp> //object for fishing mortality
#include <Simulation.hpp>

// R Simulation Function //

// [[Rcpp::export]]
List Simulate(  List OM, 
				List MP, 
				int sim_no, 
				int y_sim,
				String SimType = "Age-based") {
					
	//Surplus Production Model
	if(	str_is( SimType, "SPM" ) ){	
	
		return SPMSim::SurplusModel<double> ( OM, 
											  MP, 
											  sim_no, 
											  y_sim );
								
	}else if( str_is( SimType, "Age-based" ) ){
		
		return PopSim::Simulation<double>().AgeModel(OM,
												     MP, 
												     sim_no,
												     y_sim);
		
	}				

	return List::create(R_NilValue);
							
}