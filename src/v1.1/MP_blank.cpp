#include <math.h> //this is for std::pow, std::exp, and std::log
#include <stdio.h> //this is for std::printf
#include <string> //this is for std::string in 'str_is'

// #include <Rcpp.h> //not needed if including Armadillo
#include <RcppArmadillo.h> //this is for use of Rcpp with mvnorm
#include <mvnorm.h> // this is a sub-header within Arma
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
//This ^ is required for Arma, and is essential for rmvnorm

using Rcpp::_;
using Rcpp::max;
using Rcpp::sum;

using Rcpp::Named;
using Rcpp::String;
using Rcpp::NumericVector;
using Rcpp::IntegerVector;
using Rcpp::StringVector;
using Rcpp::LogicalVector;
using Rcpp::NumericMatrix;
using Rcpp::List;
using Rcpp::DataFrame;


#include <mseliteutils.cpp> //simple functions for conversions and containers
using namespace mseutils;


#include <Mortality.cpp> //class object for natural mortality
#include <Ncalc.cpp> //objects and functions for abundances, indices, lengths, and data simulation
#include <Recruitment.cpp> //objects and functions for recruitment
#include <HCR.cpp> //objects and classes for MPs
// #include <Rec_SPRf.cpp> //initializes Rec functions that require SPR
#include <Fishing.cpp> //object for fishing mortality

#include <simulation.cpp>

/* 
								Instructions

Provided are blank templates for Management Procedures. Any numer of MPs may be
defined within these functions. Either or both provided MP types may be defined. 
A secondary input (i.e. 'name') is provided for each function, and decides which 
MP is active for a given String input.

Ex.			if(name == "MP1") { value = 0.5 * SSB; }
			else if(name == "MP2") { value = 0.4 * SSB; } etc.
			
Values for 'name' may shift between MP types if multiple MP types are used within 
a model. For example, TAC_MP may have 'name == "MP1"' and 'name == "MP3"', while F_MP
may have 'name == "MP2"'.

If only one MP is provided, no 'if' statement or 'name' value is required (although the 
argument 'String name' is still required in the function definition).

Each MP type is a void function, and thus does not return a value. Values instead are defined 
within their internal object using local containers F & TAC. These are accessed as this->F and
this->TAC, respectively.

Ex.			if(name == "MP5"){ this->TAC = (SSB - 20) * 0.1); }
			if(name == "MP6"){ this->F = (SSB/Blim * 0.1); }
			
##Both TAC & F can be defined within a single MP using a baranov catch or inverse baranov catch.

Regardless of the function blocks, each MP function must end with 'return;'

This should be Saved As a separate .cpp file and called using Rcpp::sourcecpp( /filename/ ) in R.
Doing so will create a function in R called Simulate, which will take OM and MP inputs, along with
number of simulations and number of years, to run the simulation. By default this function will run
a Age-based model, but it can also run a SPM.

*/

namespace PopSim{

	// template<class T>
	// // HCR<T>::F_MP(T SSB, String name){ return; }

	// template<class T>
	// HCR<T>::SSB_MP(T SSB, String name){ return; }

	// template<class T>
	// HCR<T>::Survey_MP(NumericVector survey_indices, String name){ return; }

	// template<class T>
	// HCR<T>::F_MP(NumericVector &item, String name, String call_item = "SSB")

	// Custom MP Function //
	
	template<class T>
	T MP_decision(T SSB, NumericVector parms){ return 0.; }

	// Custom Recruitment Function //

	template<class T>
	T SR<T>::CustomRecruitment(T SSB, NumericVector &parms){ return 0.; }

}


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
		
		return PopSim::AgeModel<double> ( OM, 
										  MP, 
										  sim_no, 
										  y_sim );
		
	}				

	return List::create(R_NilValue);
							
}