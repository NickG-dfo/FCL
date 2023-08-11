#include <math.h> //this is for std::pow, std::exp, and std::log
#include <stdio.h> //this is for std::printf
#include <string> //these are for std::string in 'str_is'
#include <cstring>

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
#include <Fishing.cpp> //object for fishing mortality

#include <simulation.cpp>


namespace PopSim{

	// Custom MP Function //
	
	template<class T>
	T Simulation<T>::MP_decision(Simulation<T>* psim, int n)
	{ 
		// Possible calls within MP_decision are SSB, Biomass, Biomass index, and Abundance index
		// Calls are all matrix of dimensions Years x Age, i.e. Matrix(Y, A);
		// argument 'n' allows multiple possible MPs in one model, 
		//		and can be changed via OM$MP_n in R
		
		NumericMatrix SSB = psim -> MP_ssb,
					  Biomass = psim -> MP_biomass,
					  BIndex = psim -> MP_bindex,
					  NIndex = psim -> MP_nindex;
		int A = SSB.cols(),
			Y = SSB.rows();
		
		T output = 0.;
					
		if(n == 0){
			output = sum(SSB(SSB.rows(), _)) * 0.1;
		}
		
		return output; 
	}

	// Custom Recruitment Function //

	template<class T>
	T SR<T>::CustomRecruitment(T SSB, NumericVector &parms)
	{
		T Rec = SSB * parms(1) + parms(0);		
		return Rec; 
	}

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
		
		return PopSim::Simulation<double>().AgeModel(OM,
												   MP, 
												   sim_no,
												   y_sim);
		
	}				

	return List::create(R_NilValue);
							
}