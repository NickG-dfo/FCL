#include <math.h> //this is for std::pow, std::exp, and std::log
//#include <stdio.h> //this is for std::printf
#include <iostream>
#include <string> //this is for std::string in 'str_is'

// #include <Rcpp.h> //not needed if including Armadillo
#include <RcppArmadillo.h> //this is for use of Rcpp with mvnorm
#include <mvnorm.h> // this is a sub-header within Arma
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
//This ^ is required for Arma, and is essential for rmvnorm
using namespace Rcpp;

#include <SimSource.hpp>

//Save MP_blank.cpp as a separate file, then source it through R//
//This will create a function called 'Simulate' in R//

namespace PopSim{

	// MP Function //
	
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
		int Ages = SSB.cols(),
			Years = SSB.rows(),
			Surveys = NIndex.rows();
		
		T output = 0.;
					
		if(n == 0){
			output = sum( SSB(SSB.nrow()-1, _) ) * 0.1;
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

