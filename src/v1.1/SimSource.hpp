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
											  // MP, 
											  sim_no, 
											  y_sim );
								
	}else if( str_is( SimType, "Age-based" ) ){
		
		return PopSim::Simulation<double>().AgeModel(OM,
												     MP, 
												     sim_no,
												     y_sim);
		
	}				

	std::cout << "'SimType' does not match any available value in Simulate. Returning NULL.";
	return List::create(R_NilValue);
							
}