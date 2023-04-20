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

// namespace PopSim{

	// // template<class T>
	// // // HCR<T>::F_MP(T SSB, String name){ return; }

	// // template<class T>
	// // HCR<T>::SSB_MP(T SSB, String name){ return; }

	// // template<class T>
	// // HCR<T>::Survey_MP(NumericVector survey_indices, String name){ return; }

	// // template<class T>
	// // HCR<T>::F_MP(NumericVector &item, String name, String call_item = "SSB")

	// // Custom MP Function //
	
	// template<class T>
	// SR::CustomMP = 

	// // Custom Recruitment Function //

	// template<class T>
	// SR::CustomRecruitment(T SSB, const NumericVector &parms){ return 0.; }

// };


// R Simulation Function //

// [[Rcpp::export]]
List Simulate(  List OM, 
				List MP, 
				int sim_no, 
				int y_sim,
				String SimType = "Age-based") {
					
	//Surplus Production Model
	if(	str_is( SimType, "SPM" ) ){	
	
		return SurplusModel<double> ( OM, 
									  MP, 
									  sim_no, 
									  y_sim );
								
	}else if( str_is( SimType, "Age-based" ) ){
		
		return AgeModel<double> ( OM, 
								  MP, 
								  sim_no, 
								  y_sim );
		
	}						
								
}