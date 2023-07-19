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
#include <Mortality.cpp> //class object for natural mortality
#include <Ncalc.cpp> //objects and functions for abundances, indices, lengths, and data simulation
#include <Recruitment.cpp> //objects and functions for recruitment
#include <HCR.cpp> //objects and classes for MPs
// #include <Rec_SPRf.cpp> //initializes Rec functions that require SPR
#include <Fishing.cpp> //object for fishing mortality


// Age-based Simulation Function //

template<class T>
List AgeModel(  List OM, 
				List MP, 
				int sim_no, 
				int y_sim
				){
					
	using namespace PopSim;

	// Call OM and MP elements //
	
	int r_lag = OM["SR_lag"], // Age of 'recruitment' (or minimum in age model, against which an SRR is modeled)
		y0 = MP["y0"], // Number of years with pre-defined TAC/F
		ymp = MP["delay"]; // this is delay on MP implementation -- either 0 or 1 year
	
	T terminalYield = OM["terminalYield"], // Yield from last year of assessment (or most recent Yield)
	  R_std = OM["Rstd"], // Constant uncertainty in recruitment
	  
	  Cmin = MP["Min_Catch"], // Minimum allowable catch for all MPs
	  TAC_error = MP["TAC_Error"]; // Allows discrpepancy between TAC and Yield (implementation error)
	  
	String OM_name = OM["Model_Name"],  //These may be unnecessary!
		   M_name = OM["Mortality_Name"],
 		   R_name = OM["Recruitment_Name"],
		   MP_name = MP["Name"];
		   
	NumericVector Catch_wt = OM["Catch_wt"],
				  Stock_wt = OM["Stock_wt"],
				  maturity = OM["Maturity"],
				  
				  N1 = OM["N1"], // First year abundance, for projecting
				  Rec_parms = OM["Rparameters"], // Parameters needed for SRR
				  
				  F_mu = OM["F_Mu"], // Means for error in F
				  F_sel = OM["Selectivity"], // Selectivity vector
				  F_terminal = OM["terminalF"], // F-at-age for the last year of the assessment
				  
				  M_terminal = OM["terminalM"], // M-at-age for the last year of the assessment
				  // M_est = OM["Mpar_est"], // change these
				  // M_std = OM["Mpar_sd"], // 		these
				  
				  PE_Std = OM["PE_Std"], // Constant uncertainty in process error for each age
				  
				  TAC0 = MP["TAC0"]; // Fixed TAC for start of simulation (this should be the same size as y0)
				  
	IntegerVector //M_ind = OM["M_index"], // CHECK WHAT THIS IS FOR
				  // Findex = OM["keyF"],
				  Fbar_ind = OM["Fbar_ind"], // Ages (indices) over which to average popwt F; vector of length 2
				  Mbar_ind = OM["Mbar_ind"]; // Ages (indices) over which to average popwt M; vector of length 2
				  // list_of_MPs = OM["MP_list"]; // This should be a vector of Strings named with MP names
				  
	NumericMatrix Rec_CoV = OM["RCoV"], // (Optional) Covariance matrix for SRR parameters
				  Rec_Kernel = OM["RKernel"], // (Optional) Kernel matrix for SRR parameters
				  N_rec = OM["Nrec"]; // This is a matrix of abundance for recruitment -- nrow() depending on SR lag-time
				  
	List Minfo = OM["Minfo"]; // This is just for input into the M object, not for local calls
	
	T max_rec = OM["Max_Rec"];
	String Rparm_type = OM["Rparm_type"];
	
	// Ftype F_std = OM["F_sigma"]; 
	//if F_type works properly, this should allow F_sigma to be either a matrix OR vector
	F_type F_std = OM["F_sigma"]; //Uncertainty matrix/vector for error in F
	
	RFunType Custom_SR = OM["SR_function"];
	List R_info = List::create(Named("Name") = R_name, Named("Function") = Custom_SR);
	
	//			-------------------			  //
			
	// Swtiches //		 
	LogicalVector R_settings = OM["Rec_Settings"],
				  F_settings = OM["F_settings"],
				  M_settings = OM["Mort_Settings"];
	bool Rparm_error = R_settings["Rparm_error"],
		 R_error = R_settings["R_error"],
		 bias_correction = R_settings["bias_correction"],
		 RKseed = R_settings["Kernel_seed"]; // Randomize Kernel index for fixed parm recruitment (only works for use_kernel)
	bool F_error = F_settings["F_error"],
		 F_exp = F_settings["F_experror"];
	bool M_cond_on = M_settings["condition"],
		 M_age_on = M_settings["age"],
		 M_year_on = M_settings["year"],
		 M_AR_on = M_settings["AR"];
	bool MP_Type = MP["Type"];		
	
	const int A = OM["A"]; // A isn't needed but makes some implementation easier (length of Age vector, not max age)
	int y, a, sim; //for loop values
	
	IntegerVector years = OM["years"]; //all years in assessment (may remove and just read last year)
	int years1 = years(years.size()-1); //last year of assessment = first year of the MSE
	
	
	// Create objects for sim //	
	
	T ssby, LogN,
	  logrec, std_log_r, // Recruitment containers
	  Fbar_top, Fbar_bot, // Fbar_popwt
	  Mbar_top, Mbar_bot; // Mbar_popwt
	
	NumericVector Nstock(A); //Just for ssby calculations
				  // Merr(A), Ferr(A); //vectors for ages
	
	 //vectors for years
	NumericVector TAC(y_sim), Fbar(y_sim), 
				  Fbar_popwt(y_sim), Mbar_popwt(y_sim), 
				  Ydiff(y_sim), Yield(y_sim), 
				  Abundance(y_sim), Biomass(y_sim), SSB(y_sim);
				  
	//matrices for years and ages
	NumericMatrix //logNage(y_sim, A), Zage(y_sim, A),
				  Cage(y_sim, A), Yage(y_sim, A), 
				  Bage(y_sim, A), SSBage(y_sim, A), 
				  Nage(y_sim, A), PE(y_sim, A);
								
	//matrices for years and sims
	NumericMatrix all_TAC(sim_no, y_sim), all_Fbar(sim_no, y_sim), //all_Mbar(sim_no, y_sim), 
				  all_Fbar_popwt(sim_no, y_sim), all_Mbar_popwt(sim_no, y_sim),
				  all_Ydiff(sim_no, y_sim), all_Yield(sim_no, y_sim),
				  all_Biomass(sim_no, y_sim), all_SSB(sim_no, y_sim),
				  all_Rec(sim_no, y_sim), all_Abundance(sim_no, y_sim);
	
	NumericMatrix3D YaA(sim_no, y_sim, A),
					BaA(sim_no, y_sim, A),
					SaA(sim_no, y_sim, A);


	// -- Create class objects -- //
		
	//Natural Mortality Object//
	M_obj<T> Mya_obj( Minfo, y_sim );
	
	//List for Rec_obj & RP_obj
	List LHC_list = List::create(Named("A") = A, //Max age
								 Named("Mort") = M_terminal, //Assumes M & Sel are constant
								 Named("Sel") = F_sel,
								 Named("Mat") = maturity,
								 Named("Sw") = Stock_wt,
								 Named("Cw") = Catch_wt
								);
								
	//Recrutiment Stuff//
	SR<T> Rec_obj(R_info, R_std, 
				  LHC_list, //Just for SPR0
				  y_sim, max_rec, 
				  bias_correction, RKseed,
				  Rparm_error, R_error);
	
	if( str_is(Rparm_type, "use_CoV") ){
		Rec_obj.set_Parms(Rec_parms);
		Rec_obj.set_CoV(Rec_CoV);
	}else if( str_is(Rparm_type, "use_Kernel") ){
		Rec_obj.set_Kernel(Rec_Kernel);
	}else if( str_is(Rparm_type, "find_CoV") ){
		Rec_obj.find_CoV(Rec_Kernel);
	}
	
	RP<T> RP_obj( LHC_list, Rec_obj );	//Object for Reference Points	
	
	//MP Stuff//
	HCR<T> MP_obj; //Object for HCR outputs/functions
	Rcpp::RObject MP_inputs = MP["Function"]
	
	//Fishing Mortality Stuff//
	F_obj<T> Fya_obj( F_mu, F_std, F_error, F_experror ); //Hopefully this reads F_type for F_std
	Fya_obj.Sel = F_sel; //Set selectivity
	
		
	// -- Start simulation -- //
	
	for(sim = 0; sim < sim_no; sim++){
				
		std::printf("Thread: %d / %d \r", sim+1, sim_no);
		
		// set as zero each iter
		Yield.fill(0.);
		Ydiff.fill(0.);
		Fbar_popwt.fill(0.);
		Mbar_popwt.fill(0.);
		Abundance.fill(0.);
		Biomass.fill(0.);
		SSB.fill(0.);
		
		// avgC_short_val = 0.;
		// avgC_long_val = 0.;
		// aav_val = 0.;
		
		Fbar_top = 0.;
		Fbar_bot = 0.;
		Mbar_top = 0.;
		Mbar_bot = 0.;
		
		y = 0;
					
		//first year of simulation is same as terminal year from assessment	
		PE(y, _) = Rcpp::rnorm(A, rep(0., A), PE_Std);
		logN = std::log(N1(a)) + PE(y, _); //Start using N1 with PE
		Nage(y, _) = Rcpp::exp(logN_a);
				
		Bage(y, _) = Nage(y, _) * Stock_wt;
		SSBage(y, _) = Bage(y, _) * maturity;

		Fya_obj.Fya(y, _) = Fterminal;
		Mya_obj.Mya(y, _) = Mterminal;
		Zage(y, _) = Fterminal + Mterminal;
		Cage(y, _) = baranov_catch(Nage(y, _), Fterminal, Mterminal); // Nage(y, _) * (1. - Rcpp::exp(-Zage(y, _))) * Fterminal / Zage(y, _);
		Yage(y, _) = Cage(y, _) * Catch_wt; // Cage(y, _) * Catch_wt(_);
		
		for(a = 0; a < A; a++){
			Fbar_top += ( (a > Fbar_ind(0)) & (a < Fbar_ind(1)) ) * Nage(y, a) * Fya_obj.Fya(y, a);
			Fbar_bot += ( (a > Fbar_ind(0)) & (a < Fbar_ind(1)) ) * Nage(y, a);
			Mbar_top += ( (a > Mbar_ind(0)) & (a < Mbar_ind(1)) ) * Nage(y, a) * Mya_obj.Mya(y, a);
			Mbar_bot += ( (a > Mbar_ind(0)) & (a < Mbar_ind(1)) ) * Nage(y, a);
		}
		
		Ydiff(y) = sum(Yage(y, _)) - terminalYield;
		Yield(y) = sum(Yage(y, _));
		TAC(y) = TAC0(y);
		Fbar(y) = solveF(TAC(y), Fya_obj.F_y(y), Nage(y, _), Mort_obj.M_ya(y, _), Catch_wt);
		//Fya_obj( Fbar(y) );
		Abundance(y) = sum(Nage(y, _));
		Biomass(y) = sum(Bage(y, _));
		SSB(y) = sum(SSBage(y, _));
		
		//3D matrices with age components, e.g. for yield5-8
		YaA(a, sim, y) = Yage(y, a);
		BaA(a, sim, y) = SSBage(y, a);
		SaA(a, sim, y) = Bage(y, a);
		
		Fbar_popwt(y) = Fbar_top / Fbar_bot;
		Mbar_popwt(y) = Mbar_top / Mbar_bot;
		
		//start of projection at year 2022
		for(y = 1; y < y_sim; y++){
			
			if(y < r_lag){
				Nstock = N_rec(y, _) +  Rcpp::rnorm(A, rep(0., A), PE_Std);
			}else{
				Nstock = Nage(y-r_lag, _);
			}
			ssby = sum( Nstock * Stock_wt * maturity );
			
			// Recruitment
			logrec = Rec_obj.logRec(ssby);

			// Abundance
			PE(y, _) = Rcpp::rnorm(A, rep(0., A), PE_Std); 
			logN = logN_a( logNage(y-1, _), 
							Mort_obj.M_y(y-1) + Fya_obj.F_y(y-1), 
							PE(y, _), 
							logrec );
			Nage(y, _) = Rcpp::exp(logN(y, _));
    		
			// Natural Mortality
			if( M_cond_on ){
				Mya_obj.MCondition(y);
			}else if( M_age_on ){
				Mya_obj.MAgeEffect(y);
			}else if( M_year_on ){
				Mya_obj.MYearEffect(y);
			}else if( M_AR_on ){
				Mya_obj.MAR(y);
			}
						
			if(y < y0){ //Fixed TAC/F years
			
				TAC(y) = TAC0(y);
				Fbar(y) = solveF( TAC(y), Fya_obj.Fya(y, _), Nage(y, _), Mya_obj.Mya(y, _), Catch_wt );

			}else{ //Years for Variable MPs	
			
				MP_obj.Calc(MP_inputs);
				Fbar(y) = MP_obj.F;
				TAC(y) = MP_obj.TAC;

				// MP_obj.F_MP( SSB(y - ymp), list_of_MPs.offset( MP_name ) );
				// MP_obj.TAC_MP( SSB(y - ymp), list_of_MPs.offset( MP_name ) );
				
				// if( !MP_obj.F ){ // if F=0, therefore TAC MP
					// MP_obj.F = solveF( MP_obj.TAC, Fya_obj.F_y(y), Nage(y, _), Mort_obj.M_y(y), Catch_wt );
				// }else if( !MP_obj.TAC ){ // if TAC=0, therefore F MP
					// MP_obj.TAC = yield<T>( baranov_catch( Nage(y, _), Fya_obj.F_y(y), Mort_obj.M_y(y) ), Catch_wt );
				// }
				
				// TAC(y) = MP_obj.TAC + R::rnorm(0., TAC_error);
				// Fbar(y) = MP_obj.F;
				
			}	
			// If TAC is too low, set minimum to Cmin
			if( !str_is(mp, "No fishing") | (TAC(y) < Cmin) ) {
				TAC(y) = Cmin;
				Fbar(y) = solveF( TAC(y), Fya_obj.Fya(y, _), Nage(y, _), Mya_obj.Mya(y, _), Catch_wt );
			}			
						
			// Fishing Mortality
			Fya_obj(Fbar(y));	

			// Fill containers
			Bage(y, _) = Nage(y, _) * Stock_wt;
			SSBage(y, _) = Bage(y, _) * maturity;
			
			// Fage(y, _) = Fbar(y) * Fya_obj.Sel;
			Zage(y, _) = Fya_obj.Fya(y, _) + Mya_obj.Mya(y, _);
			Cage(y, _) = baranov_catch(Nage(y, _), Fya_obj.Fya(y, _), Mya_obj.Mya(y, _)) // Nage(y, _) * (1. - std::exp(-Zage(y, _))) * Fage(y, _) / Zage(y, a);
			Yage(y, _) = Cage(y, _) * Catch_wt; // Cage(y, a) * Catch_wt(a);			
						
			Fbar_top = 0.;
			Fbar_bot = 0.;
			Mbar_top = 0.;
			Mbar_bot = 0.;
			
			for(a = 0; a < A; a++){
				
				Fbar_top += ( (a > Fbar_ind(0)) & (a < Fbar_ind(1)) ) * Nage(y, a) * Fage(y, a);
				Fbar_bot += ( (a > Fbar_ind(0)) & (a < Fbar_ind(1)) ) * Nage(y, a);
				Mbar_top += ( (a > Mbar_ind(0)) & (a < Mbar_ind(1)) ) * Nage(y, a) * Mage(y, a);
				Mbar_bot += ( (a > Mbar_ind(0)) & (a < Mbar_ind(1)) ) * Nage(y, a);
				
			}
								
			Abundance(y) = sum(Nage(y, _));
			Ydiff(y) = sum(Yage(y, _) - Yage(y-1, _));
			Yield(y) = sumYage(y, _));
			Biomass(y) = sum(Bage(y, _));
			SSB(y) = sum(SSBage(y, _));
			
			//3D matrices with age components, e.g. for yield3-8
			YaA(sim, y, a) = Yage(y, a);
			BaA(sim, y, a) = SSBage(y, a);
			SaA(sim, y, a) = Bage(y, a);
			
			Fbar_popwt(y) = Fbar_top / Fbar_bot;
			Mbar_popwt(y) = Mbar_top / Mbar_bot;
			
		} // end projection
		
		//Annual values
		all_TAC(sim, _) = clone(TAC);
		all_Fbar(sim, _) = clone(Fbar);
		all_Fbar_popwt(sim, _) = clone(Fbar_popwt);
		all_Mbar_popwt(sim, _) = clone(Mbar_popwt);
		all_Yield(sim, _) = clone(Yield);
		all_Ydiff(sim, _) = clone(Ydiff);
		all_Abundance(sim, _) = clone(Abundance);
		all_Biomass(sim, _) = clone(Biomass);
		all_SSB(sim, _) = clone(SSB);
		all_Rec(sim, _) = clone( Rcpp::exp(Rec_obj.logRy) );
		
	} // end simulation
	
	return  List::create(
			Named("TAC") = all_TAC,
			Named("Fbar") = all_Fbar,
			Named("Fbar_popwt") = all_Fbar_popwt,
			Named("Yield") = all_Yield,
			Named("Abundance") = all_Abundance,
			Named("Biomass") = all_Biomass,
			Named("SSB") = all_SSB,
			Named("Rec") = all_Rec,
			
			//These will need fixing to allow output in R
			Named("YaA") = YaA,
			Named("BaA") = BaA,
			Named("SaA") = SaA
			);
	
				
} // end function




// Surplus Production Model Simulation Function //

List SPM( List OM, 
		  List MP, 
		  int sim_no, 
		  int y_sim ){
			  
	using namespace SPMSim;
	
	NumericMatrix Abundance(sim_no, y_sim),
				  Biomass(sim_no, y_sim);
	
	int sim, y;
	
	for(sim = 0; sim < sim_no; sim++){
		for(y = 0; y < y_sim; y++){
		
			;
		
		}
	}
			
	return List::Create(
		   Named("X") = 0.		   
		   );
			
}




