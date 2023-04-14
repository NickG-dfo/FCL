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
	
	T std_pe = OM["PE_std"], // Constant uncertainty in process error for all ages and years (maybe allow variability in some dimension?)
	  terminalYield = OM["terminalYield"], // Yield from last year of assessment (or most recent Yield)
	  // R_std = OM["Rstd"], // Constant uncertainty in recruitment
	  
	  Cmin = MP["Min_Catch"], // Minimum allowable catch for all MPs
	  TAC_error = MP["TAC_Error"]; // Allows discrpepancy between TAC and Yield (implementation error)
	  
	  // SSB_decline = OM["SSB_decline"], // Threshold to calculate risk/decline PM
	  // TAC_decline = OM["TAC_decline"]; // Threshold to calculate prob. of achieving TAC
	  
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
				  M_est = OM["Mpar_est"], // change these
				  M_std = OM["Mpar_sd"], // 		these
				  
				  Std_PE = OM["PE_Std"], // Constant uncertainty in process error for each age
				  
				  TAC0 = MP["TAC0"]; // Fixed TAC for start of simulation (this should be the same size as y0)
				  
	IntegerVector //M_ind = OM["M_index"], // CHECK WHAT THIS IS FOR
				  // Findex = OM["keyF"],
				  Fbar_ind = OM["Fbar_ind"], // Ages (indices) over which to average popwt F; vector of length 2
				  Mbar_ind = OM["Mbar_ind"], // Ages (indices) over which to average popwt M; vector of length 2
				  // list_of_MPs = OM["MP_list"]; // This should be a vector of Strings named with MP names
				  
	NumericMatrix //Fmat = OM["F_err"],  // CHECK WHAT THIS IS FOR
				  Rec_CoV = OM["RCoV"], // (Optional) Covariance matrix for SRR parameters
				  Rec_Kernel = OM["RKernel"], // (Optional) Kernel matrix for SRR parameters
				  N_rec = OM["Nrec"]; // This is a matrix of abundance for recruitment -- nrow() depending on SR lag-time
				  
	List Minfo = OM["Minfo"]; // This is just for input into the M object, not for local calls
	
	// Ftype F_std = OM["F_sigma"]; 
	//if F_type works properly, this should allow F_sigma to be either a matrix OR vector
	F_type F_std = OM["F_sigma"]; //Uncertainty matrix/vector for error in F
			
			
	// Swtiches //		 
	LogicalVector R_settings = OM["Rec_Settings"],
				  F_settings = OM["F_settings"],
				  M_settings = OM["Mort_Settings"];
	T max_rec = OM["Max_Rec"];
	String Rparm_type = OM["Rparm_type"];
	bool Rparm_error = R_settings["Rparm_error"],
		 R_error = R_settings["R_error"],
		 bias_correction = R_settings["bias_correction"],
		 RKseed = R_settings["Rec_Kernel_seed"]; // Randomize Kernel index for fixed parm recruitment (only works for use_kernel)
	bool F_error = F_settings["F_error"],
		 F_exp = F_settings["F_experror"];
	bool M_cond_on = M_settings["condition"],
		 M_age_on = M_settings["age"],
		 M_year_on = M_settings["year"],
		 M_AR_on = M_settings["AR"];
		
	
	const int A = OM["A"]; // A isn't needed but makes some implementation easier (length of Age vector, not max age)
	int y, a, sim;
	
	IntegerVector years = OM["years"]; //all years in assessment (may remove and just read last year)
	int years1 = years(years.size()-1); //last year of assessment = first year of the MSE
	
	
	// Create objects for sim //	
	
	T ssby, logrec, std_log_r, // Recruitment containers
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
	NumericMatrix logNage(y_sim, A), Zage(y_sim, A),
				  Cage(y_sim, A), Yage(y_sim, A), 
				  Bage(y_sim, A), SSBage(y_sim, A), 
				  Nage(y_sim, A), PE(y_sim, A);
								
	//matrices for years and sims
	NumericMatrix all_TAC(sim_no, y_sim), all_Fbar(sim_no, y_sim), //all_Mbar(sim_no, y_sim), 
				  all_Fbar_popwt(sim_no, y_sim), all_Mbar_popwt(sim_no, y_sim),
				  all_Ydiff(sim_no, y_sim), all_Yield(sim_no, y_sim),
				  all_Biomass(sim_no, y_sim), all_SSB(sim_no, y_sim),
				  all_Rec(sim_no, y_sim), all_Abundance(sim_no, y_sim);
	
	NumericMatrix3D YaA(A, sim_no, y_sim),
					BaA(A, sim_no, y_sim),
					SaA(A, sim_no, y_sim);
	
	// Performance Metrics //
	
	//timeframes
	// int short_time, //short-term timeframe
		// medium_time, //medium-term timeframe
		// long_time, //long-term timeframe
		// aav_time; //time-frame to calculate AAV

	// //annual values
	// T avgC_short_val = 0., //Short-term average catches
	  // avgC_medium_val = 0., //Medium-term average catches
	  // avgC_long_val = 0., //Long-term average catches
	  // medC_short_val = 0., //Short-term median catches
	  // medC_medium_val = 0., //Medium-term median catches
	  // medC_long_val = 0., //Long-term median catches
	  // aav_val = 0., // Average (inter)annual variability in catches
	  // aav_sum = 0.; // Sum of catches over aav_time for AAV
	  
	// //simulation values
	// NumericVector aav_all(sim_no), 
				  // medC_long_all(sim_no),
				  // medC_medium_all(sim_no),
				  // medC_short_all(sim_no),
				  // avgC_long_all(sim_no),
				  // avgC_medium_all(sim_no),
				  // avgC_short_all(sim_no);
	
	// //Probability-by-year values
	// NumericVector Prob_Decline(y_sim), //Probability of SSB being below pre-determined value (also called risk)
				  // Prob_lowTAC(y_sim), //Probability of TAC being below pre-determined value
				  // Prob_underTAC(y_sim); //Probability of Yield being below TAC (Only for non-zero implementation error)


		
	//Natural Mortality Object//
	M_obj<T> Mort_obj( Minfo, y_sim );
	
	//Recrutiment Stuff//
	SR<T> Rec_obj(R_name, R_std, 
				  y_sim, max_rec, 
				  bias_correction, RKseed);
	
	if( str_is(Rparm_type, "use_CoV") ){
		Rec_obj.set_Parms(Rec_parms);
		Rec_obj.set_CoV(Rec_CoV);
	}else if( str_is(Rparm_type, "use_Kernel") ){
		Rec_obj.set_Kernel(Rec_Kernel);
	}else if( str_is(Rparm_type, "find_CoV") ){
		Rec_obj.find_CoV(Kernel);
	}
	
	//MP Stuff//
	HCR<T> MP_obj; //Object for HCR outputs/functions
	//List for RP_obj
	List LHC_list = List::create(Named("A") = A, //Max age
								 Named("Mort") = Mterminal, //Assumes M & Sel are constant
								 Named("Sel") = F_sel,
								 Named("Mat") = maturity,
								 Named("Sw") = Stock_wt,
								 Named("Cw") = Catch_wt
								);
	RP<T> RP_obj( LHC_List, Rec_obj );	//Object for Refernce Points
	
	//Fishing Mortality Stuff//
	F_obj Fya_obj( F_mu, F_std ); //Hopefully this reads F_type for F_std
	Fya_obj.Sel = F_sel; //Set selectivity
		
		
	// Start simulation //
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
		for(a = 0; a < A; a++){
			
			PE(y, a) = (T)R::rnorm(0., std_pe); //R::rnorm returns 1 value
			logNage(y, a) = std::log(N1(a)) + PE(y, a); //Start using N1 with PE
			Nage(y, a) = std::exp( logNage(y, a) );
			
			Bage(y, a) = Nage(y, a) * Stock_wt(a);
			SSBage(y, a) = Bage(y, a) * maturity(a);

			Fage(y, a) = Fterminal(a);
			Mage(y, a) = Mterminal(a);
			Zage(y, a) = Fage(y, a) + Mage(y, a);
			Cage(y, a) = Nage(y, a) * (1. - std::exp(-Zage(y, a))) * Fage(y, a) / Zage(y, a);
			Yage(y, a) = Cage(y, a) * Catch_wt(a);
			
			Fbar_top += ( (a > Fbar_ind(0)) & (a < Fbar_ind(1)) ) * Nage(y, a) * Fage(y, a);
			Fbar_bot += ( (a > Fbar_ind(0)) & (a < Fbar_ind(1)) ) * Nage(y, a);
			Mbar_top += ( (a > Mbar_ind(0)) & (a < Mbar_ind(1)) ) * Nage(y, a) * Mage(y, a);
			Mbar_bot += ( (a > Mbar_ind(0)) & (a < Mbar_ind(1)) ) * Nage(y, a);
			
			Ydiff(y) += Yage(y, a);
			Yield(y) += Yage(y, a);
			Abundance(y) += Nage(y, a);
			Biomass(y) += Bage(y, a);
			SSB(y) += SSBage(y, a);
			
			//3D matrices with age components, e.g. for yield5-8
			YaA(a, sim, y) = Yage(y, a);
			BaA(a, sim, y) = SSBage(y, a);
			SaA(a, sim, y) = Bage(y, a);
			
		}
		
		Ydiff(y) -= terminalYield;
		TAC(y) = TAC0(y);
		Fbar(y) = solveF(TAC(y), Fya_obj.F_y(y), Nage(y, _), Mort_obj.M_ya(y, _), Catch_wt);
		Fya_obj( Fbar(y) );
		
		Fbar_popwt(y) = Fbar_top / Fbar_bot;
		Mbar_popwt(y) = Mbar_top / Mbar_bot;
		
		//start of projection at year 2022
		for(y = 1; y < y_sim; y++){
			
			if(y < r_lag){
				Nstock = N0 + Rcpp::rnorm(A, 0., std_pe);
			}else{
				Nstock = Nage(y-r_lag, _);
			}
			ssby = sum( Nstock * Stock_wt * maturity );
			
			// Recruitment
			logrec = Rec_obj.logRec(ssby, Rparm_error, R_error);

			// Abundance
			PE(y, _) = Rcpp::rnorm(A, 0., std_pe); 
			logNage(y, _) = logN_a( logNage(y-1, _), 
									Mort_obj.M_y(y-1) + Fya_obj.F_y(y-1), 
									PE(y, _), 
									logrec );
			Nage(y, _) = Rcpp::exp(logNage(y, _));
    		
			// Natural Mortality
			if( M_cond_on ){
				Mort_obj.MCondition(y);
			}else if( M_age_on ){
				Mort_obj.MAgeEffect(y);
			}else if( M_year_on ){
				Mort_obj.MYearEffect(y);
			}else if( M_AR_on ){
				Mort_obj.MAR(y);
			}
						
			// Fishing Mortality
			Fya_obj(F_error, F_experror); //Calculates error in Sel each year
						
			if(y < y0){ //Fixed TAC/F years
			
				TAC(y) = TAC0(y);
				Fbar(y) = solveF( TAC(y), Fya_obj.F_y(y), Nage(y, _), Mort_obj.M_y(y), Catch_wt );

			}else{ //Years for Variable MPs

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
			if( !str_is(mp, "No fishing") & (TAC(y) < Cmin) ) {
				TAC(y) = Cmin;
				Fbar(y) = solveF( Cmin, Fya_obj.F_y(y), Nage(y, _), Mort_obj.M_y(y), Catch_wt );
			}			
			
			Fya_obj.F( Fbar(y) );					
						
			Fbar_top = 0.;
			Fbar_bot = 0.;
			Mbar_top = 0.;
			Mbar_bot = 0.;
			
			for(a = 0; a < A; a++){
			
				Bage(y, a) = Nage(y, a) * Stock_wt(a);
				SSBage(y, a) = Bage(y, a) * maturity(a);
				
				Fage(y, a) = Fbar(y) * sel(y, a);
				Zage(y, a) = Fage(y, a) + Mage(y, a);
				Cage(y, a) = Nage(y, a) * (1. - std::exp(-Zage(y, a))) * Fage(y, a) / Zage(y, a);
				Yage(y, a) = Cage(y, a) * Catch_wt(a);
				
				Fbar_top += ( (a > Fbar_ind(0)) & (a < Fbar_ind(1)) ) * Nage(y, a) * Fage(y, a);
				Fbar_bot += ( (a > Fbar_ind(0)) & (a < Fbar_ind(1)) ) * Nage(y, a);
				Mbar_top += ( (a > Mbar_ind(0)) & (a < Mbar_ind(1)) ) * Nage(y, a) * Mage(y, a);
				Mbar_bot += ( (a > Mbar_ind(0)) & (a < Mbar_ind(1)) ) * Nage(y, a);
								
				Abundance(y) += Nage(y, a);
				Ydiff(y) += Yage(y, a) - Yage(y-1, a);
				Yield(y) += Yage(y, a);
				Biomass(y) += Bage(y, a);
				SSB(y) += SSBage(y, a);
				
				//3D matrices with age components, e.g. for yield3-8
				YaA(a, sim, y) = Yage(y, a);
				BaA(a, sim, y) = SSBage(y, a);
				SaA(a, sim, y) = Bage(y, a);
				
			}
			
			Fbar_popwt(y) = Fbar_top / Fbar_bot;
			Mbar_popwt(y) = Mbar_top / Mbar_bot;
			
			//For PMs
			// Prob_Decline(y) += SSB(y) < SSB_decline;
			// Prob_lowTAC(y) += TAC(y) < TAC_decline;
			// Prob_underTAC(y) += Yield(y) < TAC(y);
			
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
				
		//PM stuff//
		// for(y = 2; y < 27; y++){ //y = 2 is 2023
			
			// if(y < 7){
				// avgC_short_val += Yield(y);
			// }
			// avgC_long_val += Yield(y);
			// aav_val += std::abs(Yield(y) - Yield(y-1));
			
		// }
						
		
		// aav_val /= avgC_long_val;
		// avgC_long_val /= 25.;
		// avgC_short_val /= 5.;
		
		// aav_all(sim) = aav_val;
		// avgC_long_all(sim) = avgC_long_val;
		// avgC_short_all(sim) = avgC_short_val;
		
	} // end simulation
	
	//Preventable Decline
	// NumericVector PD_y = apply(PD, sum, 0);
	// Decline = Decline / (T)sim_no;
	
	return  List::create(
			Named("TAC") = all_TAC,
			Named("Fbar") = all_Fbar,
			Named("Fbar_popwt") = all_Fbar_popwt,
			Named("Yield") = all_Yield,
			Named("Abundance") = all_Abundance,
			Named("Biomass") = all_Biomass,
			Named("SSB") = all_SSB,
			Named("Rec") = all_Rec
			
			//These will need fixing to allow output in R
			Named("YaA") = YaA,
			Named("BaA") = BaA,
			Named("SaA") = SaA
			
			// Named("AAV") = aav_all,
			// Named("L_AvgCatch") = avgC_long_all,
			// Named("M_AvgCatch") = avgC_medium_all,
			// Named("S_AvgCatch") = avgC_short_all,
			// Named("L_MedCatch") = medC_long_all,
			// Named("M_MedCatch") = medC_medium_all,
			// Named("S_MedCatch") = medC_short_all,
			
			// Named("P_Decline") = Prob_Decline,
			// Named("P_lowTAC") = Prob_lowTAC,
			// Named("P_underTAC") = Prob_underTAC
			);
	
				
} // end function




// Surplus Production Model Simulation Function //

List SPM( List OM, 
		  List MP, 
		  int sim_no, 
		  int y_sim ){
			  
	using namespace SPMSim;
			
	return List::Create(
		   Named("X") = 0.		   
		   );
			
}




