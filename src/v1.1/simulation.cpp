
// Age-based Simulation Function //
namespace PopSim{

	template<class T>
	struct Simulation{
		
	NumericMatrix MP_biomass, MP_ssb,
				  MP_bindex, MP_nindex;
				  
	T MP_decision(Simulation<T>*, int);
				  
	List AgeModel(  List OM, 
					List MP, 
					int sim_no, 
					int y_sim
					){					

		// Call OM and MP elements //
		
		int r_lag = OM["SR_lag"], // Age of 'recruitment' (or minimum in age model, against which an SRR is modeled)
			MP_n = MP["n"];
		
		T terminalYield = OM["terminalYield"], // Yield from last year of assessment (or most recent Yield)
		  R_std = OM["Rstd"], // Constant uncertainty in recruitment
		  
		  Cmin = MP["Min_Catch"], // Minimum allowable catch for all MPs
		  TAC_Error = MP["TAC_Error"]; // Allows discrpepancy between TAC and Yield (implementation error)
		  
		String OM_name = OM["Model_Name"],  //These may be unnecessary!
			   M_name = OM["Mortality_Name"],
			   R_name = OM["Recruitment_Name"],
			   MP_name = MP["Name"];
			   
		NumericVector Catch_wt = OM["Catch_wt"],
					  Stock_wt = OM["Stock_wt"],
					  Maturity = OM["Maturity"],
					  
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
		
		//if F_type works properly, this should allow F_sigma to be either a matrix OR vector
		Rcpp::RObject F_temp = OM["F_sigma"];
		F_type F_std(F_temp); //Uncertainty matrix/vector for error in F
		
		// String Custom_SR = OM["SR_function"];
		// List R_info = List::create(Named("Name") = R_name, Named("Function") = Custom_SR);
		
		bool MP_Type = MP["F_MP"];
		// Rcpp::RObject MP_inputs = MP["MP_inputs"];
		// String MP_Function = MP["Function"];
		// Rcpp::Function MP_decision(MP_Function); //defining MP function
		
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
		
		int A = OM["A"]; // A isn't needed but makes some implementation easier (length of Age vector, not max age)
		int y, a, sim; //for loop values
		int y0 = TAC0.size();
		
		IntegerVector years = OM["years"]; //all years in assessment (may remove and just read last year)
		// int years1 = years(years.size()-1); //last year of assessment = first year of the MSE
		
		
		// Create objects for sim //	
		
		T ssby,
		  logrec, // std_log_r, // Recruitment containers
		  Fbar_top, Fbar_bot, // Fbar_popwt
		  Mbar_top, Mbar_bot; // Mbar_popwt
		
		NumericVector Nstock(A), //Just for ssby calculations
					  logN(A);
					  // Merr(A), Ferr(A); //vectors for ages
		
		 //vectors for years
		NumericVector TAC(y_sim), TAC_disp(y_sim),
					  Fbar(y_sim), Fbar_popwt(y_sim), 
					  Mbar_popwt(y_sim), 
					  Ydiff(y_sim), Yield(y_sim), 
					  Abundance(y_sim), Biomass(y_sim), SSB(y_sim);
					  
		//matrices for years and ages
		NumericMatrix //logNage(y_sim, A), 
					  Zage(y_sim, A),
					  Cage(y_sim, A), Yage(y_sim, A), 
					  Bage(y_sim, A), SSBage(y_sim, A), 
					  Nage(y_sim, A), PE(y_sim, A);
									
		//matrices for years and sims
		NumericMatrix all_TAC(sim_no, y_sim), all_TAC_disp(sim_no, y_sim),
					  all_Fbar(sim_no, y_sim), //all_Mbar(sim_no, y_sim), 
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
									 Named("Mat") = Maturity,
									 Named("Sw") = Stock_wt,
									 Named("Cw") = Catch_wt
									);
									
		//Recrutiment Stuff//
		SR<T> Rec_obj(R_name,
					  LHC_list, //Just for SPR0
					  y_sim, R_std, max_rec, 
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
		
		//Fishing Mortality Stuff//
		F_obj<T> Fya_obj( F_mu, F_std, F_error, F_exp ); //Hopefully this reads F_type for F_std
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
			
			Rec_obj.clear(); // reset all Rec values
			
			Fbar_top = 0.;
			Fbar_bot = 0.;
			Mbar_top = 0.;
			Mbar_bot = 0.;
			
			y = 0;
						
			//first year of simulation is same as terminal year from assessment	
			PE(y, _) = rnorm(rep(0., A), PE_Std);
			Nage(y, _) = Rcpp::exp( Rcpp::log(N1) + PE(y, _) ); //Start using N1 with PE
			Bage(y, _) = Nage(y, _) * Stock_wt;
			SSBage(y, _) = Bage(y, _) * Maturity;

			Fya_obj.Fya(y, _) = F_terminal;
			Mya_obj.Mya(y, _) = M_terminal;
			Zage(y, _) = F_terminal + M_terminal;
			Cage(y, _) = baranov_catch(Nage(y, _), F_terminal, M_terminal); // Nage(y, _) * (1. - Rcpp::exp(-Zage(y, _))) * Fterminal / Zage(y, _);
			Yage(y, _) = Cage(y, _) * Catch_wt; // Cage(y, _) * Catch_wt(_);
			
			for(a = 0; a < A; a++){
				//pop-wt Fbar and Mbar
				Fbar_top += ( (a > Fbar_ind(0)) & (a < Fbar_ind(1)) ) * Nage(y, a) * Fya_obj.Fya(y, a);
				Fbar_bot += ( (a > Fbar_ind(0)) & (a < Fbar_ind(1)) ) * Nage(y, a);
				Mbar_top += ( (a > Mbar_ind(0)) & (a < Mbar_ind(1)) ) * Nage(y, a) * Mya_obj.Mya(y, a);
				Mbar_bot += ( (a > Mbar_ind(0)) & (a < Mbar_ind(1)) ) * Nage(y, a);
				
				//3D matrices with age components, e.g. for yield5-8
				YaA(a, sim, y) = Yage(y, a);
				BaA(a, sim, y) = SSBage(y, a);
				SaA(a, sim, y) = Bage(y, a);
			}
					
			Fbar_popwt(y) = Fbar_top / Fbar_bot;
			Mbar_popwt(y) = Mbar_top / Mbar_bot;
			
			Ydiff(y) = sum(Yage(y, _)) - terminalYield;
			Yield(y) = sum(Yage(y, _));
			TAC(y) = TAC0(y);
			TAC_disp(y) = TAC(y) - Yield(y);
			Fbar(y) = solveF(TAC(y), 0.1, Nage(y, _), Mya_obj.Mya(y, _), Catch_wt); //Fya_obj.F_(y)
			Abundance(y) = sum(Nage(y, _));
			Biomass(y) = sum(Bage(y, _));
			SSB(y) = sum(SSBage(y, _));
			
			//start of projection
			for(y = 1; y < y_sim; y++){
				
				if(y < r_lag){
					Nstock = N_rec(y, _) + rnorm(rep(0., A), PE_Std);
				}else{
					Nstock = Nage(y-r_lag, _);
				}
				ssby = sum( Nstock * Stock_wt * Maturity );
				
				// Recruitment
				logrec = Rec_obj.logRec(ssby);

				// Abundance
				PE(y, _) = rnorm(rep(0., A), PE_Std); 
				logN = logN_a(  Rcpp::log(Nage(y-1, _)), 
								Mya_obj.Mya(y-1, _) + Fya_obj.Fya(y-1, _), 
								PE(y, _), 
								logrec );
				Nage(y, _) = Rcpp::exp(logN);
				
				// Natural Mortality
				if( M_cond_on ){
					Mya_obj.MCondition(y);
				}
				if( M_age_on ){
					Mya_obj.MAgeEffect(y);
				}
				if( M_year_on ){
					Mya_obj.MYearEffect(y);
				}
				if( M_AR_on ){
					Mya_obj.MAR(y);
				}
							
				if(y < y0){ //Fixed TAC/F years
				
					TAC(y) = TAC0(y);			// init F -v-
					Fbar(y) = solveF( TAC(y), 0.1, Nage(y, _), Mya_obj.Mya(y, _), Catch_wt );
					
				}else{ //Years for Variable MPs	
				
					MP_biomass = Bage;
					MP_ssb = SSBage;
				
					if(MP_Type){
						
						Fbar(y) = MP_decision(this, MP_n);
						TAC(y) = yield<T>( baranov_catch(Nage(y, _), Fbar(y), Mya_obj.Mya(y, _)), Catch_wt );
						
					}else if(!MP_Type){					
					
						Fbar(y) = MP_decision(this, MP_n);
						
					}
					
				}	
				TAC(y) = TAC(y) < Cmin ? Cmin : TAC(y);
				TAC(y) *= std::exp( R::rnorm(0., TAC_Error) );
				Fbar(y) = solveF( TAC(y), .1, Nage(y, _), Mya_obj.Mya(y, _), Catch_wt );
				// If TAC is too low, set minimum to Cmin
				// if( !str_is(mp, "No fishing") | (TAC(y) < Cmin) ) {
					// TAC(y) = Cmin;
					// Fbar(y) = solveF( TAC(y), Fya_obj.Fya(y, _), Nage(y, _), Mya_obj.Mya(y, _), Catch_wt );
				// }			
							
				// Fishing Mortality
				Fya_obj(Fbar(y));

				// Fill containers
				Bage(y, _) = Nage(y, _) * Stock_wt;
				SSBage(y, _) = Bage(y, _) * Maturity;
				
				// Fage(y, _) = Fbar(y) * Fya_obj.Sel;
				Zage(y, _) = Fya_obj.Fya(y, _) + Mya_obj.Mya(y, _);
				Cage(y, _) = baranov_catch(Nage(y, _), Fya_obj.Fya(y, _), Mya_obj.Mya(y, _)); // Nage(y, _) * (1. - std::exp(-Zage(y, _))) * Fage(y, _) / Zage(y, a);
				Yage(y, _) = Cage(y, _) * Catch_wt; // Cage(y, a) * Catch_wt(a);
							
				Fbar_top = 0.;
				Fbar_bot = 0.;
				Mbar_top = 0.;
				Mbar_bot = 0.;
				
				for(a = 0; a < A; a++){
					//pop-wt Fbar & Mbar
					Fbar_top += ( (a > Fbar_ind(0)) & (a < Fbar_ind(1)) ) * Nage(y, a) * Fya_obj.Fya(y, a);
					Fbar_bot += ( (a > Fbar_ind(0)) & (a < Fbar_ind(1)) ) * Nage(y, a);
					Mbar_top += ( (a > Mbar_ind(0)) & (a < Mbar_ind(1)) ) * Nage(y, a) * Mya_obj.Mya(y, a);
					Mbar_bot += ( (a > Mbar_ind(0)) & (a < Mbar_ind(1)) ) * Nage(y, a);
					
					//3D matrices with age components, e.g. for yield3-8
					YaA(sim, y, a) = Yage(y, a);
					BaA(sim, y, a) = SSBage(y, a);
					SaA(sim, y, a) = Bage(y, a);
				}
				
				Fbar_popwt(y) = Fbar_top / Fbar_bot;
				Mbar_popwt(y) = Mbar_top / Mbar_bot;					
						
				Abundance(y) = sum(Nage(y, _));
				Ydiff(y) = sum(Yage(y, _) - Yage(y-1, _));
				Yield(y) = sum(Yage(y, _));
				TAC_disp(y) = TAC(y) - Yield(y);
				Biomass(y) = sum(Bage(y, _));
				SSB(y) = sum(SSBage(y, _));
				
			} // end projection
			
			//Annual values
			all_TAC(sim, _) = clone(TAC);
			all_TAC_disp(sim, _) = clone(TAC_disp);
			all_Fbar(sim, _) = clone(Fbar);
			all_Fbar_popwt(sim, _) = clone(Fbar_popwt);
			all_Mbar_popwt(sim, _) = clone(Mbar_popwt);
			all_Yield(sim, _) = clone(Yield);
			all_Ydiff(sim, _) = clone(Ydiff);
			all_Abundance(sim, _) = clone(Abundance);
			all_Biomass(sim, _) = clone(Biomass);
			all_SSB(sim, _) = clone(SSB);
			all_Rec(sim, _) = clone( Rcpp::exp(Rec_obj.getRec()) );
			
		} // end simulation
		
		return  List::create(
				Named("OM") = OM_name,
				Named("Mort") = M_name,
				Named("Rec") = R_name,
				Named("MP") = MP_name,
		
				Named("TAC") = all_TAC,
				Named("TAC_disp") = TAC_disp,
				Named("Fbar") = all_Fbar,
				Named("Fbar_popwt") = all_Fbar_popwt,
				Named("Mbar_popwt") = all_Mbar_popwt,
				Named("Yield") = all_Yield,
				Named("Ydiff") = all_Ydiff,
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
	
	};
	
};




// Surplus Production Model Simulation Function //
namespace SPMSim{

	template<class T>
	List SurplusModel( List OM, 
			  List MP, 
			  int sim_no, 
			  int y_sim ){
				  
		NumericMatrix Abundance(sim_no, y_sim),
					  Biomass(sim_no, y_sim);
		
		int sim, y;
		
		for(sim = 0; sim < sim_no; sim++){
			for(y = 0; y < y_sim; y++){
			
				;
			
			}
		}
				
		return List::create(
			   Named("X") = 0.		   
			   );
				
	}
	
};

