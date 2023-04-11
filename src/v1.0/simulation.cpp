#include <math.h> //this is for std::pow, std::exp, and std::log
#include <stdio.h> //this is for std::printf
#include <string> //this is for std::string in 'str_is'

// #include <Rcpp.h> //not needed if including Armadillo
#include <RcppArmadillo.h>
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
using Rcpp::NumericMatrix;
using Rcpp::List;
using Rcpp::DataFrame;


#include <mseliteutils.cpp>
#include <Mortality.cpp>
#include <Ncalc.cpp>
#include <Recruitment.cpp>
#include <HCR.cpp>

//	Simulation Function //

template<class T>
List simulation(List OM, 
				List MP, 
				int sim_no, 
				int y_sim, 
				T Cmin ){			

	// Call OM and MP elements //
	
	int y_hcr = MP["y_hcr"];
	T std_pe = OM["std_pe"],
	  terminalYield = OM["yieldterminal"],
   	  max_S = OM["max_S"];
	String om = OM["name"], 
		   mproj = OM["mproj"],
 		   rproj = OM["rproj"],
		   mp = MP["name"];
	NumericVector Catch_wt = OM["catch_wt"],
				  Stock_wt = OM["stock_wt"],
				  maturity = OM["maturity"],
				  Favg = OM["Favg"],
				  mscale = OM["mscale"],
				  mscale_low = OM["mscale_low"],
				  mscale_high = OM["mscale_high"],
				  N1 = OM["N1"],
				  Fterminal = OM["Fterminal"],
				  Mterminal = OM["Mterminal"],
				  N2020 = OM["N2020"],
				  m_est = OM["mpar_est"],
				  mpar_sd = OM["mpar_sd"],
				  init_TAC = MP["TAC"],
				  ssbr_bracket = MP["ssbrbracket"], 
				  fbracket = MP["fbracket"];
				  // Fbar_ind = OM["Fbar_ind"] //ages over which to average popwt F; vector of length 2
	IntegerVector m_ind = OM["m_ind"],
				  Findex = OM["keyF"];
	NumericMatrix Fmat = OM["F_err"];
		
	int A = OM["A"];
	int y, a, sim;
	
	IntegerVector years = OM["years"];
	int y0 = years(years.size()-1); //last year of assessment = first year of the MSE
					
	// Create objects for sim //
	
	T ssby, logrec, std_log_r, Fbar_top, Fbar_bot;
	
	NumericVector Nstock(A), Merr(A), Ferr(A); //vectors for ages
	
	 //vectors for years
	NumericVector TAC(y_sim), Fbar(y_sim), Fbar_popwt(y_sim), 
				  Ydiff(y_sim), Yield(y_sim), Abundance(y_sim),
				  Biomass(y_sim), SSB(y_sim);
				  
	//matrices for years and ages
	NumericMatrix logNage(y_sim, A), Zage(y_sim, A), Fage(y_sim, A), 
				  sel(y_sim, A), Mage(y_sim, A), Mpar(y_sim, A),
				  Cage(y_sim, A), Yage(y_sim, A), Bage(y_sim, A), 
				  PE(y_sim, A), Nage(y_sim, A), SSBage(y_sim, A);
								
	//matrices for years and sims -- for Mains and All_summary values
	NumericMatrix all_Yield(sim_no, y_sim), all_TAC(sim_no, y_sim), all_Fbar(sim_no, y_sim),
				  all_Fbar_popwt(sim_no, y_sim), all_Ydiff(sim_no, y_sim),
				  all_Biomass(sim_no, y_sim), all_SSB(sim_no, y_sim),
				  all_Rec(sim_no, y_sim), all_Abundance(sim_no, y_sim);
	
	//catch PMs for dt_avgs
	T avgC_short_val = 0., avgC_long_val = 0., aav_val = 0.; //annual values
	NumericVector aav_all(sim_no), avgC_long_all(sim_no), avgC_short_all(sim_no);  //all sim values
	
	//Preventable Decline container
	// NumericMatrix PD(sim_no, y_sim);
	NumericVector PD(y_sim);
		
	HCR<T> MP_vals; //HCR is a struct defined in utils used for MP
	
	T max_F; //for selectivity
	NumericVector Ferr_vec(8); //for err in F
	NumericVector Fvec(8, 0.); //a vector of zeroes for mvnorm
	NumericVector sel_vec(A); //for selectivity
	NumericVector scaledMindex(A); //for M cond
	NumericVector sim_rec(y_sim); //temporary container for all_Rec
	
	// if(OM["mproj"] == "avg3yr") { vector<T> scaledMindex = apply(OM["scaledMindex"], 1, mean); }
	if( str_is(mproj, "termyr")) { scaledMindex = mscale; }
	if( str_is(mproj, "min") ) { scaledMindex = mscale_low; } //apply<T>(mscale_low, 0, mean)
	if( str_is(mproj, "max") ) { scaledMindex = mscale_high; } //apply<T>(mscale_high, 0, mean)
	
	//Recrutiment Stuff
	NumericMatrix srr_obj(2, 4);
		
	// Start simulation //
	for(sim = 0; sim < sim_no; sim++){
		
		//SRR is defined in mseliteutils
		srr_obj = SRR<T>(OM);
		//Do this here so parms stay consistence within each projection
				
		std::printf("Thread: %d \n", sim+1);
		
		// set as zero each iter
		Yield.fill(0.);
		Ydiff.fill(0.);
		Fbar_popwt.fill(0.);
		Abundance.fill(0.);
		Biomass.fill(0.);
		SSB.fill(0.);
		
		avgC_short_val = 0.;
		avgC_long_val = 0.;
		aav_val = 0.;
		
		Fbar_top = 0.;
		Fbar_bot = 0.;
		
		y = 0;
		
		//first year of simulation is same as terminal year from assessment
		for(a = 0; a < A; a++){
			
			PE(y, a) = (T)R::rnorm(0., std_pe); //R::rnorm returns 1 value
			logNage(y, a) = std::log(N1(a)) + PE(y, a); //sample of terminal N
			Nage(y, a) = std::exp(logNage(y, a));
			
			Bage(y, a) = Nage(y, a) * Stock_wt(a);
			SSBage(y, a) = Bage(y, a) * maturity(a);

			Fage(y, a) = Fterminal(a);
			Mage(y, a) = Mterminal(a);
			Zage(y, a) = Fage(y, a) + Mage(y, a);
			Cage(y, a) = Nage(y, a) * (1. - std::exp(-Zage(y, a))) * Fage(y, a) / Zage(y, a);
			Yage(y, a) = Cage(y, a) * Catch_wt(a);
			
			Fbar_top += ( (a > 2) & (a < 7) ) * Nage(y, a) * Fage(y, a);
			Fbar_bot += ( (a > 2) & (a < 7) ) * Nage(y, a);
			Ydiff(y) += Yage(y, a);
			Yield(y) += Yage(y, a);
			Abundance(y) += Nage(y, a);
			Biomass(y) += Bage(y, a);
			SSB(y) += SSBage(y, a);
			
		}
		sel(y, _) = scale(Favg, 1./max(Favg));
		
		Ydiff(y) -= terminalYield;
		TAC(y) = init_TAC(y);
		Fbar(y) = solveF(TAC(y), sel(y, _), Nage(y, _), Mage(y, _), Catch_wt);
		
		Fbar_popwt(y) = Fbar_top / Fbar_bot;
		
		//start of projection at year 2022
		for(y = 1; y < y_sim; y++){
			
			ssby = 0.;
			if(y == 1){
				for(a = 0; a < A; a++){
					Nstock(a) = N2020(a) + (T)R::rnorm(0., std_pe);
				}
			}else{
				Nstock = Nage(y-2, _);
			}
			for(a = 0; a < A; a++){
				ssby += Nstock(a) * Stock_wt(a) * maturity(a);
			}
			
			// Recruitment Here //
			if( str_is(rproj, "BH") ){
				logrec = BH(ssby, srr_obj(0, _));
				std_log_r = OM["bherr"];
			}else if( str_is(rproj, "SigBH") ){
				logrec = SigBH(ssby, max_S, srr_obj(1, _));
				std_log_r = OM["bhsigerr"];
			}else{
				logrec = OM["logR"];
				std_log_r = OM["log_std_log_r"];
				std_log_r = std::exp( std_log_r );
			}
			logrec = logrec + R::rnorm(0., std_log_r);

			PE(y, _) = Rcpp::rnorm(A, 0., std_pe);
			logNage(y, _) = calc_logN(A, logNage(y-1, _), Mage(y-1, _), Fage(y-1, _), PE(y, _), logrec);
			Nage(y, _) = Rcpp::exp(logNage(y, _));
    		
			Mage(y, _) = rep(0.3, A);
			if( !str_is(mproj, "base") ){
				
				for(a = 0; a < A; a++){
					Merr(a) = R::rnorm(0., mpar_sd(a));
					Mpar(y, a) = m_est(m_ind(a)-1) + Merr(a);
					Mage(y, a) = std::exp( std::log(Mage(y, a)) + Mpar(y, a) * scaledMindex(a) );
				}
				
			}
			
			Ferr_vec = mvrnorm(Fvec, Fmat);
			for(a = 0; a < A; a++){
				Ferr(a) = Ferr_vec(Findex(a));
				sel_vec(a) = std::exp( std::log(Favg(a)) + Ferr(a) );
			}
			max_F = max(sel_vec);
			sel(y, _) = scale(sel_vec, 1./max_F);
			
			if(y < y_hcr){ //2021 and 2022
			
				TAC(y) = init_TAC(y);
				Fbar(y) = solveF(TAC(y), sel(y, _), Nage(y, _), Mage(y, _), Catch_wt);

			}else{ //2023 and higher

			    // TAC is defined by MP_decision function
				MP_vals = MP_decision<T>( //Ubound, Fbound,
										  // harvestrate, threshold,
										  // ssbr_bracket, fbracket,
										  mp,
										  Nage(y, _), Mage(y, _), sel(y, _), 
										  maturity, Stock_wt, Catch_wt);

				TAC(y) = MP_vals.TAC;
				Fbar(y) = MP_vals.F;
				
			}			  
			// std::printf("F = %f, ", MP_vals.F);
				// If TAC is too low, set minimum to Cmin
				if( !str_is(mp, "No fishing") & (TAC(y) < Cmin) ) {
					TAC(y) = Cmin;
					Fbar(y) = solveF(Cmin, sel(y, _), Nage(y, _), Mage(y, _), Catch_wt);
				}
						
			Fbar_top = 0.;
			Fbar_bot = 0.;
			for(a = 0; a < A; a++){
			
				Bage(y, a) = Nage(y, a) * Stock_wt(a);
				SSBage(y, a) = Bage(y, a) * maturity(a);
				
				Fage(y, a) = Fbar(y) * sel(y, a);
				Zage(y, a) = Fage(y, a) + Mage(y, a);
				Cage(y, a) = Nage(y, a) * (1. - std::exp(-Zage(y, a))) * Fage(y, a) / Zage(y, a);
				Yage(y, a) = Cage(y, a) * Catch_wt(a);
				
				Fbar_top += ( (a > 2) & (a < 7) ) * Nage(y, a) * Fage(y, a);
				Fbar_bot += ( (a > 2) & (a < 7) ) * Nage(y, a);
				Ydiff(y) += Yage(y, a) - Yage(y-1, a);
				Yield(y) += Yage(y, a);
				Abundance(y) += Nage(y, a);
				Biomass(y) += Bage(y, a);
				SSB(y) += SSBage(y, a);
				
			}
			
			Fbar_popwt(y) = Fbar_top / Fbar_bot;
			
			//Preventable Decline
			PD(y) += SSB(y) < SSB(0);
			
		} // end projection
		
		//this is for Mains
		all_TAC(sim, _) = clone(TAC);
		all_Fbar(sim, _) = clone(Fbar);
		all_Fbar_popwt(sim, _) = clone(Fbar_popwt);
		all_Yield(sim, _) = clone(Yield);
		all_Ydiff(sim, _) = clone(Ydiff);
		all_Abundance(sim, _) = clone(Abundance);
		all_Biomass(sim, _) = clone(Biomass);
		all_SSB(sim, _) = clone(SSB);
		sim_rec = Nage(_, 0);
		all_Rec(sim, _) = clone(sim_rec);
				
		//catch PM stuff
		// aav_val = std::abs(Yield(0) - terminalYield);
		for(y = 2; y < 27; y++){ //y = 2 is 2023
			
			if(y < 7){
				avgC_short_val += Yield(y);
			}
			avgC_long_val += Yield(y);
			aav_val += std::abs(Yield(y) - Yield(y-1));
						
		}
		aav_val /= avgC_long_val;
		avgC_long_val /= 25.;
		avgC_short_val /= 5.;
		
		aav_all(sim) = aav_val;
		avgC_long_all(sim) = avgC_long_val;
		avgC_short_all(sim) = avgC_short_val;
		
	} // end simulation
	
	//Preventable Decline
	// NumericVector PD_y = apply(PD, sum, 0);
	PD = scale(PD, 1./sim_no);
	
	IntegerVector sim_count = rep(seq(1, sim_no), y_sim),
				  main_years = rep(seq(y0, y0+y_sim-1), sim_no, true),
				  sim_seq = seq(1, sim_no),
				  summ_years = rep(seq(y0, y0+y_sim-1), 8);
	 
	//Mains is outputs for each year of each thread
	DataFrame Mains = DataFrame::create(Named("year") = main_years, // Named("") could also be _[""]
									  Named("TAC") = clone(as_vector(all_TAC)),
									  Named("Fbar") = clone(as_vector(all_Fbar)),
									  Named("Fbar_popwt") = clone(as_vector(all_Fbar_popwt)),
									  Named("Yield") = clone(as_vector(all_Yield)),
									  Named("Ydiff") = clone(as_vector(all_Ydiff)),
									  Named("Abundance") = clone(as_vector(all_Abundance)),
									  Named("Biomass") = clone(as_vector(all_Biomass)),
									  Named("SSB") = clone(as_vector(all_SSB)),
									  Named("Rec") = clone(as_vector(all_Rec)),
									  Named("sim") = sim_count,
									  Named("om") = rep(om, y_sim*sim_no),
									  Named("mp") = rep(mp, y_sim*sim_no),
									  Named("mproj") = rep(mproj, y_sim*sim_no),
									  Named("rproj") = rep(rproj, y_sim*sim_no));
									  
	//Avgs is PMs for each thread
	DataFrame Avgs = DataFrame::create(Named("avgC_short") = clone(avgC_short_all),
									 Named("avgC_long ") = clone(avgC_long_all),
									 Named("aav") = clone(aav_all),
									 Named("sim") = sim_seq,
									 Named("om") = rep(om, sim_no),
									 Named("mp") = rep(mp, sim_no),
									 Named("mproj") = rep(mproj, sim_no),
									 Named("rproj") = rep(rproj, sim_no));

	//means of matrices
	// NumericVector TAC_mean = Mmean(all_TAC, 1),
				  // Fbar_mean = Mmean(all_Fbar, 1),
				  // Fbar_popwt_mean = Mmean(all_Fbar_popwt, 1),
				  // Yield_mean = Mmean(all_Yield, 1),
				  // Ydiff_mean = Mmean(all_Ydiff, 1),
				  // Biomass_mean = Mmean(all_Biomass, 1),
				  // SSB_mean = Mmean(all_SSB, 1),
				  // Rec_mean = Mmean(all_Rec, 1);
	// //std of matrices
	// NumericVector TAC_std = Mstd(all_TAC, 1),
				  // Fbar_std = Mstd(all_Fbar, 1),
				  // Fbar_popwt_std = Mstd(all_Fbar_popwt, 1),
				  // Yield_std = Mstd(all_Yield, 1),
				  // Ydiff_std = Mstd(all_Ydiff, 1),
				  // Biomass_std = Mstd(all_Biomass, 1),
				  // SSB_std = Mstd(all_SSB, 1),
				  // Rec_std = Mstd(all_Rec, 1);

	//quantile matrices
	NumericMatrix p05(y_sim, 8), p10(y_sim, 8), p50(y_sim, 8), p75(y_sim, 8), p90(y_sim, 8), p95(y_sim, 8);

	StringVector var_names =  {"TAC", "Fbar", "Fbar_popwt", "Yield",
							   "Ydiff", "SSB", "Biomass", "Rec"};
	StringVector all_vars = rep( var_names, y_sim, true ); // 'at = true'
												
	for(int item = 0; item < var_names.size(); item++){
		
		// NumericVector mean_vec(y_sim), std_vec(y_sim);
		
		for(y = 0; y < y_sim; y++){
		
			switch(item){
				case 0:
					p05(y, item) = quantile(all_TAC(_, y), .05);
					p10(y, item) = quantile(all_TAC(_, y), .10);
					p50(y, item) = quantile(all_TAC(_, y), .50);
					p75(y, item) = quantile(all_TAC(_, y), .75);
					p90(y, item) = quantile(all_TAC(_, y), .90);
					p95(y, item) = quantile(all_TAC(_, y), .95);
					// mean_vec = TAC_mean;
					// std_vec = TAC_std;
					break;
				case 1: 
					p05(y, item) = quantile(all_Fbar(_, y), .05);
					p10(y, item) = quantile(all_Fbar(_, y), .10);
					p50(y, item) = quantile(all_Fbar(_, y), .50);
					p75(y, item) = quantile(all_Fbar(_, y), .75);
					p90(y, item) = quantile(all_Fbar(_, y), .90);
					p95(y, item) = quantile(all_Fbar(_, y), .95);
					// mean_vec = Fbar_mean;
					// std_vec = Fbar_std;
					break;
				case 2: 
					p05(y, item) = quantile(all_Fbar_popwt(_, y), .05);
					p10(y, item) = quantile(all_Fbar_popwt(_, y), .10);
					p50(y, item) = quantile(all_Fbar_popwt(_, y), .50);
					p75(y, item) = quantile(all_Fbar_popwt(_, y), .75);
					p90(y, item) = quantile(all_Fbar_popwt(_, y), .90);
					p95(y, item) = quantile(all_Fbar_popwt(_, y), .95);
					// mean_vec = Fbar_popwt_mean;
					// std_vec = Fbar_popwt_std;
					break;
				case 3: 
					p05(y, item) = quantile(all_Yield(_, y), .05);
					p10(y, item) = quantile(all_Yield(_, y), .10);
					p50(y, item) = quantile(all_Yield(_, y), .50);
					p75(y, item) = quantile(all_Yield(_, y), .75);
					p90(y, item) = quantile(all_Yield(_, y), .90);
					p95(y, item) = quantile(all_Yield(_, y), .95);
					// mean_vec = Yield_mean;
					// std_vec = Yield_std;
					break;
				case 4: 
					p05(y, item) = quantile(all_Ydiff(_, y), .05);
					p10(y, item) = quantile(all_Ydiff(_, y), .10);
					p50(y, item) = quantile(all_Ydiff(_, y), .50);
					p75(y, item) = quantile(all_Ydiff(_, y), .75);
					p90(y, item) = quantile(all_Ydiff(_, y), .90);
					p95(y, item) = quantile(all_Ydiff(_, y), .95);
					// mean_vec = Ydiff_mean;
					// std_vec = Ydiff_std;
					break;
				case 5: 
					p05(y, item) = quantile(all_SSB(_, y), .05);
					p10(y, item) = quantile(all_SSB(_, y), .10);
					p50(y, item) = quantile(all_SSB(_, y), .50);
					p75(y, item) = quantile(all_SSB(_, y), .75);
					p90(y, item) = quantile(all_SSB(_, y), .90);
					p95(y, item) = quantile(all_SSB(_, y), .95);
					// mean_vec = SSB_mean;
					// std_vec = SSB_std;
					break;
				case 6: 
					p05(y, item) = quantile(all_Biomass(_, y), .05);
					p10(y, item) = quantile(all_Biomass(_, y), .10);
					p50(y, item) = quantile(all_Biomass(_, y), .50);
					p75(y, item) = quantile(all_Biomass(_, y), .75);
					p90(y, item) = quantile(all_Biomass(_, y), .90);
					p95(y, item) = quantile(all_Biomass(_, y), .95);
					// mean_vec = Biomass_mean;
					// std_vec = Biomass_std;
					break;
				case 7: 
					p05(y, item) = quantile(all_Rec(_, y), .05);
					p10(y, item) = quantile(all_Rec(_, y), .10);
					p50(y, item) = quantile(all_Rec(_, y), .50);
					p75(y, item) = quantile(all_Rec(_, y), .75);
					p90(y, item) = quantile(all_Rec(_, y), .90);
					p95(y, item) = quantile(all_Rec(_, y), .95);
					// mean_vec = Rec_mean;
					// std_vec = Rec_std;
					break;
			}			
			
		}
		
	}

	//All summary is qunatiles over all sims for each year
	DataFrame All_summary = DataFrame::create(Named("p05") = clone(as_vector(p05)),
												Named("p10") = clone(as_vector(p10)),
												Named("p50") = clone(as_vector(p50)),
												Named("p75") = clone(as_vector(p75)),
												Named("p90") = clone(as_vector(p90)),
												Named("p95") = clone(as_vector(p95)),
												Named("year") = summ_years,
												Named("var") = all_vars,
												Named("om") = rep(om, y_sim*8),
												Named("mp") = rep(mp, y_sim*8),
												Named("mproj") = rep(mproj, y_sim*8),
												Named("rproj") = rep(rproj, y_sim*8)
												);
												
	DataFrame PD_Data = DataFrame::create(Named("Year") = seq(y0, y0+y_sim-1),
										  Named("PD") = PD
										  );
	
	List MSEout = List::create( Named("Avgs") = Avgs,
								Named("Mains") = Mains,
								Named("All_summary") = All_summary, 
								Named("Preventable_Decline") = PD_Data
								);
	
	return MSEout;
				
} // end function



// [[Rcpp::export]]
List Simulate(  List OM, 
				List MP, 
				int sim_no, 
				int y_sim, 
				double Cmin = 0.1) {
					
	return simulation<double> ( OM, 
								MP, 
								sim_no, 
								y_sim, 
								Cmin );
								
}