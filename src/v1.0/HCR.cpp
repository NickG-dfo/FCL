
// Define MP_decision function //

//Create class  for reference point functions
template<class T>
class RP{
	
		const int A;
		NumericVector Mage,
					  Sage, 
					  Mat,
					  Stock_wt,
					  Catch_wt;
		SR::R_obj<T> r_obj;
	
		RP(List &LHC, SR::R_obj<T> &Rec_object): r_obj(Rec_object) {
			A = LHC["A"];
			Mage = LHC["M"];
			Sage = LHC["S"];
			Mat = LHC["Mt"];
			Stock_wt = LHC["Sw"];
			Catch_wt = LHC["Cw"];
		}
	
	public:
		
		// template<class T>
		//No typedef because its decalred in Recruitment.cpp
		SPR(T F){
			NumericVector M = cat(0., Mage), 
						  S = cat(0., Sage),
						  cZ = cumsum(M - S * F);
			T spr = sum( (1. - Rcpp::exp(-cZ)) * Mat * Stock_wt );
			return spr;
		}
		
		// template<class T>
		T YPR(T F){
			NumericVector M = cat(0., Mage), 
						  S = cat(0., Sage),
						  Z = (M - S * F),
						  cZ = cumsum(Z);
			T ypr = sum( (1. - Rcpp::exp(-cZ)) * (F/Z) * Catch_wt );
			return ypr;
		}
		
		// template<class T>
		T SSBeq(T F){ // Newton-Raphson root
			if( str_is(r_obj, "Const") ){
				std::printf("Constant recruitment does not allow equilbirum states.");
				return 0.;
			}
		
			T delta = 1.,
			  Si = 1.;
			T y, dy;
			do{
				y = (r_obj.Rec(Si)*SPR(F) - Si);
				dy = (r_obj.Rec(Si+.01)*SPR(F) - r_obj.Rec(Si-.01)*SPR(F))/.02 - 1.;
				delta = dy/y;
				Si -= delta;
			} while(delta > .001);
			return Si;
		}
		
		// template<class T>
		T Yeq(T F){
			return SSBeq(F) * SPR(F);
		}
		
		// template<class T>
		T F_SPRx(T percent){
			T delta = 1.,
			  Fi = 1.;
			T y, dy;
			do{
				y = (percent*SPR(0.) - SPR(Fi));
				dy = (SPR(Fi+.0001) - SPR(Fi-.0001))/.0002;
				delta = -dy/y;
				Fi -= delta;
			} while(delta > .001);
			return Fi;
		}
		
		// template<class T>
		T FMSY(void){
			T delta = 1.,
			  Fi = 1.;
			T y, dy;
			do{
				y = -std::pow(Yeq(Fi), 2.);
				dy = (Yeq(Fi+.01) - Yeq(Fi-.01))/.02 - 1.;
				delta = dy/y;
				Fi -= delta;
			} while(delta > .001);
			return Fi;
		}
		
		// template<class T>
		T BMSY(void){
			return SSBeq( FMSY() );
		}
		
		// template<class T>
		T MSY(void){
			return Yeq( FMSY() );
		}
	
};

//Define struct object called HCR for MP functions and values
template<class T>
struct HCR{
	
	T TAC;
	T F;
	HCR& operator= (const HCR tf){
		this->TAC = tf.TAC;
		this->F = tf.F;
		return *this;
	}
	
	// Like TMB, this is a function that will be defined externally //
	
	void SSB_MP(T, int);
	// Define as -
	//		template<class T>
	// 		HCR<T>::SSB_MP(T S, int n){ *MP block* }
	// S is SSB (for SSB-based MP), n is MP selection number (if multiple are defined)
	
	void Selctivity_MP(NumericVector, int);
	// Define as -
	//		template<class T>
	// 		HCR<T>::Selectivity_MP(NumericVector Sel, int n){ *MP block* }
	// Sel is selectivity-at-age (for Selectivity-based MP), n is MP selection number (if multiple are defined)
	
	void Closure_MP(IntegerVector, int);
	// Define as -
	//		template<class T>
	// 		HCR<T>::Closure_MP(IntegerVector strata, int n){ *MP block* }
	// strata is strata number (ints) (for Closure-based MP), n is MP selection number (if multiple are defined)
	
	// All MP functions should be defined, but may just 'return' if unneeded //
	// For example, Closure_funtion(IntegerVector strata, int n) { return; }
	// 'int n' may also be used to choose between different types of MPs
	//		i.e. n = 0 for SSB, n = 1 for Selectivity, and n = 2 for Closure
	
};
	
//Define function for MP HCRs (for 3Ps RP)
template<class T>
HCR<T> MP_decision(//T F, T Catch, 
				// const T Ubound, const T Fbound,
				// const T harvestrate, const T threshold,
				// NumericVector ssbr_bracket, const NumericVector f_bracket,
				// NumericVector fbreak,
				const String MP_name,
				const NumericVector &Na, const NumericVector &Ma, const NumericVector &sel,
				const NumericVector &mat, const NumericVector &Sw, const NumericVector &Cw){
	
	// declare values for MPs
	HCR<T> MP_out;
	short A = Na.size();
	// T Lbound = 66.;
	NumericVector temp_catch(A), temp_f(A);
	
	T ssb = 0.;	
	for(short i = 0; i < A; i++){
		ssb += Na(i) * mat(i) * Sw(i);
	}
	
	//define TAC/F based on MP_name
	if(MP_name == "No fishing"){
		MP_out.TAC = 0.;
		MP_out.F = 0.001;
	}
	else if(MP_name == "Critical 100"){
		MP_out.F = 0.001 +
					( (ssb > 66.) & (ssb <= 1.75*66.) ) * (ssb - 66.)/(.75*66.) * .1 +
					(ssb > 1.75*66.) * .1;
		
		temp_f = scale(sel, MP_out.F);
		temp_catch = baranov_catch(Na, temp_f, Ma);
		MP_out.TAC = yield<T>(temp_catch, Cw);
	}
	else if(MP_name == "Elbow C"){
		T m1 = (.035-.025)/(.6*66.);
		T b1 = .035-m1*66.;
		T m2 = (.1-.035)/(66.);
		T b2 = .1-m2*(2.*66.);
		MP_out.F = 0.001 +
				  ( (ssb > .4*66.) & (ssb <= 66.) ) * (m1*ssb+b1) +
				  ( (ssb > 66.) & (ssb <= 2.*66.) ) * (m2*ssb+b2) +
				  (ssb > 2.*66.) * .1;
		temp_catch = baranov_catch(Na, scale(sel, MP_out.F), Ma);
		MP_out.TAC = yield<T>(temp_catch, Cw);
	}
	else if(MP_name == "AGC"){
		T m1 = (.2-.05)/(66.);
		T b1 = .2-m1*(2.*66.);
		MP_out.F = 0.001 +
				  ( (ssb > .4*66.) & (ssb <= 66.) ) * .05 +
				  ( (ssb > 66.) & (ssb <= 2.*66.) ) * (m1*ssb+b1) +
				  (ssb > 2.*66.) * .2;
		temp_catch = baranov_catch(Na, scale(sel, MP_out.F), Ma);
		MP_out.TAC = yield<T>(temp_catch, Cw);
	}
	else if(MP_name == "FFAW"){
		T m1 = (3.-1.)/(.4*66.);
		T b1 = 3.-m1*1.2*66.;
		T m2 = (20.-3.)/(.8*66.);
		T b2 = 20.-m2*(2.*66.);
		MP_out.TAC = 0.001 +
				  ( (ssb > .4*66.) & (ssb <= .8*66.) ) * 1. +
				  ( (ssb > .8*66.) & (ssb <= 1.2*66.) ) * (m1*ssb+b1) +
				  ( (ssb > 1.2*66.) & (ssb <= 2.*66.) ) * (m2*ssb+b2) +
				  (ssb > 2.*66.) * 20.;
		MP_out.F = solveF(MP_out.TAC, sel, Na, Ma, Cw);
	}
	
	// else if(MP_name == "harvestrate"){
		// T surplus = (ssb - threshold) * harvestrate;
		// MP_out.TAC = surplus > 0. ? surplus : 0.;
		// MP_out.F = solveF(MP_out.TAC, sel, Ma, Na, Cw);
	// }
	// else if(MP_name == "f_bracket"){
		// bool which;
		// short i = 0;
		// ssbr_bracket = cat( cat(0., ssbr_bracket), 100.);
		// for(; i < ssbr_bracket.size(); i++){
			// which = ( (ssb/66. >= ssbr_bracket(i)) & (ssb/66. < ssbr_bracket(i+1)) );
			// if(which) { break; }
		// }
		// MP_out.F = f_bracket(i);
		
		// temp_catch = baranov_catch(Na, scale(sel, MP_out.F), Ma);
		// MP_out.TAC = yield<T>(temp_catch, Cw);
	// }
	// else if(MP_name == "Elbow_PA"){
		// T m1 = (0.05 - 0.025)/(0.6 * 66.);
		// T b1 = 0.05 - m1 * 66.;
		// T m2 = (0.26 - 0.05)/66.;
		// T b2 = 0.26 - m2 * (2*66.);
		// MP_out.F =  ( (ssb > 0.4 * 66.) & (ssb <= 66.) ) * (m1 * ssb + b1) +
					// ( (ssb > 66.) & (ssb <= 2 * 66.) ) * (m2 * ssb + b2) +
					// (ssb > 2 * 66.) * 0.26;
					
		// temp_catch = baranov_catch(Na, scale(sel, MP_out.F), Ma);
		// MP_out.TAC = yield<T>(temp_catch, Cw);
	// }
	// else if(MP_name == "Elbow_Spectrum1"){
		// T m1 = (0.045 - 0.055)/(0.6 * 66.);
		// T b1 = 0.055 - m1 * 66.;
		// T m2 = (0.2 - 0.055)/66.;
		// T b2 = 0.2 - m2 * (2*66.);
		// MP_out.F =  ( (ssb > 0.4 * 66.) & (ssb <= 66.) ) * (m1 * ssb + b1) +
					// ( (ssb > 66.) & (ssb <= 2 * 66.) ) * (m2 * ssb + b2) +
					// (ssb > 2 * 66.) * 0.2;
					
		// temp_catch = baranov_catch(Na, scale(sel, MP_out.F), Ma);
		// MP_out.TAC = yield<T>(temp_catch, Cw);
	// }
	// else if(MP_name == "Elbow_Spectrum2"){
		// T m1 = (0.07 - 0.03)/(0.6 * 66.);
		// T b1 = 0.07 - m1 * 66.;
		// T m2 = (0.2 - 0.07)/66.;
		// T b2 = 0.2 - m2 * (2*66.);
		// MP_out.F =  ( (ssb > 0.4 * 66.) & (ssb <= 66.) ) * (m1 * ssb + b1) +
					// ( (ssb > 66.) & (ssb <= 2 * 66.) ) * (m2 * ssb + b2) +
					// (ssb > 2 * 66.) * 0.2;
					
		// temp_catch = baranov_catch(Na, scale(sel, MP_out.F), Ma);
		// MP_out.TAC = yield<T>(temp_catch, Cw);
	// }
	
	//return struct object with TAC and F values
	
	return MP_out;
	
}