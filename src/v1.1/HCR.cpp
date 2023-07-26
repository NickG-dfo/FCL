
namespace PopSim{

	//Create class  for reference point functions
	template<class T>
	class RP{
		
			int A;
			NumericVector Mage,
						  Sage, 
						  Mat,
						  Stock_wt,
						  Catch_wt;
			SR<T> Rec;
		
		public:	
			
			RP(List &LHC, SR<T> Recruitment): Rec(Recruitment){				
				A = LHC["A"];
				Mage = LHC["Mort"];
				Sage = LHC["Sel"];
				Mat = LHC["Mat"];
				Stock_wt = LHC["Sw"];
				Catch_wt = LHC["Cw"];
			}
						
			SPR(T F){
				NumericVector M = cat(0., Mage), 
							  S = cat(0., Sage),
							  cZ = cumsum(M - S * F);
				T spr = sum( (1. - Rcpp::exp(-cZ)) * Mat * Stock_wt );
				return spr;
			}
			
			T YPR(T F){
				NumericVector M = cat(0., Mage), 
							  S = cat(0., Sage),
							  Z = (M - S * F),
							  cZ = cumsum(Z);
				T ypr = sum( (1. - Rcpp::exp(-cZ)) * (F/Z) * Catch_wt );
				return ypr;
			}
			
			T SSBeq(T F){
				if( str_is( Rec.Rec_Type, "Const" ) ){
					std::printf("Constant recruitment does not allow equilbirum states.");
					return 0.;
				}
				 // Newton-Raphson root
				T delta = 1.,
				  Si = 1.;
				T y, dy;
				do{
					y = ((*Rec->recruitment)(Si, 0, 0)*SPR(F) - Si);
					dy = ((*Rec->recruitment)(Si+.01, 0, 0)*SPR(F) - (*Rec->recruitment)(Si-.01, 0, 0)*SPR(F))/.02 - 1.;
					delta = dy/y;
					Si -= delta;
				} while(delta > .001);
				return Si;
			}
			
			T Yeq(T F){
				return SSBeq(F) * SPR(F);
			}
			
			T F_SPRx(T percent){ // Newton-Raphson root
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
			
			T FMSY(void){ // Newton-Raphson root
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

	template<class T>
	T MP_decision(T, NumericVector);

	//Define struct object called HCR for MP functions and values
	// template<class T>
	// struct HCR{
		
		// //These should initialize as 0.
		// T TAC; 
		// T F;
		// Rcpp::Function MPf;
		// HCR& operator= (const HCR tf){
			// this->TAC = tf.TAC;
			// this->F = tf.F;
			// return *this;
		// }
		
		// HCR(String MPname){
			// MPf = Rcpp::Function(MPname);
		// }
		
		// //But how to change which value input becomes?
		// T Calc(Rcpp::RObject inputs){
			// return Rcpp::as<T>( MPf(inputs) );
		// }
		
		// Like TMB, this is a function that will be defined externally //
		
		// void F_MP(T, int);
		// void F_MP(T, String);
		
		// // void TAC_MP(T, int);
		// void TAC_MP(T, String);
		
	/// Old MP functions ///	
		
		// void Selctivity_MP(NumericVector, int);
		// Define as -
		//		template<class T>
		// 		HCR<T>::Selectivity_MP(NumericVector Sel, int n){ *MP block* }
		// Sel is selectivity-at-age (for Selectivity-based MP), n is MP selection number (if multiple are defined)
		
		// void Closure_MP(IntegerVector, int);
		// Define as -
		//		template<class T>
		// 		HCR<T>::Closure_MP(IntegerVector strata, int n){ *MP block* }
		// strata is strata number (ints) (for Closure-based MP), n is MP selection number (if multiple are defined)
		
		
		// All MP functions should be defined, but may just 'return' if unneeded //
		// For example, Closure_funtion(IntegerVector strata, int n) { return; }
		// 'int n' may also be used to choose between different types of MPs
		//		i.e. n = 0 for SSB, n = 1 for Selectivity, and n = 2 for Closure
		
		
	// };

}