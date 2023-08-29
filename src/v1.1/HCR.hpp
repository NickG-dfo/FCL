
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
			T Feps, Seps;
		
		public:	
			
			RP(List &LHC, SR<T> Recruitment): Rec(Recruitment){				
				A = LHC["A"];
				Mage = LHC["Mort"];
				Sage = LHC["Sel"];
				Mat = LHC["Mat"];
				Stock_wt = LHC["Sw"];
				Catch_wt = LHC["Cw"];
				
				Feps = 1E-6;
				Seps = 1E-3;
			}
						
			SPR(T F){
				NumericVector M = cat(0., Mage), 
							  S = cat(0., Sage),
							  cZ = cumsum(M - S * F);
				T spr = sum( (Rcpp::exp(-cZ)) * Mat * Stock_wt );
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
					y = ((*Rec->recruitment)(Si, Rec.parms) * SPR(F) - Si);
					dy = ((*Rec->recruitment)(Si+Seps, Rec.parms)*SPR(F) - (*Rec->recruitment)(Si-Seps, Rec.parms) * SPR(F))/(2*Seps) - 1.;
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
					y = SPR(Fi) - percent*SPR(0.);
					dy = (SPR(Fi+Feps) - SPR(Fi-Feps))/(2*Feps);
					delta = dy/y;
					Fi -= delta;
				} while(delta > .001);
				return Fi;
			}
			
						
			F_BMSYx(T percent)
			{
				T delta = 1.,
				  Fi = 1.;
				T y, dy;
				do{
					y = SSBeq(Fi) - percent*SSBeq(0);
					dy = (SSBeq(Fi+Feps) - SSBeq(Fi-Feps))/(2*Feps);
					delta = dy/y;
					Fi -= delta;
				} while(delta > .001);
				return Fi;
			}
			
			T FMSY(void){ // Newton-Raphson root
				T delta = 1.,
				  Fi = 1.;
				T y, dy;
				do{
					y = std::pow(Yeq(Fi), 2.); //y and dy are already simplified for delta
					dy = 2 * (Yeq(Fi+Feps) - Yeq(Fi-Feps))/(2*Feps);
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

}