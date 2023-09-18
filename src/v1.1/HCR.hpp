
namespace PopSim{

	//Create class for reference point functions
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
			
			RP(List &LHC, SR<T> Recruitment): Rec(Recruitment)
			{				
				A = LHC["A"];
				Mage = LHC["Mort"];
				Sage = LHC["Sel"];
				Mat = LHC["Mat"];
				Stock_wt = LHC["Sw"];
				Catch_wt = LHC["Cw"];
				
				Feps = 1E-6;
				Seps = 1E-3;
			}
						
			T SPR(T F)
			{
				NumericVector M = Mage, 
							  S = Sage;
				M.push_front(0.); M.erase(A);
				S.push_front(0.); M.erase(A);
				NumericVector cZ = cumsum(M+S*F);
				T spr = sum( (exp(-cZ) ) * Mat * Stock_wt );
				return spr;
			}
			
			T YPR(T F)
			{
				NumericVector M = Mage, 
							  S = Sage;
				M.push_front(0.); M.erase(A);
				S.push_front(0.); M.erase(A);
				NumericVector Z = M+F*S,
							  cZ = cumsum(Z);
				
				T ypr = sum( (1. - exp(-cZ)) * (F/Z) * Catch_wt );
				return ypr;
			}
			
			T SSBeq(T F)
			{
				if( str_is( Rec.Rec_Type, "Const" ) ){
					std::printf("Constant recruitment does not allow equilbirum states.");
					return Rec->recruitment((T)0., Rec.parms);
				}
				 // Newton-Raphson root
				T delta = 1.,
				  Si = 1.;
				T y, dy;
				do{
					y = ((Rec->recruitment)(Si, Rec.parms) * SPR(F) - Si);
					dy = ((Rec->recruitment)(Si+Seps, Rec.parms)*SPR(F) - (Rec->recruitment)(Si-Seps, Rec.parms)) * SPR(F)/(2*Seps) - 1.;
					delta = dy/y;
					Si -= delta;
				} while(delta > Seps);
				return Si;
			}
			
			T Yeq(T F)
			{
				return SSBeq(F) * SPR(F);
			}
			
			T F_SPRx(T percent)
			{ // Newton-Raphson root
				T delta = 1.,
				  Fi = 1.;
				T y, dy;
				do{
					y = SPR(Fi) - percent*SPR(0.);
					dy = (SPR(Fi+Feps) - SPR(Fi-Feps))/(2*Feps);
					delta = dy/y;
					Fi -= delta;
				} while(delta > Feps);
				return Fi;
			}			
						
			F_B0x(T percent)
			{
				T delta = 1.,
				  Fi = 1.;
				T y, dy;
				do{
					y = SSBeq(Fi) - percent*SSBeq(0);
					dy = (SSBeq(Fi+Feps) - SSBeq(Fi-Feps))/(2*Feps);
					delta = dy/y;
					Fi -= delta;
				} while(delta > Feps);
				return Fi;
			}
			
			T FMSY(void)
			{ // Gradient Descent
				T delta = 1.,
				  gamma = .05,
				  Fi = 1.;
				T y, dy;
				do{
					// y = Yeq(Fi);
					dy = (Yeq(Fi+Feps) - Yeq(Fi-Feps))/(2*Feps);
					// ddy = (Yeq(log_F+2.*Feps, Recruitment, parms, Sel, Mort, Sw, Cw, Mat) - 
							 // 2.*Y +
							 // Yeq(log_F-2.*Feps, Recruitment, parms, Sel, Mort, Sw, Cw, Mat) 
							// )/(4.*Feps*Feps);
					delta = gamma*dy;
					Fi += delta;
				} while(delta > Feps);
				return Fi;
			}
			
			T BMSY(void)
			{
				return SSBeq( FMSY() );
			}
			
			T MSY(void)
			{
				return Yeq( FMSY() );
			}
		
	};

}