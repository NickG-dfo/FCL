 
template<class T>
class M_obj{
	
		int a, y;
		const int Y, A, pn;
		const IntegerVector Ages;
				
		T Mbase;
		NumericVector Mstd, // single std for Merr (might allow multiple uncertainties for each component)
							//If so, name the vectors for easy access to each component
					  Mcond, // year effect
					  Maeffect, // loglinear age affects
					  Myeffect; // loglinear year affects
		NumericMatrix Mcorr, // matrix of AR parms
					  Mage,
					  Merr;
		// Mcorr uses first value for AR1d and uses full matrix for AR2d
			
		const String Mtype;
		const String Error; // "Cond", "Temp", "Age", "Year", "Base" //Consider making this a vector for each component
		
	public:
		
		M_obj(const Rcpp::List &MInfo, const int Y_in, const int pn_in): 
			Y(Y_in), pn(pn_in) {
						
			A = MInfo["A"];
			Ages = MInfo["Ages"];
			Mbase = MInfo["base"];
			Mstd = MInfo["std"]; 
			Mcond = MInfo["condition"];
			Maeffect = MInfo["age_effect"];
			Myeffect = MInfo["year_effect"];
			Mcorr = MInfo["correlates"];
			Mtype = MInfo["mtype"];
			Error = MInfo["error_type"];
			
			for(y = 0; y < Y; y++){
				for(a = 0; a < A; a++){
					Merr(y, a) = R::rnorm(0., (T)Mstd["base"]);
					Mage(y, a) = std::log(Mbase) + Merr(y, a) * ( str_is(Error, "Base") );
				}
			}
			
		}
		
		void MCondition(int y){
			
			// Remember to define Maeffect as NumericVector Mest(Mind) in List! //
			// And remember to call MAgeEffect before condition, since condition applies to age effect! //
			Maeffect = Maeffect * ( Mcond(y) + Merr(y, _) * ( str_is(Error, "Cond") ) );
			
		}
		
		void MAgeEffect(int y){
			
			Mage(y, _) = Mage(y, _) + Maeffect + Merr(y, _) * ( str_is(Error, "Age") );
			
		}
		
		void MYearEffect(int y){
			
			Mage(y, _) = Mage(y, _) + Myeffect(y) + Merr(y, _) * ( str_is(Error, "Year") );
		
		}
		
		void MAR1(int y, const String corrtype = "year"){ // corrtype = "year", "age", "yearage"
			
			AR1d<T> AR(1);
			// AR2d ARya(1);
			// If AR1 process is ON, might have to do it first before other M components
			if( str_is(corrtype, "age") ){
				for(a = 1; a < A; a++){
					Mage(y, a) = AR( Mage(y, a-1), Mcorr(0, _) );
				}
			}else if( str_is(corrtype, "year") ){
				for(y = 1; y < Y; y++){
					Mage(y, a) = AR( Mage(y-1, a), Mcorr(0, _) );
				}
			// }else if( str_is(corrtype, "yearage") ){
				// for(y = 1; y < Y; y++){
					// for(a = 0; a < A; a++){
						// Mage(y, a) = ARya( Mage, Mcorr );
					// }
				// }
			}
			
			// AR1 process requires a prior value to simulate?
			// Maybe function takes 1 value and propogates from there? Ask Divya
			
			// Also maybe construct AR1/M values at start of sim and pull from matrix throughout sim
		
		}
		
		void MError(int y){
			
			//if ( Error == R_NilValue )
			for(a = 0; a < A; a++){
				Mage(y, a) += Merr(y, a);
			}
			
		}
	
		NumericMatrix logM(void){
			return Mage;
		}
	
		T M_ya(int y, int a){
			return std::exp( Mage(y, a) );
		}
		
		NumericVector M_y(int y){
			return Rcpp::exp( Mage(y, _) );
		}
	
};
