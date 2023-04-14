 
 namespace PopSim{
 
	template<class T>
	class M_obj{
		
			int a, y;
			const int Y, A;
			const IntegerVector Ages;
					
			//Error containers (assuming constant error throughout sim)
			T err_base,
			  err_cond,
			  err_age,
			  err_year;
			NumericVector Mstd, // Multiple uncertainties for each M component (can set to 0. if uneeded, or set Error value to 0)
						  Mcond, // loglinear year effect applied to all ages
						  //Age and Year effects are numeric vectors of length A and Y, respectively
						  //but if they are only length 1 they will be the same for each age/year
						  Maeffect, // loglinear age affects (length A)
						  Myeffect, // loglinear year affects (length Y)
						  //Error container for ARp
						  err_corr;
			NumericMatrix Mcorr; // 3x3 matrix of AR parms
								 // Mcorr uses (0, 0) value for AR1d
								 // Mcorr uses upper-left matrix for AR2d
								 // Mcorr uses full matrix for AR3d
								 
						  Mage; // Internal storage for M values
				
			const String corrtype; //this will be either: a, y, r, ya, ra, ry, rya
			
			const LogicalVector Error; // Named Binary vector: "condition", "age_effect", "year_effect", "base", "correlates"
			
		public:
			
			M_obj(const Rcpp::List &MInfo, const int Y_in): 
				Y(Y_in){
							
				A = MInfo["A"];
				Ages = MInfo["Ages"]; //This might be unnecessary
				Mbase = MInfo["base"]; //This should be a vector of length A
				Mstd = MInfo["std"]; //This should be a vector  with same length/names as Error
				Mcond = MInfo["condition"];
				Maeffect = MInfo["age_effect"];
				Myeffect = MInfo["year_effect"];
				Mcorr = MInfo["correlates"];
				Error = MInfo['errors']
				// Mtype = MInfo["mtype"];
				
				Error = MInfo["error_type"]; //This should be a Named bool vector
				//Error assignment
				err_base = Error["base"];
				err_cond = Error["condition"];
				err_age = Error["age_effect"];
				err_year = Error["year_effect"];
				err_corr = Error["correlates"];
				
				for(y = 0; y < Y; y++){
					for(a = 0; a < A; a++){
						//Mbase should be input as it standard value, but internally its stored as it log value
						Mage(y, a) = Rcpp::log(Mbase) + Rcpp::rnorm( A, 0., Mstd["base"] ) * err_base;
					}
				}
				
			}
			
			void MCondition(int y){
				
				// Remember to define Maeffect as NumericVector Mest(Mind) in List! //
				// And remember to call MAgeEffect after condition, since condition applies to age effect! //
				Maeffect = Maeffect * ( Mcond(y) +  R::rnorm( 0., Mstd["condition"] ) * err_cond );
				
			}
			
			void MAgeEffect(int y){
				
				Mage(y, _) = Mage(y, _) + Maeffect +  R::rnorm( 0., Mstd["age_effect"] ) * err_age;
				
			}
			
			void MYearEffect(int y){
				
				Mage(y, _) = Mage(y, _) + Myeffect(y) +  R::rnorm( 0., Mstd["year_effect"] ) * err_year;
			
			}
			
			void MAR(int n = 1, const String corrtype = "year"){
				
				AR1d<T> AR(n);
				AR2d<T> AR2(n);
				// AR3d<T> ARrya(3)
				T Sigma = Mstd["correlates"];
				
				T Mcorr_a = Mcorr(0, 0),
				  Mcorr_y = Mcorr(1, 1),
				  Mcorr_r = Mcorr(2, 2);
				NumericMatrix Mcorr_ya = sub( Mcorr, {0, 1}, {0, 1} ),
							  Mcorr_ry = sub( Mcorr, {1, 2}, {1, 2} ),
							  Mcorr_ra = sub( Mcorr, {0, 2}, {0, 2} ),
							  Mcorr_rya = Mcorr;
				
				// If AR1d process is ON, might have to do it first before other M components
				if( str_is(corrtype, "a") ){
					for(a = 1; a < A; a++){
						
						Mage(y, a) = AR( Mage(y, _), Mcorr_a, Sigma * err_corr )
					}
				}else if( str_is(corrtype, "y") ){
					for(y = 1; y < Y; y++){
						Mage(y, a) = AR( Mage(_, a), Mcorr_y, Sigma * err_corr );
					}
				}else if( str_is(corrtype, "ya") ){
					for(y = 1; y < Y; y++){
						for(a = 0; a < A; a++){
							Mage(y, a) = ARya( Mage, Mcorr_ya, Sigma * err_corr );
						}
					}
				}
				// Also maybe construct AR1/M values at start of sim and pull from matrix throughout sim
				
				//AR3d process (or anything with spatial component) might not be necessary
				//Do we model different M values in space? Probably not
			
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

}