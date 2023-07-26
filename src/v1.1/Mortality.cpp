 
 namespace PopSim{
 
	template<class T>
	class M_obj{
		
			int a, y;
			int Y, A;
			IntegerVector Ages;
					
			//Error containers (assuming constant error throughout sim)
			// bool err_base,
			  // err_cond,
			  // err_age,
			  // err_year;
			NumericVector Mbase,
						  Mcond, // loglinear year effect applied to all ages
						  //Age and Year effects are numeric vectors of length A and Y, respectively
						  //but if they are only length 1 they will be the same for each age/year
						  Maeffect, // loglinear age affects (length A)
						  Myeffect; // loglinear year affects (length Y)
						  //Error container for ARp
						  // err_corr;
			NumericMatrix Mcorr; // 3x3 matrix of AR parms
								 // Mcorr uses (0, 0) value for AR1d
								 // Mcorr uses upper-left matrix for AR2d
								 // Mcorr uses full matrix for AR3d
			List Mstd;  //Error effects -- list of vectors for each Error type across ages
			NumericVector Mbase_std,
						  Mcond_std,
						  Maeffect_std,
						  Myeffect_std;
			NumericMatrix Mcorr_std;
						  //(can set to 0. if uneeded, or set Error value to 0)
				
			String corrtype; //this will be either: a, y, r, ya, ra, ry, rya
			
			LogicalVector Error; // Named Binary vector: "condition", "age_effect", "year_effect", "base", "correlates"
			
		public:
		
			NumericMatrix Mya; // Internal storage for M values
			
			M_obj(const Rcpp::List &MInfo, const int Y_in): 
				Y(Y_in){
							
				A = MInfo["A"];
				Ages = MInfo["Ages"]; //This might be unnecessary
				Mbase = MInfo["base"]; //This should be a vector of length A
				
				Mstd = MInfo["std"]; //This should be a vector  with same length/names as Error
				Mbase_std = Mstd["base"];
				Mcond_std = Mstd["condition"];
				Maeffect_std = Mstd["age_effect"];
				Myeffect_std = Mstd["Myear_effect"];
				
				NumericMatrix tempM1 = Mstd["Mcorr_std"],
							  tempM2 = MInfo["correlates"];
				Mcorr_std = tempM1;
				Mcorr = tempM2;
				
				Mcond = MInfo["condition"];
				Maeffect = MInfo["age_effect"];
				Myeffect = MInfo["year_effect"];
				
				String tempS = MInfo["correlation_type"];
				corrtype = tempS;
				
				// Error = MInfo["error_type"]; //This should be a Named bool vector
				// //Error assignment
				// err_base = Error["base"];
				// err_cond = Error["condition"];
				// err_age = Error["age_effect"];
				// err_year = Error["year_effect"];
				// err_corr = Error["correlates"];
				
				Mya = NumericMatrix(Y, A);
				for(y = 0; y < Y; y++){
					// for(a = 0; a < A; a++){
						//Mbase should be input as a standard value, but internally its stored as log value
						Mya(y, _) = Rcpp::log(Mbase) + rnorm( rep(0., A), Mbase_std );
					// }
				}
				
			}
			
			void MCondition(int y){
				
				// Call MAgeEffect after condition, since condition applies to age effect //
				Maeffect = Maeffect * ( Mcond + rnorm( rep(0., A), Mcond_std ) );
				
			}
			
			void MAgeEffect(int y){
				
				Mya(y, _) = Mya(y, _) + Maeffect +  rnorm( rep(0., A), Maeffect_std );
				
			}
			
			void MYearEffect(int y){
				
				Mya(y, _) = Mya(y, _) + Myeffect(y) +  R::rnorm( 0., Myeffect_std(y) );
			
			}
			
			void MAR(int n = 1){
				
				AR1d AR(n);
				// AR2d AR2(n);
				// AR3d ARrya(3)
				NumericMatrix Sigma = Mcorr_std;
				
				T Mcorr_a = Mcorr(0, 0),
				  Mcorr_y = Mcorr(1, 1),
				  Mcorr_r = Mcorr(2, 2);
				NumericMatrix Mcorr_ya = sub( Mcorr, {0, 1}, {0, 1} ),
							  Mcorr_ry = sub( Mcorr, {1, 2}, {1, 2} ),
							  Mcorr_ra = sub( Mcorr, {0, 2}, {0, 2} ),
							  Mcorr_rya = Mcorr;
				
				// If AR1d process is ON, might have to do it first before other M components
				if( str_is(corrtype, "a") ){
					// for(a = 1; a < A; a++){						
						Mya(y, _) = AR( Mya(y, _), rep((T)Mcorr_a, A), (T)Sigma(0,0) );
					// }
				}else if( str_is(corrtype, "y") ){
					// for(y = 1; y < Y; y++){
						Mya(_, a) = AR( Mya(_, a), rep((T)Mcorr_y, A), (T)Sigma(1,1) );
					// }
				}else if( str_is(corrtype, "r") ){
						Mya(_, a) = AR( Mya(_, a), rep((T)Mcorr_r, A), (T)Sigma(2,2) );
				}
				// }else if( str_is(corrtype, "ya") ){
					// for(y = 1; y < Y; y++){
						// for(a = 0; a < A; a++){
							// Mya(y, a) = ARya( Mya, Mcorr_ya, Sigma(1,0) * err_corr );
						// }
					// }
				// }
				
				// Also maybe construct AR1/M values at start of sim and pull from matrix throughout sim
				
				//AR3d process (or anything with spatial component) might not be necessary
				//Do we model different M values in space? Probably not
			
			}
		
	};

}