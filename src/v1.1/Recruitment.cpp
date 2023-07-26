
namespace PopSim{
	
	//Used for calling function depending on NULL type or Function type
	// union RFunType{
		
		// bool Nil;
		// Rcpp::Function Fun;
		
		// RFunType(void) { Nil = 0; }
		
		// RFunType& operator= (Rcpp::Function& call_input){
			// Fun = call_input;
			// return *this;
		// }
		// RFunType& operator= (Rcpp::RObject call_input){
			// if(call_input == R_NilValue){
				// Nil = 1;
			// }else{
				// std::printf("Unusable value for RFunType.");
			// }
			// return *this;
		// }
		
		// double operator() (double SSB, NumericVector& parms){
				// if(Nil){ std::printf("Warning: No custom function available. Returning 0."); return 0.; }
				// return Fun( SSB, parms );
				// // return Fun( (Rcpp::RObject)SSB, (Rcpp::RObject)parms );
				// //May not use named elements for function call
		// }
		
	// };

	template<class T>
	class RP;

	template<class T>
	class SR{
		
		T SPR0;

		//For steepness formulations, parms(0) is 'h' and parms(1) is 'S0'
		//For shepherd and SigmoidBH, parms(2) is 'c'
		
		T Constant(T SSB, NumericVector &parms){
			return parms(0);
		}
		
		T BevertonHolt1(T SSB, NumericVector &parms){
			return SSB * parms(0) / (parms(1) + SSB);
		}
		
		T BevertonHolt2(T SSB, NumericVector &parms){
			return SSB * parms(0) / (1. + parms(1) * SSB);
		}
		
		T BevertonHolth(T SSB, NumericVector &parms){
			T h = parms(0),
			  R0 = parms(1);
			return .8*R0*h*SSB / (.2*SPR0*R0*(1.-h) + (h-.2)*SSB);
		}
		
		T BevertonHoltComp(T SSB, NumericVector &parms){
			T h = parms(1) / (4. + parms(1));
			NumericVector tempparm = {h, parms(1)};
			return BevertonHolth(SSB, tempparm);
		}
		
		T Ricker(T SSB, NumericVector &parms){
			return parms(0) * SSB * std::exp( -parms(1) * SSB );
		}
		
		T Rickerh(T SSB, NumericVector &parms){
			T b = std::log(5. * parms(0)) / (.8 * parms(1)),
			  a = std::exp(b * parms(1) / SPR0 );
			NumericVector tempparm = {a, b};
			return Ricker(SSB, tempparm);
		}
		
		T RickerComp(T SSB, NumericVector &parms){
			T h = std::pow(.2 * parms(0), .8);
			NumericVector tempparm = {h, parms(1)};
			return Rickerh(SSB, tempparm);
		}
		
		T SigmoidBevertonHolt(T SSB, NumericVector &parms){
			return parms(0) * SSB / (1. + std::pow(parms(1) / SSB), parms(2));
		}
		
		T Shepherd(T SSB, NumericVector &parms){
			return parms(0) * SSB / (1. + std::pow(SSB / parms(1)), parms(2));
		}
		
		T Shepherdh(T SSB, NumericVector &parms){
			T b = parms(1) * ((.2 - parms(0)) / (parms(0) * std::pow(.2, parms(2)) - .2)) * (-1./parms(2)),
			  a = (1. + std::pow(parms(1)/b), parms(2)) / SPR0;
			NumericVector tempparm = {a, b};
			return Shepherd(SSB, tempparm); // a * SSB / (1. + std::pow(SSB / b), parms(2));
		}
		
		T Schaefer(T SSB, NumericVector &parms){
			return parms(0) * SSB * (1. - parms(1) * SSB);
		}
		
		T DerisoShnute(T SSB, NumericVector &parms){
			return parms(0) * SSB * std::pow(1. - parms(2)*parms(1)*SSB, 1./parms(2));
		}
		
		T CustomRecruitment(T, NumericVector&);
		
		
		//				---				//
							
		bool kernel_call; //whether or not to use a kernel for parms
		int yi, //year, max of Y
			kernel_i; //kernel index being used
		NumericMatrix cov, kernel; //cov and kernel matrices (const)
		NumericVector parms; //vector of parms values for SRR
		
		T recruitment(T, NumericVector&); //function
		
		//Constructor args//
		String rec_opt;
		int Y; //max number of years in sim
		T std_log_r, max_rec;		
		bool bias_correction, //whether or not b-c is needed
			 rkseed, //whether or not to randomize seed each iter
			 perror, //whether or not to include parameter uncertainty
			 rerror; //whether or not to include error in mean recruitment estimates

		NumericVector logRy; // This is recruitment values for all years
		
		NumericMatrix derive_VCoV(NumericMatrix &Kernel){
			
			// Kernel dimensions should be: parms as columns & samples as rows
			
			int i, j, k;
			const int ncol = Kernel.ncol(),
					  nrow = Kernel.nrow();
			NumericVector Mu(ncol);
			NumericMatrix inner(nrow, ncol),
						  vcov(ncol, ncol);
						
			for(j = 0; j < ncol; ++j){
				Mu(j) = mean(Kernel(j, _));
				for(i = 0; i < nrow; ++i){
					inner(i, j) = Kernel(i, j) - Mu(j);
				}
			}
			for(k = 0; k < ncol; ++k){
				for(j = 0; j < ncol; ++j){
					if(k > j) {
						vcov(k, j) = vcov(j, k);
						break;
					}
					for(i = 0; i < nrow; ++i){						
						vcov(k, j) += inner(i, j) * Rcpp::transpose(inner)(k, i);
					}
					vcov(k, j) /= nrow;
				}
			}
			
			return vcov;
			
		}

		public:
			
			SR(String Rec_option, 
			   List &LHC, // Just for SPR0
			   int y, T R_std, T Max_R,
			   bool bc, bool RKseed,
			   bool p_err = 1, bool r_err = 1): 
				rec_opt(Rec_option), 
				Y(y), std_log_r(R_std), max_rec(Max_R),
				bias_correction(bc), rkseed(RKseed),
				perror(p_err), rerror(r_err)
				{
					
				//Define function and settings
				SPR0 = RP<T>(LHC, *this).SPR(0.); // Might work
				
				kernel_i = 0;
				if(rkseed) { kernel_i = (int)R::runif(0, kernel.rows()); }
				
				yi = 0;
				logRy = rep(0., Y);

					  if( str_is(Rec_option, "Const") ){ 
					  T recruitment(T SSB, NumericVector &parms){ return Constant(SSB, parms); }
				}else if( str_is(Rec_option, "BH1") ){ 
					T recruitment(T SSB, NumericVector &parms){ return BevertonHolt1(SSB, parms); }
				}else if( str_is(Rec_option, "BH2") ){ 
					T recruitment(T SSB, NumericVector &parms){ return BevertonHolt2(SSB, parms); }
				}else if( str_is(Rec_option, "RK") ){ 
					T recruitment(T SSB, NumericVector &parms){ return Ricker(SSB, parms); }
				}else if( str_is(Rec_option, "BHh") ){ 
					T recruitment(T SSB, NumericVector &parms){ return BevertonHolth(SSB, parms); }
				}else if( str_is(Rec_option, "RKh") ){ 
					T recruitment(T SSB, NumericVector &parms){ return Rickerh(SSB, parms); }
				}else if( str_is(Rec_option, "BHComp") ){ 
					T recruitment(T SSB, NumericVector &parms){ return BevertonHoltComp(SSB, parms); }
				}else if( str_is(Rec_option, "RKComp") ){ 
					T recruitment(T SSB, NumericVector &parms){ return RickerComp(SSB, parms); }
				}else if( str_is(Rec_option, "SigBH") ){ 
					T recruitment(T SSB, NumericVector &parms){ return SigmoidBevertonHolt(SSB, parms); }
				}else if( str_is(Rec_option, "Shepherd") ){ 
					T recruitment(T SSB, NumericVector &parms){ return Shepherd(SSB, parms); }
				}else if( str_is(Rec_option, "Shepherdh") ){ 
					T recruitment(T SSB, NumericVector &parms){ return Shepherdh(SSB, parms); }
				}else if( str_is(Rec_option, "Schaefer") ){ 
					T recruitment(T SSB, NumericVector &parms){ return Schaefer(SSB, parms); }
				}else if( str_is(Rec_option, "DerisoShnute") ){ 
					T recruitment(T SSB, NumericVector &parms){ return DerisoShnute(SSB, parms); }
				
				}else if( str_is(Rec_option, "Custom") ){  
				
					T recruitment(T SSB, NumericVector &parms){ return CustomRecruitment(SSB, parms); } 
					
				}else{
					std::printf("No valid recruitment option selected. Setting to default (constant) recruitment.");
					rec_opt = "Const";
					T recruitment(T SSB, NumericVector &parms){ return Constant(SSB, parms); }
				}
				
			}
		
			NumericVector get_Parms(void){ return parms; }
			
			NumericMatrix get_CoV(void){ return cov; }
			
			NumericMatrix get_Kernerl(void){ return kernel; }
			
			void set_Parms(NumericVector &Parms){
				parms = Parms;
				kernel_call = 0;
			}
			
			void set_CoV(NumericMatrix &CoV){ cov = CoV; }
			
			void find_CoV(NumericMatrix &Kernel){
				cov = derive_VCoV(Kernel);
				NumericVector temp_parms(Kernel.ncol());
				for(int i = 0; i < Kernel.ncol(); i++){
					temp_parms(i) = Rcpp::median(Kernel(_, i)); 
					//what if we want to use mean
				}
				parms = temp_parms;
				kernel_call = 0;
			}
			
			void set_Kernel(NumericMatrix &Kernel){
				kernel = Kernel;
				kernel_call = 1;
			}
			
			T logRec(T S){
				
				T rec;
				
				//Parameter Error		
				NumericVector inp_parms( parms.size() );
				if(perror){
					if(kernel_call){ 
						inp_parms = kernel( (int)R::runif(0, kernel.rows()), _ );		
					}else{
						inp_parms = mvrnorm(parms, cov);
					}
				}else{
					if(kernel_call){ 
						inp_parms = kernel( kernel_i, _ );	// kernel with p_err = F will always use first row in kernel
					}else{
						inp_parms = parms;
					}
				}
				
				//Recruitment Estimate
				// if( is_in(rec_opt, SRRs_require_SPR0) ){ //check is rec_opt is in StringVector of SRR names
					// rec = std::log( recruitment(S, inp_parms, spr0) ) // this spr0 is a nuisance
					// //Find a better way to do this...
				// }else{
				rec = std::log( recruitment(S, inp_parms) );
				// }
				
				//Recruitment Error
				rec += R::rnorm(0., std_log_r) * rerror;
				
				if(rec > max_rec) { rec = max_rec; }
				if(bias_correction) { rec -= std_log_r*std_log_r / 2.; }
		
				logRy(yi) = rec;
				yi++;				
				
				return rec;
				
			}
			
			void clear(void) { yi = 0; }
			
			String Rec_Type(void) { return rec_opt; }
			
			NumericVector getRec(void){ return logRy; }
			
			// Could add other loglinear/AR1 effects to rec, but for now parm error is fine
		
	}; //End Class

}
