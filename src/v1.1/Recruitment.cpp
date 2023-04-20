
namespace PopSim{
	
	//Used for calling function depending on NULL type or Function type
	template<class T>
	union RFuncType{
		
		bool Nil;
		Rcpp::Function Fun;
		
		RFuncType(void) { Nil = 0; }
		
		RFunType& operator= (Rcpp::R_NilValue& call_input){
			Nil = 1;
		}
		RFunType& operator= (Rcpp::Function& call_input){
			Fun = call_input;
		}
		
		T operator() (T SSB, NumericVector& parms){
			if(Nil){ std::printf("Warning: No custom function available. Returning 0."); return 0.; }
			return Fun(SSB, parms);
			//May not use named elements for function call
		}
		
	} CallCustomFun ;

	template<class T>
	class RP;

	template<class T>
	class SR{
		
		T SPR0;

		//For steepness formulations, parms(0) is 'h' and parms(1) is 'S0'
		//For shepherd and SigmoidBH, parms(2) is 'c'
		
		T Constant(T SSB, const NumericVector &parms){
			return parms(0);
		}
		
		T BevertonHolt1(T SSB, const NumericVector &parms){
			return SSB * parms(0) / (parms(1) + SSB);
		}
		
		T BevertonHolt2(T SSB, const NumericVector &parms){
			return SSB * parms(0) / (1. + parms(1) * SSB);
		}
		
		T BevertonHolth(T SSB, const NumericVector &parms){
			T h = parms(0),
			  R0 = parms(1);
			return .8*R0*h*SSB / (.2*SPR0*R0*(1.-h) + (h-.2)*SSB);
		}
		
		T BevertonHolth(T, NumericVector);
		
		T BevertonHoltComp(T SSB, const NumericVector &parms){
			T h = parms(1) / (4. + parms(1));
			NumericVector tempparm = {h, parms(1)};
			return BevertonHolth<T>(SSB, tempparm);
		}
		
		T Ricker(T SSB, const NumericVector &parms){
			return parms(0) * SSB * std::exp( -parms(1) * SSB );
		}
		
		T Rickerh(T SSB, const NumericVector &parms){
			T b = std::log(5. * parms(0)) / (.8 * parms(1)),
			  a = std::exp(b * parms(1) / SPR0 );
			NumericVector tempparm = {a, b};
			return Ricker(SSB, tempparm);
		}
		
		T Rickerh(T, NumericVector);
		
		T RickerComp(T SSB, const NumericVector &parms){
			T h = std::pow(.2 * parms(0), .8);
			NumericVector tempparm = {h, parms(1)};
			return Rickerh<T>(SSB, tempparm);
		}
		
		T SigmoidBevertonHolt(T SSB, const NumericVector &parms){
			return parms(0) * SSB / (1. + std::pow(parms(1) / SSB), parms(2));
		}
		
		T Shepherd(T SSB, const NumericVector &parms){
			return parms(0) * SSB / (1. + std::pow(SSB / parms(1)), parms(2));
		}
		
		T Shepherdh(T SSB, const NumericVector &parms){
			T b = parms(1) * ((.2 - parms(0)) / (parms(0) * std::pow(.2, parms(2)) - .2)) * (-1./parms(2)),
			  a = (1. + std::pow(parms(1)/b), parms(2)) / SPR0;
			NumericVector tempparm = {a, b};
			return Shepherd(SSB, tempparm); // a * SSB / (1. + std::pow(SSB / b), parms(2));
		}
		
		T Shepherdh(T, NumericVector);
		
		T Schaefer(T SSB, const NumericVector &parms){
			return parms(0) * SSB * (1. - parms(1) * SSB);
		}
		
		T DerisoShnute(T SSB, const NumericVector &parms){
			return parms(0) * SSB * std::pow(1. - parms(2)*parms(1)*SSB, 1./parms(2));
		}
		
		CustomRecruitment(T SSB, NumericVector &parms) { return 0.; }
		
		
		//				---				//
						
		String rec_opt;
		bool kernel_call, //whether or not to use a kernel for parms
			 bias_correction, //whether or no b-c is needed
			 rkseed; //whether or not to randomize seed each iter
		T recruitment; //function
		
		const int Y; //number of sim years
		int yi, //year, max of Y
			kernel_i; //kernel index being used
		const NumericMatrix cov, kernel; //cov and kernel matrices (const)
		NumericVector parms; //vector of parms values for SRR
		const T std_log_r, max_rec;

		NumericVector logRy; // This is recruitment values for all years
		
		NumericMatrix derive_VCoV(const NumericMatrix &Kernel){
			
			// Kernel dimensions should be parms as columns & samples as rows
			
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
			
			SR(const String &Rec_option, const T &R_std, 
			   const List &LHC, // Just for SPR0
			   const int y, T Max_R, bool bc, bool RKseed): 
				rec_opt(Rec_option), std_log_r(R_std), Y(y),
				max_rec(Max_R), bc(bias_correction), rkseed(RKseed)
				{
					
				//Check for 'Function' input in Rec_option
				// if( Rec_option["Function"] == R_NilValue ){
					
				// }else if( Rf_isFunction(Rec_option["Function"]) ){
					// //this may not work because Rec_option["Function"] isn't assigned to a container
					// CustomRecruitment = Rec_option["Function"];
				// }		
				//OR... if prior assignment is needed
				CallCustomFun = Rec_option["Function"];
				//This should automatically detect and store either a NULL of Function type
				if( !CallCustomFun.Nil ){
					CustomRecruitment = CallCustomFun.Fun;
				}
				
					
				//Define function and settings
				SPR0 = RP<T>(LHC, this->recruitment).SPR(0.); // Might work
				
				kernel_i = 0;
				if(rkseed) { kernel_i = (int)R::runif(0, kernel.rows()) }
				
				yi = 0;
				logRy = rep(0., Y); //Might not work because logRy is in public below
				
					 if( str_is(Rec_option, "Const") ){ this->recruitment = Constant; }
				else if( str_is(Rec_option, "BH1") ){ this->recruitment = BevertonHolt1; }
				else if( str_is(Rec_option, "BH2") ){ this->recruitment = BevertonHolt2; }
				else if( str_is(Rec_option, "RK") ){ this->recruitment = Ricker; }
				else if( str_is(Rec_option, "BHh") ){ this->recruitment = BevertonHolth; }
				else if( str_is(Rec_option, "RKh") ){ this->recruitment = Rickerh; }
				else if( str_is(Rec_option, "BHComp") ){ this->recruitment = BevertonHoltComp; }
				else if( str_is(Rec_option, "RKComp") ){ this->recruitment = RickerComp; }
				else if( str_is(Rec_option, "SigBH") ){ this->recruitment = SigmoidBevertonHolt; }
				else if( str_is(Rec_option, "Shepherd") ){ this->recruitment = Shepherd; }
				else if( str_is(Rec_option, "Shepherdh") ){ this->recruitment = Shepherdh; }
				else if( str_is(Rec_option, "Schaefer") ){ this->recruitment = Schaefer; }
				else if( str_is(Rec_option, "DerisoShnute") ){ this->recruitment = DerisoShnute; }
				
				else if( str_is(Rec_option, "Custom") ){ this->recruitment = CustomRecruitment; }
				
				else{
					std::printf("No valid recruitment option selected. Setting to default (constant) recruitment.");
					this->rec_opt = "Const";
					this->recruitment = Constant; 
				}
				
			}
		
			NumericVector get_Parms(void){ return parms; }
			
			NumericMatrix get_CoV(void){ return cov; }
			
			NumericMatrix get_Kernerl(void){ return kernel; }
			
			void set_Parms(const NumericVector &Parms){
				parms = Parms;
				kernel_call = 0;
			}
			
			void set_CoV(const NumericMatrix &CoV){ cov = CoV; }
			
			void find_CoV(const NumericMatrix &Kernel){
				cov = derive_VCoV(Kernel);
				NumericVector temp_parms(Kernel.ncol());
				for(int i = 0; i < Kernel.ncol(); i++){
					temp_parms(i) = Rcpp::median(Kernel(_, i)); 
					//what if we want to use mean
				}
				parms = temp_parms;
				kernel_call = 0;
			}
			
			void set_Kernel(const NumericMatrix &Kernel){
				kernel = Kernel;
				kernel_call = 1;
			}
			
			// T logRec(T, bool, bool);
			T logRec(T S, bool p_err = 1, bool r_err = 1){
				
				T rec;
				
				//Parameter Error		
				NumericVector inp_parms( parms.size() );
				if(p_err){
					if(kernel_call){ 
						inp_parms = kernel( (int)R::runif(0, kernel.rows()), _ );		
					}else{
						inp_parms = mvrnorm(parms, cov);
					}
				}else{
					if(kernel_call){ 
						inp_parms = kernel( kernel_i, _ );	// kernel with no parm error will randomize parm
					}else{
						inp_parms = parms;
					}
				}
				
				//Recruitment Estimate
				if( is_in(rec_opt, SRRs_require_SPR0) ){ //check is rec_opt is in StringVector of SRR names
					rec = std::log( recruitment(S, inp_parms, spr0) ) // this spr0 is a nuisance
					//Find a better way to do this...
				}else{
					rec = std::log( recruitment(S, inp_parms) );
				}
				
				//Recruitment Error
				rec += R::rnorm(0., std_log_r) * r_err;
				
				if(rec > max_rec) { rec = max_rec; }
				if(bc) { rec -= std::pow(b_sigma, 2.) / 2.; }
		
				logRy(yi) = rec;
				++yi;				
				
				return rec;
				
			}
			
			T Rec(T S, bool p_err = 1, bool r_err = 1){
				
				return std::exp( logRec(S, p_err, r_err) );
				
			}
			
			void clear(void) { yi = 0; }
			
			String Rec_Type(void) { return rec_opt; }
			
			// Could add other loglinear/AR1 effects to rec, but for now parm error is fine
		
	}; //End Class

}
