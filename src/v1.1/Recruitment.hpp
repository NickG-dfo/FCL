
namespace PopSim{

	template<class T>
	class RP;

	template<class T>
	class SR{
		
		T SPR0;	
		
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
		
			NumericVector parms; //vector of parameter values for SRR
			T (*recruitment)(T, NumericVector&); //function
			void Rassign(T (*f)(T, NumericVector&))
			{
				this->recruitment = f;
			}
			
			SR(String& Rec_option, 
			   List &LHC, // Just for SPR0
			   int y, T& R_std, T& Max_R,
			   bool& bc, bool& RKseed,
			   bool& p_err = 1, bool& r_err = 1): 
				rec_opt(Rec_option), 
				Y(y), std_log_r(R_std), max_rec(Max_R),
				bias_correction(bc), rkseed(RKseed),
				perror(p_err), rerror(r_err)
				{
					
				//Define function and settings
				SPR0 = RP<T>(LHC, *this).SPR(0.); // Might work
				
				kernel_i = (int)R::runif(0, kernel.rows());
				if(rkseed) { kernel_i = (int)R::runif(0, kernel.rows()); }
				
				yi = 0;
				logRy = rep(0., Y);

					  if( str_is(Rec_option, "Const") ){ 
					// recruitment = this->Constant;
					Rassign(Constant);
				}else if( str_is(Rec_option, "BH1") ){ 
					recruitment = this->BevertonHolt1;
				}else if( str_is(Rec_option, "BH2") ){ 
					recruitment = this->BevertonHolt2;
				}else if( str_is(Rec_option, "RK") ){ 
					recruitment = this->Ricker;
				}else if( str_is(Rec_option, "BHh") ){ 
					recruitment = this->BevertonHolth;
				}else if( str_is(Rec_option, "RKh") ){ 
					recruitment = this->Rickerh;
				}else if( str_is(Rec_option, "BHComp") ){ 
					recruitment = this->BevertonHoltComp;
				}else if( str_is(Rec_option, "RKComp") ){ 
					recruitment = this->RickerComp;
				}else if( str_is(Rec_option, "SigBH") ){ 
					recruitment = this->SigmoidBevertonHolt;
				}else if( str_is(Rec_option, "Shepherd") ){ 
					recruitment = this->Shepherd;
				}else if( str_is(Rec_option, "Shepherdh") ){ 
					recruitment = this->Shepherdh;
				}else if( str_is(Rec_option, "Schaefer") ){ 
					recruitment = this->Schaefer;
				}else if( str_is(Rec_option, "DerisoShnute") ){ 
					recruitment = this->DerisoShnute;
				
				}else if( str_is(Rec_option, "Custom") ){
					recruitment = this->CustomRecruitment; 
				}else{
					std::printf("No valid recruitment option selected. Setting to default (constant) recruitment.");
					rec_opt = "Const";
					recruitment = this->Constant;
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
					temp_parms(i) = median( (NumericVector)Kernel(_, i) );
					//what if we want to use mean
				}
				parms = temp_parms;
				kernel_call = 0;
			}
			
			void set_Kernel(NumericMatrix &Kernel){
				kernel = Kernel;
				kernel_call = 1;
			}
			
			T logRec(T S)
			{
				
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
						inp_parms = kernel( kernel_i, _ );	// kernel with p_err = F will always use same row of kernel
					}else{
						inp_parms = parms;
					}
				}
				
				//Recruitment Estimate
				rec = std::log( recruitment(S, inp_parms) );
				
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
