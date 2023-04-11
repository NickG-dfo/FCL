
namespace SR{

	//For steepness formulations, parms(0) is 'h' and parms(1) is 'S0'
	//For shepherd and SigmoidBH, parms(2) is 'c'
	
	template<class T>
	T Constant(T SSB, const NumericVector &parms){
		return parms(0);
	}
	
	template<class T>
	T BevertonHolt1(T SSB, const NumericVector &parms){
		return SSB * parms(0) / (parms(1) + SSB);
	}
	
	template<class T>
	T BevertonHolt2(T SSB, const NumericVector &parms){
		return SSB * parms(0) / (1. + parms(1) * SSB);
	}
	
	template<class T>
	T BevertonHolth(T SSB, const NumericVector &parms, T SPR0){
		T h = parms(0),
		  R0 = parms(1);
		return .8*R0*h*SSB / (.2*SPR0*R0*(1.-h) + (h-.2)*SSB);
	}
	
	template<class T>
	T BevertonHoltComp(T SSB, const NumericVector &parms, T SPR0){
		T h = parms(1) / (4. + parms(1));
		NumericVector tempparm = {h, parms(1)};
		return BevertonHolth<T>(SSB, tempparm, SPR0);
	}
	
	template<class T>
	T Ricker(T SSB, const NumericVector &parms){
		return parms(0) * SSB * std::exp( -parms(1) * SSB );
	}
	
	template<class T>
	T Rickerh(T SSB, const NumericVector &parms, T SPR0){
		T b = std::log(5. * parms(0)) / (.8 * parms(1)),
		  a = std::exp(b * parms(1) / SPR0 );
		NumericVector tempparm = {a, b};
		return Ricker<T>(SSB, tempparm);
	}
	
	template<class T>
	T RickerComp(T SSB, const NumericVector &parms, T SPR0){
		T h = std::pow(.2 * parms(0), .8);
		NumericVector tempparm = {h, parms(1)};
		return Rickerh<T>(SSB, tempparm, SPR0);
	}
	
	template<class T>
	T SigmoidBevertonHolt(T SSB, const NumericVector &parms){
		return parms(0) * SSB / (1. + std::pow(parms(1) / SSB), parms(2));
	}
	
	template<class T>
	T Shepherd(T SSB, const NumericVector &parms){
		return parms(0) * SSB / (1. + std::pow(SSB / parms(1)), parms(2));
	}
	
	template<class T>
	T Shepherdh(T SSB, const NumericVector &parms, T SPR0){
		T b = parms(1) * ((.2 - parms(0)) / (parms(0) * std::pow(.2, parms(2)) - .2)) * (-1./parms(2)),
		  a = (1. + std::pow(parms(1)/b), parms(2)) / SPR0;
		return a * SSB / (1. + std::pow(SSB / b), parms(2));
	}
	
	template<class T>
	T Schaefer(T SSB, const NumericVector &parms){
		return parms(0) * SSB * (1. - parms(1) * SSB);
	}
	
	template<class T>
	T DerisoShnute(T SSB, const NumericVector &parms){
		return parms(0) * SSB * std::pow(1. - parms(2)*parms(1)*SSB, 1./parms(2));
	}
	
// }

	template<class T>
	class R_obj{
					
			String rec_opt;
			bool kernel_call;
			T recruitment;// (const &T, const &NumericVector);
			
			int yi;
			NumericVector logRy;
			NumericMatrix cov, kernel;
			NumericVector parms;
			T std_log_r;
			
			R_obj(const String &Rec_option, 
				  const NumericVector &Parameters, 
				  const T &R_std, const int Y): rec_opt(Rec_option) {
				
				yi = 0;
				std_log_r = R_std;
				
					 if( str_is(Rec_option, "Const") ){ this->recruitment = SR::Constant<T>; }
				else if( str_is(Rec_option, "BH1") ){ this->recruitment = SR::BevertonHolt1<T>; }
				else if( str_is(Rec_option, "BH2") ){ this->recruitment = SR::BevertonHolt2<T>; }
				else if( str_is(Rec_option, "RK") ){ this->recruitment = SR::Ricker<T>; }
				else if( str_is(Rec_option, "BHh") ){ this->recruitment = SR::BevertonHolth<T>; }
				else if( str_is(Rec_option, "RKh") ){ this->recruitment = SR::Rickerh<T>; }
				else if( str_is(Rec_option, "BHComp") ){ this->recruitment = SR::BevertonHoltComp<T>; }
				else if( str_is(Rec_option, "RKComp") ){ this->recruitment = SR::RickerComp<T>; }
				else if( str_is(Rec_option, "SigBH") ){ this->recruitment = SR::SigmoidBevertonHolt<T>; }
				else if( str_is(Rec_option, "Shepherd") ){ this->recruitment = SR::Shepherd<T>; }
				else if( str_is(Rec_option, "Shepherdh") ){ this->recruitment = SR::Shepherdh<T>; }
				else if( str_is(Rec_option, "Schaefer") ){ this->recruitment = SR::Schaefer<T>; }
				else if( str_is(Rec_option, "DerisoShnute") ){ this->recruitment = SR::DerisoShnute<T>; }
				
				else{
					std::printf("No valid recruitment option selected. Setting to default (constant) recruitment.");
					this->rec_opt = "Const";
					this->recruitment = SR::Constant<T>; 
				}
				
			}
			
			NumericMatrix derive_VCoV(const NumericMatrix &Kernel){
				
				// Kernel dimensions should be parms as columns & samples as rows
				
				int i, j, k;
				const int ncol = Kernel.ncol(),
						  nrow = Kernel.nrow();
				NumericVector Mu(ncol);
				NumericMatrix inner(nrow, ncol),
							  vcov(ncol, ncol);
							
				for(j = 0; j < ncol; j++){
					Mu(j) = mean(Kernel(j, _));
					for(i = 0; i < nrow; i++){
						inner(i, j) = Kernel(i, j) - Mu(j);
					}
				}
				for(k = 0; k < ncol; k++){
					for(j = 0; j < ncol; j++){
						if(k > j) {
							vcov(k, j) = vcov(j, k);
							break;
						}
						for(i = 0; i < nrow; i++){						
							vcov(k, j) += inner(i, j) * Rcpp::transpose(inner)(k, i);
						}
						vcov(k, j) /= nrow;
					}
				}
				
				return vcov;
				
			}
		
		public:
		
			NumericVector get_Parms(void){ return parms; }
			
			NumericMatrix get_CoV(void){ return cov; }
			
			NumericMatrix get_Kernerl(void){ return kernel; }
			
			void set_CoV(const NumericMatrix &CoV){
				if( (CoV.nrow() != cov.nrow()) | (CoV.ncol() != cov.ncol()) ){
					std::printf("Input CoV dimensions do not match Recruitment parameter vector-size.");
					return;
				}
				
				cov = CoV;
				kernel_call = 0;
			}
			
			void find_CoV(const NumericMatrix &Kernel){
				if( Kernel.ncol() != cov.ncol() ){
					std::printf("Input kernel row-size does not match Recruitment parameter vector-size.");
					return;
				}
				
				cov = derive_VCoV(Kernel);
				kernel_call = 0;
			}
			
			void use_Kernel(const NumericMatrix &Kernel){
				if( Kernel.ncol() != parms.size() ){
					std::printf("Input kernel row-size does not match Recruitment parameter vector-size.");
					return;
				}		
				
				kernel_call = 1;
				kernel = Kernel;
			}
			
			T logRec(T S, bool p_err = 1, bool r_err = 1){
				
				//parameter error		
				NumericVector inp_parms( parms.size() );
				if(p_err){
					if(kernel_call){ 
						inp_parms = kernel( (int)R::runif(0, kernel.rows()), _ );		
					}else{
						inp_parms = mvrnorm(parms, cov);
					}
				}else{
					inp_parms = parms;
				}
				//recruitment
				T rec = std::log( recruitment(S, inp_parms) );
				//normal error
				rec += R::rnorm(0., std_log_r) * r_err;
				
				logRy(yi) = rec; 
				yi++;
				
				return rec;
				
			}
			
			T Rec(T S, bool p_err = 1, bool r_err = 1){
				
				return std::exp( logRec(S, p_err, r_err) );
				
			}
			
			// Could add other loglinear/AR1 effects to rec, but for now parm error is fine
		
	};
	
} //End Namespace


// ------------- Old Rec Functions ------------- //

// Recruitment functions //
//functions return logRec because it makes implementation easier
//make sure all parms are defined in OM for recruitment!//

template<class T>
T BH(T SSB, const NumericVector &parms){
	
	return parms(0) + std::log(SSB) - std::log(1. + std::exp(parms(1)) * SSB);
	
}

template<class T>
T SigBH(T SSB, const T &max_S, const NumericVector &parms){
		
	//&parms = c(Rinf, S50, c), which are extracted from kernel in SRR f(x)
	T max_R = parms(0) / ( 1. + std::pow(parms(1)/max_S, parms(2)) );
	max_R = std::exp( std::log(max_R) - std::pow(parms(3), 2.)/2. ); // bias-correction
	
	T Rec = parms(0) / ( 1. + std::pow(parms(1)/SSB, parms(2)) ) * (SSB < max_S) + max_R * (SSB >= max_S);
	Rec = std::log(Rec) - std::pow(parms(3), 2.)/2.; // bias-correction
	
	return Rec;
	
}

NumericVector call_kernel(const NumericMatrix &kernel){
	
	// default_random_engine gen;
	// uniform_int_distribution<int> index(0, kernel.rows()-1);
	
	// return kernel(index(gen), _); 
	
	int ind = (int)R::runif(0, kernel.rows());
	return kernel(ind, _);
	
}

template<class T>
NumericMatrix SRR(const List &OM){
		
	NumericMatrix out(2, 4);
	
	NumericVector bh_parms = OM["BH_parms"];
	NumericMatrix bh_cov = OM["BH_cov"];
	NumericVector bh_est = mvrnorm(bh_parms, bh_cov);
	
	NumericMatrix kernel = OM["sigmoid_kernel"];
	T bhsigerr = OM["bhsigerr"];
	NumericVector sig_est = call_kernel(kernel);
	
	out(0, 0) = bh_est(0); out(0, 1) = bh_est(1); //BH parms
	out(1, 0) = sig_est(0); out(1, 1) = sig_est(1); out(1, 2) = sig_est(2); out(1, 3) = bhsigerr; //sigBH parms
	
	return out;
	
}