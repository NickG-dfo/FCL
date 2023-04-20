
namespace PopSim{

	struct F_type{
		
		NumericVector _V;
		NumericMatrix _M;
			
		F_type& operator= (NumericVector& vec){
			this->_V = vec;
			return *this;
		}
		F_type& operator= (NumericMatrix& mat){
			this->_M = mat;
			return *this;
		}
		
		F_type operator= (F_type ft){
			ft._M = this->_M;
			ft._V = this->_V;
			return ft;
		}
		
		int size(void) { return _V.size(); } //Will return 0 if matrix
		
	};


	//---------------------------------------------------------------------------//


	//This objective is effectively selectivity (f-at-age and error)
	template<class T>
	class F_obj{
		
		int y;
		const int _A;
		NumericVector Mu,
					  N_sigma
					  //Error+Sel values (replaced each year)
					  Fa,
					  //Annual F values (from HCR)
					  Fy;
		NumericMatrix MVN_Sigma;
		bool mvn;
		 
		public:
		
			NumericVector Sel,
			NumericMatrix Fya;
			
			F_obj(NumericVector Mu, F_type Std):
				Mu(Mu){
				mvn = Std.size(); //size() will be 0 if a matrix
				if(mvn){
					N_sigma = Std._V;
				}else{
					mvn = 1;
					MVN_Sigma = Std._M;
				}
			}
			
			void operator() (bool Error, bool Exp_err){
							
				if(mvn){
					Fa = mvrnorm(MNV_Mu, MNV_Sigma);
				}else if(!mvn){
					Fa = Rcpp::rnorm(_A, N_mu, N_sigma);
				}
				if(Exp_err) { Fa = Rcpp::exp( Fa ); }
				
				if( max(Sel) > 1. ){ Sel / max(Sel); } //Set max sel to 1
				
				if(Error) { Fa *= Sel; }
				else	  { Fa = Sel; }
				
			}
			
			void F(T Fi){
				
				Fy.push_back(Fi);
				Fya.push_bacK( Fi * Fa ); //May not work
				
			}
			
			NumericVector F_y(int y){ return Fya(y, _); }
		
	};

}
