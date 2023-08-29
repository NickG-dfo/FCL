
namespace PopSim{

	struct F_type{
		
		NumericVector _V;
		NumericMatrix _M;
		
		F_type(Rcpp::RObject O){
			if(Rf_isVector(O))
			{ _V = Rcpp::as<NumericVector>(O); }
			else if(Rf_isMatrix(O))
			{ _M = Rcpp::as<NumericMatrix>(O); }
			else
			{ std::printf("RObject not compatible with 'F_Type'."); }
		}
			
		F_type& operator= (NumericVector& vec){
			this->_V = vec;
			return *this;
		}
		F_type& operator= (NumericMatrix& mat){
			this->_M = mat;
			return *this;
		}
		
		F_type& operator= (Rcpp::RObject O){
			if(Rf_isVector(O)){ this->_V = Rcpp::as<NumericVector>(O); }
			if(Rf_isMatrix(O)){ this->_M = Rcpp::as<NumericMatrix>(O); }
			return *this;
		}
		
		F_type& operator= (F_type ft){
			ft._M = this->_M;
			ft._V = this->_V;
			return *this;
		}
		
		int size(void) { return _V.size(); } //Will return 0 if matrix
		
	};


	//---------------------------------------------------------------------------//


	//This objective is effectively selectivity (f-at-age and error)
	template<class T>
	class F_obj{
		
		int y;
		int _A;
		NumericVector Mu,
					  N_Sigma,
					  //Error+Sel values (replaced each year)
					  Fa,
					  //Annual F values (from HCR)
					  Fy;
		NumericMatrix MVN_Sigma;
		bool mvn, e, ee;
		 
		public:
		
			NumericVector Sel;
			NumericMatrix Fya;
			
			F_obj(NumericVector Mu, F_type Std, 
				  bool Error, bool exp_Error):
				Mu(Mu), e(Error), ee(exp_Error){
				_A = Mu.size();
				mvn = Std.size(); //size() will be 0 if a matrix
				N_Sigma = Std._V;
				MVN_Sigma = Std._M;
			}
			
			void operator() (T Fi){
							
				if(mvn){
					Fa = mvrnorm(Mu, MVN_Sigma);
				}else if(!mvn){
					Fa = rnorm(Mu, N_Sigma);
				}
				if(ee) { Fa = Rcpp::exp( Fa ); }
				
				if( max(Sel) > 1. ){ Sel / max(Sel); } //Set max sel to 1				
				if(e) { Fa = Fa * Sel; }
				else  { Fa = Sel; }
				
				Fy.push_back(Fi);
				Fya.push_back(Fi * Fa);
				
			}
		
	};

}
