
// Simple Functions //

bool str_is(String String1, String String2){
	
	std::string string1 = String1.get_cstring();
	std::string string2 = String2.get_cstring();
	if( string1.size() != string2.size() ) { return false; }
	for(unsigned int i = 0; i < string1.size(); i++){
		if( string1[i] != string2[i] ) { return false; }
	}
	return true;
	
}

int which(const StringVector vec, const String val){
	
	int i = 0;
	for(; i < vec.size(); i++){
		if( vec(i) == val ) { break; }
	}
	return i;
	
}

template<class T>
NumericVector scale(NumericVector vec, T k){
	
	for(int i = 0; i < vec.size(); i++){
		vec(i) *= k;
	}
	
	return vec;
	
}

IntegerVector seq(int a, int b, int d = 1){

	int range = (b - a)/d;
	IntegerVector vec(range+1);
	for(int i = 0; i <= range; i++){
		vec(i) = a + i*d;
	}
	
	return vec;
	
}

template<class T>
NumericVector seq(T a, T b, T d = 1.){

	int range = (int)(b - a)/d;
	NumericVector vec(range);
	for(int i = 0; i < range; i++){
		vec(i) = a + i*d;
	}
	
	return vec;
	
}

template<class T>
NumericVector cat(const T x, NumericVector vec){
	
	vec.push_front(x);
	return vec;
	
}

IntegerVector cat(const int x, const IntegerVector vec){
	
	IntegerVector nvec(vec.size()+1);
	for(int i = 1; i < vec.size()+1; i++){
		nvec(i) = vec(i-1);
	}
	// vec.push_front(x);	
	nvec(0) = x;
	return nvec;
	
}

StringVector cat(const String x, const StringVector vec){
	
	StringVector nvec(vec.size()+1);
	for(int i = 1; i < vec.size()+1; i++){
		nvec(i) = vec(i-1);
	}
	// vec.push_front(x);	
	nvec(0) = x;
	return nvec;
	
}

template<class T>
NumericVector cat(NumericVector vec, const T x){
	
	vec.push_back(x);	
	return vec;
	
}

IntegerVector cat(const IntegerVector vec, const int x){
	
	int len = vec.size();
	IntegerVector nvec(len+1);
	for(int i = 0; i < len; i++){
		nvec(i) = vec(i);
	}
	// vec.push_front(x);	
	nvec(len) = x;
	return nvec;
	
}

StringVector cat(const StringVector vec, const String x){
	
	int len = vec.size();
	StringVector nvec(len+1);
	for(int i = 0; i < len; i++){
		nvec(i) = vec(i);
	}
	// vec.push_front(x);	
	nvec(len) = x;
	return nvec;
	
}

NumericVector cat(const NumericVector &vec1, const NumericVector &vec2){
	
	int len1 = vec1.size(), len2 = vec2.size();
	NumericVector vecn(len1 + len2);
	
	for(int i = 0; i < len1+len2; i++){
		if(i < len1){
			vecn(i) = vec1(i);
		}else{
			vecn(i) = vec2(i-len1);
		}		
	}
	
	return vecn;
	
}

IntegerVector cat(const IntegerVector &vec1, const IntegerVector &vec2){
	
	int len1 = vec1.size(), len2 = vec2.size();
	IntegerVector vecn(len1 + len2);
	
	for(int i = 0; i < len1+len2; i++){
		if(i < len1){
			vecn(i) = vec1(i);
		}else{
			vecn(i) = vec2(i-len1);
		}		
	}
	
	return vecn;
	
}

StringVector cat(const StringVector &vec1, const StringVector &vec2){
	
	int len1 = vec1.size(), len2 = vec2.size();
	StringVector vecn(len1 + len2);
	
	for(int i = 0; i < len1+len2; i++){
		if(i < len1){
			vecn(i) = vec1(i);
		}else{
			vecn(i) = vec2(i-len1);
		}		
	}
	
	return vecn;
	
}

template<class T>
NumericVector rep(const T x, const int n){
	
	NumericVector vec(n);
	for(int i = 0; i < n; i++){
		vec(i) = x;
	}
	
	return vec;
	
}

IntegerVector rep(const int x, const int n){
	
	IntegerVector vec(n);
	for(int i = 0; i < n; i++){
		vec(i) = x;
	}
	
	return vec;
	
}

StringVector rep(const String x, const int n){
	
	StringVector vec(n);
	for(int i = 0; i < n; i++){
		vec(i) = x;
	}
	
	return vec;
	
}

NumericVector rep(const NumericVector x, int n, bool at = false){
	
	NumericVector vec(0);
	if(!at){
		for(int i = 0; i < x.size(); i++){ //this is like each = n, repeating each val at a time
			vec = cat(vec, rep(x(i), n));
			//this is a nested loop -> O(n^2) time; optimize??
		}	
	}else if(at){
		for(int i = 0; i < x.size(); i++){ // this is repeating all values in order they're given
			vec = cat(vec, x);
			//this is a nested loop -> O(n^2) time; optimize??
		}
	}
	
	return vec;
	
}

IntegerVector rep(const IntegerVector x, int n, bool at = false){
	
	IntegerVector vec(0);
	if(!at){
		for(int i = 0; i < x.size(); i++){ //this is like each = n, repeating each val at a time
			vec = cat(vec, rep((int)x(i), n));
			//this is a nested loop -> O(n^2) time; optimize??
		}	
	}else if(at){
		for(int i = 0; i < x.size(); i++){ // this is repeating all values in order they're given
			vec = cat(vec, x);
			//this is a nested loop -> O(n^2) time; optimize??
		}
	}
	
	return vec;
	
}

StringVector rep(const StringVector x, int n, bool at = false){
	
	StringVector vec(0);
	if(!at){
		for(int i = 0; i < x.size(); i++){ //this is like each = n, repeating each val at a time
			vec = cat(vec, rep((String)x(i), n));
			//this is a nested loop -> O(n^2) time; optimize??
		}
	}else if(at){
		for(int i = 0; i < n; i++){
			vec = cat(vec, x); // this is repeating all values in order they're given
		}
	}
	
	return vec;
	
}

NumericVector as_vector(const NumericMatrix &M){
	
	int J = M.cols(), I = M.rows();
	NumericVector V(J*I);
	for(int i = 0; i < I; i++){
		for(int j = 0; j < J; j++){
			V(j+i*J) = M(i, j);
		}
	}
	
	return V;
	
}


// -------------------------------------------------------------------------------------------------//


// MSE-Lite Functions //

NumericVector baranov_catch(const NumericVector &Na, const NumericVector &Fa, const NumericVector &Ma){
	
	short A = Na.size();
	NumericVector Za(A), EC(A);
	
	for(int i = 0; i < A; i++){
		
		Za(i) = Fa(i) + Ma(i);
		EC(i) = Na(i) * (1. - std::exp(-Za(i))) * Fa(i)/Za(i);
		
	}

	return EC;
	
}

template<class T>
T yield(NumericVector catches, const NumericVector &Cw){
	
	for(int i = 0; i < catches.size(); i++){
		catches(i) *= Cw(i);
	}
	
	return sum(catches);
	
}

template<class T>
T solveF(T TAC, const NumericVector &sel, const NumericVector &Na, const NumericVector &Ma, const NumericVector &Cw){
	
	short niter = 7, A = Na.size();
	T F = 0.1, Catch, d_Catch;
	NumericVector Za(A), Ca(A), d_Ca(A);
	
	for(short n = 0; n < niter; n++){
			
		for(short i = 0; i < A; i++){
			
			Za(i) = F * sel(i) + Ma(i);
			Ca(i) = F * sel(i) * (Na(i) * (1. - std::exp(-Za(i))) * Cw(i) / Za(i)); 
					// Catch weight at age
			d_Ca(i) = Na(i) * sel(i) * Cw(i) * ( (1. - std::exp(-Za(i))) * Ma(i) / Za(i) + std::exp(-Za(i)) * F * sel(i) ) / Za(i); 
					// Derivative of catch weight at age
			
		}
		
		Catch = sum(Ca);
		d_Catch = sum(d_Ca);
		
		// if(d_Catch == 0.){
			// F = 0.001;
		// }else{
			F += (TAC - Catch) / d_Catch;
		// }
		
	}
	
	if(F > 1.){ F = 1.; }
	if(F < 0.001){ F = 0.001; }
	
	return F;
	
}


// Armadillo Stuff //

//This class is to convert arma 3Dcube to output in Rcpp
class NumericMatrix3D{
	
		const int x, y, z;
		arma::cube Storage3D;
		// Note: arma::cube is a form of the generic type arma::Cube<Type> where <Type> is <double>
		
	public:
			
		NumericMatrix3D(const int ix, const int iy, const int iz): x(ix), y(iy), z(iz) {
			Storage3D.reshape(x, y, z);
		}
	
	
		IntegerVector dims(void){ return {x, y, z}; }
	
		NumericMatrix slice_z(int depth){
			NumericMatrix out = Rcpp::wrap( Storage3D.slice(depth) );
			return out;			
		}
		NumericMatrix slice_y(int height){
			NumericMatrix out = Rcpp::wrap( Storage3D.col(height) );
			return out;			
		}
		NumericMatrix slice_x(int width){
			NumericMatrix out = Rcpp::wrap( Storage3D.row(width) );
			return out;			
		}
		
		NumericVector trim_xy(int width, int height){
			return slice_x(width)(height, _);
		}
		NumericVector trim_yx(int height, int width){
			return slice_x(width)(height, _);
		}
		NumericVector file_xz(int width, int depth){
			return slice_x(width)(_, depth);
		}
		NumericVector file_zx(int depth, int width){
			return slice_x(width)(_, depth);
		}
		NumericVector rank_yz(int height, int depth){
			return slice_z(depth)(height, _);
		}
		NumericVector rank_zy(int depth, int height){
			return slice_z(depth)(height, _);
		}
		
		double& operator() (int width, int height, int depth){
			double *ret = &( Storage3D(width, height, depth) );
			return *ret;
		}
		
		NumericMatrix3D& operator= (NumericMatrix3D M){
			IntegerVector Dims = M.dims();
			for(int i = 0; i < Dims(0); i++){
				for(int j = 0; j < Dims(1); j++){
					for(int k = 0; k < Dims(2); k++){
						(*this)(i, j, k) = M(i, j, k);
					}
				}
			}
			return *this;
		}
		
};

//This function converts arma rmvnorm to a Rcpp equivalent
NumericVector mvrnorm(const NumericVector &mu, const NumericMatrix &S){

	const unsigned int n = 1;
	arma::vec amu = Rcpp::as<arma::vec>(mu);
	arma::mat aS = Rcpp::as<arma::mat>(S);
	arma::mat x = rmvnorm(n, amu, aS);
	NumericMatrix X = Rcpp::wrap(x);
	NumericVector out = X(0, _);
	return out;

}



// Auto-correlation class/functions //

template<class T>
T randWalk (T X, T &Std){ // Random Walk - takes two doubles
	return X + R::rnorm(0., Std); // X = X_n-1, returns X_n
}

NumericVector randWalk (NumericVector X, NumericVector &Std){ // Random Walk - takes two vectors of equal length
	for(int i = 0; i < X.size(); i++){
		X(i) += R::rnorm(0., Std(i));
	}
	return X; // X = X_n-1, returns X_n
}

template<class T>
struct AR1d{
	
	//Does AR(1) have a constant coefficient??//
		const int n;
		AR1d(const int size): n(size) {
			if(n < 1) { n = 1; }
		}
		
		int nforward;
		NumericVector Y;
		void forward(int __n){
			nforward += __n;
			Y.resize(nforward);
		}
		
		//Numeric Vector should be length of n, so AR(1) X is a double in vector form
		double operator() (NumericVector &X, T &Phi){
			
			if(Parms.size() > n) { std::printf("Number or parameters exceed requirement for AR(%d) process; some parameters will be ignored.", n); }
			if(Parms.size() < n) {
				std::printf("Required number of parameters for AR(%d) process exceed number provided; using '0' as default(s).", n);
				Parms = cat( Parms, (NumericVector)rep(0., n-Parms.size()) );
			}
			
			if(X.size() > n) { std::printf("Number or prior values exceed requirement for AR(%d) process; some inputs will be ignored.", n); }
			if(X.size() < n) { 
				std::printf("Required number of prior values for AR(%d) process exceed number provided; using '0' as default(s).", n);
				X = cat( X, (NumericVector)rep(0., n-X.size()) );
			}
		
			for(int i = 0; i < nforward; i++){
				if(i < n){
					Y(i) = X(i);
				}else{
					Y(i) += Phi * Y(i-n); //No error in AR process	
				}
			}
			return Y;
		}
		
		//Maybe remove error conmponent internally and just add externally if Error_type == 'AR'

};

// struct AR2d{
	
		// const int n, m;
	
		// AR2d(const int rows, const int cols): n(rows), m(cols) {;};
		
		// //Numeric Vector should be length of n, so AR(1) X is a double in a vector object
		// double operator() (NumericMatrix &X, NumericMatrix &Parms, const double intercept = 0.){
			
			// NumericMatrix Parmss(n+1, m+1);
			// if(Parms.nrow() > n+1 || Parms.ncol() > m+1) {
				// std::printf("Number or parameters exceed requirement for AR(%d) process; some parameters will be ignored.", n);
			// }
			// if(Parms.nrow() < n+1 || Parms.ncol() < m+1) {
				// std::printf("Required number of parameters for AR(%d) process exceed number provided; using '0' as default(s).", n);
				// // Parms = cat( Parms, (NumericVector)rep(0., n+1-Parms.size()) );
			// }
			
			// if(X.nrow() > n || X.ncol() > m) { 
				// std::printf("Number or prior values exceed requirement for AR(%d) process; some inputs will be ignored.", n);
				// }
			// if(X.nrow() < n || X.ncol() < m) { 
				// std::printf("Required number of prior values for AR(%d) process exceed number provided; using '0' as default(s).", n);
				// // X = cat( X, (NumericVector)rep(0., n-X.size()) );
			// }
		
			// NumericMatrix Y(n, m);// = mvrnorm( rep(0., Parms.nrows()), Parms );
			// for(int i = 0; i < n; i++){
					// Y += Parms(i) * X(i);
			// }
			// return Y + intercept;
		// }

// };