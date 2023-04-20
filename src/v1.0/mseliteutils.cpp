
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



//quantile function

template<class T>
T quantile(NumericVector X, T prob){
	
	X = X.sort();
	
	T h = ((T)X.size() + .25)*prob + 3./8.;
	NumericVector hb = Rcpp::floor( NumericVector {h} ),
				  ht = Rcpp::ceil( NumericVector {h} );
	T Qp = X(hb(0)) + (h - (T)hb(0)) * (X(ht(0)) - X(hb(0)));
		
	return Qp;
	// NumericVector Y = Rcpp::qnorm();
	// return (T)Y(0);
}



// -------------------------------------------------------------------------------------------------//


// Armadillo Stuff //

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


// -------------------------------------------------------------------------------------------------//


// -- // MSELite Functions // -- //

// Catch & F Functions //

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


// Abundance Function //

template<class T>
NumericVector calc_logN(const int A, const NumericVector &logN, const NumericVector &Ma, 
						const NumericVector &Fa, const NumericVector &PE, T logR){

	NumericVector new_logN(A);
	NumericVector Za(A), logNplus(A);
	
	for(int i = 0; i < A; i++){
	
		Za(i) = Ma(i) + Fa(i);
		logNplus(i) = logN(i) - Za(i) + PE(i);
		
		if( (i != 0) & (i != A-1) ){
			new_logN(i) = logNplus(i-1);
		}
		
	}
	
	new_logN(0) = logR;
	new_logN(A-1) = std::log( std::exp(logNplus(A-2)) + std::exp(logNplus(A-1)) );
	
	return new_logN;

}


// Define Function for MP HCRs (for 3Ps RP) //

template<class T>
struct HCR{
	
	T TAC;
	T F;
	HCR& operator= (const HCR tf){
		this->TAC = tf.TAC;
		this->F = tf.F;
		return *this;
	}
	
};

template<class T>
HCR<T> MP_decision(//T F, T Catch, 
				// const T Ubound, const T Fbound,
				// const T harvestrate, const T threshold,
				// NumericVector ssbr_bracket, const NumericVector f_bracket,
				// NumericVector fbreak,
				const String MP_name,
				const NumericVector &Na, const NumericVector &Ma, const NumericVector &sel,
				const NumericVector &mat, const NumericVector &Sw, const NumericVector &Cw){
	
	// declare values for MPs
	HCR<T> MP_out;
	short A = Na.size();
	// T Lbound = 66.;
	NumericVector temp_catch(A), temp_f(A);
	
	T ssb = 0.;	
	for(short i = 0; i < A; i++){
		ssb += Na(i) * mat(i) * Sw(i);
	}
	
	//define TAC/F based on MP_name
	if(MP_name == "No fishing"){
		MP_out.TAC = 0.; //Fix this
		MP_out.F = 0.001;
	}
	else if(MP_name == "Critical 100"){
		MP_out.F = 0.001 +
					( (ssb > 66.) & (ssb <= 1.75*66.) ) * (ssb - 66.)/(.75*66.) * .1 +
					(ssb > 1.75*66.) * .1;
		
		temp_f = scale(sel, MP_out.F);
		temp_catch = baranov_catch(Na, temp_f, Ma);
		MP_out.TAC = yield<T>(temp_catch, Cw);
	}
	else if(MP_name == "Elbow C"){
		T m1 = (.035-.025)/(.6*66.);
		T b1 = .035-m1*66.;
		T m2 = (.1-.035)/(66.);
		T b2 = .1-m2*(2.*66.);
		MP_out.F = 0.001 +
				  ( (ssb > .4*66.) & (ssb <= 66.) ) * (m1*ssb+b1) +
				  ( (ssb > 66.) & (ssb <= 2.*66.) ) * (m2*ssb+b2) +
				  (ssb > 2.*66.) * .1;
		temp_catch = baranov_catch(Na, scale(sel, MP_out.F), Ma);
		MP_out.TAC = yield<T>(temp_catch, Cw);
	}
	else if(MP_name == "AGC"){
		T m1 = (.2-.05)/(66.);
		T b1 = .2-m1*(2.*66.);
		MP_out.F = 0.001 +
				  ( (ssb > .4*66.) & (ssb <= 66.) ) * .05 +
				  ( (ssb > 66.) & (ssb <= 2.*66.) ) * (m1*ssb+b1) +
				  (ssb > 2.*66.) * .2;
		temp_catch = baranov_catch(Na, scale(sel, MP_out.F), Ma);
		MP_out.TAC = yield<T>(temp_catch, Cw);
	}
	else if(MP_name == "FFAW"){
		T m1 = (3.-1.)/(.4*66.);
		T b1 = 3.-m1*1.2*66.;
		T m2 = (20.-3.)/(.8*66.);
		T b2 = 20.-m2*(2.*66.);
		MP_out.TAC = 0.001 +
				  ( (ssb > .4*66.) & (ssb <= .8*66.) ) * 1. +
				  ( (ssb > .8*66.) & (ssb <= 1.2*66.) ) * (m1*ssb+b1) +
				  ( (ssb > 1.2*66.) & (ssb <= 2.*66.) ) * (m2*ssb+b2) +
				  (ssb > 2.*66.) * 20.;
		MP_out.F = solveF(MP_out.TAC, sel, Na, Ma, Cw);
	}
	
	// else if(MP_name == "harvestrate"){
		// T surplus = (ssb - threshold) * harvestrate;
		// MP_out.TAC = surplus > 0. ? surplus : 0.;
		// MP_out.F = solveF(MP_out.TAC, sel, Ma, Na, Cw);
	// }
	// else if(MP_name == "f_bracket"){
		// bool which;
		// short i = 0;
		// ssbr_bracket = cat( cat(0., ssbr_bracket), 100.);
		// for(; i < ssbr_bracket.size(); i++){
			// which = ( (ssb/66. >= ssbr_bracket(i)) & (ssb/66. < ssbr_bracket(i+1)) );
			// if(which) { break; }
		// }
		// MP_out.F = f_bracket(i);
		
		// temp_catch = baranov_catch(Na, scale(sel, MP_out.F), Ma);
		// MP_out.TAC = yield<T>(temp_catch, Cw);
	// }
	// else if(MP_name == "Elbow_PA"){
		// T m1 = (0.05 - 0.025)/(0.6 * 66.);
		// T b1 = 0.05 - m1 * 66.;
		// T m2 = (0.26 - 0.05)/66.;
		// T b2 = 0.26 - m2 * (2*66.);
		// MP_out.F =  ( (ssb > 0.4 * 66.) & (ssb <= 66.) ) * (m1 * ssb + b1) +
					// ( (ssb > 66.) & (ssb <= 2 * 66.) ) * (m2 * ssb + b2) +
					// (ssb > 2 * 66.) * 0.26;
					
		// temp_catch = baranov_catch(Na, scale(sel, MP_out.F), Ma);
		// MP_out.TAC = yield<T>(temp_catch, Cw);
	// }
	// else if(MP_name == "Elbow_Spectrum1"){
		// T m1 = (0.045 - 0.055)/(0.6 * 66.);
		// T b1 = 0.055 - m1 * 66.;
		// T m2 = (0.2 - 0.055)/66.;
		// T b2 = 0.2 - m2 * (2*66.);
		// MP_out.F =  ( (ssb > 0.4 * 66.) & (ssb <= 66.) ) * (m1 * ssb + b1) +
					// ( (ssb > 66.) & (ssb <= 2 * 66.) ) * (m2 * ssb + b2) +
					// (ssb > 2 * 66.) * 0.2;
					
		// temp_catch = baranov_catch(Na, scale(sel, MP_out.F), Ma);
		// MP_out.TAC = yield<T>(temp_catch, Cw);
	// }
	// else if(MP_name == "Elbow_Spectrum2"){
		// T m1 = (0.07 - 0.03)/(0.6 * 66.);
		// T b1 = 0.07 - m1 * 66.;
		// T m2 = (0.2 - 0.07)/66.;
		// T b2 = 0.2 - m2 * (2*66.);
		// MP_out.F =  ( (ssb > 0.4 * 66.) & (ssb <= 66.) ) * (m1 * ssb + b1) +
					// ( (ssb > 66.) & (ssb <= 2 * 66.) ) * (m2 * ssb + b2) +
					// (ssb > 2 * 66.) * 0.2;
					
		// temp_catch = baranov_catch(Na, scale(sel, MP_out.F), Ma);
		// MP_out.TAC = yield<T>(temp_catch, Cw);
	// }
	
	//return struct object with TAC and F values
	
	return MP_out;
	
}
