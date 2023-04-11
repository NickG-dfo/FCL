
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


namespace LaA{
	
	//Standard VonBertalanffy
	template<class T>
	T VonB(const T age, const T Linf, const T k,
			const T age0 = 0., const T L0 = 0.){
		return Linf * (1. - std::exp(-k * age - age0)) + L0;
	}
	
	template<class T>
	NumericVector VonB(const NumericVector ages, const T Linf, const T k, 
						const T age0 = 0., const T L0 = 0.){
		return Linf * (1. - Rcpp::exp(-k * ages - age0)) + L0;
	}
	
	template<class T>
	T Inv_VonB(const T Length, const T Linf, const T k, 
				const T age0 = 0., const T L0 = 0.){
		return (std::log(1. - (Length - L0)/Linf) + age0) / (-k);
	}
	
	template<class T>
	NumericVector Inv_VonB(const NumericVector Lengths, const T Linf, const T k,
							const T age0 = 0., const T L0 = 0.){
		return (Rcpp::log(1. - (Lengths - L0)/Linf) + age0) / (-k);
	}
	
	//Temperature VonBertalanffy
	
	// template<class T>
	// T VonB(const int age, const T Linf, const T k, const T Temp, 
			// bool plusGr = 1, const T age0 = 0., const T L0 = 0.){
		// return Linf * (1. - std::exp(-k * age - age0) + L0;
	// }
	
	// NumericVector VonB(const NumericVector ages, const T Linf, const T k, const T Temp, 
						// bool plusGr = 1, const T age0 = 0., const T L0 = 0.){
		// return Linf * (1. - Rcpp::exp(-k * ages - age0) + L0;
	// }
	
	// template<class T>
	// T Inv_VonB(const T Length, const T Linf, const T k, const T Temp,
				// bool plusGr = 1, const T age0 = 0., const T L0 = 0.){
		// return (std::log(1. - (Length - L0)/Linf) + age0) / (-k);
	// }
	
	// NumericVector Inv_VonB(const NumericVector Lengths, const T Linf, const T k, const T Temp,
							// bool plusGr = 1, const T age0 = 0., const T L0 = 0.){
		// return (Rcpp::log(1. - (Lengths - L0)/Linf) + age0) / (-k);
	// }
	
}


namespace SurvSim{
	
	template<class T>
	NumericVector Indices_A(const NumericVector &Ny, const NumericVector &qy){
		
		// const int A = Ny.size();
		// NumericVector Indices(A);
		
		// for(int a = 0; a < A; a++){
			// Indices(a) = Ny(a) * qy(a);
		// }
		return Ny * qy;
		
	}
	
	template<class T>
	NumericVector Indices_tA(const NumericVector &Ny, const NumericVector &qy,
							 const T &t_, const NumericVector &Zy){
		
		// const int A = Ny.size();
		// NumericVector Indices(A);
		
		// for(int a = 0; a < A; a++){
			// Indices(a) = Ny(a) * qy(a) * std::exp(-t_ * Zy(a));
		// }
		return Ny * qy * Rcpp::exp(-t_ * Zy);
		
	}
	
	template<class T>
	NumericMatrix Indices_YA(const NumericMatrix &Nya, const NumericMatrix &qya){
		
		const int A = Nya.ncol(),
				  Y = Nya.nrow();
		NumericMatrix Indices(Y, A);
		
		for(int y = 0; y < Y; y++){
			for(int a = 0; a < A; a++){
				Indices(y, a) = Nya(y, a) * qya(y, a);
			}
		}
		
		return Indices;
		
	}
	
	template<class T>
	NumericMatrix Indices_tYA(const NumericMatrix &Nya, const NumericMatrix &qya,
							  const NumericVector &t_y, const NumericMatrix &Zya){
		
		const int A = Nya.ncol(),
				  Y = Nya.nrow();
		NumericMatrix Indices(Y, A);
		
		for(int y = 0; y < Y; y++){
			for(int a = 0; a < A; a++){
				Indices(y, a) = Nya(y, a) * qya(y, a) * std::exp(-t_y(y) * Zya(y, a));
			}
		}
		
		return Indices;
		
	}
	
	template<class T>
	NumericMatrix3D Indices_RYA(NumericMatrix3D &Nrya, NumericMatrix3D &qrya){
		
		IntegerVector dims = Nrya.dims();
		const int A = dims(0),
				  Y = dims(1),
				  R = dims(2);
		NumericMatrix3D Indices(R, Y, A);
		
		for(int r = 0; r < R; r++){
			for(int y = 0; y < Y; y++){
				for(int a = 0; a < A; a++){
					Indices(r, y, a) = Nrya(r, y, a) * qrya(r, y, a);
				}
			}
		}
		
		return Indices;
		
	}
	
	template<class T>
	NumericMatrix3D Indices_tRYA(NumericMatrix3D &Nrya, NumericMatrix3D &qrya,
								 NumericMatrix &t_ry, NumericMatrix3D &Zrya){
		
		IntegerVector dims = Nrya.dims();
		const int A = dims(0),
				  Y = dims(1),
				  R = dims(2);
		NumericMatrix3D Indices(R, Y, A);
		
		for(int r = 0; r < R; r++){
			for(int y = 0; y < Y; y++){
				for(int a = 0; a < A; a++){
					Indices(r, y, a) = Nrya(r, y, a) * qrya(r, y, a) * std::exp(-t_ry(r, y) * Zrya(r, y, a));
				}
			}
		}
		
		return Indices;
		
	}
	
	NumericVector Selectivity(NumericVector Sel, 
							  NumericVector &F_err, 
							  NumericVector &scale){
	
		Sel = Sel + F_err;
		if(scale == R_NilValue) { scale = rep(1., Sel.size()); }
		Sel = Sel * scale;
		
		return Sel;
	
	}
	
	template<class T>
	NumericVector logN_a(//const int A, //this is the number of ages to include?
						const NumericVector &logN, 
						const NumericVector &Za,
						const NumericVector &PE, 
						const T logR
						//const String plusGr = "sum", // sum, trunc, or none
						// const String ProcessType = "expN" //expN, M, I
						){

		const int A = logN.size();

		NumericVector new_logN(A),
					  logNplus(A);
		
		for(int a = 0; a < A; a++){
		
			// Za(a) = Ma(a) + Fa(a);
			logNplus(a) = logN(a) - Za(a) + PE(a);
			
			if(a){ new_logN(a) = logNplus(a-1); }
			
		}
		
		new_logN(0) = logR;
		// if( str_is(plusGr, "sum") ){
			new_logN(A-1) = std::log( std::exp(logNplus(A-2)) + std::exp(logNplus(A-1)) );
		// }else if( str_is(plusGr, "trunc") ){
			// new_logN(A-1) = logNplus(A-2);
		// }
		
		return new_logN;

	}
	
	template<class T>
	NumericMatrix logN_ra(
						const NumericMatrix &logN, 
						// const NumericMatrix &Mra, 
						// const NumericMatrix &Fra, 
						const NumericMatrix &Zra, 
						const NumericMatrix &PE,
						const NumericVector logR
						// const String plusGr = "sum", // sum, trunc, or none
						// const String ProcessType = "expN" //expN, M, I
						){

		const int R = logN.nrow(),
				  A = logN.ncol();

		NumericMatrix new_logN(R, A),
					  logNplus(R, A);
					  
		for(int r = 0; r < R; r++){
			
			for(int a = 0; a < A; a++){
			
				// Za(r, a) = Ma(r, a) + Fa(r, a);
				logNplus(r, a) = logN(r, a) - Zra(r, a) + PE(r, a);
				
				if(a){ new_logN(r, a) = logNplus(r, a-1); }
				
			}
		
			//Add regional mixing eventually
			
			new_logN(r, 0) = logR(r);
			new_logN(r, A-1) = std::log( std::exp(logNplus(r, A-2)) + std::exp(logNplus(r, A-1)) );
			
		}
		
		return new_logN;

	}
	
}