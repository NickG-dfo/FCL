
// template<class T>
// NumericVector calc_logN(const int A, const NumericVector &logN, const NumericVector &Ma, 
						// const NumericVector &Fa, const NumericVector &PE, T logR){

	// NumericVector new_logN(A);
	// NumericVector Za(A), logNplus(A);
	
	// for(int i = 0; i < A; i++){
	
		// Za(i) = Ma(i) + Fa(i);
		// logNplus(i) = logN(i) - Za(i) + PE(i);
		
		// if( (i != 0) & (i != A-1) ){
			// new_logN(i) = logNplus(i-1);
		// }
		
	// }
	
	// new_logN(0) = logR;
	// new_logN(A-1) = std::log( std::exp(logNplus(A-2)) + std::exp(logNplus(A-1)) );
	
	// return new_logN;

// }


// SPMSim is called in SPM model ; PopSim is called in AgeModel //

namespace SPMSim{
	
	template<class T>
	struct SPMitem{
		
		//Note: Process Error will be implemented outside the item, in case unperturbed Biomass is needed
		
		int item;
		const NumericVector b, h, c, m;
		
		SPMitem(const T &B){
			b.push_back(B);
			item = 0; // item starts at 0 instead of 1
		}
		
		void Discrete(const NumericVector &Parms){
			
			T H = Parms["HarvestRate"],
			  M = Parms["MortalityRate"],
			  R = Parms["Recruitment"],
			  G = Parms["GrowthRate"];
			
			item++;
			h.push_back( H );
			m.push_back( M );
			c.push_back( H * b(item-1) );
			b.push_back( b(item-1) + R + G - C - M );
			
		}
		
		void Continuous(const NumericVector &Parms, int years = 1){
			
			T C = Parms["Catch"],
			  r = Parms["GrowthRate"],
			  K = Parms["CarryingCapacity"];
			
			for(int i = 0; i < years; i++){
				item++;
				c.push_back( C );
				b.push_back( r * b(item-1) * (1. - b(item-1)/K ) );
				h.push_back( C / b(item-1) );
				m.push_back( r/K * b(item-1) );
			}
			
		}
		
		int size(void){ return item + 1; }
		
		NumericVector Biomass(void){ return b; }
		NumericVector Catch(void){ return c; }
		NumericVector HarvestRate(void){ return h; }
		NumericVector Mortality(void){ return m; }

	};
	
}


namespace PopSim{
	
	//Survey Indices//
	
	template<class T>
	NumericVector Indices_A(const NumericVector &Ny, const NumericVector &qy){
		return Ny * qy;		
	}
	
	template<class T>
	NumericVector Indices_tA(const NumericVector &Ny, const NumericVector &qy,
							 const T &t_, const NumericVector &Zy){
		return Ny * qy * Rcpp::exp(-t_ * Zy);		
	}
	
	template<class T>
	NumericMatrix Indices_YA(const NumericMatrix &Nya, const NumericMatrix &qya){
		
		const int A = Nya.ncol(),
				  Y = Nya.nrow();
		NumericMatrix Indices(Y, A);
		
		for(int y = 0; y < Y; y++){
			Indices(y, _) = Nya(y, _) * qya(y, _);
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
			Indices(y, _) = Nya(y, _) * qya(y, _) * std::exp(-t_y(y) * Zya(y, _));
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
	
	//Abundances//
	
	template<class T>
	NumericVector logN_a(
						const NumericVector &logN, 
						const NumericVector &Za,
						const T &PE, 
						const T logR
						){

		const int A = logN.size();

		NumericVector new_logN(A);
					  
		new_logN = logN - Za + PE;
		
		new_logN(A-2) += new_logN(A-1);
		new_logN.push_front(logR);
		new_logN.erase(A);
		
		return new_logN;

	}
	
	template<class T>
	NumericVector logN_a(
						const NumericVector &logN, 
						const NumericVector &Za,
						const NumericVector &PE, 
						const T logR
						){

		const int A = logN.size();

		NumericVector new_logN(A);
		
		new_logN = logN - Za + PE;;
		
		new_logN(A-1) += new_logN(A-2);
		new_logN.push_front(logR);
		new_logN.erase(A);
		
		return new_logN;

	}
	
	template<class T>
	NumericMatrix logN_ra(
						const NumericMatrix &logN,
						const NumericMatrix &Zra, 
						const NumericVector &PE,
						const NumericVector logR
						){

		const int R = logN.nrow(),
				  A = logN.ncol();

		NumericMatrix new_logN(R, A);
					  
		for(int r = 0; r < R; r++){
			
			new_logN(r, _) = logN(r, _) - Zra(r, _) + PE(r);
			
			new_logN(r, A-1) += new_logN(r, A-2);
			new_logN(r, _).push_front(logR(r));
			new_logN(r, _).erase(A);
			
			//Add regional mixing eventually
			
		}
		
		return new_logN;

	}
	
	template<class T>
	NumericMatrix logN_ra(
						const NumericMatrix &logN,
						const NumericMatrix &Zra, 
						const NumericMatrix &PE,
						const NumericVector logR
						){

		const int R = logN.nrow(),
				  A = logN.ncol();
					  
		NumericMatrix new_logN(R, A);
					  
		for(int r = 0; r < R; r++){
			
			new_logN(r, _) = logN(r, _) - Zra(r, _) + PE(r, _);
			
			new_logN(r, A-1) += new_logN(r, A-2);
			new_logN(r, _).push_front(logR(r));
			new_logN(r, _).erase(A);
			
			//Add regional mixing eventually
			
		}
		
		return new_logN;

	}
	
	
	//Standard VonBertalanffy
	template<class T>
	class LaA{
		
		//const _A, _Y;
		int _y;
		T _Linf, _k, _age0, _L0;
				
		// T Inv_VonB(const T Length, const T Linf, const T k, 
					// const T age0 = 0., const T L0 = 0.){
			// return (std::log(1. - (Length - L0)/Linf) + age0) / (-k);
		// }
		
		// NumericVector Inv_VonB(const NumericVector Lengths, const T Linf, const T k,
								// const T age0 = 0., const T L0 = 0.){
			// return (Rcpp::log(1. - (Lengths - L0)/Linf) + age0) / (-k);
		// }
		
		public:
		
			NumericMatrix Lengths;
			
			LaA(//const int A, const int Y,
				const T Linf, const T k,
				const T L0, const T ag0):
				//_A(A), _Y(Y), 
				_Linf(Linf), _k(k), _L0(L0), _age0(age0){
			}
			
			T VonB(T &age, int y = _y){
				T L = Linf * (1. - std::exp(-k * age - age0)) + L0;
				if( y == Lengths.size() ){
					Lengths.push_back( rep(0., _A) );
				}else if( y > Lengths.size() ){
					std::printf("Problem pushing back Length Matrix; year index incorrect.")
					return 0.
				}
				Length(y, _).push_back( L ); //Might not work
				return L;
			}
			
			NumericVector VonB(NumericVector &ages, int y = _y){
				NumericVector L = Linf * (1. - Rcpp::exp(-k * ages - age0)) + L0;
				Length(y, _) = L;
				return L;
			}
		
	};
	
}