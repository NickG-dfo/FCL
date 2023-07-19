
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
			  C = Parms["Catch"],
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
			Indices(y, _) = Nya(y, _) * qya(y, _) * Rcpp::exp(-t_y(y) * Zya(y, _));
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
		
		new_logN(A-2) += new_logN(A-1); //plus group
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
		NumericVector tempN(A+1);
					  
		for(int r = 0; r < R; ++r){
			
			tempN = logN(r, _) - Zra(r, _) + PE(r);
			
			tempN(A-2) += tempN(r, A-1); //plus group
			tempN.push_front(logR(r));
			tempN.erase(A);
			
			new_logN(r, _) = tempN;
			
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
		NumericVector tempN(A+1);
					  
		for(int r = 0; r < R; ++r){
			
			tempN = logN(r, _) - Zra(r, _) + PE(r, _);
			
			tempN(A-2) += tempN(r, A-1); //plus group
			tempN.push_front(logR(r));
			tempN.erase(A);
			
			new_logN(r, _) = tempN;
			
			//Add regional mixing eventually
			
		}
		
		return new_logN;

	}
	
	
	//Standard VonBertalanffy
	// template<class T>
	// class LaA{
		
		// // const _A, _Y;
		// // int y;
		// const T Linf, k, age0, L0;
		
		// public:
		
			// NumericVector Lengths;
			
			// LaA(const T _Linf, const T _k,
				// const T _L0, const T _age0):
				// Linf(_Linf), k(_k), L0(_L0), age0(_age0){
			// }
			
			// T VonB(T &age, int y){
				// T L = Linf * (1. - std::exp(-k * age - age0)) + L0;
				// Lengths.push_back(L);
				// return L;
			// }
			
			// NumericVector VonB(NumericVector &ages, int y){
				// NumericVector L = Linf * (1. - Rcpp::exp(-k * ages - age0)) + L0;
				// Lengths = cat( Lengths, L );
				// return L;
			// }
		
	// };
	
}