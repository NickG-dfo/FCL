
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
			
			// for(int i = 0; i < years; i++){
			item++;
			c.push_back( C );
			b.push_back( r * b(item-1) * (1. - b(item-1)/K ) );
			h.push_back( C / b(item-1) );
			m.push_back( r/K * b(item-1) );
			// }
			
		}
		
		int size(void){ return item + 1; }
		
		NumericVector Biomass(void){ return b; }
		NumericVector Catch(void){ return c; }
		NumericVector HarvestRate(void){ return h; }
		NumericVector Mortality(void){ return m; }

	};
	
}


namespace PopSim{
	
	//Length Estimation//
	
	//Standard VonBertalanffy
	// template<class T>
	// class LaA{
		
		// int y;
		// const T Linf, k, age0, L0;
		
	// public:
	
		// NumericMatrix Lengths;
		
		// LaA(int Y, int A){				
			// Lengths = NumericMatrix(Y, A);
			// y = 0;
		// }
		
		// void parametrize(const T _Linf, const T _k,
						 // const T _L0, const T _age0)
		// {
			// Linf = _Linf;
			// k = _k;
			// L0 = _L0;
			// age0 = _age0;
		// }
		
		// NumericVector VonB(NumericVector ages){
			// NumericVector L = Linf * (1. - exp(-k * ages - age0)) + L0;
			// Lengths(y, _) = L;
			// y++;
			// return L;
		// }
		
		// void clear(void){
			// Lengths.fill(0.);
			// y = 0;
		// }
		
	// };
	
	//Survey Indices//
	
	NumericVector Indices_A(const NumericVector Ny, const NumericVector &qy){
		return Ny * qy;
	}
	
	NumericVector Indices_tA(const NumericVector Ny, const NumericVector &qa,
							 const NumericVector Zy, const double &tf){		
		return Ny * qa * Rcpp::exp(-tf * Zy);
	}
	
	NumericMatrix Indices_tSA(const NumericVector Ny, const NumericMatrix &qsa,
							  const NumericVector Zy, const NumericVector &tf)
	{
		int S = qsa.rows();
		NumericMatrix Indices(S, qsa.cols());
		
		for(int s = 0; s < S; ++s){
			Indices(s, _) = Ny * qsa(s, _) * exp(-tf(s) * Zy);
		}
		
		return Indices;		
	}
	
	//Abundances//
	
	template<class T>
	NumericVector logN_a(const NumericVector &logN, 
						 const NumericVector &Z,
						 const T &PE, 
						 const T logR,
						 bool plus_group = 1
						 ){

		const int A = logN.size();

		NumericVector new_logN(A);
					  
		new_logN = logN - Z + PE;		
		if(plus_group){
			new_logN(A-2) += new_logN(A-1);
		}
		new_logN.push_front(logR);
		new_logN.erase(A-1);
		
		return new_logN;

	}
	
	template<class T>
	NumericVector logN_a(const NumericVector &logN, 
						 const NumericVector &Za,
						 const NumericVector &PE, 
						 const T logR,
						 bool plus_group = 1
						 ){

		const int A = logN.size();

		NumericVector new_logN(A);
		
		new_logN = logN - Za + PE;
		if(plus_group){
			new_logN(A-2) += new_logN(A-1);
		}
		new_logN.push_front(logR);
		new_logN.erase(A);
		
		return new_logN;

	}
	
	template<class T>
	NumericMatrix logN_ra(const NumericMatrix &logN,
						  const NumericMatrix &Zra, 
						  const double &PE,
						  const NumericVector logR,
						  bool plus_group = 1
						  ){

		const int R = logN.nrow(),
				  A = logN.ncol();

		NumericMatrix new_logN(R, A);
		NumericVector* tempN = new NumericVector(new_logN);
		// *tempN = NumericVector(R, A);
					  
		for(int r = 0; r < R; ++r){
			
			*tempN = logN(r, _) - Zra(r, _) + PE;
			if(plus_group){
				(*tempN)(A-2) += (*tempN)(A-1);
			}
			tempN->push_front(logR(r));
			tempN->erase(A);
			
			new_logN(r, _) = *tempN;
			
			//Add regional mixing eventually
			
		}
		delete tempN;
		return new_logN;
	}
	
	template<class T>
	NumericMatrix logN_ra(const NumericMatrix &logN,
						  const NumericMatrix &Zra, 
						  const NumericVector &PE,
						  const NumericVector logR,
						  bool plus_group = 1
						  ){

		const int R = logN.nrow(),
				  A = logN.ncol();

		NumericMatrix new_logN(R, A);
		NumericVector* tempN = new NumericVector(new_logN);
		// *tempN = NumericVector(R, A);
					  
		for(int r = 0; r < R; ++r){
			
			*tempN = logN(r, _) - Zra(r, _) + PE(r);
			if(plus_group){
				(*tempN)(A-2) += (*tempN)(A-1);
			}
			tempN->push_front(logR(r));
			tempN->erase(A);
			
			new_logN(r, _) = *tempN;
			
			//Add regional mixing eventually
			
		}
		delete tempN;
		return new_logN;
	}
	
	template<class T>
	NumericMatrix logN_ra(const NumericMatrix &logN,
						  const NumericMatrix &Zra, 
						  const NumericMatrix &PE,
						  const NumericVector logR,
						  bool plus_group = 1
						  ){

		const int R = logN.nrow(),
				  A = logN.ncol();

		NumericMatrix new_logN(R, A);
		NumericVector* tempN = new NumericVector(new_logN);
		// *tempN = NumericVector(R, A);
					  
		for(int r = 0; r < R; ++r){
			
			*tempN = logN(r, _) - Zra(r, _) + PE(r, _);
			if(plus_group){
				(*tempN)(A-2) += (*tempN)(A-1);
			}
			tempN->push_front(logR(r));
			tempN->erase(A);
			
			new_logN(r, _) = *tempN;
			
			//Add regional mixing eventually
			
		}
		delete tempN;
		return new_logN;
	}
	
}