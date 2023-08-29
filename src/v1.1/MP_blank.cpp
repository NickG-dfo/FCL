#include <SimSource.hpp>

//Save MP_blank.cpp as a separate file, the source this file through R//
//This will create a function called 'Simulate' in R//

namespace PopSim{

	// Custom MP Function //
	
	template<class T>
	T Simulation<T>::MP_decision(Simulation<T>* psim, int n)
	{ 
		// Possible calls within MP_decision are SSB, Biomass, Biomass index, and Abundance index
		// Calls are all matrix of dimensions Years x Age, i.e. Matrix(Y, A);
		// argument 'n' allows multiple possible MPs in one model, 
		//		and can be changed via OM$MP_n in R
		
		NumericMatrix SSB = psim -> MP_ssb,
					  Biomass = psim -> MP_biomass,
					  BIndex = psim -> MP_bindex,
					  NIndex = psim -> MP_nindex;
		int A = SSB.cols(),
			Y = SSB.rows();
		
		T output = 0.;
					
		if(n == 0){
			output = sum( SSB(SSB.nrow()-1, _) ) * 0.1;
		}
		
		return output; 
	}

	// Custom Recruitment Function //

	template<class T>
	T SR<T>::CustomRecruitment(T SSB, NumericVector &parms)
	{
		T Rec = SSB * parms(1) + parms(0);		
		return Rec; 
	}

}

