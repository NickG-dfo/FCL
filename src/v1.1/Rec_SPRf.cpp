
// namespace SR{

	// template<class T>
	// T R_obj<T>::BevertonHolth(T SSB, const NumericVector &parms){
		// T h = parms(0),
		  // R0 = parms(1);
		// return .8*R0*h*SSB / (.2* rp.SPR(0.) *R0*(1.-h) + (h-.2)*SSB);
	// }

	// template<class T>
	// T R_obj<T>::Rickerh(T SSB, const NumericVector &parms){
		// T b = std::log(5. * parms(0)) / (.8 * parms(1)),
		  // a = std::exp(b * parms(1) / rp.SPR(0.) );
		// NumericVector tempparm = {a, b};
		// return (this->Ricker)<T>(SSB, tempparm); //Might nee SR:: here
	// }

	// template<class T>
	// T R_obj<T>::Shepherdh(T SSB, const NumericVector &parms){
		// T b = parms(1) * ((.2 - parms(0)) / (parms(0) * std::pow(.2, parms(2)) - .2)) * (-1./parms(2)),
		  // a = (1. + std::pow(parms(1)/b), parms(2)) / rp.SPR(0.);
		// NumericVector tempparm = {a, b};
		// return (this->Shepherd)<T>(SSB, tempparm); //Might nee SR:: here
	// }

// }



// template<class T>
// SR::R_obj<T>::T logRec(T S, bool p_err = 1, bool r_err = 1){
	
	// T rec;
	
	// //Parameter Error		
	// NumericVector inp_parms( (this->parms).size() );
	// int k_index;
	
	// if(p_err){
		// if( this->kernel_call) { 
			// k_index = R::runif(0, (this->kernel).rows());
			// inp_parms = (this->kernel)( k_index, _ );		
		// }else{
			// inp_parms = mvrnorm(this->parms, this->cov);
		// }
	// }else{
		// if( this->kernel_call ){ 
			// k_index = this->kernel_i;
			// inp_parms = (this->kernel)( k_index, _ );	// kernel with no parm error will randomize parm
		// }else{
			// inp_parms = this->parms;
		// }
	// }
	
	// //Recruitment Estimate
	// rec = std::log( (this->recruitment)(S, inp_parms, spr0) )
	
	// //Recruitment Error
	// rec += R::rnorm(0., this->std_log_r) * this->r_err;
	
	// if(rec > this->max_rec) { rec = this->max_rec; }
	// if(bc) { rec -= std::pow(b_sigma, 2.) / 2.; } //b_sigma?
	
	// if( (this->logR).size() < this->Y ){
		// (this->logR).push_back(rec);
	// }else{ 
		// (this->logRy)(yi) = rec;
	// }
	// (this->y)++;				
	
	// return rec;
	
// }