#include"sam_ippm.h"

sam_fits::sam_fits(){};
sam_fits::~sam_fits(){};

void sam_fits::initialise( const int &nObs, const int &nG, const int &nS, const int &nP, const int &NAnum){
	
	// array for catching Mu	
	allMus.resize(nObs*nS*nG, NAnum);
	
	// array for catching species loglike contribution to the model.
	log_like_species_contrib.resize(nS, NAnum);

	// array for catching species/group loglike contribution to the model.
	log_like_species_group_contrib.resize(nS*nG, NAnum);
	
	// array for catching dMu	
	all_deriv_mus.resize(nObs*nS*nG, NAnum);
	
	// array for catching for dlogdalpha
	dlogdalpha.resize(nS*nG, NAnum);
	
	// array for catching for dlogdbeta
	dlogdbeta.resize(nG*nP*nS, NAnum); // This is actually nG*nP*nSP
	
	// array for catching for dlogdpi
	dlogdpi.resize(nG, NAnum);
	
	// array for catching for dlogddispersion
	dlogdpi.resize(nS, NAnum);
		
};

void sam_fits::zero(const int &NAnum){
	
	allMus.assign(allMus.size(), NAnum);
	log_like_species_contrib.assign(log_like_species_contrib.size(), NAnum);
	log_like_species_group_contrib.assign(log_like_species_group_contrib.size(), NAnum);
	all_deriv_mus.assign(all_deriv_mus.size(), NAnum);
	dlogdalpha.assign(dlogdalpha.size(), NAnum);
	dlogdbeta.assign(dlogdbeta.size(), NAnum);
	dlogdpi.assign(dlogdpi.size(), NAnum);
	dlogddispersion.assign(dlogddispersion.size(), NAnum);
			
};

