/*
 *  fox.cpp
 *  
 *  Created by Daniele Linaro on 9/2/10.
 *  daniele.linaro@unige.it
 *
 */

#include "fox.h"

int fox::simulate() {
	double t;
	
	FILE *fid;
	char fname[30];

	sprintf(fname, "fox_V=%.1f.dat", v);
	fid = fopen(fname,"w");
	if(fid == NULL)
		return 0;
	fprintf(stdout, "Saving to %s\n", fname);

	tau_n = 1.0/(an+bn);
	n_inf = an*tau_n;
	n = n_inf;
	var_n = (1.0 / N) * an*bn / ((an+bn)*(an+bn));
	
	#ifdef EXACT
	mu_n = exp(-dt/tau_n);
	#endif
	
	t = 0.;
	while(t < tend) {
		consistency();
		fprintf(fid, "%e %e\n", t, n*n*n*n);

		rates();

		/* exact solution taken from "Gillespie, Physical Review E, 1994". */
		#ifdef EXACT
		n_aux = n_aux*mu_n + noise_n;
		n = n_inf + n_aux;
		#endif
		
		/* Euler-Maruyama */
		#ifdef EULER
		n = n - dt * n/tau_n + noise_n;
		#endif
		
		t += dt;
	}
	
	fclose(fid);

	return 1;
}

void fox::consistency() {
	// n
	if(n < 0) {
		n = 0;
	}
	if(n > 1) {
		n = 1;
	}
}

void fox::rates() {
	#ifdef EULER
	noise_n = sqrt(2*var_n*dt/tau_n) * norm->deviate();
	#endif
	#ifdef EXACT
	/***
	 * The solution, as it appears in the paper by Gillespie,
	 * should be as follows:
	 *
	 *	noise = sqrt((var*tau_n/2)*(1-mu*mu)) * norm->deviate();
	 *	noise4 = sqrt((var4*tau_n4/2)*(1-mu4*mu4)) * norm->deviate();
	 *
	 * I have taken away the coefficient (tau_n/2) in order to have a
	 * resulting process whose variance is equal to var, which is what we
	 * want in order to have results which are coherent with the effective
	 * model.
	 ***/
	noise_n = sqrt(var_n*(1-mu_n*mu_n)) * norm->deviate();
	#endif
}

double fox::vtrap(double x, double y) const {
	if (fabs(x/y) < 1e-6)
		return y*(1. - x/y/2.);
	return x/(exp(x/y) - 1.);
}

double fox::alphan() const {
	return 0.01*vtrap(-(v+55.),10.);
}
	
double fox::betan() const {
	return 0.125*exp(-(v+65.)/80.);
}
	

