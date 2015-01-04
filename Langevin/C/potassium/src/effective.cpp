/*
 *  effective.cpp
 *  
 *
 *  Created by Daniele Linaro on 5/2/10.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "effective.h"

int effective::simulate() {
	double t, n4_inf, one_minus_n, Y;
	int i;
	
	FILE *fid;
	char fname[30];

	sprintf(fname, "effective_V=%.1f.dat", v);
	fid = fopen(fname,"w");
	if(fid == NULL)
		return 0;
	printf("Saving to %s\n", fname);

	n = 0.0;
	n_inf = an/(an+bn);
	tau_n = 1.0/(an+bn);
	n4_inf = n_inf*n_inf*n_inf*n_inf;
	
	Y = 0.0;
	for(i=0; i<ngates; i++) {
		y[i] = 0.0;
		Y += y[i];
		tau_y[i] = tau_n/(i+1);
#ifdef EXACT
		mu_y[i] = exp(-dt/tau_y[i]);
#endif
	}
	one_minus_n = 1-n_inf;
	var_y[0] = (1.0/N) * 4*n4_inf*n_inf*n_inf*n_inf * one_minus_n;
	var_y[1] = (1.0/N) * 6*n4_inf*n_inf*n_inf * one_minus_n*one_minus_n;
	var_y[2] = (1.0/N) * 4*n4_inf*n_inf * one_minus_n*one_minus_n*one_minus_n;
	var_y[3] = (1.0/N) * n4_inf * one_minus_n*one_minus_n*one_minus_n*one_minus_n;
	
	t = 0.;
	while(t < tend) {
		consistency();
		fprintf(fid, "%e %e\n", t, n*n*n*n+Y);

		rates();

		/* forward Euler for n (it is deterministic) */
		n = n + dt * (n_inf-n)/tau_n;

		Y = 0.0;
		for(i=0; i<ngates; i++) {
		/* exact solution taken from "Gillespie, Physical Review E, 1994". */
#ifdef EXACT
			y[i] = y[i]*mu_y[i] + noise_y[i];
#endif
		/* Euler-Maruyama */
#ifdef EULER
			y[i] = y[i] - dt * y[i]/tau_y[i] + noise_y[i];
#endif
			Y += y[i];
		}
		
		t += dt;
	}
	
	fclose(fid);

	return 1;
}

void effective::consistency() {
	// n
	if(n < 0) {
		n = 0;
	}
	if(n > 1) {
		n = 1;
	}
}

void effective::rates() {
	#ifdef EULER
	for(int i=0; i<ngates; i++)
		noise_y[i] = sqrt(2*var_y[i]*dt/tau_y[i]) * norm->deviate();
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
	for(int i=0; i<ngates; i++)
		noise_y[i] = sqrt(var_y[i]*(1-mu_y[i]*mu_y[i])) * norm->deviate();
	#endif
}

double effective::vtrap(double x, double y) const {
	if (fabs(x/y) < 1e-6)
		return y*(1. - x/y/2.);
	return x/(exp(x/y) - 1.);
}

double effective::alphan() const {
	return 0.01*vtrap(-(v+55.),10.);
}
	
double effective::betan() const {
	return 0.125*exp(-(v+65.)/80.);
}
	

