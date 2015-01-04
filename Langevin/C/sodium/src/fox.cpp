/*
 *  fox.cpp
 *  
 *
 *  Created by Daniele Linaro on 9/2/10.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
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

	// m
	tau_m = 1.0/(am+bm);
	m_inf = am*tau_m;
	m = m_inf;
	var_m = (1.0/N) * am*bm / ((am+bm)*(am+bm));
	// h
	tau_h = 1.0/(ah+bh);
	h_inf = ah*tau_h;
	h = h_inf;
	var_h = (1.0/N) * ah*bh / ((ah+bh)*(ah+bh));
	
	#ifdef EXACT
	mu_m = exp(-dt/tau_m);
	mu_h = exp(-dt/tau_h);
	#endif
	
	t = 0.;
	while(t < tend) {
		consistency();
		fprintf(fid, "%e %e\n", t, m*m*m*h);

		rates();
				
		/* exact solution taken from "Gillespie, Physical Review E, 1994". */
		#ifdef EXACT
		m_aux = m_aux*mu_m + noise_m;
		h_aux = h_aux*mu_h + noise_h;
		m = m_inf + m_aux;
		h = h_inf + h_aux;
		#endif
		
		/* Euler-Maruyama */
		#ifdef EULER
		m = m - dt * m/tau_m + noise_m;
		h = h - dt * h/tau_h + noise_h;
		#endif
		
		t += dt;
	}
	
	fclose(fid);

	return 1;
}

void fox::consistency() {
	// m
	if(m < 0) {
		m = 0;
	}
	if(m > 1) {
		m = 1;
	}
	// h
	if(h < 0) {
		h = 0;
	}
	if(h > 1) {
		h = 1;
	}
}

void fox::rates() {
	#ifdef EULER
	noise_m = sqrt(2*var_m*dt/tau_m) * norm->deviate();
	noise_h = sqrt(2*var_h*dt/tau_h) * norm->deviate();
	#endif
	#ifdef EXACT
	noise_m = sqrt(var_m*(1-mu_m*mu_m)) * norm->deviate();
	noise_h = sqrt(var_h*(1-mu_h*mu_h)) * norm->deviate();
	#endif
}

double fox::vtrap(double x, double y) const {
	if (fabs(x/y) < 1e-6)
		return y*(1. - x/y/2.);
	return x/(exp(x/y) - 1.);
}


double fox::alpham() const {
	return 0.1 * vtrap(-(v+40.),10.);
}

double fox::betam() const {
	return 4. * exp(-(v+65.)/18.);
}

double fox::alphah() const {
	return 0.07 * exp(-(v+65.)/20.);
}

double fox::betah() const {
	return 1.0 / (exp(-(v+35.)/10.) + 1.);
}

