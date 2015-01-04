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
	double t, m3_inf, one_minus_m, one_minus_h, Z;
	int i;
	
	FILE *fid;
	char fname[30];

	sprintf(fname, "effective_V=%.1f.dat", v);
	fid = fopen(fname,"w");
	if(fid == NULL)
		return 0;
	printf("Saving to %s\n", fname);

	// m
	m = 0.0;
	m_inf = am/(am+bm);
	tau_m = 1.0/(am+bm);
	m3_inf = m_inf*m_inf*m_inf;

	// h
	h = 0.0;
	h_inf = ah/(ah+bh);
	tau_h = 1.0/(ah+bh);

	// auxiliary variables
	one_minus_m = 1.-m_inf;
	one_minus_h = 1.-h_inf;
	// tau's
	tau_z[0] = tau_h;
	tau_z[1] = tau_m;
	tau_z[2] = tau_m/2.;
	tau_z[3] = tau_m/3.;
	tau_z[4] = tau_m*tau_h/(tau_m+tau_h);
	tau_z[5] = tau_m*tau_h/(tau_m+2*tau_h);
	tau_z[6] = tau_m*tau_h/(tau_m+3*tau_h);
	// variances
	var_z[0] = (1.0/N) * m3_inf*m3_inf*h_inf * one_minus_h;
	var_z[1] = (1.0/N) * 3*m3_inf*m_inf*m_inf*h_inf*h_inf * one_minus_m;
	var_z[2] = (1.0/N) * 3*m3_inf*m_inf*h_inf*h_inf * one_minus_m*one_minus_m;
	var_z[3] = (1.0/N) * m3_inf*h_inf*h_inf * one_minus_m*one_minus_m*one_minus_m;
	var_z[4] = (1.0/N) * 3*m3_inf*m_inf*m_inf*h_inf * one_minus_m*one_minus_h;
	var_z[5] = (1.0/N) * 3*m3_inf*m_inf*h_inf * one_minus_m*one_minus_m*one_minus_h;
	var_z[6] = (1.0/N) * m3_inf*h_inf * one_minus_m*one_minus_m*one_minus_m*one_minus_h;	
	// z
	Z = 0.0;
	for(i=0; i<nstates-1; i++) {
		z[i] = 0.0;
		Z += z[i];
#ifdef EXACT
		mu_z[i] = exp(-dt/tau_z[i]);
#endif
	}
	
	t = 0.;
	while(t < tend) {
		consistency();
		fprintf(fid, "%e %e\n", t, m*m*m*h+Z);

		rates();
		
		/* forward Euler for m and h (they are deterministic) */
		m = m + dt * (m_inf-m)/tau_m;
		h = h + dt * (h_inf-h)/tau_h;
		
		Z = 0.0;
		/* exact solution taken from "Gillespie, Physical Review E, 1994". */
		#ifdef EXACT
		for(i=0; i<nstates-1; i++) {
			z[i] = z[i]*mu_z[i] + noise_z[i];
			Z += z[i];
		}
		#endif
		
		/* Euler-Maruyama */
		#ifdef EULER
		for(i=0; i<nstates-1; i++) {
			z[i] = z[i] - dt * z[i]/tau_z[i] + noise_z[i];
			Z += z[i];
		}
		#endif
		
		t += dt;
	}
	
	fclose(fid);

	return 1;
}

void effective::consistency() {
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

void effective::rates() {
	for(int i=0; i<nstates-1; i++) {
	#ifdef EULER
		noise_z[i] = sqrt(2*var_z[i]*dt/tau_z[i]) * norm->deviate();
	#endif
	#ifdef EXACT
		noise_z[i] = sqrt(var_z[i]*(1-mu_z[i]*mu_z[i])) * norm->deviate();
	#endif
	}
}

double effective::vtrap(double x, double y) const {
	if (fabs(x/y) < 1e-6)
		return y*(1. - x/y/2.);
	return x/(exp(x/y) - 1.);
}


double effective::alpham() const {
	return 0.1 * vtrap(-(v+40.),10.);
}

double effective::betam() const {
	return 4. * exp(-(v+65.)/18.);
}

double effective::alphah() const {
	return 0.07 * exp(-(v+65.)/20.);
}

double effective::betah() const {
	return 1.0 / (exp(-(v+35.)/10.) + 1.);
}

