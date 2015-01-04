/*
 *  microscopic.cpp
 *  
 *
 *  Created by Daniele Linaro on 5/2/10.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "microscopic.h"

int microscopic::simulate() {
	double t;
	
	// data logging
	FILE *fid;
	char fname[40];
	sprintf(fname, "microscopic_V=%.1f.dat", v);
	fid = fopen(fname,"w");
	if(fid == NULL)
		return 0;
	printf("Saving to %s\n", fname);

	t = 0.0;
	while(t < tend) {
		fprintf(fid, "%e %e\n", t, (double) m3h1/N);
		update_states();
		t += dt;
	}
	
	fclose(fid);	
	return 1;
}

void microscopic::update_states() {
	double p, gamma, rnd;
	int i;
	
	setlimits();
	
	// m0h0
	gamma = 3*am+ah;
    p = 1-exp(-dt*gamma);
    for(i=0; i<m0h0now; i++) {
		if (r->doub() <= p)	{
			m0h0--;
			rnd = r->doub();
			if (rnd <= (3*am)/gamma)
				m1h0++;
			else
				m0h1++;
		}
    }

	// m1h0
	gamma = 2*am+ah+bm;
    p = 1-exp(-dt*gamma);
    for(i=0; i<m1h0now; i++) {
		if (r->doub() <= p)	{
			m1h0--;
			rnd = r->doub();
			if (rnd <= (2*am)/gamma)
				m2h0++;
			else if (rnd <= (2*am+ah)/gamma)
				m1h1++;
			else
				m0h0++;
		}
    }

	// m2h0
	gamma = am+ah+2*bm;
    p = 1-exp(-dt*gamma);
    for(i=0; i<m2h0now; i++) {
		if (r->doub() <= p)	{
			m2h0--;
			rnd = r->doub();
			if (rnd <= (am)/gamma)
				m3h0++;
			else if (rnd <= (am+ah)/gamma)
				m2h1++;
			else
				m1h0++;
		}
    }

	// m3h0
	gamma = ah+3*bm;
    p = 1-exp(-dt*gamma);
    for(i=0; i<m3h0now; i++) {
		if (r->doub() <= p)	{
			m3h0--;
			rnd = r->doub();
			if (rnd <= (ah)/gamma)
				m3h1++;
			else
				m2h0++;
		}
    }

	// m0h1
	gamma = 3*am+bh;
    p = 1-exp(-dt*gamma);
    for(i=0; i<m0h1now; i++) {
		if (r->doub() <= p)	{
			m0h1--;
			rnd = r->doub();
			if (rnd <= (3*am)/gamma)
				m1h1++;
			else
				m0h0++;
		}
    }

	// m1h1
	gamma = 2*am+bh+bm;
    p = 1-exp(-dt*gamma);
    for(i=0; i<m1h1now; i++) {
		if (r->doub() <= p)	{
			m1h1--;
			rnd = r->doub();
			if (rnd <= (2*am)/gamma)
				m2h1++;
			else if (rnd <= (2*am+bh)/gamma)
				m1h0++;
			else
				m0h1++;
		}
    }

	// m2h1
	gamma = am+bh+2*bm;
    p = 1-exp(-dt*gamma);
    for(i=0; i<m2h1now; i++) {
		if (r->doub() <= p)	{
			m2h1--;
			rnd = r->doub();
			if (rnd <= (am)/gamma)
				m3h1++;
			else if (rnd <= (am+bh)/gamma)
				m2h0++;
			else
				m1h1++;
		}
    }

	// m3h1
	gamma = bh+3*bm;
    p = 1-exp(-dt*gamma);
    for(i=0; i<m3h1now; i++) {
		if (r->doub() <= p)	{
			m3h1--;
			rnd = r->doub();
			if (rnd <= (bh)/gamma)
				m3h0++;
			else
				m2h1++;
		}
    }
}

double microscopic::vtrap(double x, double y) const {
	if (fabs(x/y) < 1e-6)
		return y*(1. - x/y/2.);
	return x/(exp(x/y) - 1.);
}

double microscopic::alpham() const {
	return 0.1 * vtrap(-(v+40.),10.);
}

double microscopic::betam() const {
	return 4. * exp(-(v+65.)/18.);
}

double microscopic::alphah() const {
	return 0.07 * exp(-(v+65.)/20.);
}

double microscopic::betah() const {
	return 1.0 / (exp(-(v+35.)/10.) + 1.);
}

void microscopic::setlimits() {
	m0h0now = m0h0;
	m1h0now = m1h0;
	m2h0now = m2h0;
	m3h0now = m3h0;
	m0h1now = m0h1;
	m1h1now = m1h1;
	m2h1now = m2h1;
	m3h1now = m3h1;
}

