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
		fprintf(fid, "%e %e\n", t, (double) n4/N);
		update_states();
		t += dt;
	}
	
	fclose(fid);	
	return 1;
}

void microscopic::update_states() {
	double p;
	int i;
	
	setlimits();

	// n0
    p = 1-exp(-dt*(4*an));
    for(i=0; i<n0now; i++) {
		if (r->doub() <= p)	{
			n0--;
			n1++;
		}
    }
	
	// n1
    p = 1-exp(-dt*(bn+3*an));
    for(i=0; i<n1now; i++) {
		if (r->doub() <= p)	{
			n1--;
			if (r->doub() <= bn/(bn+3*an))
				n0++;
			else
				n2++;
		}
    }

	// n2
    p = 1-exp(-dt*(2*bn+2*an));
    for(i=0; i<n2now; i++) {
		if (r->doub() <= p)	{
			n2--;
			if (r->doub() <= (2*bn)/(2*bn+2*an))
				n1++;
			else
				n3++;
		}
    }

	// n3
    p = 1-exp(-dt*(3*bn+an));
    for(i=0; i<n3now; i++) {
		if (r->doub() <= p)	{
			n3--;
			if (r->doub() <= (3*bn)/(3*bn+an))
				n2++;
			else
				n4++;
		}
    }

	// n4
    p = 1-exp(-dt*(4*bn));
    for(i=0; i<n4now; i++) {
		if (r->doub() <= p)	{
			n4--;
			n3++;
		}
    }
}

double microscopic::vtrap(double x, double y) const {
	if (fabs(x/y) < 1e-6)
		return y*(1. - x/y/2.);
	return x/(exp(x/y) - 1.);
}

double microscopic::alphan() const {
	return 0.01*vtrap(-(v+55.),10.);
}

double microscopic::betan() const {
	return 0.125*exp(-(v+65.)/80.);
}

void microscopic::setlimits() {
	n0now = n0;
	n1now = n1;
	n2now = n2;
	n3now = n3;
	n4now = n4;
}

