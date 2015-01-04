/*
 *  effective.h
 *  
 *
 *  Created by Daniele Linaro on 5/2/10.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef EFFECTIVE
#define EFFECTIVE

#include "randlib.h"
#include <cstdlib>
#include <cstdio>
#include <ctime>

//#define EULER
#define EXACT

const int ngates = 4;	// number of gates in each channel

class effective {
	public:
		effective(long seed, double area, double gamma, double gbar) {
			/*
			 * area [um^2] - gamma [pS] - gbar [S/cm^2]
			 */
			N = (int) ceil(10000 * (area*gbar/gamma));
			printf("EFFECTIVE model>> the number of channels is %d.\n", N);
			init(seed);
		}
	
		~effective() {
			delete norm;
		}
		
		void init(long seed) {
			// random number generator
			norm = new NormalBM(0., 1., (ullong) seed);
			
			setV(-65);
			tend = 10;
			dt = 0.001;
		}
	
		void setV(double Vm) {
			v = Vm; 
			an = alphan();
			bn = betan();
		}

		void setTend(double Tend) { tend = Tend; }
		void setDt(double Dt) { dt = Dt; }
		
		double getVm() const { return v; }
		double getTend() const { return tend; }
		double getDt() const { return dt; }

		int simulate();
		
	private:
		
		void rates();
		void consistency();
		
		double vtrap(double x, double y) const;
		double alphan() const;
		double betan() const;
		
	private:
		int N;	// number of channels
	
		double n;	// gating variable
		double n_inf;	// steady state
		double tau_n;	// time constants
		double an, bn;
		
		double v;	// membrane potential
		double tend;	// total simulation time
		double dt;	// current time and integration step

		NormalBM *norm;

		double y[ngates], tau_y[ngates], noise_y[ngates], var_y[ngates];	// noise
		#ifdef EXACT
		double mu_y[ngates];
		#endif
};


#endif
