/*
 *  fox.h
 *  
 *
 *  Created by Daniele Linaro on 5/2/10.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef FOX
#define FOX

#include "randlib.h"
#include <cstdlib>
#include <cstdio>
#include <ctime>

//#define EULER
#define EXACT

class fox {
	public:
		fox(long seed, int NN = 10000) : N(NN) {
			init(seed);
		}
			
		fox(long seed, double area, double gamma, double gbar) {
			/*
			 * area [um^2] - gamma [pS] - gbar [S/cm^2]
			 */
			N = ceil(10000 * (area*gbar/gamma));
			printf("FOX EFFECTIVE model>> the number of channels is %d.\n", N);
			init(seed);
		}
	
		~fox() {
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
			am = alpham();
			bm = betam();
			ah = alphah();
			bh = betah();		
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
		double alpham() const;
		double betam() const;
		double alphah() const;
		double betah() const;
		
	private:
		int N;	// number of channels

		double m, h;	// gating variables
		double m_inf, h_inf;	// steady states
		double tau_m, tau_h;	// time constants
		double noise_m, noise_h;// noisy terms
		double var_m, var_h;	// variance
		double am, bm, ah, bh;
		
		double v;	// membrane potential
		double tend;	// total simulation time
		double dt;	// current time and integration step

		NormalBM *norm;

		#ifdef EXACT
		double mu_m, mu_h;
		double m_aux, h_aux;
		#endif
};


#endif
