/*
 *  fox.h
 *  
 *
 *  Created by Daniele Linaro on 9/2/10.
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
		double noise_n;	// noisy term
		double var_n;	// variance
		double an, bn;
		
		double v;	// membrane potential
		double tend;	// total simulation time
		double dt;	// current time and integration step

		NormalBM *norm;

		#ifdef EXACT
		double mu_n;
		double n_aux;
		#endif
};


#endif
