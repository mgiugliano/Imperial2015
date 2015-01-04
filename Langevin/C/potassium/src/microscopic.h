/*
 *  microscopic.h
 *  
 *
 *  Created by Daniele Linaro on 5/2/10.
 *
 */

#ifndef MICROSCOPIC
#define MICROSCOPIC

#include "randlib.h"
#include <ctime>
#include <cstdio>

class microscopic {
	public:
		microscopic(long seed, double area, double gamma, double gbar) {
			/*
			 * area [um^2] - gamma [pS] - gbar [S/cm^2]
			 */
			N = (int) ceil(10000 * (area*gbar/gamma));
			printf("MICROSCOPIC model>> the number of channels is %d.\n", N);
			init(seed);
		}
		
		~microscopic() {
			delete r;
		}
		
		void init(long seed) {
			// random number generator
			r = new UniformRandom(seed);
			
			setV(-65);
			tend = 10;
			dt = 0.001;

			n0 = N;
			n1 = 0;
			n2 = 0;
			n3 = 0;
			n4 = 0;
		}
	
		void setV(double Vm) {
			v = Vm;
			an = alphan();
			bn = betan();
		}
		
		void setTend(double Tend) { tend = Tend; }
		void setDt(double Dt) { dt = Dt; }
		
		double getV() const { return v; }
		double getTend() const { return tend; }
		double getDt() const { return dt; }
		
		int simulate();
		
	private:
		void update_states();
		
		double vtrap(double x, double y) const;
		double alphan() const;
		double betan() const;
		
		void setlimits();
		
	private:
		int N;	// number of channels
		int n0, n1, n2, n3, n4;
		int n0now, n1now, n2now, n3now, n4now;
		
		double an, bn;

		double v;	// membrane potential
		double tend;	// total simulation time
		double dt;	// integration step
		
		UniformRandom * r;
};


#endif

