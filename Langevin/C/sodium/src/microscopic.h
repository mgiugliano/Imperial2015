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
#include <cmath>

class microscopic {
	public:
		microscopic(long seed, double area, double gamma, double gbar) {
			/*
			 * area [um^2] - gamma [pS] - gbar [S/cm^2]
			 */
			N = (int) ceil(10000.0 * (area*gbar/gamma));
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
			
			m0h0 = 0;
			m1h0 = 0;
			m2h0 = 0;
			m3h0 = 0;
			m0h1 = N;
			m1h1 = 0;
			m2h1 = 0;
			m3h1 = 0;
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
		
		double getV() const { return v; }
		double getTend() const { return tend; }
		double getDt() const { return dt; }
		
		int simulate();
		
	private:
		void update_states();
		
		double vtrap(double x, double y) const;
		double alpham() const;
		double betam() const;
		double alphah() const;
		double betah() const;
		
		void setlimits();
		
	private:
		int N;	// number of channels
		int m0h0, m1h0, m2h0, m3h0, m0h1, m1h1, m2h1, m3h1;
		int m0h0now, m1h0now, m2h0now, m3h0now, m0h1now, m1h1now, m2h1now, m3h1now;
		
		double am, bm, ah, bh;

		double v;	// membrane potential
		double tend;	// total simulation time
		double dt;	// integration step
				
		UniformRandom * r;
};


#endif

