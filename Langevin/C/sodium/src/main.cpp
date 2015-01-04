/*
 *  main.cpp
 *  
 *
 *  Created by Daniele Linaro on 5/2/10.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "randlib.h"
#include "microscopic.h"
#include "effective.h"
#include "fox.h"

// SIMULATIONS PROPERTIES
const long seed = 5061983L;	// random number generator's seed
const double tend = 20000.0;// [ms]
const double dt = 0.01;	// [ms]

// NEURON PROPERTIES
const double area = 10.0;	// [um^2]
const double gammana = 10.0;	// [pS]
const double gnabar = 0.12;		// [S/cm^2]


typedef enum { 
	MODMICRO, MODEFF, MODFOX
} model;


typedef struct {
	double V;
	model type;
} options;


void printHelp(const char *progname);
int parseArgs(int argc, char *argv[], options *opt);

void printHelp(const char *progname) {
	fprintf(stderr, "\nUsage: %s [-t <type] [-V <voltage>]\n", progname);
	fprintf(stderr, "\twhere <type> is one of 'microscopic', 'effective' or 'fox'.\n");
}

int parseArgs(int argc, char *argv[], options *opt) {
	opt->V = -65;
	opt->type = MODEFF;
	
	for(int i=1; i<argc; i++) {
		if(strcmp(argv[i],"-V") == 0) {
			i++;
			if(i<argc)
				opt->V = atof(argv[i]);
			else
				return 1;
		}
		else if(strcmp(argv[i],"-t") == 0) {
			i++;
			if(i<argc) {
				if(strcmp(argv[i],"microscopic") == 0 || strcmp(argv[i],"micro") == 0)
					opt->type = MODMICRO;
				else if(strcmp(argv[i],"effective") == 0 || strcmp(argv[i],"eff") == 0)
					opt->type = MODEFF;
				else if(strcmp(argv[i],"fox") == 0)
					opt->type = MODFOX;
				else
					return 1;
			}
			else
				return 1;
		}
	}

	return 0;
}

int main(int argc, char *argv[]) {
	
	fprintf(stdout, "\nSimulator of the stochastic open-close kinetics of potassium\n");
	fprintf(stdout, "ion channels in the HH neuron model.\n");
	fprintf(stdout, "\nAuthor: Daniele Linaro - daniele.linaro@unige.it\n\n");
	
	options opt;
	microscopic micro(seed,area,gammana,gnabar);
	effective eff(seed,area,gammana,gnabar);
	fox fx(seed,area,gammana,gnabar);
	
	if(parseArgs(argc,argv,&opt) != 0) {
		printHelp(argv[0]);
		exit(1);
	}
	
	switch(opt.type) {
		case MODMICRO:
			micro.setV(opt.V);
			micro.setDt(dt);
			micro.setTend(tend);
			printf("Simulating the microscopic model...\n");
			micro.simulate();
			break;
		case MODEFF:
			eff.setV(opt.V);
			eff.setDt(dt);
			eff.setTend(tend);
			fprintf(stderr, "Simulating the effective model...\n");
			eff.simulate();
			break;
		case MODFOX:
			fx.setV(opt.V);
			fx.setDt(dt);
			fx.setTend(tend);
			fprintf(stderr, "Simulating Fox's model...\n");
			fx.simulate();
			break;
	}
	
	return 0;
}


