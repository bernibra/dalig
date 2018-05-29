
// autotools
//#include <config.h>

// c++ header files
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

// gsl header files
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_siman.h>
#include <gsl/gsl_ieee_utils.h>

// for option parsing
#include <unistd.h>

// my header files
#include <common.hpp>
#include <alignment.hpp>
#include <network.hpp>
#include <simulated_annealing.hpp>

// namespaces
using namespace std;

// the networks are stored as global variables
Network n1;
Network n2;
double nullcost=1;


void help(){
	cerr << "Incorrect usage. Please RTFM.\n";
	exit(1);
}

int main(int argc, char *argv[])
{
    // relevant parameters for simulated annealing
    long degree = 0;
    int cost_function = 2;
    int normalization = 1;
    bool pali = false, pnei = false;


    // set the above parameters with command line options
    int opt;

	while((opt = getopt(argc, argv, "ask:l:n:")) != -1) {
		switch (opt) {
	    		case 'a':
	    			pali = true;
	    			break;
		   	case 's':
		        	pnei = true;
		       		break;
			case 'k':
				if(optarg)
					degree = strtoul(optarg, NULL, 0);
				else
					help();
				break;
			case 'l':
				if(optarg)
					cost_function = strtoul(optarg, NULL, 0);
				else
					help();
				break;
			case 'n':
				if(optarg)
					normalization = strtoul(optarg, NULL, 0);
				else
					help();
				break;
			default: // '?' //
				help();
		}
	}

	// set up the random number generator
	gsl_rng_env_setup();
	gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);

	// read in two files of networks
	Alignment * alignment = read_alignment_data(' ',n1,n2);

	alignment->degree = degree;

	if(cost_function==0){
		alignment->dfunc = &role_euclidean_distance;
	}else{
		if(cost_function==1){
			alignment->dfunc = &role_correlation;
		}else{
			alignment->dfunc = &role_chisquared;
		}
	}

	print_energy(alignment, cost_function, degree);

	print_normalized_energy(alignment, cost_function, degree);		
	print_nullmalized_energy(alignment, cost_function, degree);

	if (pali){
		overlap_pairs(alignment, true, 0);
	}

	if (pnei){
		print_neighbours(alignment);
	}

	// free allocated memory
	alignment_free(alignment);
	gsl_rng_free(r);

	return 0;
}
