
#include "ConfigFile.h"

#include "nr3.h"
#include "ran.h"

#include "xyz.h"
#include "bfield.h"
#include "system.h"
#include "bd.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <vector>
#include <string>


using namespace std;



int main()
{

	// read input into config
	ConfigFile config("input.txt");

	// random number generator
	Ranq2 ranNR(config.read<unsigned int>("seed"));

	// read integration parameters
	Integration int_params(config);
	// read system parameters
	System system(config);



	// start with random uniform dist. in r
	// start with ~ Boltzmann in v
	// start with random unit sphere in p
	system.init_random(ranNR);

	STEPPER step;


	// integrate Nt_init time steps
	double t_init = 0;
	unsigned int ti;
	cout << "Starting with equilibration ...\n";
	for( ti = 0; ti < int_params.Nt_init; ++ ti) {
		if( (ti%int_params.print_freq) == 0 ) {
			cout << (int_params.Nt_init + int_params.Nt) << '\t';
			cout << ti << endl;
		}
		step(int_params.dt,t_init, system,ranNR);		
	}
	cout << "Ended equilibration. Starting sampling ... \n";


	//double t = 0;
	for(; ti < (int_params.Nt+int_params.Nt_init); ++ti) {
		if( (ti%int_params.print_freq) == 0 ) {
			cout << (int_params.Nt_init + int_params.Nt) << '\t';
			cout << ti << endl;
		}

		step(int_params.dt,t_init, system,ranNR);		

		if( (ti%int_params.sample_freq) == 0 ) {
			// sample
		}
	}



	return 0;
}








