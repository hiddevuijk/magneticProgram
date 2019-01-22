
#include "ConfigFile.h"

#include "nr3.h"
#include "ran.h"

#include "xyz.h"
#include "bfield.h"
#include "system.h"
#include "bd.h"

#include <iostream>
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
	cout << system.r[0].x << endl;

	// integrate Nt_init time steps
	double t_init = 0;
	for(unsigned int ti = 0; ti < int_params.Nt_init; ++ ti) {
		step(int_params.dt,t_init, system,ranNR);		
	}



	double t = 0;
	for(unsigned int ti = 0; ti < int_params.Nt; ++ti) {
		step(int_params.dt,t_init, system,ranNR);		

	}
	cout << system.r[0].x << endl;


	return 0;
}








