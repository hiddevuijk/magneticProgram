
#include "ConfigFile.h"

#include "nr3.h"
#include "ran.h"

#include "xyz.h"
#include "bfield.h"
#include "system.h"
#include "bd.h"
#include "density.h"
#include "orientation.h"
#include "flux.h"

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


	// increments coordinates one time step
	//STEPPER step;

	Density_xy density(int_params.bs,system.L);
	Orientation_xy orientation(int_params.bs,system.L);
	Flux_xy flux(int_params.bs,system.L);

	// integrate Nt_init time steps
	double t_init = 0;
	unsigned int ti;
	cout << "Starting with equilibration ...\n";
	for( ti = 0; ti < int_params.Nt_init; ++ ti) {
		if( (ti%int_params.print_freq) == 0 ) {
			cout << (int_params.Nt_init + int_params.Nt) << '\t';
			cout << ti << endl;
		}
		stepper::step(system,int_params.dt,t_init,ranNR);		
	}
	cout << "Ended equilibration. Starting sampling ... \n";


	double t = t_init;
	for(; ti < (int_params.Nt+int_params.Nt_init); ++ti) {
		if( (ti%int_params.print_freq) == 0 ) {
			cout << (int_params.Nt_init + int_params.Nt) << '\t';
			cout << ti << endl;
		}

		stepper::step(system,int_params.dt,t,ranNR);		

		if( (ti%int_params.sample_freq) == 0 ) {
			density.sample(system);
			orientation.sample(system);
			flux.sample(system);
		}
	}

	cout << "Done simulation. Normalizing and writing results ..." << endl;
	// normalize and save density
	density.normalize(system);
	density.write("rho.dat");
	density.write_bins("rho_bins.dat");

	// normalize and save orientation
	orientation.normalize();
	orientation.writeX("px.dat");
	orientation.writeY("py.dat");
	orientation.write_bins("p_bins.dat");


	// normalize and save flux
	flux.normalize(system);
	flux.writeX("fx.dat");
	flux.writeY("fy.dat");
	flux.writeZ("fz.dat");
	flux.write_bins("f_bins.dat");

	return 0;
}



