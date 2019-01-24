
#include "ConfigFile.h"

#include "nr3.h"
#include "ran.h"

#include "xyz.h"
#include "bfield.h"
#include "system.h"
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

	// read integration parameters
	Integration int_params(config);

	// read system parameters
	System system(config);

	// random number generator
	Ranq2 ranNR(config.read<unsigned int>("seed"));

	// start with random config. 
	system.init_random(ranNR);

	// objects to sample density, orientation and flux
	Density_xy density(int_params.bs,system.L);
	Orientation_xy orientation(int_params.bs,system.L);
	Flux_xy flux(int_params.bs,system.L);

	// integrate Nt_init time steps
	unsigned int ti;
	cout << "Starting with equilibration ...\n";
	for( ti = 0; ti < int_params.Nt_init; ++ ti) {
		if( (ti%int_params.print_freq) == 0 ) {
			cout << (int_params.Nt_init + int_params.Nt) << '\t';
			cout << ti << endl;
		}
		for(int tti=0; tti < int_params.t_unit; ++tti)		
			system.step(ranNR);		
	}
	cout << "Ended equilibration. Starting sampling ... \n";


	for(; ti < (int_params.Nt+int_params.Nt_init); ++ti) {
		if( (ti%int_params.print_freq) == 0 ) {
			cout << (int_params.Nt_init + int_params.Nt) << '\t';
			cout << ti << endl;
		}
		for(int tti = 0;tti<int_params.t_unit;++tti)
			system.step(ranNR);		

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
	orientation.writeZ("pz.dat");
	orientation.write_bins("p_bins.dat");

	// normalize and save flux
	flux.normalize(system);
	flux.writeX("fx.dat");
	flux.writeY("fy.dat");
	flux.writeZ("fz.dat");
	flux.write_bins("f_bins.dat");

	// write final configuration
	system.write("final_config.dat");
	return 0;
}



