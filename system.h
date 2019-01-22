#ifndef GUARD_SYSTEM_H
#define GUARD_SYSTEM_H



#include "ConfigFile.h"

#include "xyz.h"
#include "bfield.h"


#include <vector>

struct System {
public:
	// initialize from ConfigFile object
	System(ConfigFile config);

	unsigned int N;
	double L;
	double m;
	double v0;
	double Dr;
	double sqrt_2Dr;


	std::vector<XYZ> r;
	std::vector<XYZ> dr;
	std::vector<XYZ> v;
	std::vector<XYZ> p;
	
	
	BsinY fieldSineY;

	BNone noField;


	Bfield *bfield_ptr;

	// initialize with random coordinates.
	void init_random(Ranq2 &ranNR);

};


System::System(ConfigFile config)
: fieldSineY(config.read<double>("B"), config.read<double>("w") )
{
	N = config.read<unsigned int>("N");
	L = config.read<double>("L");
	m = config.read<double>("m");
	v0 = config.read<double>("v0");
	Dr = config.read<double>("Dr");
	sqrt_2Dr = std::sqrt(2*Dr);

	r = std::vector<XYZ>(N);
	dr = std::vector<XYZ>(N);
	v = std::vector<XYZ>(N);
	p = std::vector<XYZ>(N);
	
	std::string BType = config.read<std::string>("BType");

		
	if( BType == "none") {
		bfield_ptr = &noField;	
	} else if( BType == "sineY" ) {
		bfield_ptr = &fieldSineY;
	}
	

}


void System::init_random(Ranq2 &ranNR)
{
	XYZ zeta;
	double d;
	for(unsigned int i=0;i<N; ++i) {
		r[i].x = ranNR.doub()*L;
		r[i].y = ranNR.doub()*L;
		r[i].z = ranNR.doub()*L;

		// check!!!!!!	
		v[i].x = (ranNR.doub()-1.)/std::sqrt(m);
		v[i].y = (ranNR.doub()-1.)/std::sqrt(m);
		v[i].z = (ranNR.doub()-1.)/std::sqrt(m);

		do {
			zeta.x = 2*ranNR.doub() - 1.;
			zeta.y = 2*ranNR.doub() - 1.;
			zeta.z = 2*ranNR.doub() - 1.;
			d = zeta.length_sq();	
		} while (d > 1.);
		zeta.normalize();
		p[i] = zeta;
	}


}



class Integration {
public:
	Integration( ConfigFile config)
	{
		Nt_init = config.read<unsigned int>("Nt_init");
		Nt = config.read<unsigned int>("Nt");
		sample_freq = config.read<unsigned int>("sample_freq");
		print_freq = config.read<unsigned int>("print_freq");
		dt = config.read<double>("dt");
	}

	unsigned int Nt_init;
	unsigned int Nt;
	unsigned int sample_freq;
	unsigned int print_freq;
	double dt;
};


#endif
