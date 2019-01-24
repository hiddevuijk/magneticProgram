#ifndef GUARD_SYSTEM_H
#define GUARD_SYSTEM_H



#include "ConfigFile.h"

#include "xyz.h"
#include "bfield.h"
#include "box_muller.h"

#include <vector>

namespace system_func {
inline void xyz_random_normal(XYZ &r, Ranq2 &ranNR) {
	r.x = ndist(ranNR);
	r.y = ndist(ranNR);
	r.z = ndist(ranNR);
}
};

struct System {
public:
	// initialize from ConfigFile object
	System(ConfigFile config);

	unsigned int N;
	double L;
	double m;
	double v0;
	double Dr;
	double dt;
	double dt2;

	double sqrt_dt;
	double sqrt_2dt;
	double sqrt_dt_Dr;
	double sqrt_2dt_Dr;

	double t;
	std::vector<XYZ> r;
	std::vector<XYZ> dr;
	std::vector<XYZ> v;
	std::vector<XYZ> p;
	
	
	BNone noField;
	BsinY fieldSineY;

	Bfield *bfield_ptr;

	// initialize with random coordinates.
	void init_random(Ranq2 &ranNR);

	// increment time
	void step(Ranq2 &ranNR );


	void write(const char* outname);
	
	// temporary containers
	double Bri;
	XYZ xi,eta,dp;	


};

void System::write(const char* outname)
{
	std::ofstream out;
	out.open(outname);
	XYZ temp;
	for(unsigned int i=0;i<N;++i) {
		temp = r[i];
		temp.pbc(L);	
		out << temp.x << '\t';
		out << temp.y << '\t';
		out << temp.z;
		if(i<(N-1)) out << '\n';
	}

	out.close();

}
void System::step(Ranq2 &ranNR )
{
	for(unsigned int i=0;i<N;++i) {

		// make half a time step with
		// deterministic forces onlyl
		Bri = bfield_ptr->f(r[i]);

		v[i].x += ( Bri*v[i].y*dt2 + v0*p[i].x*dt2)/m;
		v[i].y += (-Bri*v[i].x*dt2 + v0*p[i].x*dt2)/m;
		v[i].z += (v0*p[i].z*dt2)/m;

		r[i] += dt2*v[i]; 
			// half a time step in p space
		if( v0 > 0) {
			system_func::xyz_random_normal(eta,ranNR);
			eta *= sqrt_dt_Dr;
			dp = xyz::cross(eta,p[i]);
			p[i] += dp;
			p[i].normalize();
		}

		// make a full time step with stochastic
		// and friction forces
		system_func::xyz_random_normal(xi,ranNR);
		eta = std::exp(-dt/m)*v[i] + 
			std::sqrt( (1-std::exp(-2*dt/m))/m )*xi;	


		// make half a time step with eta 
		// in stead of v
		r[i] += dt2*eta;

			// half a time step in p space
		if( v0 > 0) {
			system_func::xyz_random_normal(eta,ranNR);
			eta *= sqrt_dt_Dr;
			dp = xyz::cross(eta,p[i]);
			p[i] += dp;
			p[i].normalize();
		}

		Bri = bfield_ptr->f(r[i]);
		v[i].x = eta.x +( Bri*eta.y*dt2 + v0*p[i].x*dt2)/m;
		v[i].y = eta.y +(-Bri*eta.x*dt2 + v0*p[i].x*dt2)/m;
		v[i].z = eta.z + (v0*p[i].z*dt2)/m;


	}

	t += dt;
}



System::System(ConfigFile config)
: noField(), fieldSineY(config.read<double>("B"), config.read<double>("w") )
{
	XYZ rr;
	N = config.read<unsigned int>("N");
	L = config.read<double>("L");
	m = config.read<double>("m");
	v0 = config.read<double>("v0");
	Dr = config.read<double>("Dr");
	dt = config.read<double>("dt");
	dt2 = 0.5*dt;
	sqrt_dt = std::sqrt(dt);
	sqrt_2dt = std::sqrt(2*dt);
	sqrt_dt_Dr = std::sqrt(dt*Dr);
	sqrt_2dt_Dr = std::sqrt(2*dt*Dr);

	t = 0.0;
	r = std::vector<XYZ>(N);
	dr = std::vector<XYZ>(N);
	v = std::vector<XYZ>(N);
	p = std::vector<XYZ>(N);
	
	std::string BType = config.read<std::string>("BType");

	if( BType == "none") {
		bfield_ptr = &noField;	
	} else if( BType == "sineY" ) {
		bfield_ptr = &fieldSineY;
	} else {
		cerr << "ERROR: " << BType << " is not a valid option.";
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
		t_unit = config.read<unsigned int>("t_unit");
		dt = config.read<double>("dt");
		bs = config.read<double>("bs");
	}

	unsigned int Nt_init;
	unsigned int Nt;
	unsigned int sample_freq;
	unsigned int print_freq;
	unsigned int t_unit;
	double dt;
	double bs;
};


#endif
