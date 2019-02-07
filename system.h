#ifndef GUARD_SYSTEM_H
#define GUARD_SYSTEM_H



#include "ConfigFile.h"

#include "xyz.h"
#include "bfield.h"
#include "walls.h"

#include <iostream>
#include <vector>
#include <boost/random.hpp>

namespace system_func {
template<class RNDIST> 
inline void xyz_random_normal(XYZ &r, RNDIST &rndist) {
	r.x = rndist();
	r.y = rndist();
	r.z = rndist();
}
};

struct System {
public:
	// initialize from ConfigFile object
	System(ConfigFile config);

	const boost::normal_distribution<double> ndist;
	const boost::uniform_real<double> udist;

	int seed;
	boost::mt19937 rng;		
	boost::variate_generator<boost::mt19937&,
		boost::normal_distribution<double> > rndist;
	boost::variate_generator<boost::mt19937&,
		boost::uniform_real<double> > rudist;

	unsigned int N;
	double L;
	double m;
	double v0;
	double Dr;
	double dt;


	double sqrt_2dt;
	double sqrt_2dt_Dr;

	double t;
	std::vector<XYZ> r;
	std::vector<XYZ> dr;
	std::vector<XYZ> v;
	std::vector<XYZ> p;
	
	
	BNone noField;
	Bsin fieldSine;
	Bhill fieldHill;

	Bfield *bfield_ptr;

	std::string wallType;
	NoWall nowall;
	SquareWall sqwall;

	Wall *wall_ptr;


	// initialize with random coordinates.
	void init_random();
	void init_random(double);

	// increment time
	void step();

	void write(const char* outname);
	
	// temporary containers
	double Bri;
	XYZ xi,eta,dp,dv;	

	// force on particle i due to the walls
	std::vector<XYZ> Fwall;

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



void System::step()
{
	for(unsigned int i=0;i<N;++i) {

		Bri = bfield_ptr->f(r[i]);
		
		system_func::xyz_random_normal(xi,rndist); xi *= sqrt_2dt;

		r[i] += v[i]*dt;

		dv.x = ( Bri*v[i].y*dt - v[i].x*dt + 
				Fwall[i].x*dt + v0*p[i].x*dt + xi.x)/m;	
		dv.y = (-Bri*v[i].x*dt - v[i].y*dt +
				Fwall[i].y*dt + v0*p[i].y*dt + xi.y)/m;	
		dv.z = (-v[i].z*dt + Fwall[i].z*dt + v0*p[i].z*dt + xi.z)/m;

		v[i] += dv;
	
		Fwall[i] = wall_ptr->f(r[i]);

		if( v0 > 0) {
			
			system_func::xyz_random_normal(eta,rndist);
			eta *= sqrt_2dt_Dr;

			dp = xyz::cross(eta,p[i]);

			p[i] += dp;
			p[i].normalize();
		}
	}

	t += dt;
}



System::System(ConfigFile config)
:
	ndist(0.,1.),udist(0,1),
	seed(config.read<unsigned int>("seed")),
	rng(seed), rndist(rng,ndist), rudist(rng,udist),
	noField(), fieldSine(config.read<double>("B"),config.read<double>("w") ),
	fieldHill(config.read<double>("B"),config.read<double>("L"))
{
	XYZ rr;
	N = config.read<unsigned int>("N");
	L = config.read<double>("L");
	m = config.read<double>("m");
	v0 = config.read<double>("v0");
	Dr = config.read<double>("Dr");
	dt = config.read<double>("dt");
	sqrt_2dt = std::sqrt(2*dt);
	sqrt_2dt_Dr = std::sqrt(2*dt*Dr);

	t = 0.0;
	r = std::vector<XYZ>(N);
	dr = std::vector<XYZ>(N);
	v = std::vector<XYZ>(N);
	p = std::vector<XYZ>(N);
	
	std::string BType = config.read<std::string>("BType");

	if( BType == "none") {
		bfield_ptr = &noField;	
	} else if( BType == "sine" ) {
		bfield_ptr = &fieldSine;
	} else if( BType == "hill" ) {
		bfield_ptr = &fieldHill;
	}

	wallType = config.read<std::string>("WallType");

	nowall = NoWall();
	sqwall = SquareWall(config.read<double>("sigma"),
			config.read<double>("epsilon"),
			config.read<double>("L") );

	if( wallType == "none") {
		wall_ptr = &nowall;
	} else if( wallType == "square") {
		wall_ptr = &sqwall;
	}

	Fwall = std::vector<XYZ>(N,XYZ(0,0,0));
}

void System::init_random()
{
	double l = 0.;
	if(wallType == "square") {
		l = wall_ptr->get_sigma()*pow(2.,1./6.);
	}
	XYZ zeta;
	double d;
	for(unsigned int i=0;i<N; ++i) {
		r[i].x = l+rudist()*(L-2*l);
		r[i].y = l+rudist()*(L-2*l);
		r[i].z = l+rudist()*(L-2*l);

		// check!!!!!!	
		v[i].x = (rudist()-1.)/std::sqrt(m);
		v[i].y = (rudist()-1.)/std::sqrt(m);
		v[i].z = (rudist()-1.)/std::sqrt(m);

		do {
			zeta.x = 2*rudist() - 1.;
			zeta.y = 2*rudist() - 1.;
			zeta.z = 2*rudist() - 1.;
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
