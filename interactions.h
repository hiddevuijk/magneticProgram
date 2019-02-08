#ifndef GUARD_INTERACTIONS_H
#define GUARD_INTERACTIONS_H

#include "xyz.h"

#include <math.h>

class Interactions {
public:
	Interactions() : epsilon(0.), sigma(0.), rco(0.) {}

	Interactions(double eps, double sig)
		: epsilon(eps), sigma(sig), rco(sig*pow(2.,1./6.)){}


	// calculate forces between particles in 
	// the vector r, and store in matrix F
	void get_forces(std::vector<XYZ>& F,
		const std::vector<XYZ>& r);

	double get_epsilon() const { return epsilon; };
private:

	// force between r1 and r2	
	XYZ force(const XYZ& r1,const XYZ& r2);


	double epsilon;
	double sigma;
	double rco;


	// obj. used in funcs.
	XYZ d;
	double dist, d6, f;
};


XYZ Interactions::force(const XYZ& r1,const XYZ& r2)
{
	d = r1 - r2;
	dist = sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
	if( dist < rco) {
		d6 = sigma/dist;
		d6 = d6*d6*d6*d6*d6*d6;
		f = 48*epsilon*d6*(d6-.5)/(dist*dist);
		return f*d;
	}
	return 0*d;
}	


void Interactions::get_forces(
	std::vector<XYZ>& F, const std::vector<XYZ>& r)
{
	XYZ f;
	std::fill(F.begin(),F.end(),XYZ(0.,0.,0.));

	for(unsigned int i=0; i<(r.size()-1); ++i) {
		for(unsigned int j=i+1; j<r.size(); ++j) {
			f = force(r[i],r[j]);
			F[i] += f;
			F[j] -= f;
		}
	}
}




#endif	// GUARD_INTERACTIONS_H
