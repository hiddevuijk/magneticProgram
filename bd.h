#ifndef GUARD_BD_H
#define GUARD_BD_H

#include "xyz.h"
#include "system.h"
#include "box_muller.h"


#include <cmath>
#include <vector>

namespace stepper {
// make r to a random normal distributed vector
inline void xyz_random_normal(XYZ &r, Ranq2 &ranNR) {
	r.x = ndist(ranNR);
	r.y = ndist(ranNR);
	r.z = ndist(ranNR);
}
const long double sqrt2 = std::sqrt(2.);


void step(System &system, double dt, double &t, Ranq2 &ranNR )
{
	double sqrt_dt = std::sqrt(dt);
	double Bri;
	XYZ xi,eta,dp;	



	for(unsigned int i=0;i<system.N;++i) {

		Bri = system.bfield_ptr->f(system.r[i]);

		xyz_random_normal(xi,ranNR);
		xi *= sqrt_dt*stepper::sqrt2;

		system.v[i].x += (Bri*system.v[i].y*dt - 
				system.v[i].x*dt + system.v0*system.p[i].x*dt +
				xi.x)/system.m;	
		system.v[i].y += (-Bri*system.v[i].x*dt - 
				system.v[i].y*dt + system.v0*system.p[i].y*dt +
				xi.y)/system.m;	
		system.v[i].z += (-system.v[i].z*dt + 
				system.v0*system.p[i].z*dt + xi.z)/system.m;


		system.dr[i] = system.v[i]*dt;
		system.r[i] += system.dr[i];

		if(system.v0 > 0) {
			
			xyz_random_normal(eta,ranNR);
			eta *= sqrt_dt*system.sqrt_2Dr;

			dp = xyz::cross(eta,system.p[i]);

			system.p[i] += dp;
			system.p[i].normalize();
		}
	}

	t += dt;

}



};

#endif


