#ifndef GUARD_BD_H
#define GUARD_BD_H

#include "xyz.h"
#include "system.h"
#include "box_muller.h"


#include <cmath>
#include <vector>

struct STEPPER {
public:
	STEPPER(){	sqrt2 = std::sqrt(2.); }

	// evolves the state of the system to t+dt
	void  operator() ( double dt,double &t, System &system, Ranq2 &ranNR);

private:
	double sqrt2;
};



void STEPPER::operator() ( double dt, double &t, System &system, Ranq2 &ranNR )
{
	double sqrt_dt = std::sqrt(dt);

	// random numbers for r, v increment
	XYZ xi;	
	// random numbers for p increment
	XYZ eta;
	// p increment
	XYZ dp;

	double Bri; // magnetic field at position ri

	
	
	for(unsigned int i=0;i<system.N;++i) {

		//Bri = system.bfield_ptr->f(system.r[i],t);
		Bri = 0;

		xi.x = ndist(ranNR);
		xi.y = ndist(ranNR);
		xi.z = ndist(ranNR);

		system.v[i].x += (Bri*system.v[i].y*dt - 
				system.v[i].x*dt + system.v0*system.p[i].x*dt +
				xi.x*sqrt_dt*sqrt2)/system.m;	
		system.v[i].y += (-Bri*system.v[i].x*dt - 
				system.v[i].y*dt + system.v0*system.p[i].y*dt +
				xi.y*sqrt_dt*sqrt2)/system.m;	
		system.v[i].z += (-system.v[i].z*dt + 
				system.v0*system.p[i].z + xi.z*sqrt_dt*sqrt2)/system.m;


		system.dr[i].x = system.v[i].x*dt;
		system.dr[i].y = system.v[i].y*dt;
		system.dr[i].z = system.v[i].z*dt;

		system.r[i].x += system.dr[i].x;
		system.r[i].y += system.dr[i].y;
		system.r[i].z += system.dr[i].z;

		if(system.v0 > 0) {
			eta.x = ndist(ranNR)*sqrt_dt*system.sqrt_2Dr;
			eta.y = ndist(ranNR)*sqrt_dt*system.sqrt_2Dr;
			eta.z = ndist(ranNR)*sqrt_dt*system.sqrt_2Dr;

			dp.x = eta.y*system.p[i].z - eta.z*system.p[i].y;	
			dp.y = eta.z*system.p[i].x - eta.x*system.p[i].z;	
			dp.z = eta.x*system.p[i].y - eta.y*system.p[i].x;	
			system.p[i].x += dp.x;
			system.p[i].y += dp.y;
			system.p[i].z += dp.z;
			system.p[i].normalize();
		}
	}

	t += dt;

}





#endif


