#ifndef GUARD_WALLS_H
#define GUARD_WALLS_H


#include <math.h>
#include <iostream>

class Wall {
public:
	virtual XYZ f( const XYZ& r) = 0;


	// force from LJ potential
	double flj(double dist) {
		double dist6 = pow(sig/dist,6);
		return 24*eps*dist6*(2*dist6-1)/dist;
		}

	double get_sigma() {return sig;}
protected:
	double L;

	double sig;
	double eps;
	double rWCA;
	double L_minus_rWCA;
};

class NoWall: public Wall {
	XYZ f(const XYZ& r)
		{ XYZ Fwall(0,0,0);
			return Fwall;
		}
};

class SquareWall: public Wall {
public:
	SquareWall() {};
	SquareWall(double sigg, double epss, double l)
		 {	sig = sigg; eps = epss;
			rWCA = sigg*pow(2.,1./6.); 
			L_minus_rWCA = l-rWCA; L = l;}

	XYZ f(const XYZ& r)
	{
		XYZ Fwall(0,0,0);
		if(r.x>L_minus_rWCA) {
			Fwall.x = -flj(L-r.x);
		} else if(r.x<rWCA) {
			Fwall.x = flj(r.x);
		}


		if(r.y>L_minus_rWCA) {
			Fwall.y = -flj(L-r.y);
		} else if(r.y<rWCA) {
			Fwall.y = flj(r.y);
		} 

		return Fwall;
	}
};

#endif
