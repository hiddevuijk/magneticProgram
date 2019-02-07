#ifndef GUARD_BFIELD_H
#define GUARD_BFIELD_H

#include <vector>
#include "xyz.h"

class Bfield {
public:
	virtual double f(const XYZ& r) = 0; 

protected:
	double B;
	double B0;
	double w;
	double wt;
	double L;
};

class BNone: public Bfield {
	double f(const XYZ& r)
			{return 0;}
};


class Bsin: public Bfield {
public:
	Bsin(double BB, double ww) {
		B = BB; w=ww; }

	double f(const XYZ& r) {
			return B*std::sin(w*r.x); }

};

class Bhill: public Bfield {
public:
	Bhill(double BB, double LL) {
		B = BB; L = LL; }

	double f(const XYZ& r) {
		if( (r.x < 0.2*L) or (r.x>0.8*L) ){
			return 0.;
		} else if( r.x<0.5*L) {
			return B*(r.x-0.2*L);
		} else {
			return B*(0.8*L-r.x);
		} 
	}
};

#endif

