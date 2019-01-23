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

class Bsigmoid: public Bfield {
public:
	Bsigmoid(double b0, double b, double l) {
			B = b; B0 = b0; L=l;}

	double f(const XYZ& r) {
			return B0*tanh(B*(r.y-L/2.));}

};


class Bsaw: public Bfield {
public:
	Bsaw(double b0, double b, double l) {
			B = b; B0 = b0; L=l;}

	double f(const XYZ& r) {
			if(r.y<L/4.) {
				return B0+B*0.25*L-B*(r.y-0.25*L);	
			} else if(r.y>3*L/4.) {
				return B0+B*0.75*L-B*(r.y-0.75*L);
			} else {
				return B0+B*r.y;
			}
	}

};




class BsinY: public Bfield {
public:
	BsinY(double BB, double ww) {
		B = BB; w=ww; }

	double f(const XYZ& r) {
			return B*std::sin(w*r.y); }

};

class BlinearY: public Bfield {
public:
	BlinearY(double b0, double b) {
			B = b; B0 = b0;}

	double f(const XYZ& r) {
			return B0+B*r.y;}

};

class BlinearR: public Bfield {
public:
	BlinearR(double b0,double b,double l) {
		B = b; B0 = b0; L = l;}

	double f(const XYZ& r) {
		double d = sqrt( (r.x-L/2)*(r.x-L/2) + 
						 (r.y-L/2)*(r.y-L/2));
		return B0+B*d;
	}

};





#endif

