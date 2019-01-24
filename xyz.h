#ifndef	GUARD_XYZ_H
#define	GUARD_XYZ_H

#include <cmath>
#include <ostream>

class XYZ {
public:
	XYZ(): x(0.), y(0.), z(0.) {}
	XYZ(double xx, double yy, double zz)
		: x(xx), y(yy), z(zz) {}

	double x,y,z;

	double length()
		{ return std::sqrt(x*x + y*y + z*z); }

	double length_sq()
		{ return x*x + y*y + z*z; } 

	void normalize( double d = 1.)
	{	double len = std::sqrt(x*x+y*y+z*z);
		x *= d/len;
		y *= d/len;
		z *= d/len;
	}

	void pbc(double L) {
		x -= L*std::floor(x/L);
		y -= L*std::floor(y/L);
		z -= L*std::floor(z/L);
	}

	// addition, subtraction, multiplication and division	
	XYZ operator+=(const XYZ &r) {
		x += r.x;
		y += r.y;
		z += r.z;
		return *this;
	}

	XYZ operator+=(const double & add) {
		x += add;
		y += add;
		z += add;
		return *this;
	}

	XYZ operator-=(const XYZ &r) {
		x -= r.x;
		y -= r.y;
		z -= r.z;
		return *this;
	}

	XYZ operator-=(const double & sub) {
		x -= sub;
		y -= sub;
		z -= sub;
		return *this;
	}

	XYZ operator*=(const XYZ &r) {
		x *= r.x;
		y *= r.y;
		z *= r.z;
		return *this;
	}

	XYZ operator*=(const double & mult) {
		x *= mult;
		y *= mult;
		z *= mult;
		return *this;
	}

	XYZ operator/=(const XYZ &r) {
		x /= r.x;
		y /= r.y;
		z /= r.z;
		return *this;
	}

	XYZ operator/=(const double & div) {
		x /= div;
		y /= div;
		z /= div;
		return *this;
	}
};


std::ostream& operator<< (std::ostream &out, XYZ const& r)
{
	out << r.x << '\t' << r.y << '\t' << r.z;
	return out;
}

// arithmatic operators
XYZ operator+ (const XYZ &r1,const XYZ &r2)
{ return XYZ(r1.x+r2.x,r1.y+r2.y,r1.z+r2.z); }

XYZ operator+ (const XYZ &r1,double add)
{ return XYZ(r1.x+add,r1.y+add,r1.z+add); }

XYZ operator- (const XYZ &r1,const XYZ &r2)
{ return XYZ(r1.x-r2.x,r1.y-r2.y,r1.z-r2.z); }

XYZ operator- (const XYZ &r1,double sub)
{ return XYZ(r1.x-sub,r1.y-sub,r1.z-sub); }

XYZ operator* (const XYZ &r1,const XYZ &r2)
{ return XYZ(r1.x*r2.x,r1.y*r2.y,r1.z*r2.z); }

XYZ operator* (const XYZ &r1, double mult)
{ return XYZ(r1.x*mult,r1.y*mult,r1.z*mult); }

XYZ operator/ (const XYZ &r1,const XYZ &r2)
{ return XYZ(r1.x/r2.x,r1.y/r2.y,r1.z/r2.z); }

XYZ operator/ (const XYZ &r1, double div)
{ return XYZ(r1.x/div,r1.y/div,r1.z/div); }



// implement
namespace xyz {


	double dist(const XYZ &r1, const XYZ &r2);
	double dist2(const XYZ &r1, const XYZ &r2);
	double dist_pbc(const XYZ &r1, const XYZ &r2,double L);
	double dist2_pbc(const XYZ &r1, const XYZ &r2, double L);

	// dot product cross product
	double dot(const XYZ &r1, const XYZ &r2) {
		return r1.z*r2.x + r1.y*r2.y + r1.z*r2.z;
	}

	XYZ cross(const XYZ &r1, const XYZ &r2) {
		return XYZ(r1.y*r2.z - r1.z*r2.y,
				r1.z*r2.x - r1.x*r2.z,
				r1.x*r2.y - r1.y*r2.x);
	}
};




#endif
