#ifndef	GUARD_XYZ_H
#define	GUARD_XYZ_H


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

	void normalize( double d=1.)
	{	double len = std::sqrt(x*x+y*y+z*z);
		x *= d/len;
		y *= d/len;
		z *= d/len;
	}


	// addition, subtraction, multiplication and division	

};


// implement
namespace xyz {


	double dist(const XYZ &r1, const XYZ &r2);
	double dist2(const XYZ &r1, const XYZ &r2);
	double dist_pbc(const XYZ &r1, const XYZ &r2,double L);
	double dist2_pbc(const XYZ &r1, const XYZ &r2, double L);


};




#endif
