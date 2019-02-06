
#include "density.h"

using namespace std;

Density_xy::Density_xy(double bss,double max)
{
	bs = bss;
	Nbin = (unsigned int) ceil(max/bs);

	rho = vector<vector<double> >(Nbin,vector<double>(Nbin,0.));

	bins = vector<double>(Nbin,0.);
	for(unsigned int i=0;i<Nbin;++i)
		bins[i] = (i+0.5)*bs;
	Nsample = 0;
}	



void Density_xy::sample(const System &system)
{
	++Nsample;
	double x,y;
	unsigned int jx,jy;
	for(unsigned int i=0;i<system.N;++i ) {
		x = system.r[i].x;
		x -= system.L*floor(x/system.L);
		jx = floor(x/bs);

		y = system.r[i].y;
		y -= system.L*floor(y/system.L);	
		jy = floor(y/bs);

		if((jx<Nbin) && (jy<Nbin) ) {
			rho[jx][jy] += 1.;
		}
	}
		
}




void Density_xy::normalize(const System &system) 
{
	double norm = 1./(bs*bs*system.N*Nsample);
	for(unsigned int jx = 0;jx < Nbin; ++jx ) {
		for(unsigned int jy = 0;jy < Nbin; ++jy ) {
			rho[jx][jy] *= norm;
		}
	}
}

void Density_xy::write(ostream &out)
{
	for(unsigned int jy=0;jy<Nbin;++jy) {
		for(unsigned int jx=0;jx<Nbin;++jx) {
			out << rho[jx][jy];
			if(jx<(Nbin-1)) out << ' '; 
		}
		if(jy<(Nbin-1)) out << '\n';
	}
}

void Density_xy::write(const char* outname)
{
	ofstream out;
	out.open(outname);

	for(unsigned int jy=0;jy<Nbin;++jy) {
		for(unsigned int jx=0;jx<Nbin;++jx) {
			out << rho[jx][jy];
			if(jx<(Nbin-1)) out << ' '; 
		}
		if(jy<(Nbin-1)) out << '\n';
	}

	out.close();

}

void Density_xy::write_bins(ostream &out)
{

	for(unsigned int j=0;j<Nbin;++j) {
		out << bins[j];
		if(j<(Nbin-1)) out << ' '; 
	}

}

void Density_xy::write_bins(const char* outname)
{
	ofstream out;
	out.open(outname);

	for(unsigned int j=0;j<Nbin;++j) {
		out << bins[j];
		if(j<(Nbin-1)) out << ' '; 
	}

	out.close();

}

