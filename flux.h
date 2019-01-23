#ifndef GUARD_FLUX_H
#define GUARD_FLUX_H


class Flux_xy {

public:
	Flux_xy(double bs, double max);

	void sample(const System &system);

	void normalize(const System &system);
 
	// write to out stream
	void writeX(std::ostream &out);
	void writeY(std::ostream &out);
	void writeZ(std::ostream &out);
	void write_bins(std::ostream &out);

	// write to file named outname
	void writeX(const char* outname);
	void writeY(const char* outname);
	void writeZ(const char* outname);
	void write_bins(const char* outname);

	unsigned int get_Nsample() {return Nsample;}

private:
	double bs;
	unsigned int Nbin;

	std::vector<std::vector<XYZ> > f;
	std::vector<double> bins;

	unsigned int Nsample;
};


Flux_xy::Flux_xy(double bss,double max)
{
	bs = bss;
	Nbin = (unsigned int) std::ceil(max/bs);

	f = std::vector<std::vector<XYZ> >(Nbin,
			std::vector<XYZ>(Nbin));

	bins = std::vector<double>(Nbin,0.);
	for(unsigned int i=0;i<Nbin;++i)
		bins[i] = (i+0.5)*bs;
	Nsample = 0;
}	



void Flux_xy::sample(const System &system)
{
	++Nsample;
	double x,y;
	unsigned int jx,jy;
	for(unsigned int i=0;i<system.N;++i ) {
		x = system.r[i].x;
		x -= system.L*std::floor(x/system.L);
		jx = std::floor(x/bs);

		y = system.r[i].y;
		y -= system.L*std::floor(y/system.L);	
		jy = std::floor(y/bs);

		if((jx<Nbin) && (jy<Nbin) ) 
			f[jx][jy] += system.v[i];
	}
		
}

void Flux_xy::normalize(const System &system) 
{
	double norm = 1./(Nsample*system.N*bs*bs);
	for(unsigned int jx = 0;jx < Nbin; ++jx ) {
		for(unsigned int jy = 0;jy < Nbin; ++jy ) {
			f[jx][jy] *= norm;
		}
	}
}

void Flux_xy::writeX(std::ostream &out)
{
	for(unsigned int jy=0;jy<Nbin;++jy) {
		for(unsigned int jx=0;jx<Nbin;++jx) {
			out << f[jx][jy].x;
			if(jx<(Nbin-1)) out << ' '; 
		}
		if(jy<(Nbin-1)) out << '\n';
	}
}

void Flux_xy::writeY(std::ostream &out)
{
	for(unsigned int jy=0;jy<Nbin;++jy) {
		for(unsigned int jx=0;jx<Nbin;++jx) {
			out << f[jx][jy].y;
			if(jx<(Nbin-1)) out << ' '; 
		}
		if(jy<(Nbin-1)) out << '\n';
	}
}


void Flux_xy::writeZ(std::ostream &out)
{
	for(unsigned int jy=0;jy<Nbin;++jy) {
		for(unsigned int jx=0;jx<Nbin;++jx) {
			out << f[jx][jy].z;
			if(jx<(Nbin-1)) out << ' '; 
		}
		if(jy<(Nbin-1)) out << '\n';
	}
}


void Flux_xy::writeX(const char* outname)
{
	std::ofstream out;
	out.open(outname);

	for(unsigned int jy=0;jy<Nbin;++jy) {
		for(unsigned int jx=0;jx<Nbin;++jx) {
			out << f[jx][jy].x;
			if(jx<(Nbin-1)) out << ' '; 
		}
		if(jy<(Nbin-1)) out << '\n';
	}

	out.close();

}

void Flux_xy::writeY(const char* outname)
{
	std::ofstream out;
	out.open(outname);

	for(unsigned int jy=0;jy<Nbin;++jy) {
		for(unsigned int jx=0;jx<Nbin;++jx) {
			out << f[jx][jy].y;
			if(jx<(Nbin-1)) out << ' '; 
		}
		if(jy<(Nbin-1)) out << '\n';
	}

	out.close();

}
void Flux_xy::writeZ(const char* outname)
{
	std::ofstream out;
	out.open(outname);

	for(unsigned int jy=0;jy<Nbin;++jy) {
		for(unsigned int jx=0;jx<Nbin;++jx) {
			out << f[jx][jy].z;
			if(jx<(Nbin-1)) out << ' '; 
		}
		if(jy<(Nbin-1)) out << '\n';
	}

	out.close();

}
void Flux_xy::write_bins(std::ostream &out)
{

	for(unsigned int j=0;j<Nbin;++j) {
		out << bins[j];
		if(j<(Nbin-1)) out << ' '; 
	}

}

void Flux_xy::write_bins(const char* outname)
{
	std::ofstream out;
	out.open(outname);

	for(unsigned int j=0;j<Nbin;++j) {
		out << bins[j];
		if(j<(Nbin-1)) out << ' '; 
	}

	out.close();

}


#endif