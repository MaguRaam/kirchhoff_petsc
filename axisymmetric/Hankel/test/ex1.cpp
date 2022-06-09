// boost headers:
#include <boost/math/special_functions.hpp>

// c++ headers:
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>
#include <cstdlib>


double cylindrical_wave(double r, double t){

    double omega = 1.0;                     //frequency
    double k     = 1.0;                     //wave number
    std::complex<double> A(1.0, 0.0);       //amplitude
    std::complex<double> i(0.0,1.0);        //iota
    
    return (A*boost::math::cyl_hankel_2(0, k*r)*exp(i*omega*t)).real();
}


void write(const std::vector<double>& x, const std::vector<double>& p){

    std::ofstream File("plot.dat", std::ios::out) ;
    File.flags( std::ios::dec | std::ios::scientific );
    File.precision(16) ;
    if(!File) {std::cerr<< "Error: Output file couldnot be opened.\n";}

    for (int i = 0; i < x.size(); ++i) File << x[i] << "\t" << p[i] << "\t" << std::endl; 

    File.close();
}


int main(){

	//parameters:
	double rmin = 0.001;
	double rmax = 20.0;
	int nr = 10000;
	double dr = (rmax - rmin)/static_cast<double>(nr - 1);
	double t = 1.0;


	//create grid:
	std::vector<double> r(nr);
	for (int i = 0; i < nr; i++) 
		r[i] = rmin + i*dr;

	//compute cylindrical wave at time t:
	std::vector<double> p(nr);
	for (int i = 0; i < nr; ++i) p[i] = cylindrical_wave(r[i], t);


	//write data:
	write(r, p);


	return 0;
}
