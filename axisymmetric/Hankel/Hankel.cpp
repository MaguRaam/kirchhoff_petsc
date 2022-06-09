//Computes cylindrical wave analytical solution using the Hankel function and writes the pressure data for Kirchhoff solver

// boost headers:
#include <boost/math/special_functions.hpp>

// c++ headers:
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>
#include <cstdlib>





class CylindricalWave{

public:
    
    //constructor:
    CylindricalWave(double omega, double k, std::complex<double> A) : omega(omega), k(k), A(A), i(0.0, 1.0){}


    //pressure p
    double p(double r, double t)
    {    
        return (A*boost::math::cyl_hankel_2(0, k*r)*exp(i*omega*t)).real();
    }


    //pressure derivative dp/dr:
    double dpdr(double r, double t)
    {
        std::complex<double> hankel_derivative = k*0.5*(boost::math::cyl_hankel_2(-1, r) - boost::math::cyl_hankel_2(1, r));
        return (A*hankel_derivative*exp(i*omega*t)).real();
    }

    //pressure derivative dp/dt:
    double dpdt(double r, double t)
    {
        return (i*omega*A*boost::math::cyl_hankel_2(0, k*r)*exp(i*omega*t)).real();   
    }



private:
    const double omega;
    const double k;
    const std::complex<double> A;
    const std::complex<double> i;
};


int main(){

    //cylindrical wave parameters:
    double omega = 1.0;                     //frequency
    double k     = 1.0;                     //wave number
    double c     = omega/k;                 //wave speed
    std::complex<double> A(1.0, 0.0);       //amplitude
    CylindricalWave wave(omega, k, A);      //cylindrical wave object 


    //cylindrical Kirchhoff surface radius:
    double r = 1.0;                        
    

    //observer radius:
    double ro = 2.0;


    //start time, end time and time step:
    double t   = 0.0;
    double tf  = 10.0;
    double dt  = 0.001;



    //write pressure data for kirchhoff solver:

    std::ofstream kirchhoff("kirchhoff.dat", std::ios::out) ;
    kirchhoff.flags( std::ios::dec | std::ios::scientific );
    kirchhoff.precision(16) ;
    if(!kirchhoff) {std::cerr<< "Error: Output file couldnot be opened.\n";}


    //write pressure data at observer point ro:

    std::ofstream observer("p_exact.dat", std::ios::out) ;
    observer.flags( std::ios::dec | std::ios::scientific );
    observer.precision(16) ;
    if(!observer) {std::cerr<< "Error: Output file couldnot be opened.\n";}


    //loop over time:
    while (t < tf){


        //write p, dp/dr, dp/dt at Kirchhoff surface r:
        kirchhoff << t << "\t" << wave.p(r, t) << "\t" << wave.dpdr(r, t) << "\t" << wave.dpdt(r, t) << std::endl;

        //write pressure at observer radius ro:
        observer << t << "\t" << wave.p(ro, t) << std::endl;


        t+=dt;
    }

    kirchhoff.close();
    observer.close();

    return 0;
}
