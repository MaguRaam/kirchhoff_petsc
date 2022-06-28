// boost headers:
#include <boost/math/special_functions.hpp>
#include <boost/multi_array.hpp>

// c++ headers:
#include <array>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <string>
#include <sstream>


//pi
const long double pi = 3.141592653589793238462643L;


// 2d 2-point gauss quadrature on a reference cell [-0.5,0.5]x[-0.5,0.5]
namespace Gauss_2d_2pt{

	static constexpr int nqpts = 4;
	static constexpr double psi[nqpts] = {-0.28867513459481287, -0.28867513459481287,  0.28867513459481287, 0.28867513459481287};
	static constexpr double eta[nqpts] = {-0.28867513459481287,  0.28867513459481287, -0.28867513459481287, 0.28867513459481287};
	static constexpr double w[nqpts] = {0.25, 0.25, 0.25, 0.25};

}; 


//map from reference to real cell:
double x(double psi, double xc, double dx){
    return xc + psi*dx;
}

double y(double eta, double yc, double dy){
    return yc + eta*dy;
}


//compute area of a circle given radius and cell size 
double circle_area(double R, double dh){

    using namespace Gauss_2d_2pt;
	

    //circle domain:
    double r_min = 0.0;
	double r_max = R;
	double theta_min = 0.0;
	double theta_max = 2.0 * pi;

    //no of cells:
    int Ncells_r = (r_max - r_min)/dh;
    int Ncells_theta = (theta_max - theta_min)/dh;

    //cell size:
    double dr = (r_max - r_min) / static_cast<double>(Ncells_r);
    double dtheta = (theta_max - theta_min) / static_cast<double>(Ncells_theta);


    // compute nodal coordinates along each direction
    std::vector<double> r(Ncells_r), theta(Ncells_theta);			

    for (int i = 0; i < Ncells_r; ++i) r[i] = r_min + (i + 0.5)*dr;
	for (int j = 0; j < Ncells_theta; ++j) theta[j] = theta_min + (j + 0.5)*dtheta;

    //compute area:
    double cell_integral = 0.0;
	double global_integral = 0.0;


    //loop over cells:
	for (boost::multi_array<double, 3>::index i = 0; i < Ncells_r; ++i){
		for (boost::multi_array<double, 3>::index j = 0; j < Ncells_theta; ++j){

            
            //compute cell integral
            cell_integral = 0.0;

            //loop over quad points on the cell:
			for (boost::multi_array<double, 3>::index q = 0; q < nqpts; ++q){

                //get r and theta at quad points:
				double r_q = x(psi[q], r[i], dr);

                cell_integral += log(r_q)*r_q*dr*dtheta*w[q]; //TODO
            }

            //accumulate cell integral to global integral
            global_integral += cell_integral;
        }
    }

    
    return global_integral;

}



int main()
{
    double R = 1.0;     //Radius of the circle
    double dh = 0.1;   //grid size:

    //write area error:
    std::ofstream file("error.dat", std::ios::app);
    file.flags(std::ios::dec | std::ios::scientific);
    file.precision(16);
    if (!file)    {
    std::cerr << "Error: Output file couldnot be opened.\n";
    }

    for (int i = 1; i < 5; ++i){

        double exact_area = -1.5707963344618505e+00;

        double error = fabs(circle_area(R, dh) -  exact_area);

        file << dh << "\t" << error << std::endl;

        //reduce cell size to half
        dh = 0.5*dh;
    }

    
    return 0;
}


//0.25*R*R*(2.*log(R) - 1.0)*2.0*pi