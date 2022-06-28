// Read pressure data on a cylindrical surface and compute the far-field pressure:

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


// R3 vector:
template <typename T>
struct Vector{ T x, y, z;};


// add two Vectors:
template <typename T>
inline Vector<T> operator+(const Vector<T> &p1, const Vector<T> &p2){ 
	return {p1.x + p2.x, p1.y + p2.y, p1.z + p2.z};
}


// subtract two Vectors:
template <typename T>
inline Vector<T> operator-(const Vector<T> &p1, const Vector<T> &p2){ 
	return {p1.x - p2.x, p1.y - p2.y, p1.z - p2.z};
}


// dot product:
template <typename T>
T dot(const Vector<T> &p1, const Vector<T> &p2){ 
	return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}


// euclidean norm of a vector:
template <typename T>
T norm(const Vector<T> &p) { 
	return sqrt(dot(p, p)); 
}


// normalize vector:
template <typename T>
Vector<T> normalize(const Vector<T> &p){ 
	return {p.x / norm(p), p.y / norm(p), p.z / norm(p)}; 
}


// cos of angle between two vectors:
template <typename T>
T cos(const Vector<T> &p1, const Vector<T> &p2){ 
	return dot(normalize(p1), normalize(p2));
}


// print Vector:
template <typename T>
std::ostream &operator<<(std::ostream &os, const Vector<T> &p){ os << "(" << p.x << "," << p.y << "," << p.z << ")\n"; return os;}


// point in 3d space:
using Point = Vector<double>;


// emission time  τ = t − |x−y|/c:
inline double emission_time(const Point &x, const Point &y, double t, double c){ 
	return t - norm(x - y) / c;
}


// map r,theta,z -> x,y,z
inline Point cylinder_to_cartesian(double r, double theta, double z){ 
	return {r * cos(theta), r * sin(theta), z};
}

// map x,y,z -> r
inline double radius(double x, double y, double z){
	return sqrt(x*x + y*y);
}



// 2d 2-point gauss quadrature on a reference cell [-0.5,0.5]x[-0.5,0.5]
namespace Gauss_2d_2pt{

	static constexpr int nqpts = 4;
	static constexpr double psi[nqpts] = {-0.28867513459481287, -0.28867513459481287,  0.28867513459481287, 0.28867513459481287};
	static constexpr double eta[nqpts] = {-0.28867513459481287,  0.28867513459481287, -0.28867513459481287, 0.28867513459481287};
	static constexpr double w[nqpts] = {0.25, 0.25, 0.25, 0.25};

}; 


//linear map from reference to real cell:
inline double linear_map(double psi, double xc, double dx){
	return xc + psi*dx;
}



//Analytical solution for axisymmetric wave equation:
class CylindricalWave{

public:
    
    //constructor:
    CylindricalWave(double omega, double k, std::complex<double> A) : omega(omega), k(k), A(A), i(0.0, 1.0){}

    //pressure p
    double p(double r, double t){    
        return (A*boost::math::cyl_hankel_2(0, k*r)*exp(i*omega*t)).real();
    }

    //pressure derivative dp/dr:
    double dpdr(double r, double t){
        std::complex<double> hankel_derivative = k*0.5*(boost::math::cyl_hankel_2(-1, r) - boost::math::cyl_hankel_2(1, r));
        return (A*hankel_derivative*exp(i*omega*t)).real();
    }

    //pressure derivative dp/dt:
    double dpdt(double r, double t){
        return (i*omega*A*boost::math::cyl_hankel_2(0, k*r)*exp(i*omega*t)).real();   
    }

private:
    const double omega;
    const double k;
    const std::complex<double> A;
    const std::complex<double> i;
};





// Kirchhoff solver class:
class Kirchhoff
{

public:
	Kirchhoff(double c, double r_min, double r_max, double theta_min, double theta_max, double z_min, double z_max, double dh);
	void create_grid();
	void write_grid();
	void compute_emission_time(double, Point);
	void compute_pressure_at_emission_time(CylindricalWave);
    double compute_kirchhoff_integral();

private:
	void create_top_bottom_surface();
	void create_curved_surface();

private:
	double c;
	double r_min, r_max, theta_min, theta_max, z_min, z_max;
	double dh;
	int Ncells_r, Ncells_theta, Ncells_z;						
	double dr, dtheta, dz;	
	std::vector<double> r, theta, z;									
	boost::multi_array<double, 4> xq_top, xq_bottom, xq_curved; 
	boost::multi_array<double, 4> nq_top, nq_bottom, nq_curved; 
	boost::multi_array<double, 3> rq_top, rq_bottom, rq_curved;
	boost::multi_array<double, 3> cosq_top, cosq_bottom, cosq_curved;
	boost::multi_array<double, 3> tauq_top, tauq_bottom, tauq_curved;
	boost::multi_array<double, 3> pq_top, pq_bottom, pq_curved;
	boost::multi_array<double, 3> pnq_top, pnq_bottom, pnq_curved;
	boost::multi_array<double, 3> ptq_top, ptq_bottom, ptq_curved;

};
