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



//2d 2-point Lagrange basis function defined on [-1, 1]x[-1, 1]
namespace Lagrange_2d_2pt{

	// nodal basis functions:
	inline double N0(double psi, double eta) { return 0.25 * (1 - psi) * (1 - eta); }
	inline double N1(double psi, double eta) { return 0.25 * (1 + psi) * (1 - eta); }
	inline double N2(double psi, double eta) { return 0.25 * (1 + psi) * (1 + eta); }
	inline double N3(double psi, double eta) { return 0.25 * (1 - psi) * (1 + eta); }

	// interpolation on reference cell [-1, 1]x[-1, 1], given nodal values
	inline double interpolate(double psi, double eta, const std::array<double, 4> &p){
		return p[0] * N0(psi, eta) + p[1] * N1(psi, eta) + p[2] * N2(psi, eta) + p[3] * N3(psi, eta);
	}

};



// 2d 2-point gauss quadrature on a reference cell [-1, 1]x[-1, 1]
namespace Gauss_2d_2pt{

	static constexpr int nqpts = 4;
	static constexpr double psi[nqpts] = {-0.577350269189626, -0.577350269189626, 0.577350269189626, 0.577350269189626};
	static constexpr double eta[nqpts] = {-0.577350269189626, 0.577350269189626, -0.577350269189626, 0.577350269189626};
	static constexpr double w[nqpts] = {1.0, 1.0, 1.0, 1.0};

}; 



//1d 3-point Lagrange basis function defined on arbirtary domain
namespace Lagrange_1d_1pt{

	//nodal basis functions:
	inline double N0(double tau, const std::array<double, 3> &t){
		return ((tau - t[1]) * (tau - t[2])) / ((t[0] - t[1]) * (t[0] - t[2]));
	}

	inline double N1(double tau, const std::array<double, 3> &t){
		return ((tau - t[0]) * (tau - t[2])) / ((t[1] - t[0]) * (t[1] - t[2]));
	}

	inline double N2(double tau, const std::array<double, 3> &t){
		return ((tau - t[0]) * (tau - t[1])) / ((t[2] - t[0]) * (t[2] - t[1]));
	}

	// derivative of basis functions:
	inline double dN0(double tau, const std::array<double, 3> &t){
		return (tau - t[1]) / ((t[0] - t[1]) * (t[0] - t[2])) + (tau - t[2]) / ((t[0] - t[1]) * (t[0] - t[2]));
	}

	inline double dN1(double tau, const std::array<double, 3> &t){
		return (tau - t[0]) / ((t[1] - t[0]) * (t[1] - t[2])) + (tau - t[2]) / ((t[1] - t[0]) * (t[1] - t[2]));
	}

	inline double dN2(double tau, const std::array<double, 3> &t){
		return (tau - t[0]) / ((t[2] - t[0]) * (t[2] - t[1])) + (tau - t[1]) / ((t[2] - t[0]) * (t[2] - t[1]));
	}

	// interpolation:
	inline double interpolate(double tau, const std::array<double, 3> &t, const std::array<double, 3> &p){
		return p[0] * N0(tau, t) + p[1] * N1(tau, t) + p[2] * N2(tau, t);
	}

	// derivative of interpolation function:
	inline double derivative(double tau, const std::array<double, 3> &t, const std::array<double, 3> &p){
		return p[0] * dN0(tau, t) + p[1] * dN1(tau, t) + p[2] * dN2(tau, t);
	}

};




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
	double c = 1.0;											    		// wave speed

	// cylindrical domain:
	double r_min = 0.0;
	double r_max = 0.5;
	double theta_min = 0.0;
	double theta_max = 2.0 * M_PI;
	double z_min = -100.;
	double z_max =  100.;

	//cell size dh = dr = dtheta = dz
	double dh = 0.1;

	// no of nodes: ERROR!
	int Nnodes_r = (r_max - r_min)/dh;											// no of nodes along r direction
	int Nnodes_theta = (theta_max - theta_min)/dh;					// no of nodes along theta direction
	int Nnodes_z = (z_max - z_min)/dh;											// no of nodes along z direction

public:
	Kirchhoff();
	void create_grid();
	void write_grid();
	void compute_emission_time(double, Point);
	void compute_pressure_at_emission_time(CylindricalWave);
  double compute_kirchhoff_integral();

private:
	void create_top_bottom_surface();
	void create_curved_surface();

private:
	int Ncells_r, Ncells_theta, Ncells_z;						// no of cells along each direction
	double dr, dtheta, dz;										// cell size in each direction
	std::vector<double> r, theta, z;							// nodal coordinates along each direction
	boost::multi_array<double, 4> xq_top, xq_bottom, xq_curved; // coordinates of quad points located on top, bottom and curved surface of the cylinder
	boost::multi_array<double, 4> nq_top, nq_bottom, nq_curved; // normal vector at quad points located on top, bottom and curved surface of the cylinder
	boost::multi_array<double, 3> rq_top, rq_bottom, rq_curved;
	boost::multi_array<double, 3> cosq_top, cosq_bottom, cosq_curved;
	boost::multi_array<double, 3> tauq_top, tauq_bottom, tauq_curved;
	boost::multi_array<double, 3> pq_top, pq_bottom, pq_curved;
	boost::multi_array<double, 3> pnq_top, pnq_bottom, pnq_curved;
	boost::multi_array<double, 3> ptq_top, ptq_bottom, ptq_curved;

};
