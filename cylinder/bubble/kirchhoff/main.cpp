// Read pressure data on a cylindrical surface and compute the far-field pressure:

// boost headers:
#include <boost/multi_array.hpp>

// c++ headers:
#include <iostream>
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
#include <cassert>

// R3 vector:
template <typename T>
struct Vector
{
	T x, y, z;
};

// add two Vectors:
template <typename T>
inline Vector<T> operator+(const Vector<T> &p1, const Vector<T> &p2)
{
	return {p1.x + p2.x, p1.y + p2.y, p1.z + p2.z};
}

// subtract two Vectors:
template <typename T>
inline Vector<T> operator-(const Vector<T> &p1, const Vector<T> &p2)
{
	return {p1.x - p2.x, p1.y - p2.y, p1.z - p2.z};
}

// dot product:
template <typename T>
T dot(const Vector<T> &p1, const Vector<T> &p2)
{
	return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}

// euclidean norm of a vector:
template <typename T>
T norm(const Vector<T> &p)
{
	return sqrt(dot(p, p));
}

// normalize vector:
template <typename T>
Vector<T> normalize(const Vector<T> &p)
{
	return {p.x / norm(p), p.y / norm(p), p.z / norm(p)};
}

// cos of angle between two vectors:
template <typename T>
T cos(const Vector<T> &p1, const Vector<T> &p2)
{
	return dot(normalize(p1), normalize(p2));
}

// scalar multiplied by vector
template <typename T>
Vector<T> operator*(T alpha, const Vector<T> &p)
{
	return {p.x * alpha, p.y * alpha, p.z * alpha};
}

// print Vector:
template <typename T>
std::ostream &operator<<(std::ostream &os, const Vector<T> &p)
{
	os << "(" << p.x << "," << p.y << "," << p.z << ")\n";
	return os;
}

// point in 3d space:
using Point = Vector<double>;

// emission time  τ = t − |x−y|/c:
inline double emission_time(const Point &x, const Point &y, double t, double c)
{
	return t - norm(x - y) / c;
}

// map r,theta,z -> x,y,z
inline Point cylinder_to_cartesian(double r, double theta, double z)
{
	return {r * cos(theta), r * sin(theta), z};
}

// map x,y,z -> r
inline double radius(double x, double y, double z)
{
	return sqrt(x * x + y * y);
}

// 2d 2-point gauss quadrature on a reference cell [-0.5,0.5]x[-0.5,0.5]
namespace Gauss_2d_2pt
{

	static constexpr int nqpts = 4;
	static constexpr double psi[nqpts] = {-0.28867513459481287, 0.28867513459481287, -0.28867513459481287, 0.28867513459481287};
	static constexpr double eta[nqpts] = {-0.28867513459481287, -0.28867513459481287, 0.28867513459481287, 0.28867513459481287};
	static constexpr double w[nqpts] = {0.25, 0.25, 0.25, 0.25};

};

// linear map from reference to real cell:
inline double linear_map(double psi, double xc, double dx)
{
	return xc + psi * dx;
}

//linear time interpolation
inline double linear_interpolate(double tau, double t1, double t2, double p1, double p2)
{
	return p1 + ((p2 - p1) / (t2 - t1)) * (tau - t1);
}

//forward difference in time
inline double forward_difference(double t1, double t2, double p1, double p2)
{
	return (p2 - p1) / (t2 - t1);
}

//map 2d surface data to 1d line data:
inline int map(int i, int j, int q){
	if (q < 2)
		return 2*i + q;
	else
		return 2*(i - 1) + q ;
}


inline int map_curved(int i, int j, int q){
	if (q < 2)
		return 2*j;
	else
		return 2*j + 1;
}


//read pressure and pressure derivative from a file:
void read_pressure_data(std::string filename, std::vector<double> &p, std::vector<double> &pn)
{

	std::string line;
	std::ifstream file(filename);

	// check if file is open:
	if (!file.is_open())
	{
		std::cerr << "pressure file is not open\n"
				  << std::endl;
		exit(1);
	};

	// skip first two lines:
	getline(file, line); // Vec Object: 4 MPI processes
	getline(file, line); // type: mpi

	//read pressure vector:
	while (getline(file, line) && (line.find("Vec") == std::string::npos)){

		// if there is a process key word skip that line:
		if (line.find("Process") != std::string::npos)
			continue;

		// push back values to vec1
		std::istringstream iss(line);

		double val;
		iss >> val;
		p.push_back(val);
	}

	getline(file, line); // type: mpi

	//read pressure derivative vector
	while (getline(file, line)){

		// if there is a process key word skip that line:
		if (line.find("Process") != std::string::npos)
			continue;

		// push back values to vec1
		std::istringstream iss(line);

		double val;
		iss >> val;
		pn.push_back(val);
	}
}

// Kirchhoff solver class:
class Kirchhoff
{

public:
	Kirchhoff(double c, double r_min, double r_max, double theta_min, double theta_max, double z_min, double z_max,
			  int Ncells_r, int Ncells_theta, int Ncells_z, double dh, int write_interval);
	void create_grid();
	void compute_emission_time(double, Point);
	void compute_pressure_at_emission_time();
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
	int write_interval;
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

// constructor
Kirchhoff::Kirchhoff(double c, double r_min, double r_max, double theta_min, double theta_max, double z_min, double z_max, int Ncells_r, int Ncells_theta, int Ncells_z, double dh, int write_interval) : c(c),
																																													  r_min(r_min),
																																													  r_max(r_max),
																																													  theta_min(theta_min),
																																													  theta_max(theta_max),
																																													  z_min(z_min),
																																													  z_max(z_max),
																																													  dh(dh),
																																													  Ncells_r(Ncells_r),
																																													  Ncells_theta(Ncells_theta),
																																													  Ncells_z(Ncells_z),

																																													  r(Ncells_r),
																																													  theta(Ncells_theta),
																																													  z(Ncells_z),

																																													  write_interval(write_interval),

																																													  dr(dh),
																																													  dtheta((theta_max - theta_min) / static_cast<double>(Ncells_theta)),
																																													  dz(dh),
																																													  xq_top(boost::extents[Ncells_r][Ncells_theta][Gauss_2d_2pt::nqpts][3]),
																																													  xq_bottom(boost::extents[Ncells_r][Ncells_theta][Gauss_2d_2pt::nqpts][3]),
																																													  xq_curved(boost::extents[Ncells_theta][Ncells_z][Gauss_2d_2pt::nqpts][3]),

																																													  nq_top(boost::extents[Ncells_r][Ncells_theta][Gauss_2d_2pt::nqpts][3]),
																																													  nq_bottom(boost::extents[Ncells_r][Ncells_theta][Gauss_2d_2pt::nqpts][3]),
																																													  nq_curved(boost::extents[Ncells_theta][Ncells_z][Gauss_2d_2pt::nqpts][3]),

																																													  rq_top(boost::extents[Ncells_r][Ncells_theta][Gauss_2d_2pt::nqpts]),
																																													  rq_bottom(boost::extents[Ncells_r][Ncells_theta][Gauss_2d_2pt::nqpts]),
																																													  rq_curved(boost::extents[Ncells_theta][Ncells_z][Gauss_2d_2pt::nqpts]),

																																													  cosq_top(boost::extents[Ncells_r][Ncells_theta][Gauss_2d_2pt::nqpts]),
																																													  cosq_bottom(boost::extents[Ncells_r][Ncells_theta][Gauss_2d_2pt::nqpts]),
																																													  cosq_curved(boost::extents[Ncells_theta][Ncells_z][Gauss_2d_2pt::nqpts]),

																																													  tauq_top(boost::extents[Ncells_r][Ncells_theta][Gauss_2d_2pt::nqpts]),
																																													  tauq_bottom(boost::extents[Ncells_r][Ncells_theta][Gauss_2d_2pt::nqpts]),
																																													  tauq_curved(boost::extents[Ncells_theta][Ncells_z][Gauss_2d_2pt::nqpts]),

																																													  pq_top(boost::extents[Ncells_r][Ncells_theta][Gauss_2d_2pt::nqpts]),
																																													  pq_bottom(boost::extents[Ncells_r][Ncells_theta][Gauss_2d_2pt::nqpts]),
																																													  pq_curved(boost::extents[Ncells_theta][Ncells_z][Gauss_2d_2pt::nqpts]),

																																													  pnq_top(boost::extents[Ncells_r][Ncells_theta][Gauss_2d_2pt::nqpts]),
																																													  pnq_bottom(boost::extents[Ncells_r][Ncells_theta][Gauss_2d_2pt::nqpts]),
																																													  pnq_curved(boost::extents[Ncells_theta][Ncells_z][Gauss_2d_2pt::nqpts]),

																																													  ptq_top(boost::extents[Ncells_r][Ncells_theta][Gauss_2d_2pt::nqpts]),
																																													  ptq_bottom(boost::extents[Ncells_r][Ncells_theta][Gauss_2d_2pt::nqpts]),
																																													  ptq_curved(boost::extents[Ncells_theta][Ncells_z][Gauss_2d_2pt::nqpts])
{
	for (int i = 0; i < Ncells_r; ++i)
		r[i] = r_min + (i + 0.5) * dr;
	for (int j = 0; j < Ncells_theta; ++j)
		theta[j] = theta_min + (j + 0.5) * dtheta;
	for (int k = 0; k < Ncells_z; ++k)
		z[k] = z_min + (k + 0.5) * dz;
}

// create top surface of the cylinder:
void Kirchhoff::create_top_bottom_surface()
{

	using namespace Gauss_2d_2pt;

	// get coordinates of quad points and normal vector on top and bottom surface of the cylinder:

	// loop over cells on top and bottom surface of the cylinder
	for (boost::multi_array<double, 4>::index i = 0; i < Ncells_r; ++i)
	{
		for (boost::multi_array<double, 4>::index j = 0; j < Ncells_theta; ++j)
		{

			// loop over quad points on the cell:
			for (boost::multi_array<double, 4>::index q = 0; q < nqpts; ++q)
			{

				// get r and theta at quad points:
				double r_q = linear_map(psi[q], r[i], dr);
				double theta_q = linear_map(eta[q], theta[j], dtheta);

				// top surface:

				// get the cartesian coordinate of the quad point:
				Point point_q = cylinder_to_cartesian(r_q, theta_q, z_max);

				// get coordinates of quadrature point:
				xq_top[i][j][q][0] = point_q.x;
				xq_top[i][j][q][1] = point_q.y;
				xq_top[i][j][q][2] = point_q.z;

				// get normal vector at quadrature point:
				nq_top[i][j][q][0] = 0.0;
				nq_top[i][j][q][1] = 0.0;
				nq_top[i][j][q][2] = -1.0;

				// bottom surface:

				// get the cartesian coordinate of the quad point:
				point_q = cylinder_to_cartesian(r_q, theta_q, z_min);

				// get coordinates of quadrature point:
				xq_bottom[i][j][q][0] = point_q.x;
				xq_bottom[i][j][q][1] = point_q.y;
				xq_bottom[i][j][q][2] = point_q.z;

				// get normal vector at quadrature point:
				nq_bottom[i][j][q][0] = 0.0;
				nq_bottom[i][j][q][1] = 0.0;
				nq_bottom[i][j][q][2] = 1.0;

			} // q points loop

		} // r loop
	}	  // theta loop
}

// create bottom surface of the cylinder:
void Kirchhoff::create_curved_surface()
{

	using namespace Gauss_2d_2pt;

	// get coordinates of quad points and normal vector on curved surface of the cylinder:

	// loop over cells on curved surface of the cylinder:
	for (boost::multi_array<double, 4>::index j = 0; j < Ncells_theta; ++j)
	{
		for (boost::multi_array<double, 4>::index k = 0; k < Ncells_z; ++k)
		{

			// loop over quad points on the cell:
			for (boost::multi_array<double, 4>::index q = 0; q < nqpts; ++q)
			{

				// get theta and z at quad points:
				double theta_q = linear_map(psi[q], theta[j], dtheta);
				double z_q = linear_map(eta[q], z[k], dz);

				// get the cartesian coordinate of the quad point:
				Point point_q = cylinder_to_cartesian(r_max, theta_q, z_q);

				// get coordinates of quadrature point:
				xq_curved[j][k][q][0] = point_q.x;
				xq_curved[j][k][q][1] = point_q.y;
				xq_curved[j][k][q][2] = point_q.z;

				// get normal vector at quadrature point:
				nq_curved[j][k][q][0] = cos(theta_q);
				nq_curved[j][k][q][1] = sin(theta_q);
				nq_curved[j][k][q][2] = 0.0;
			}

		} // theta loop
	}	  // z loop
}

// create cylindrical surface grid
void Kirchhoff::create_grid()
{
	create_top_bottom_surface();
	create_curved_surface();
}


// Get emission time τ = to − r/c, distance r = |xo−y| and cos of angle between r = xo−y and n vectors at each quadpoint:
void Kirchhoff::compute_emission_time(double to, Point xo)
{

	using namespace Gauss_2d_2pt;

	// loop over cells on top and bottom surface of the cylinder
	for (boost::multi_array<double, 3>::index i = 0; i < Ncells_r; ++i)
	{
		for (boost::multi_array<double, 3>::index j = 0; j < Ncells_theta; ++j)
		{

			// loop over quad points on the cell:
			for (boost::multi_array<double, 3>::index q = 0; q < nqpts; ++q)
			{

				// top:

				// get coordinate and normal vector at quadpoint:
				Point y{xq_top[i][j][q][0], xq_top[i][j][q][1], xq_top[i][j][q][2]};
				Vector<double> n{nq_top[i][j][q][0], nq_top[i][j][q][1], nq_top[i][j][q][2]};

				// compute distance, emission time and cos:
				rq_top[i][j][q] = norm(xo - y);
				tauq_top[i][j][q] = emission_time(xo, y, to, c);
				cosq_top[i][j][q] = cos(xo - y, n);

				// bottom:

				// get coordinate and normal vector at quadpoint
				y = Point{xq_bottom[i][j][q][0], xq_bottom[i][j][q][1], xq_bottom[i][j][q][2]};
				n = Vector<double>{nq_bottom[i][j][q][0], nq_bottom[i][j][q][1], nq_bottom[i][j][q][2]};

				// compute distance, emission time and cos:
				rq_bottom[i][j][q] = norm(xo - y);
				tauq_bottom[i][j][q] = emission_time(xo, y, to, c);
				cosq_bottom[i][j][q] = cos(xo - y, n);
			}
		}
	}

	// loop over cells on curved surface of the cylinder
	for (boost::multi_array<double, 3>::index j = 0; j < Ncells_theta; ++j)
	{
		for (boost::multi_array<double, 3>::index k = 0; k < Ncells_z; ++k)
		{

			// loop over quad points on the cell:
			for (boost::multi_array<double, 3>::index q = 0; q < Gauss_2d_2pt::nqpts; ++q)
			{

				Point y{xq_curved[j][k][q][0], xq_curved[j][k][q][1], xq_curved[j][k][q][2]};
				Vector<double> n{nq_curved[j][k][q][0], nq_curved[j][k][q][1], nq_curved[j][k][q][2]};

				rq_curved[j][k][q] = norm(xo - y);
				tauq_curved[j][k][q] = emission_time(xo, y, to, c);
				cosq_curved[j][k][q] = cos(xo - y, n);
			}
		}
	}
}

// compute pressure at emission time:
void Kirchhoff::compute_pressure_at_emission_time()
{
	//open t.dat:
	std::ifstream time_file("input/t.dat");

	// throw file opening error:
	if (!time_file.is_open()){ std::cerr << "t.dat is not open\n" << std::endl;exit(1);};

	//time vector and file_index
	std::vector<double> t;
	std::vector<int> file_index;
 
	int index = 0;
	while (time_file.good()){

		double value;
		time_file >> value;
		t.push_back(value);
		file_index.push_back(index);
		index += write_interval;
	}

	//loop over files (except last files : TODO think about this)

	for (int n = 0; n < t.size() - 2; ++n){
		
		//vector to store data at current time step:
		std::vector<double> p_top1, p_curved1, p_bottom1, pz_top1, pr_curved1, pz_bottom1;

		//vector to store data at next time step:
		std::vector<double> p_top2, p_curved2, p_bottom2, pz_top2, pr_curved2, pz_bottom2;

		//char to store current filename and next filename:
		char current_file[50], next_file[50];

		//TODO the padding level in file should mathch the flow solver or will throw error:

		//read top data:
		sprintf(current_file, "input/top/p-%09d.dat", file_index[n]);
		sprintf(next_file, "input/top/p-%09d.dat", file_index[n + 1]);
		read_pressure_data(current_file, p_top1, pz_top1);
		read_pressure_data(next_file, p_top2, pz_top2);

		//read bottom data:
		sprintf(current_file, "input/bottom/p-%09d.dat", file_index[n]);
		sprintf(next_file, "input/bottom/p-%09d.dat", file_index[n + 1]);
		read_pressure_data(current_file, p_bottom1, pz_bottom1);
		read_pressure_data(next_file, p_bottom2, pz_bottom2);

		//read curved data:
		sprintf(current_file, "input/curved/p-%09d.dat", file_index[n]);
		sprintf(next_file, "input/curved/p-%09d.dat", file_index[n + 1]);
		read_pressure_data(current_file, p_curved1, pr_curved1);
		read_pressure_data(next_file, p_curved2, pr_curved2);

		//interpolate data at emission time:

		//loop over cells on top and bottom surface of the cylinder
		for (boost::multi_array<double, 3>::index i = 0; i < Ncells_r; ++i){
			for (boost::multi_array<double, 3>::index j = 0; j < Ncells_theta; ++j){

				// loop over quad points on the cell:
				for (boost::multi_array<double, 3>::index q = 0; q < Gauss_2d_2pt::nqpts; ++q){

					//find the quadrature points that has emission time between current and next time interval:
					if (tauq_top[i][j][q] > t[n] && tauq_top[i][j][q] < t[n + 1]){
						
						//map the top line data to top surface data
						pq_top[i][j][q] = linear_interpolate(tauq_top[i][j][q], t[n], t[n + 1], p_top1[map(i, j, q)], p_top2[map(i, j, q)]);
						pnq_top[i][j][q] = linear_interpolate(tauq_top[i][j][q], t[n], t[n + 1], pz_top1[map(i, j, q)], pz_top2[map(i, j, q)]);
						ptq_top[i][j][q] = forward_difference(t[n], t[n + 1], p_top1[map(i, j, q)], p_top2[map(i, j, q)] );

					}

					if (tauq_bottom[i][j][q] > t[n] && tauq_bottom[i][j][q] < t[n + 1]){
						
						//map the bottom line data to bottom surface data
						pq_bottom[i][j][q] = linear_interpolate(tauq_bottom[i][j][q], t[n], t[n + 1], p_bottom1[map(i, j, q)], p_bottom2[map(i, j, q)]);
						pnq_bottom[i][j][q] = linear_interpolate(tauq_bottom[i][j][q], t[n], t[n + 1], pz_bottom1[map(i, j, q)], pz_bottom2[map(i, j, q)]);
						ptq_bottom[i][j][q] = forward_difference(t[n], t[n + 1], p_bottom1[map(i, j, q)], p_bottom2[map(i, j, q)] );

					}



				}
			}
		}	

		// loop over cells on curved surface of the cylinder
		for (boost::multi_array<double, 3>::index i = 0; i < Ncells_theta; ++i){
			for (boost::multi_array<double, 3>::index j = 0; j < Ncells_z; ++j){
				
				// loop over quad points on the cell:
				for (boost::multi_array<double, 3>::index q = 0; q < Gauss_2d_2pt::nqpts; ++q){
					if ((tauq_curved[i][j][q] > t[n]) && (tauq_curved[i][j][q] < t[n + 1]))
					{
						//map the curved line data to curved surface data
						pq_curved[i][j][q] = linear_interpolate(tauq_curved[i][j][q], t[n], t[n + 1], p_curved1[map_curved(i, j, q)], p_curved2[map_curved(i, j, q)]);
						pnq_curved[i][j][q] = linear_interpolate(tauq_curved[i][j][q], t[n], t[n + 1], pr_curved1[map_curved(i, j, q)], pr_curved2[map_curved(i, j, q)]);
						ptq_curved[i][j][q] = forward_difference(t[n], t[n + 1], p_curved1[map_curved(i, j, q)], p_curved2[map_curved(i, j, q)] );
					}
				}
			}
		}

	}

}



// compute Kirchhoff integral:
double Kirchhoff::compute_kirchhoff_integral()
{

	using namespace Gauss_2d_2pt;

	double cell_integral = 0.0;
	double global_integral = 0.0;
	double cinv = 1.0 / c;

	// loop over cells on top and bottom surface of the cylinder
	for (boost::multi_array<double, 3>::index i = 0; i < Ncells_r; ++i)
	{
		for (boost::multi_array<double, 3>::index j = 0; j < Ncells_theta; ++j)
		{

			// compute cell integral
			cell_integral = 0.0;

			// loop over quad points on the cell:
			for (boost::multi_array<double, 3>::index q = 0; q < nqpts; ++q)
			{

				double radius_q = linear_map(psi[q], r[i], dr);

				cell_integral += (0.25 / M_PI) * ((cinv * ptq_top[i][j][q] * cosq_top[i][j][q] - pnq_top[i][j][q]) / rq_top[i][j][q] + (pq_top[i][j][q] * cosq_top[i][j][q]) / (rq_top[i][j][q] * rq_top[i][j][q])) * radius_q * dr * dtheta * w[q];

				cell_integral += (0.25 / M_PI) * ((cinv * ptq_bottom[i][j][q] * cosq_bottom[i][j][q] - pnq_bottom[i][j][q]) / rq_bottom[i][j][q] + (pq_bottom[i][j][q] * cosq_bottom[i][j][q]) / (rq_bottom[i][j][q] * rq_bottom[i][j][q])) * radius_q * dr * dtheta * w[q];
			}

			global_integral += cell_integral;
		}
	}

	// loop over cells on curved surface of the cylinder
	for (boost::multi_array<double, 3>::index j = 0; j < Ncells_theta; ++j)
	{
		for (boost::multi_array<double, 3>::index k = 0; k < Ncells_z; ++k)
		{

			// compute cell integral
			cell_integral = 0.0;

			// loop over quad points on the cell:
			for (boost::multi_array<double, 3>::index q = 0; q < nqpts; ++q)
				cell_integral += (0.25 / M_PI) * ((cinv * ptq_curved[j][k][q] * cosq_curved[j][k][q] - pnq_curved[j][k][q]) / rq_curved[j][k][q] + (pq_curved[j][k][q] * cosq_curved[j][k][q]) / (rq_curved[j][k][q] * rq_curved[j][k][q])) * w[q] * r_max * dtheta * dz;

			global_integral += cell_integral;
		}
	}

	return global_integral;
}



int main()
{

	// set simulation parameters

	// wave speed
	const double c = 66.53119569044284;

	// observer time
	double t0 = 0.004;
	const double tf = 0.007;
	const double dt = 0.00001;

	// no of cells along different direction
	const int Ncells_r = 151;
	const int Ncells_theta = 151;  //this is completely arbirtary
	const int Ncells_z = 451;

	// observer point
	Point x0{0.34276, 0.0, 0.0007600000000000384};

	// cylinder geometry
	const double r_min = 0.0;
	const double r_max = 0.22952;
	const double theta_min = 0.0;
	const double theta_max = 2.0 * M_PI;
	const double z_min = -0.342;
	const double z_max = 0.34352000000000005;

	// cell size
	const double dh = 0.00152;

	// file write interval from flow solver:
	const int write_interval = 20;

	Kirchhoff solver(c, r_min, r_max, theta_min, theta_max, z_min, z_max, Ncells_r, Ncells_theta, Ncells_z, dh, write_interval);

	solver.create_grid();

	// write Kirchhoff data:
	std::ofstream observer("p_kirchhoff.dat", std::ios::out);
	observer.flags(std::ios::dec | std::ios::scientific);
	observer.precision(16);
	if (!observer)
		std::cerr << "Error: Output file couldnot be opened.\n";

	// loop over observer time:
	while (t0 <= tf)
	{
		std::cout << "observer time t0 = " << t0 << std::endl;

		solver.compute_emission_time(t0, x0);
		solver.compute_pressure_at_emission_time();
		observer << t0 << "\t" << solver.compute_kirchhoff_integral() << std::endl;
		t0 += dt;
	}

	return 0;
}
