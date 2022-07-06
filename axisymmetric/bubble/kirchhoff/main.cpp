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
	static constexpr double psi[nqpts] = {-0.28867513459481287, -0.28867513459481287, 0.28867513459481287, 0.28867513459481287};
	static constexpr double eta[nqpts] = {-0.28867513459481287, 0.28867513459481287, -0.28867513459481287, 0.28867513459481287};
	static constexpr double w[nqpts] = {0.25, 0.25, 0.25, 0.25};

};

// linear map from reference to real cell:
inline double linear_map(double psi, double xc, double dx)
{
	return xc + psi * dx;
}

// Kirchhoff solver class:
class Kirchhoff
{

public:
	Kirchhoff(double c, double r_min, double r_max, double theta_min, double theta_max, double z_min, double z_max,
			  int Ncells_r, int Ncells_theta, int Ncells_z, double dh);
	void create_grid();
	void write_grid();
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
Kirchhoff::Kirchhoff(double c, double r_min, double r_max, double theta_min, double theta_max, double z_min, double z_max, int Ncells_r, int Ncells_theta, int Ncells_z, double dh) : c(c),
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
		r[i] = (i + 0.5) * dr;
	for (int j = 0; j < Ncells_theta; ++j)
		theta[j] = theta_min + (j + 0.5) * dtheta;
	for (int k = 0; k < Ncells_z; ++k)
		z[k] = z_min + k * dz;
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

// write quad points on the cylinder to a file to be visualized using cylinder.py
void Kirchhoff::write_grid()
{

	std::ofstream file("cylinder.dat", std::ios::out);
	file.flags(std::ios::dec | std::ios::scientific);
	file.precision(16);
	if (!file)
	{
		std::cerr << "Error: cylinder.dat could not be opened.\n";
	}

	// write top and bottom surface quad points:
	for (boost::multi_array<double, 4>::index i = 0; i < Ncells_r; ++i)
	{
		for (boost::multi_array<double, 4>::index j = 0; j < Ncells_theta; ++j)
		{
			for (boost::multi_array<double, 4>::index q = 0; q < Gauss_2d_2pt::nqpts; ++q)
			{

				file << xq_top[i][j][q][0] << "\t" << xq_top[i][j][q][1] << "\t" << xq_top[i][j][q][2] << "\n";
				file << xq_bottom[i][j][q][0] << "\t" << xq_bottom[i][j][q][1] << "\t" << xq_bottom[i][j][q][2] << "\n";
			}
		}
	}

	// write curved surface quad points:
	for (boost::multi_array<double, 4>::index j = 0; j < Ncells_theta; ++j)
		for (boost::multi_array<double, 4>::index k = 0; k < Ncells_z; ++k)
			for (boost::multi_array<double, 4>::index q = 0; q < Gauss_2d_2pt::nqpts; ++q)
				file << xq_curved[j][k][q][0] << "\t" << xq_curved[j][k][q][1] << "\t" << xq_curved[j][k][q][2] << "\n";
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

void read_pressure_file(char *filename, std::vector<double> &f, std::vector<double> &fz)
{
	std::fstream f1(filename, std::ios::in);
	int line_number = 0;
	int vec_object_counter = 0;
	std::string line;
	while (std::getline(f1, line))
	{
		std::istringstream iss(line);

		size_t found = line.find("Vec Object");
		if (found != std::string::npos)
		{
			line_number = 0;
			vec_object_counter++;
		}

		if (vec_object_counter == 1)
		{
			if (line_number > 1)
			{
				size_t found = line.find("Process");
				if (found == std::string::npos)
				{
					double tmp2;
					iss >> tmp2;
					f.push_back(tmp2);
				}
			}
		}
		if (vec_object_counter == 2)
		{
			if (line_number > 1)
			{
				size_t found = line.find("Process");
				if (found == std::string::npos)
				{
					double tmp2;
					iss >> tmp2;
					fz.push_back(tmp2);
				}
			}
		}

		line_number++;
	}
}

inline double LinearInterpolate(double tau, double t1, double t2, double p1, double p2)
{
	return p1 + ((p2 - p1) / (t2 - t1)) * (tau - t1);
}

inline double ForwardDifference(double t1, double t2, double p1, double p2)
{
	return (p2 - p1) / (t2 - t1);
}

// compute pressure at emission time:
void Kirchhoff::compute_pressure_at_emission_time()
{

	// open t.dat and dt.dat file
	std::fstream f1("input/t.dat", std::ios::in);
	std::fstream f2("input/dt.dat", std::ios::in);

	// throw file opening error:
	std::vector<double> t;
	std::vector<int> it;
	// std::vector<std::string> filename;

	double tmp1;
	int it_ = 0;

	while (f1.good())
	{
		// get t and t+dt
		f1 >> tmp1;
		t.push_back(tmp1);
		it.push_back(it_);
		it_ += 500;
	}

	for (int n = 0; n < it.size() - 2; ++n)
	{
		std::cout << "n " << n << "\n";
		std::vector<double> p_top1, p_curved1, p_bottom1;
		std::vector<double> pz_top1, pr_curved1, pz_bottom1;

		std::vector<double> p_top2, p_curved2, p_bottom2;
		std::vector<double> pz_top2, pr_curved2, pz_bottom2;

		char filename1[50], filename2[50];
		sprintf(filename1, "input/top/p-%05d.dat", it[n]);
		sprintf(filename2, "input/top/p-%05d.dat", it[n + 1]);

		read_pressure_file(filename1, p_top1, pz_top1);
		read_pressure_file(filename2, p_top2, pz_top2);

		sprintf(filename1, "input/bottom/p-%05d.dat", it[n]);
		sprintf(filename2, "input/bottom/p-%05d.dat", it[n + 1]);

		read_pressure_file(filename1, p_bottom1, pz_bottom1);
		read_pressure_file(filename2, p_bottom2, pz_bottom2);

		sprintf(filename1, "input/curved/p-%05d.dat", it[n]);
		sprintf(filename2, "input/curved/p-%05d.dat", it[n + 1]);

		read_pressure_file(filename1, p_curved1, pr_curved1);
		read_pressure_file(filename2, p_curved2, pr_curved2);

		boost::multi_array<double, 3> p1_top_2d(boost::extents[Ncells_r][Ncells_theta][Gauss_2d_2pt::nqpts]),
			p1_bottom_2d(boost::extents[Ncells_r][Ncells_theta][Gauss_2d_2pt::nqpts]),
			p1_curved_2d(boost::extents[Ncells_theta][Ncells_z][Gauss_2d_2pt::nqpts]),
			pz1_top_2d(boost::extents[Ncells_r][Ncells_theta][Gauss_2d_2pt::nqpts]),
			pz1_bottom_2d(boost::extents[Ncells_r][Ncells_theta][Gauss_2d_2pt::nqpts]),
			pr1_curved_2d(boost::extents[Ncells_theta][Ncells_z][Gauss_2d_2pt::nqpts]);

		boost::multi_array<double, 3> p2_top_2d(boost::extents[Ncells_r][Ncells_theta][Gauss_2d_2pt::nqpts]),
			p2_bottom_2d(boost::extents[Ncells_r][Ncells_theta][Gauss_2d_2pt::nqpts]),
			p2_curved_2d(boost::extents[Ncells_theta][Ncells_z][Gauss_2d_2pt::nqpts]),
			pz2_top_2d(boost::extents[Ncells_r][Ncells_theta][Gauss_2d_2pt::nqpts]),
			pz2_bottom_2d(boost::extents[Ncells_r][Ncells_theta][Gauss_2d_2pt::nqpts]),
			pr2_curved_2d(boost::extents[Ncells_theta][Ncells_z][Gauss_2d_2pt::nqpts]);

		// loop over r direction:
		for (boost::multi_array<double, 3>::index j = 0; j < Ncells_theta; ++j)
		{
			for (boost::multi_array<double, 3>::index i = 0; i < Ncells_r; ++i)
			{

				// copy top data:
				p1_top_2d[i][j][0] = p_top1[2 * i];
				p1_top_2d[i][j][1] = p_top1[2 * i + 1];
				p1_top_2d[i][j][2] = p_top1[2 * i + 1];
				p1_top_2d[i][j][3] = p_top1[2 * i];

				pz1_top_2d[i][j][0] = pz_top1[2 * i];
				pz1_top_2d[i][j][1] = pz_top1[2 * i + 1];
				pz1_top_2d[i][j][2] = pz_top1[2 * i + 1];
				pz1_top_2d[i][j][3] = pz_top1[2 * i];

				// copy top data:
				p1_bottom_2d[i][j][0] = p_bottom1[2 * i];
				p1_bottom_2d[i][j][1] = p_bottom1[2 * i + 1];
				p1_bottom_2d[i][j][2] = p_bottom1[2 * i + 1];
				p1_bottom_2d[i][j][3] = p_bottom1[2 * i];

				pz1_bottom_2d[i][j][0] = pz_bottom1[2 * i];
				pz1_bottom_2d[i][j][1] = pz_bottom1[2 * i + 1];
				pz1_bottom_2d[i][j][2] = pz_bottom1[2 * i + 1];
				pz1_bottom_2d[i][j][3] = pz_bottom1[2 * i];

				// copy top data:
				p2_top_2d[i][j][0] = p_top2[2 * i];
				p2_top_2d[i][j][1] = p_top2[2 * i + 1];
				p2_top_2d[i][j][2] = p_top2[2 * i + 1];
				p2_top_2d[i][j][3] = p_top2[2 * i];

				pz2_top_2d[i][j][0] = pz_top2[2 * i];
				pz2_top_2d[i][j][1] = pz_top2[2 * i + 1];
				pz2_top_2d[i][j][2] = pz_top2[2 * i + 1];
				pz2_top_2d[i][j][3] = pz_top2[2 * i];

				// copy top data:
				p2_bottom_2d[i][j][0] = p_bottom2[2 * i];
				p2_bottom_2d[i][j][1] = p_bottom2[2 * i + 1];
				p2_bottom_2d[i][j][2] = p_bottom2[2 * i + 1];
				p2_bottom_2d[i][j][3] = p_bottom2[2 * i];

				pz2_bottom_2d[i][j][0] = pz_bottom2[2 * i];
				pz2_bottom_2d[i][j][1] = pz_bottom2[2 * i + 1];
				pz2_bottom_2d[i][j][2] = pz_bottom2[2 * i + 1];
				pz2_bottom_2d[i][j][3] = pz_bottom2[2 * i];
			}
		}

		// loop over cells on curved surface of the cylinder:
		for (boost::multi_array<double, 3>::index j = 0; j < Ncells_theta; ++j)
		{
			for (boost::multi_array<double, 3>::index k = 0; k < Ncells_z; ++k)
			{
				// copy top data:
				p1_curved_2d[j][k][0] = p_curved1[2 * j];
				p1_curved_2d[j][k][1] = p_curved1[2 * j + 1];
				p1_curved_2d[j][k][2] = p_curved1[2 * j + 1];
				p1_curved_2d[j][k][3] = p_curved1[2 * j];

				pr1_curved_2d[j][k][0] = pr_curved1[2 * j];
				pr1_curved_2d[j][k][1] = pr_curved1[2 * j + 1];
				pr1_curved_2d[j][k][2] = pr_curved1[2 * j + 1];
				pr1_curved_2d[j][k][3] = pr_curved1[2 * j];

				p2_curved_2d[j][k][0] = p_curved2[2 * j];
				p2_curved_2d[j][k][1] = p_curved2[2 * j + 1];
				p2_curved_2d[j][k][2] = p_curved2[2 * j + 1];
				p2_curved_2d[j][k][3] = p_curved2[2 * j];

				pr2_curved_2d[j][k][0] = pr_curved2[2 * j];
				pr2_curved_2d[j][k][1] = pr_curved2[2 * j + 1];
				pr2_curved_2d[j][k][2] = pr_curved2[2 * j + 1];
				pr2_curved_2d[j][k][3] = pr_curved2[2 * j];
			}
		}

		// Read data on emission time:

		// loop over cells on top and bottom surface of the cylinder

		for (boost::multi_array<double, 3>::index i = 0; i < Ncells_r; ++i)
		{
			for (boost::multi_array<double, 3>::index j = 0; j < Ncells_theta; ++j)
			{

				// loop over quad points on the cell:
				for (boost::multi_array<double, 3>::index q = 0; q < Gauss_2d_2pt::nqpts; ++q)
				{

					if (tauq_top[i][j][q] > t[n] && tauq_top[i][j][q] < t[n + 1])
					{

						pq_top[i][j][q] = LinearInterpolate(tauq_top[i][j][q], t[n], t[n + 1], p1_top_2d[i][j][q], p2_top_2d[i][j][q]);
						pnq_top[i][j][q] = LinearInterpolate(tauq_top[i][j][q], t[n], t[n + 1], pz1_top_2d[i][j][q], pz2_top_2d[i][j][q]);
						ptq_top[i][j][q] = ForwardDifference(t[n], t[n + 1], p1_top_2d[i][j][q], p2_top_2d[i][j][q]);
					}

					if (tauq_bottom[i][j][q] > t[n] && tauq_bottom[i][j][q] < t[n + 1])
					{

						pq_bottom[i][j][q] = LinearInterpolate(tauq_bottom[i][j][q], t[n], t[n + 1], p1_bottom_2d[i][j][q], p2_bottom_2d[i][j][q]);
						pnq_bottom[i][j][q] = LinearInterpolate(tauq_bottom[i][j][q], t[n], t[n + 1], pz1_bottom_2d[i][j][q], pz2_bottom_2d[i][j][q]);
						ptq_bottom[i][j][q] = ForwardDifference(t[n], t[n + 1], p1_bottom_2d[i][j][q], p2_bottom_2d[i][j][q]);
					}
				}
			}
		}

		// loop over cells on curved surface of the cylinder
		for (boost::multi_array<double, 3>::index i = 0; i < Ncells_theta; ++i)
		{
			for (boost::multi_array<double, 3>::index j = 0; j < Ncells_z; ++j)
			{
				// loop over quad points on the cell:
				for (boost::multi_array<double, 3>::index q = 0; q < Gauss_2d_2pt::nqpts; ++q)
				{
					if ((tauq_curved[i][j][q] > t[n]) && (tauq_curved[i][j][q] < t[n + 1]))
					{
						pq_curved[i][j][q] = LinearInterpolate(tauq_curved[i][j][q], t[n], t[n + 1], p1_curved_2d[i][j][q], p2_curved_2d[i][j][q]);
						pnq_curved[i][j][q] = LinearInterpolate(tauq_curved[i][j][q], t[n], t[n + 1], pr1_curved_2d[i][j][q], pr2_curved_2d[i][j][q]);
						ptq_curved[i][j][q] = ForwardDifference(t[n], t[n + 1], p1_curved_2d[i][j][q], p2_curved_2d[i][j][q]);
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

	// wave speed
	const double c = 66.53119569044284;

	// observer point and time
	double t0 = 0.005;
	double tf = 0.007;
	double dt = 0.0001;

	// observer point
	Point x0{0.36556, 0.0, 0.000760000000000038};

	// no of cells along different direction
	const int Ncells_r = 231;
	const int Ncells_theta = 200;
	const int Ncells_z = 481;

	// cylinder geometry
	const double r_min = 0.0;
	const double r_max = 0.35036;
	const double theta_min = 0.0;
	const double theta_max = 2.0 * M_PI;
	const double z_min = -0.36404000000000003;
	const double z_max = 0.36556;

	// cell size
	const double dh = 0.00152;

	Kirchhoff solver(c, r_min, r_max, theta_min, theta_max, z_min, z_max, Ncells_r, Ncells_theta, Ncells_z, dh);

	solver.create_grid();
	solver.write_grid();

	// write Kirchhoff data:
	std::ofstream observer("p_kirchhoff.dat", std::ios::out);
	observer.flags(std::ios::dec | std::ios::scientific);
	observer.precision(16);
	if (!observer)
	{
		std::cerr << "Error: Output file couldnot be opened.\n";
	}

	// loop over observer time:
	while (t0 < tf)
	{
		std::cout << "observer time t0 = " << t0 << std::endl;

		solver.compute_emission_time(t0, x0);
		solver.compute_pressure_at_emission_time();
		observer << t0 << "\t" << solver.compute_kirchhoff_integral() << std::endl;

		t0 += dt;
	}

	return 0;
}
