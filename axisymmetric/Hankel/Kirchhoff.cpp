// load pressure data on the cylindrical Kirchhoff surface and compute the far-field pressure

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
T norm(const Vector<T> &p) { return sqrt(dot(p, p)); }

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

// Bi-linear basis function defined on [-1, 1]x[-1, 1]
namespace BilinearBasis
{

  // nodal basis functions:
  inline double N0(double psi, double eta) { return 0.25 * (1 - psi) * (1 - eta); }
  inline double N1(double psi, double eta) { return 0.25 * (1 + psi) * (1 - eta); }
  inline double N2(double psi, double eta) { return 0.25 * (1 + psi) * (1 + eta); }
  inline double N3(double psi, double eta) { return 0.25 * (1 - psi) * (1 + eta); }

  // interpolation on reference cell [-1, 1]x[-1, 1], given nodal values
  inline double interpolate(double psi, double eta, const std::array<double, 4> &p)
  {
    return p[0] * N0(psi, eta) + p[1] * N1(psi, eta) + p[2] * N2(psi, eta) + p[3] * N3(psi, eta);
  }

}; // namespace BilinearBasis

// 2 point gauss quadrature on a reference cell [-1, 1]x[-1, 1]
namespace Gauss_2d_2pt
{
  static constexpr int nqpts = 4;
  static constexpr double psi[nqpts] = {-0.577350269189626, -0.577350269189626, 0.577350269189626, 0.577350269189626};
  static constexpr double eta[nqpts] = {-0.577350269189626, 0.577350269189626, -0.577350269189626, 0.577350269189626};
  static constexpr double w[nqpts] = {1.0, 1.0, 1.0, 1.0};

}; // namespace Gauss_2d_2pt


//linear interpolation in time:
inline double linear_interpolate(double tau, double t1, double t2, double p1, double p2){
 return p1 + ((p2 - p1)/(t2 - t1))*(tau - t1);
}

//forward difference in time:
inline double forward_difference(double t1, double t2, double p1, double p2){
  return (p2 - p1)/(t2 - t1);
}



class Kirchhoff
{

public:
  Kirchhoff();
  void make_grid();
  void write_grid();
  void compute_emission_time(double);
  void compute_pressure_at_emission_time();
  double compute_kirchhoff_integral();


private:
  // observer point and time:
  Point xo{2.0, 0.0, 0.0};

  // wave speed:
  double c = 1.0;

  // cylindrical domain:
  double r = 1.0;
  double z_min = -M_PI;
  double z_max =  M_PI;
  double theta_min = 0.0;
  double theta_max = 2.0 * M_PI;

  // no of nodes:
  int Nnodes_theta = 100;
  int Nnodes_z = 100;

  // no of cells:
  int Ncells_theta;
  int Ncells_z;

  // cell size:
  double dtheta;
  double dz;

private:
  std::vector<double> z, theta;       // nodal coordinates
  boost::multi_array<double, 4> xq;   // coordinates of quadrature points
  boost::multi_array<double, 4> nq;   // normal vector at quadrature points
  boost::multi_array<double, 3> rq;   // distance r = |xo - xq| at quadpoints
  boost::multi_array<double, 3> cosq; // cos of angle between r and n at quadpoints
  boost::multi_array<double, 3> tauq; // emission time at quad points

  boost::multi_array<double, 3> pq;    // pressure at quadrature point
  boost::multi_array<double, 3> pnq;   // normal derivative of pressure at quadrature point
  boost::multi_array<double, 3> ptq;   // time derivative of pressure at quadrature point

};

// constructor
Kirchhoff::Kirchhoff() : Ncells_theta(Nnodes_theta - 1),
                         Ncells_z(Nnodes_z - 1),
                         dtheta((theta_max - theta_min) / static_cast<double>(Ncells_theta)),
                         dz((z_max - z_min) / static_cast<double>(Ncells_z)),
                         z(Nnodes_z),
                         theta(Nnodes_theta),
                         xq(boost::extents[Ncells_z][Ncells_theta][Gauss_2d_2pt::nqpts][3]),
                         nq(boost::extents[Ncells_z][Ncells_theta][Gauss_2d_2pt::nqpts][3]),
                         rq(boost::extents[Ncells_z][Ncells_theta][Gauss_2d_2pt::nqpts]),
                         cosq(boost::extents[Ncells_z][Ncells_theta][Gauss_2d_2pt::nqpts]),
                         tauq(boost::extents[Ncells_z][Ncells_theta][Gauss_2d_2pt::nqpts]),
                         pq(boost::extents[Ncells_z][Ncells_theta][Gauss_2d_2pt::nqpts]),
                         pnq(boost::extents[Ncells_z][Ncells_theta][Gauss_2d_2pt::nqpts]),
                         ptq(boost::extents[Ncells_z][Ncells_theta][Gauss_2d_2pt::nqpts])
                         
{
}

// make grid:
void Kirchhoff::make_grid()
{

  using namespace BilinearBasis;
  using namespace Gauss_2d_2pt;

  // z values at nodes:
  for (int i = 0; i < Nnodes_z; ++i)
    z[i] = z_min + i * dz;

  // theta values at nodes:
  for (int j = 0; j < Nnodes_theta; ++j)
    theta[j] = theta_min + j * dtheta;

  // get cartesian coordinates of quadrature points and normal vector:

  // loop over cells:
  for (boost::multi_array<double, 4>::index i = 0; i < Ncells_z; ++i){
    for (boost::multi_array<double, 4>::index j = 0; j < Ncells_theta; ++j){

      // get z and theta at nodes:
      std::array<double, 4> z_node{z[i], z[i], z[i + 1], z[i + 1]}, theta_node{theta[j], theta[j + 1], theta[j + 1], theta[j]};

      // loop over quadpoints and coordinates of quad points:
      for (boost::multi_array<double, 4>::index q = 0; q < nqpts; ++q)
      {

        // get theta and z at quadrature point:
        double theta_q = interpolate(psi[q], eta[q], theta_node);
        double z_q = interpolate(psi[q], eta[q], z_node);

        // convert cylindrical to cartesian:
        Point point_q = cylinder_to_cartesian(r, theta_q, z_q);

        // get coordinates of quadrature point:
        xq[i][j][q][0] = point_q.x;
        xq[i][j][q][1] = point_q.y;
        xq[i][j][q][2] = point_q.z;

        // get normal vector at quadrature point:
        nq[i][j][q][0] = cos(theta_q);
        nq[i][j][q][1] = sin(theta_q);
        nq[i][j][q][2] = 0.0;
      }
    }
  }
}


//write grid:
void Kirchhoff::write_grid()
{

  std::ofstream file("cylinder.dat", std::ios::out) ;
  file.flags( std::ios::dec | std::ios::scientific );
  file.precision(16) ;
  if(!file) {std::cerr<< "Error: Output file couldnot be opened.\n";}


  // loop over cells:
  for (boost::multi_array<double, 4>::index i = 0; i < Ncells_z; ++i)
    for (boost::multi_array<double, 4>::index j = 0; j < Ncells_theta; ++j)
     for (boost::multi_array<double, 4>::index q = 0; q < Gauss_2d_2pt::nqpts; ++q)
      file << xq[i][j][q][0] << "\t" << xq[i][j][q][1] << "\t" << xq[i][j][q][2] << "\n";

  file.close();
}




// Get emission time τ = to − r/c, distance r = |xo−y| and cos of angle between r = xo−y and n vectors at each quadpoint:
void Kirchhoff::compute_emission_time(double to)
{

  // loop over cells:
  for (boost::multi_array<double, 3>::index i = 0; i < Ncells_z; ++i){
    for (boost::multi_array<double, 3>::index j = 0; j < Ncells_theta; ++j){

      // loop over quad points:
      for (boost::multi_array<double, 3>::index q = 0; q < Gauss_2d_2pt::nqpts; ++q){

        // get cartesian coordinates of quadrature point:
        Point y{xq[i][j][q][0], xq[i][j][q][1], xq[i][j][q][2]};

        // get normal vector at quadrature point:
        Vector<double> n{nq[i][j][q][0], nq[i][j][q][1], nq[i][j][q][2]};

        // compute r = |xo−y|
        rq[i][j][q] = norm(xo - y);

        // compute τ = to − r/c
        tauq[i][j][q] = emission_time(xo, y, to, c);

        // compute cos of angle between r and n vector:
        cosq[i][j][q] = cos(xo - y, n);
      }
    }
  }
}


//Compute P, Pn, Pt at emission time τ:
void Kirchhoff::compute_pressure_at_emission_time(){

  //get emission begin and end time
  double emission_begin = *std::min_element(tauq.origin(), tauq.origin() + tauq.num_elements());
  double emission_end = *std::max_element(tauq.origin(), tauq.origin() + tauq.num_elements());

  std::cout << "emission begin = "<< emission_begin << "\t" << "emission end = "<< emission_end << std::endl;


  //load t, p, dp/dr, dp/dt vector from kirchhoff.dat file:
  std::vector<double> t, p, dp_dr, dp_dt;

  //open kirchhoff.dat:
  std::ifstream file("kirchhoff.dat", std::ios::in);
  if(!file) {std::cerr<< "Error: Input file couldnot be opened.\n";}

  //read vectors from file:
  for (double t_, p_, dp_dr_, dp_dt_; file >> t_ >> p_ >> dp_dr_ >> dp_dt_;){
   t.push_back(t_);
   p.push_back(p_);
   dp_dr.push_back(dp_dr_);
   dp_dt.push_back(dp_dt_);
  }

  //loop over time and interpolate p, pn, pt at emission time:
  for (int n = 0; n < t.size(); ++n){


    //loop over quad points and interpolate p, pt, pn at emission time:

    // loop over cells:
    for (boost::multi_array<double, 3>::index i = 0; i < Ncells_z; ++i){
      for (boost::multi_array<double, 3>::index j = 0; j < Ncells_theta; ++j){

        // loop over quad points:
        for (boost::multi_array<double, 3>::index q = 0; q < Gauss_2d_2pt::nqpts; ++q){


         //check if the emission time at a quad point is in the range (t[n], t[n+1])
         if ( (tauq[i][j][q] > t[n]) && (tauq[i][j][q] < t[n+1]) ){

          pq[i][j][q] = linear_interpolate(tauq[i][j][q], t[n], t[n+1], p[n], p[n+1]);
          pnq[i][j][q] = linear_interpolate(tauq[i][j][q], t[n], t[n+1], dp_dr[n], dp_dr[n+1]);
          ptq[i][j][q] = forward_difference(t[n], t[n+1], p[n], p[n+1]);

         }


        }
      }
     }

  } //end of time loop

}



//Compute Kirchhoff Integral:
double Kirchhoff::compute_kirchhoff_integral(){

  double cell_integral = 0.0;
  double global_integral = 0.0;

  double cinv = 1.0/c;

  // loop over cells:
  for (boost::multi_array<double, 3>::index i = 0; i < Ncells_z; ++i){
    for (boost::multi_array<double, 3>::index j = 0; j < Ncells_theta; ++j){

      //compute cell integral
      cell_integral = 0.0;

      // loop over quad points:
      for (boost::multi_array<double, 3>::index q = 0; q < Gauss_2d_2pt::nqpts; ++q)
        cell_integral += (0.25/M_PI)*( (cinv*ptq[i][j][q]*cosq[i][j][q] - pnq[i][j][q])/rq[i][j][q] + (pq[i][j][q]*cosq[i][j][q])/(rq[i][j][q]*rq[i][j][q])  )*Gauss_2d_2pt::w[q]*r*dtheta*dz*0.25;

      global_integral += cell_integral;

    }
  }

  return global_integral;

}




int main()
{

  Kirchhoff solver;

  solver.make_grid();
  solver.write_grid();

  //observer time:
  double t0 = 6.0;
  double tf = 10.0;
  double dt = 0.1;


  //write Kirchhoff data:
  std::ofstream observer("p_kirchhoff.dat", std::ios::out) ;
  observer.flags( std::ios::dec | std::ios::scientific );
  observer.precision(16) ;
  if(!observer) {std::cerr<< "Error: Output file couldnot be opened.\n";}

  //loop over observer time:
  while (t0 < tf){

    solver.compute_emission_time(t0);
    solver.compute_pressure_at_emission_time();
    observer << t0 << "\t" << solver.compute_kirchhoff_integral() << std::endl;

    t0+=dt;

  }


  return 0;
}

