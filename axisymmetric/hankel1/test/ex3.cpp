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
  boost::multi_array<double, 3> xq;   // coordinates of quadrature points
  boost::multi_array<double, 3> nq;   // normal vector at quadrature points
  boost::multi_array<double, 2> rq;   // distance r = |xo - xq| at quadpoints
  boost::multi_array<double, 2> cosq; // cos of angle between r and n at quadpoints
  boost::multi_array<double, 2> tauq; // emission time at quad points

  boost::multi_array<double, 2> pq;    // pressure at quadrature point
  boost::multi_array<double, 2> pnq;   // normal derivative of pressure at quadrature point
  boost::multi_array<double, 2> ptq;   // time derivative of pressure at quadrature point

};


// constructor
Kirchhoff::Kirchhoff() : Ncells_theta(Nnodes_theta - 1),
                         Ncells_z(Nnodes_z - 1),
                         dtheta((theta_max - theta_min) / static_cast<double>(Ncells_theta)),
                         dz((z_max - z_min) / static_cast<double>(Ncells_z)),
                         xq(boost::extents[Ncells_z][Ncells_theta][3]),
                         nq(boost::extents[Ncells_z][Ncells_theta][3]),
                         rq(boost::extents[Ncells_z][Ncells_theta]),
                         cosq(boost::extents[Ncells_z][Ncells_theta]),
                         tauq(boost::extents[Ncells_z][Ncells_theta]),
                         pq(boost::extents[Ncells_z][Ncells_theta]),
                         pnq(boost::extents[Ncells_z][Ncells_theta]),
                         ptq(boost::extents[Ncells_z][Ncells_theta])
                         
{
}

// make grid:
void Kirchhoff::make_grid()
{

  // get cartesian coordinates of quadrature point and normal vector:

  // loop over cells:
  for (boost::multi_array<double, 3>::index i = 0; i < Ncells_z; ++i){
    for (boost::multi_array<double, 3>::index j = 0; j < Ncells_theta; ++j){

        // get theta and z at quadrature point:
        double theta_q = theta_min + (j+0.5)*dtheta;
        double z_q = z_min + (i + 0.5)*dz;

        // convert cylindrical to cartesian:
        Point point_q = cylinder_to_cartesian(r, theta_q, z_q);

        // get coordinates of quadrature point:
        xq[i][j][0] = point_q.x;
        xq[i][j][1] = point_q.y;
        xq[i][j][2] = point_q.z;

        // get normal vector at quadrature point:
        nq[i][j][0] = cos(theta_q);
        nq[i][j][1] = sin(theta_q);
        nq[i][j][2] = 0.0;
      
    }
  } //end of cell loop

}


//write grid:
void Kirchhoff::write_grid()
{

  std::ofstream file("cylinder.dat", std::ios::out) ;
  file.flags( std::ios::dec | std::ios::scientific );
  file.precision(16) ;
  if(!file) {std::cerr<< "Error: Output file couldnot be opened.\n";}


  // loop over cells:
  for (boost::multi_array<double, 3>::index i = 0; i < Ncells_z; ++i)
    for (boost::multi_array<double, 3>::index j = 0; j < Ncells_theta; ++j)
      file << xq[i][j][0] << "\t" << xq[i][j][1] << "\t" << xq[i][j][2] << "\n";

  file.close();
}


// Get emission time τ = to − r/c, distance r = |xo−y| and cos of angle between r = xo−y and n vectors at each quadpoint:
void Kirchhoff::compute_emission_time(double to)
{

  // loop over cells:
  for (boost::multi_array<double, 2>::index i = 0; i < Ncells_z; ++i){
    for (boost::multi_array<double, 2>::index j = 0; j < Ncells_theta; ++j){

        // get cartesian coordinates of quadrature point:
        Point y{xq[i][j][0], xq[i][j][1], xq[i][j][2]};

        // get normal vector at quadrature point:
        Vector<double> n{nq[i][j][0], nq[i][j][1], nq[i][j][2]};

        // compute r = |xo−y|
        rq[i][j] = norm(xo - y);

        // compute τ = to − r/c
        tauq[i][j] = emission_time(xo, y, to, c);

        // compute cos of angle between r and n vector:
        cosq[i][j] = cos(xo - y, n);
      
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
    for (boost::multi_array<double, 2>::index i = 0; i < Ncells_z; ++i){
      for (boost::multi_array<double, 2>::index j = 0; j < Ncells_theta; ++j){

        //check if the emission time at a quad point is in the range (t[n], t[n+1])
         if ( (tauq[i][j] > t[n]) && (tauq[i][j] < t[n+1]) ){

          pq[i][j] = linear_interpolate(tauq[i][j], t[n], t[n+1], p[n], p[n+1]);
          pnq[i][j] = linear_interpolate(tauq[i][j], t[n], t[n+1], dp_dr[n], dp_dr[n+1]);
          ptq[i][j] = forward_difference(t[n], t[n+1], p[n], p[n+1]);

         }

      }
    }

  }

}


//Compute Kirchhoff Integral:
double Kirchhoff::compute_kirchhoff_integral(){

  double integral = 0.0;
  double cinv = 1.0/c;

  // loop over cells:
  for (boost::multi_array<double, 2>::index i = 0; i < Ncells_z; ++i)
    for (boost::multi_array<double, 2>::index j = 0; j < Ncells_theta; ++j)
      integral += (0.25/M_PI)*( (cinv*ptq[i][j]*cosq[i][j] - pnq[i][j])/rq[i][j] + (pq[i][j]*cosq[i][j])/(rq[i][j]*rq[i][j])  )*r*dtheta*dz;


  return integral;

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
