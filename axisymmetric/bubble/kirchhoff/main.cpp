#include "Kirchhoff.hpp"


int main()
{
  //wave speed:
  double c = 66.53119569044284;
	
  // observer time:
  double t0 = 0.0;
  double tf = 0.1;
  double dt = 0.001;

  //observer point
  Point x0{3.0, 3.0, 3.0};

  // cylindrical domain:
  double r_min = 0.0;
  double r_max = 1.0;
  double theta_min = 0.0;
  double theta_max = 2.0 * M_PI;
  double z_min = -10.;
  double z_max =  10.;

  double dh = 0.01;


  //exact solution
  SphericalWave wave(c, 100.);

	//Kirchhoff solver
  Kirchhoff solver(c, r_min, r_max, theta_min, theta_max, z_min, z_max, dh);
  
  solver.create_grid();
  solver.write_grid();

  // write Kirchhoff data:
  std::ofstream observer("pressure.dat", std::ios::out);
  observer.flags(std::ios::dec | std::ios::scientific);
  observer.precision(16);
  if (!observer)
  {
    std::cerr << "Error: Output file couldnot be opened.\n";
  }

  // loop over observer time:
  while (t0 < tf)
  {
    std::cout << "observer time t0 = "<< t0 << std::endl;

    solver.compute_emission_time(t0, x0);
    solver.compute_pressure_at_emission_time(wave);
    observer << t0 << "\t" << solver.compute_kirchhoff_integral() << "\t"<< wave.p(norm(x0), t0) << std::endl;

    t0 += dt;
  }

  return 0;
}