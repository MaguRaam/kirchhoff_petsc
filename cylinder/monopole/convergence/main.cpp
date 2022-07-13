#include "Kirchhoff.hpp"


int main()
{
  //wave speed:
  double c = 250.0;
	
  // observer time:
  double t0 = 0.05;

  //observer point
  Point x0{3.0, 3.0, 3.0};

  // cylindrical domain:
  double r_min = 0.0;
  double r_max = 0.1;
  double theta_min = 0.0;
  double theta_max = 2.0 * M_PI;
  double z_min = -0.1;
  double z_max =  0.1;

  double dh = 0.01;


  //exact solution
  SphericalWave wave(c, 100.);

	
  // write Kirchhoff data:
  std::ofstream file("error.dat", std::ios::app);
  file.flags(std::ios::dec | std::ios::scientific);
  file.precision(16);
  if (!file)
  {
    std::cerr << "Error: Output file couldnot be opened.\n";
  }
  
  for (int i = 0; i < 5; ++i){
  	
  	//Kirchhoff solver
		Kirchhoff solver(c, r_min, r_max, theta_min, theta_max, z_min, z_max, dh);

		solver.create_grid();
		solver.compute_emission_time(t0, x0);
		solver.compute_pressure_at_emission_time(wave);

		double error = abs(solver.compute_kirchhoff_integral() - wave.p(norm(x0), t0));

		file << dh << "\t" << error << std::endl;
		
		dh = 0.5*dh;  
  }

  return 0;
}
