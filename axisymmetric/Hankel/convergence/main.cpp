#include "Kirchhoff.hpp"


int main()
{

	
  // observer time:
  double t0 = 20.0;
  
  //observer point
  Point x0{5.0, 0.0, 0.0};

	//exact solution
	CylindricalWave wave(1, 1, {1.0, 0.0});

	//Kirchhoff solver
  Kirchhoff solver;
  
  solver.create_grid();
  solver.write_grid();

  // write Kirchhoff data:
  std::ofstream file("error.dat", std::ios::app);
  file.flags(std::ios::dec | std::ios::scientific);
  file.precision(16);
  if (!file)
  {
    std::cerr << "Error: Output file couldnot be opened.\n";
  }
  
  solver.compute_emission_time(t0, x0);
  solver.compute_pressure_at_emission_time(wave);
  
  double error = abs(solver.compute_kirchhoff_integral() - wave.p(radius(x0.x, x0.y, x0.z), t0));
  
  file << dh << "\t" << error << std::endl;



  return 0;
}
