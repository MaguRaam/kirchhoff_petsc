#include "Kirchhoff.hpp"


int main()
{

	
  // observer time:
  double t0 = 0.0;
  double tf = 20.0;
  double dt = 0.1;

  //observer point
  Point x0{5.0, 0.0, 0.0};

	//exact solution
	CylindricalWave wave(1, 1, {1.0, 0.0});

	//Kirchhoff solver
  Kirchhoff solver;
  
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
    observer << t0 << "\t" << solver.compute_kirchhoff_integral() << "\t"<< wave.p(radius(x0.x, x0.y, x0.z), t0) << std::endl;

    t0 += dt;
  }

  return 0;
}
