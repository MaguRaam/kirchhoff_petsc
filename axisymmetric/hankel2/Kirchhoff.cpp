#include "Kirchhoff.hpp"

// constructor
Kirchhoff::Kirchhoff(double c, double r_min, double r_max, double theta_min, double theta_max, double z_min, double z_max, double dh) : 
						 c(c),
						 r_min(r_min),
						 r_max(r_max),
						 theta_min(theta_min),
						 theta_max(theta_max),
						 z_min(z_min),
						 z_max(z_max),
						 dh(dh),
						 Ncells_r((r_max - r_min)/dh),
						 Ncells_theta((theta_max - theta_min)/dh),
						 Ncells_z((z_max - z_min)/dh),
						 dr((r_max - r_min) / static_cast<double>(Ncells_r)),
						 dtheta((theta_max - theta_min) / static_cast<double>(Ncells_theta)),
						 dz((z_max - z_min) / static_cast<double>(Ncells_z)),
						 
						 r(Ncells_r),
						 theta(Ncells_theta),
						 z(Ncells_z),

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
	for (int i = 0; i < Ncells_r; ++i) r[i] = r_min + (i + 0.5)*dr;
	for (int j = 0; j < Ncells_theta; ++j) theta[j] = theta_min + (j + 0.5)*dtheta;
	for (int k = 0; k < Ncells_z; ++k) z[k] = z_min + (k + 0.5)*dz;

}

//create top surface of the cylinder:
void Kirchhoff::create_top_bottom_surface(){

	using namespace Gauss_2d_2pt;

	//get coordinates of quad points and normal vector on top and bottom surface of the cylinder:

	//loop over cells on top and bottom surface of the cylinder
	for (boost::multi_array<double, 4>::index i = 0; i < Ncells_r; ++i){
		for (boost::multi_array<double, 4>::index j = 0; j < Ncells_theta; ++j){

			
			//loop over quad points on the cell:
			for (boost::multi_array<double, 4>::index q = 0; q < nqpts; ++q){

				//get r and theta at quad points:
				double r_q = linear_map(psi[q], r[i], dr);
				double theta_q = linear_map(eta[q], theta[j], dtheta);


				//top surface:

				//get the cartesian coordinate of the quad point:
				Point point_q = cylinder_to_cartesian(r_q, theta_q, z_max);

				//get coordinates of quadrature point:
        		xq_top[i][j][q][0] = point_q.x;
        		xq_top[i][j][q][1] = point_q.y;
        		xq_top[i][j][q][2] = point_q.z;

				// get normal vector at quadrature point:
        		nq_top[i][j][q][0] = 0.0;
        		nq_top[i][j][q][1] = 0.0;
        		nq_top[i][j][q][2] = 1.0;


        		//bottom surface:

        		//get the cartesian coordinate of the quad point:
				point_q = cylinder_to_cartesian(r_q, theta_q, z_min);

				//get coordinates of quadrature point:
        		xq_bottom[i][j][q][0] = point_q.x;
        		xq_bottom[i][j][q][1] = point_q.y;
        		xq_bottom[i][j][q][2] = point_q.z;

				// get normal vector at quadrature point:
        		nq_bottom[i][j][q][0] = 0.0;
        		nq_bottom[i][j][q][1] = 0.0;
        		nq_bottom[i][j][q][2] = -1.0;

			} //q points loop


		} //r loop
	}	// theta loop

}



//create bottom surface of the cylinder:
void Kirchhoff::create_curved_surface(){

	using namespace Gauss_2d_2pt;

	//get coordinates of quad points and normal vector on curved surface of the cylinder:

	//loop over cells on curved surface of the cylinder:
	for (boost::multi_array<double, 4>::index j = 0; j < Ncells_theta; ++j){
		for (boost::multi_array<double, 4>::index k = 0; k < Ncells_z; ++k){

			
			//loop over quad points on the cell:
			for (boost::multi_array<double, 4>::index q = 0; q < nqpts; ++q){
				
				//get theta and z at quad points:
				double theta_q = linear_map(psi[q], theta[j], dtheta);
				double z_q = linear_map(eta[q], z[k], dz);

				//get the cartesian coordinate of the quad point:
				Point point_q = cylinder_to_cartesian(r_max, theta_q, z_q);

				//get coordinates of quadrature point:
        		xq_curved[j][k][q][0] = point_q.x;
        		xq_curved[j][k][q][1] = point_q.y;
        		xq_curved[j][k][q][2] = point_q.z;

				// get normal vector at quadrature point:
        		nq_curved[j][k][q][0] = cos(theta_q);
        		nq_curved[j][k][q][1] = sin(theta_q);
        		nq_curved[j][k][q][2] = 0.0;

			}


		} //theta loop
	}//z loop

}



//create cylindrical surface grid
void Kirchhoff::create_grid(){
	create_top_bottom_surface();
	create_curved_surface();
}



//write quad points on the cylinder to a file to be visualized using cylinder.py
void Kirchhoff::write_grid(){

	std::ofstream file("cylinder.dat", std::ios::out);
	file.flags(std::ios::dec | std::ios::scientific);
	file.precision(16);
	if (!file){std::cerr << "Error: cylinder.dat could not be opened.\n";}

	// write top and bottom surface quad points:
	for (boost::multi_array<double, 4>::index i = 0; i < Ncells_r; ++i){
		for (boost::multi_array<double, 4>::index j = 0; j < Ncells_theta; ++j){
			for (boost::multi_array<double, 4>::index q = 0; q < Gauss_2d_2pt::nqpts; ++q){

				file << xq_top[i][j][q][0] << "\t" << xq_top[i][j][q][1] << "\t" << xq_top[i][j][q][2] << "\n";
				file << xq_bottom[i][j][q][0] << "\t" << xq_bottom[i][j][q][1] << "\t" << xq_bottom[i][j][q][2] << "\n";
			}
		}
	}

	//write curved surface quad points:
	for (boost::multi_array<double, 4>::index j = 0; j < Ncells_theta; ++j)
		for (boost::multi_array<double, 4>::index k = 0; k < Ncells_z; ++k)
			for (boost::multi_array<double, 4>::index q = 0; q < Gauss_2d_2pt::nqpts; ++q)
				file << xq_curved[j][k][q][0] << "\t" << xq_curved[j][k][q][1] << "\t" << xq_curved[j][k][q][2] << "\n";

}


// Get emission time τ = to − r/c, distance r = |xo−y| and cos of angle between r = xo−y and n vectors at each quadpoint:
void Kirchhoff::compute_emission_time(double to, Point xo){

	using namespace Gauss_2d_2pt;


	//loop over cells on top and bottom surface of the cylinder
	for (boost::multi_array<double, 3>::index i = 0; i < Ncells_r; ++i){
		for (boost::multi_array<double, 3>::index j = 0; j < Ncells_theta; ++j){
			
			//loop over quad points on the cell:
			for (boost::multi_array<double, 3>::index q = 0; q < nqpts; ++q){

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


	//loop over cells on curved surface of the cylinder
	for (boost::multi_array<double, 3>::index j = 0; j < Ncells_theta; ++j){
		for (boost::multi_array<double, 3>::index k = 0; k < Ncells_z; ++k){

			//loop over quad points on the cell:
			for (boost::multi_array<double, 3>::index q = 0; q < Gauss_2d_2pt::nqpts; ++q){

        		Point y{xq_curved[j][k][q][0], xq_curved[j][k][q][1], xq_curved[j][k][q][2]};			
        		Vector<double> n{nq_curved[j][k][q][0], nq_curved[j][k][q][1], nq_curved[j][k][q][2]};

        		rq_curved[j][k][q] = norm(xo - y);
        		tauq_curved[j][k][q] = emission_time(xo, y, to, c);
        		cosq_curved[j][k][q] = cos(xo - y, n);	
        		
			}
		}
	}

}



//compute pressure at emission time:
void Kirchhoff::compute_pressure_at_emission_time(CylindricalWave wave){

	using namespace Gauss_2d_2pt;

	//loop over cells on top and bottom surface of the cylinder
	for (boost::multi_array<double, 3>::index i = 0; i < Ncells_r; ++i){
		for (boost::multi_array<double, 3>::index j = 0; j < Ncells_theta; ++j){
			
			//loop over quad points on the cell:
			for (boost::multi_array<double, 3>::index q = 0; q < nqpts; ++q){

				double r_q = linear_map(psi[q], r[i], dr);
		
				// top:
				pq_top[i][j][q] = wave.p(r_q, tauq_top[i][j][q]);
				pnq_top[i][j][q] = 0.0;
				ptq_top[i][j][q] = wave.dpdt(r_q, tauq_top[i][j][q]);

				// bottom:
				pq_bottom[i][j][q] = wave.p(r_q, tauq_bottom[i][j][q]);
				pnq_bottom[i][j][q] = 0.0;
				ptq_bottom[i][j][q] = wave.dpdt(r_q, tauq_bottom[i][j][q]);

			}
		}
	}


	//loop over cells on curved surface of the cylinder
	for (boost::multi_array<double, 3>::index j = 0; j < Ncells_theta; ++j){
		for (boost::multi_array<double, 3>::index k = 0; k < Ncells_z; ++k){

			//loop over quad points on the cell:
			for (boost::multi_array<double, 3>::index q = 0; q < nqpts; ++q){


				pq_curved[j][k][q] = wave.p(r_max, tauq_curved[j][k][q]);
				pnq_curved[j][k][q] = wave.dpdr(r_max, tauq_curved[j][k][q]);
				ptq_curved[j][k][q] = wave.dpdt(r_max, tauq_curved[j][k][q]);
			
			}
		}
	}


}


//compute Kirchhoff integral:
double Kirchhoff::compute_kirchhoff_integral(){

	using namespace Gauss_2d_2pt;

	double cell_integral = 0.0;
	double global_integral = 0.0;
	double cinv = 1.0 / c;



	//loop over cells on top and bottom surface of the cylinder
	for (boost::multi_array<double, 3>::index i = 0; i < Ncells_r; ++i){
		for (boost::multi_array<double, 3>::index j = 0; j < Ncells_theta; ++j){
			
			// compute cell integral
      		cell_integral = 0.0;

			//loop over quad points on the cell:
			for (boost::multi_array<double, 3>::index q = 0; q < nqpts; ++q){

				double radius_q = linear_map(psi[q], r[i], dr);

				cell_integral += (0.25 / M_PI) * ((cinv * ptq_top[i][j][q] * cosq_top[i][j][q] - pnq_top[i][j][q]) / rq_top[i][j][q] + (pq_top[i][j][q] * cosq_top[i][j][q]) / (rq_top[i][j][q] * rq_top[i][j][q])) * radius_q * dr * dtheta * w[q];
				
				cell_integral += (0.25 / M_PI) * ((cinv * ptq_bottom[i][j][q] * cosq_bottom[i][j][q] - pnq_bottom[i][j][q]) / rq_bottom[i][j][q] + (pq_bottom[i][j][q] * cosq_bottom[i][j][q]) / (rq_bottom[i][j][q] * rq_bottom[i][j][q])) * radius_q * dr * dtheta * w[q];
			}

			global_integral += cell_integral;
		}
	}


	//loop over cells on curved surface of the cylinder
	for (boost::multi_array<double, 3>::index j = 0; j < Ncells_theta; ++j){
		for (boost::multi_array<double, 3>::index k = 0; k < Ncells_z; ++k){

			// compute cell integral
      		cell_integral = 0.0;

			//loop over quad points on the cell:
			for (boost::multi_array<double, 3>::index q = 0; q < nqpts; ++q)
				cell_integral += (0.25 / M_PI) * ((cinv * ptq_curved[j][k][q] * cosq_curved[j][k][q] - pnq_curved[j][k][q]) / rq_curved[j][k][q] + (pq_curved[j][k][q] * cosq_curved[j][k][q]) / (rq_curved[j][k][q] * rq_curved[j][k][q])) * w[q] * r_max * dtheta * dz;


			global_integral += cell_integral;
		}
	}




	return global_integral;
}
