// boost headers:
#include <boost/math/special_functions.hpp>
#include <boost/multi_array.hpp>

// c++ headers:
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>
#include <cstdlib>


void write_tecplot(int step, double time, const boost::multi_array<double, 2> &p, const std::vector<double> &x, const std::vector<double> &y);
void write_vtk(int step, double time, const boost::multi_array<double, 2> &p, const std::vector<double> &x, const std::vector<double> &y);


double cylindrical_wave(double r, double t){

    double omega = 1.0;                     //frequency
    double k     = 1.0;                     //wave number
    std::complex<double> i(0.0,1.0);        //iota
    std::complex<double> A(1.0, 0.0);       //real amplitude

    return (A*boost::math::cyl_hankel_2(0, k*r)*exp(i*omega*t)).real();
}


int main()
{

    using array = boost::multi_array<double, 2>;
    using index = boost::multi_array<double, 2>::index;
    using boost::extents;


    // simulation parameters:
    double zmin = -10;
    double zmax = 10;
    double rmin = 0.0001;
    double rmax = 10;
    int nz = 500;
    int nr = 250;
    double dz = (zmax - zmin) / static_cast<double>(nz - 1);
    double dr = (rmax - rmin) / static_cast<double>(nr - 1);
    double t = 0.0;
    double tf = 10.0;
    double dt = 0.4;
    int step = 0;
    

    // create grid:
    std::vector<double> z(nz), r(nr);
    for (int i = 0; i < nz; ++i)
        z[i] = zmin + i * dz;

    for (int j = 0; j < nr; ++j)
        r[j] = rmin + j * dr;

    // allocate array:
    array p(extents[nz][nr]);


    // evolve in time:
    while (t < tf)
    {

        //compute analytical solution:
        for (int i = 0; i < nz; ++i)
            for (int j = 0; j < nr; ++j)
                p[i][j] = cylindrical_wave(r[j], t);


        // write data:
        if (step % 1 == 0)
            write_vtk(step, t, p, z, r);


        // update time and step
        t += dt;
        step++;
    }




    return 0;
}


void write_tecplot(int step, double time, const boost::multi_array<double, 2> &p, const std::vector<double> &x, const std::vector<double> &y)
{

    char filename[40];
    sprintf(filename, "plot/sol-%08d.dat", step);

    std::ofstream tpl;
    tpl.open(filename);
    tpl.flags(std::ios::dec | std::ios::scientific);
    tpl.precision(16);

    int nx = x.size();
    int ny = y.size();

    tpl << "TITLE = \"Wave solution\" " << std::endl
        << "VARIABLES = \"x\", \"y\", \"p\" " << std::endl;
    tpl << "Zone I = " << ny << " J = " << nx << std::endl;
    tpl << "SOLUTIONTIME = " << time << std::endl;

    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j)
            tpl << x[i] << "\t" << y[j] << "\t" << p[i][j] << std::endl;

    tpl.close();
}

void write_vtk(int step, double time, const boost::multi_array<double, 2> &p, const std::vector<double> &x, const std::vector<double> &y)
{

    // write vtk data:
    char filename[40];
    sprintf(filename, "plot/sol-%08d.vtk", step);

    std::ofstream vtk(filename, std::ios::out);
    vtk.flags(std::ios::dec | std::ios::scientific);
    vtk.precision(16);

    if (!vtk)
    {
        std::cerr << "Error: Output file couldnot be opened.\n";
        exit(1);
    }

    int nx = x.size();
    int ny = y.size();

    vtk << "# vtk DataFile Version 2.0"
        << "\n";
    vtk << "axisymmetric wave"
        << "\n";
    vtk << "ASCII"
        << "\n";
    vtk << "\nDATASET STRUCTURED_GRID"
        << "\n";
    vtk << "\nFIELD FieldData 1"
        << "\n";
    vtk << "TIME 1 1 double"
        << "\n";
    vtk << time << "\n";
    vtk << "\nDIMENSIONS " << ny << " " << nx << " " << 1 << "\n";
    vtk << "POINTS " << nx * ny << " double"
        << "\n";
    vtk << "\n";

    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            vtk << x[i] << " " << y[j] << " " << 0.0 << "\n";

    vtk << "\nPOINT_DATA " << nx * ny << "\n";

    vtk << "\nSCALARS p double 1"
        << "\n";
    vtk << "LOOKUP_TABLE default"
        << "\n";
    vtk << "\n";

    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            vtk << p[i][j] << "\n";

    vtk.close();
}


