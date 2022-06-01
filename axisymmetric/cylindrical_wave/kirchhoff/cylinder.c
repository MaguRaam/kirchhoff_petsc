//c++ headers
#include <stdio.h>
#include <stdlib.h>
#include <math.h>




//----------------------------------------------------------------------------
// Structure to define 3D point or vector
//----------------------------------------------------------------------------
typedef struct{ 
  double x, y, z;
} Point, Vector;


//----------------------------------------------------------------------------
// Distance vector R = X - Y
//----------------------------------------------------------------------------
static inline void DistanceVector(const Point* X, const Point* Y, Vector* R){

    R->x = X->x - Y->x;
    R->y = X->y - Y->y;
    R->z = X->z - Y->z;
}


//----------------------------------------------------------------------------
// Magnitude of point or vector
//----------------------------------------------------------------------------
static inline double Magnitude(const Point* R){
  return sqrt(R->x*R->x + R->y*R->y + R->z*R->z);
}


//----------------------------------------------------------------------------
// Distance between two points
//----------------------------------------------------------------------------
static inline double Distance(const Point* X, const Point* Y){

  Vector R;
  DistanceVector(X, Y, &R);
  return Magnitude(&R);
}



//----------------------------------------------------------------------------
// Normalize vector
//----------------------------------------------------------------------------
static inline void NormalizeVector(const Vector* R, Vector* r){

  double magR = Magnitude(R);

  r->x = (R->x)/magR;
  r->y = (R->y)/magR;
  r->z = (R->z)/magR;

}


//----------------------------------------------------------------------------
// Dot product <X, Y>
//----------------------------------------------------------------------------
static inline double Dot(const Vector* X, const Vector* Y){
  return X->x*Y->x + X->y*Y->y + X->z*Y->z;
}



// --------------------------------------------
// Emission time τ = t − r/c. where r = |x−y|
//---------------------------------------------
static inline double EmissionTime(double t, const Point *X, const Point *Y, double c){

  double r = Distance(X, Y);
  return t - r/c;
}



// --------------------------------------------
// Cos of angle between r = x−y and n vectors
//---------------------------------------------
static inline double CosOfNormalAndDistance(const Point *X, const Point *Y, const Point *n){

  Vector R, r;
  DistanceVector(X, Y, &R);
  NormalizeVector(&R, &r);
  return Dot(&r, n);
}





//----------------------------------------------------------------------------
// Bi-linear basis function defined on [-1, 1]x[-1, 1]
//----------------------------------------------------------------------------

double N0(double psi, double eta){
    return 0.25*(1 - psi)*(1 - eta);
}

double N1(double psi, double eta){
    return 0.25*(1 + psi)*(1 - eta);
}

double N2(double psi, double eta){
    return 0.25*(1 + psi)*(1 + eta);
}

double N3(double psi, double eta){
    return 0.25*(1 - psi)*(1 + eta);
}




//-------------------------------------------------------------------------------
// Bi-linear interpolation on reference cell [-1, 1]x[-1, 1], given nodal values
//-------------------------------------------------------------------------------

double BilinearInterpolate(double psi, double eta, double p0, double p1, double p2, double p3){

    return p0*N0(psi, eta) + p1*N1(psi, eta) + p2*N2(psi, eta) + p3*N3(psi, eta);

}


//-------------------------------------------------------------------------------
// 2 point gauss quadrature on a reference cell [-1, 1]x[-1, 1]
//-------------------------------------------------------------------------------

static const int N_gp2d = 4;
static const double psi_gp2d[] = {-0.577350269189626, -0.577350269189626,  0.577350269189626, 0.577350269189626};
static const double eta_gp2d[] = {-0.577350269189626,  0.577350269189626, -0.577350269189626, 0.577350269189626};
static const double w_gp2d[]   = {1.0,  1.0, 1.0, 1.0};



//map r,theta,z -> x,y,z
static inline Point CylindricalToCartesian(double r, double theta, double z){

    Point p;
    p.x = r*cos(theta);
    p.y = r*sin(theta);
    p.z = z;

    return p;
}










int main(){

    

    //cylindrical domain:
    double theta_min 	= 	 0.0;
    double theta_max 	= 	 2.0*M_PI;
    double z_min 		= 	-1.0;     
    double z_max		=	 1.0;
    double r            =    1.0;


    //grid:
    int Nnodes_theta	=	 100;
    int Nnodes_z		=	 100;
    int Ncells_theta	=	 Nnodes_theta - 1;
    int Ncells_z		=	 Nnodes_z - 1;
    
    //cell size:
    double dtheta 		=    (theta_max - theta_min)/(double)Ncells_theta;
    double dz  			=	 (z_max	- z_min)/(double)Ncells_z;



    //create 2d array storing pressure at nodes:
    double p[Nnodes_z][Nnodes_theta];


    //initialize pressure data:
    for (int i = 0; i < Nnodes_z; ++i){
        for (int j = 0; j < Nnodes_theta; ++j){

            p[i][j] = 1.0;

        }
    }



    //loop over cells to compute the integral:
    double integral = 0.0, cell_integral = 0.0;


     for (int i = 0; i < Ncells_z; ++i){
        for (int j = 0; j < Ncells_theta; ++j){

            //get pressure values at cell vertex:
            double p0 = p[i][j];
            double p1 = p[i][j+1];
            double p2 = p[i+1][j+1];
            double p3 = p[i+1][j];


            //compute integral on each cell by looping over quadrature points:
            cell_integral = 0.0;

            for (int q = 0; q < N_gp2d; q++){

                cell_integral += BilinearInterpolate(psi_gp2d[q], eta_gp2d[q],p0,p1,p2,p3)*w_gp2d[q]*(dtheta*dz)*0.25;

            }

            //accumulate integrals:
            integral += cell_integral;

        }
    }

    //scale it by R to compute integral on a cylinder:
    integral = r*integral;


    printf("integral = %.12g\n", integral);


    return 0;
}
