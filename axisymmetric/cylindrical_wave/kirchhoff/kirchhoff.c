/*
 * kirchhoff.c
 *      Author: Magu
 */ 
 
static char help[] = "Kirchhoff Solver\n";

//----------------------------------------------------------------------------
// Petsc headers files 
//----------------------------------------------------------------------------

#include <petscvec.h>
#include <petscmath.h> 
#include <petscdmda.h>




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


//----------------------------------------------------------------------------
// Structure defining various critical parameters controlling the simulation 
//----------------------------------------------------------------------------

typedef struct {

  //cylindrical surface:
  double theta_min;       
  double theta_max;
  double z_min;     
  double z_max;
  double r;                       //radius of the cylinder


  //grid:
  int Nnodes_theta;
  int Nnodes_z;
  int Ncells_theta;
  int Ncells_z;


  //cell size:
  double dtheta;                  //grid size along theta direction
  double dz;                      //grid size along z direction


  Point xo;                       // observer point in cartesian coordinate
  PetscReal c;                    // wave speed


  PetscReal to;                   // observer start time
  PetscReal tf;                   // observer end time
  PetscReal dt;                   // observer time step


} AppCtx;




//--------------------------------------------------------------------------------------------------------------
// Get emission time τ = t − r/c, distance r = |x−y| and cos of angle between r = x−y and n vectors at each node
//--------------------------------------------------------------------------------------------------------------
PetscErrorCode ComputeEmissionTime(DM da, Vec Tau, Vec R, Vec Cos, const AppCtx* ctx){

  PetscReal       **r_, **cos_, **tau_;
  DMDALocalInfo   info;
  PetscReal       z, theta;
  PetscInt        i, j;
  PetscErrorCode  ierr;
  Point           xs;
  Point           xo = ctx->xo;
  PetscReal       to = ctx->to;
  PetscReal       c  = ctx->c;
  Vector          normal;


  ierr = DMDAVecGetArray(da, R, &r_); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da, Cos, &cos_); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da, Tau, &tau_); CHKERRQ(ierr);

  ierr = DMDAGetLocalInfo(da, &info); CHKERRQ(ierr);


  for (j = info.ys; j < info.ys + info.ym; ++j){
    theta = j*ctx->dtheta;

    for (i = info.xs; i < info.xs + info.xm; ++i){
      z = i*ctx->dz;

      //get cartesian coordinates of the node:
      xs = CylindricalToCartesian(ctx->r, theta, z);

      //get unit normal at the node:
      normal.x = cos(theta);
      normal.y = sin(theta);
      normal.z = 0.0;


      //compute distance, emission time and cos at nodes:
      r_[j][i]   = Distance(&xo, &xs);
      tau_[j][i] = EmissionTime(to, &xo, &xs, c);
      cos_[j][i] = CosOfNormalAndDistance(&xo, &xs, &normal);

    }

  }

  ierr = DMDAVecRestoreArray(da, R, &r_); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da, Cos, &cos_); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da, Tau, &tau_); CHKERRQ(ierr);

  return 0;
}



// --------------------------------------------
// Compute Kirchhoff Integral:
//---------------------------------------------
PetscReal ComputeKirchhoffIntegral(DM da, Vec P, Vec Pn, Vec Pt, Vec Tau, Vec R, Vec Cos, const AppCtx* ctx){

  PetscErrorCode    ierr;
  const PetscReal   **p_, **pn_, **pt_, **tau_, **r_, **cos_;
  DMDALocalInfo     info;
  PetscInt          i, j;
  PetscReal         local_integral = 0.0;
  PetscReal         global_integral = 0.0;
  PetscReal         cinv = 1.0/(ctx->c);


  //get array from vector:
  ierr = DMDAVecGetArrayRead(da, P, &p_); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da, Pn, &pn_); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da, Pt, &pt_); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da, Tau, &tau_); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da, R, &r_); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da, Cos, &cos_); CHKERRQ(ierr);


  //loop over cells:
  
  for (j = info.ys; j < info.ys + info.ym; ++j){
    for (i = info.xs; i < info.xs + info.xm; ++i){

      //compute integral on each cell by looping over quadrature points:
      PetscReal cell_integral = 0.0;

      for (int q = 0; q < N_gp2d; q++){

        cell_integral += ( (cinv*ptq*cosq - pnq)/rq + (pq*cosq)/(rq*rq) )*w_gp2d[q]*(ctx->dtheta*ctx->dz)*0.25;


      }

      //accumulate integrals:
      local_integral += cell_integral;

    }
  }


  //sum kirchhoff integral from all the process:
  ierr = MPI_Allreduce(&local_integral, &global_integral,1,MPIU_REAL,MPIU_SUM, PETSC_COMM_WORLD);CHKERRQ(ierr);


  //Restore array from vector
  ierr = DMDAVecRestoreArrayRead(da, P, &p_); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da, Pn, &pn_); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da, Pt, &pt_); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da, Tau, &tau_); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da, R, &r_); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da, Cos, &cos_); CHKERRQ(ierr);


  return global_integral*ctx->r;
}







// --------------------------------------------
// Main program
//---------------------------------------------
int main(int argc, char **argv){


  // --------------------------------------------
  // Initialize MPI 
  //---------------------------------------------

  PetscErrorCode ierr;                                    
  ierr = PetscInitialize(&argc, &argv, NULL, help);CHKERRQ(ierr);

  // --------------------------------------------
  // Set important ctx defined parameters
  //---------------------------------------------
  
  AppCtx  ctx;


  ctx.theta_min 	    = 	 0.0;
  ctx.theta_max 	    = 	 2.0*M_PI;
  ctx.z_min 		      = 	 -1.0;     
  ctx.z_max		        =	   1.0;
  ctx.r               =    1.0;
  ctx.Nnodes_theta	  =	   100;
  ctx.Nnodes_z		    =	   100;
  ctx.Ncells_theta	  =	   ctx.Nnodes_theta - 1;
  ctx.Ncells_z		    =	   ctx.Nnodes_z - 1;
  ctx.dtheta 		      =    (ctx.theta_max - ctx.theta_min)/(double)ctx.Ncells_theta;
  ctx.dz  			      =	   (ctx.z_max	- ctx.z_min)/(double)ctx.Ncells_z;
  ctx.xo.x            =    2.0;
  ctx.xo.y            =    3.0;
  ctx.xo.z            =    3.0;
  ctx.c               =    1.0;
  ctx.to              =    0.0;                                
  ctx.tf              =    0.001;                              
  ctx.dt              =    0.001;                              




  // --------------------------------------------
  // Data members  
  //---------------------------------------------  
  
  DM                      da;                                 // Grid object
  Vec                     Tau, R, Cos;                        //emission time, distance and    
  Vec                     P, Pn, Pt;
  PetscMPIInt             MyPID;                              // Rank of the current processor 
  PetscMPIInt             numProcs;                           // Size of the communicator


  // --------------------------------------------
  // Obtain the rank of the process and size of 
  // the communicator 
  //---------------------------------------------

  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&numProcs);CHKERRQ(ierr); 
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&MyPID);CHKERRQ(ierr); 
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Code running with %d processes\n", numProcs);CHKERRQ(ierr); 
  

  // --------------------------------------------
  // Create cylindrical surface mesh
  //---------------------------------------------
  
  ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX, ctx.Nnodes_z, ctx.Nnodes_theta, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &da); CHKERRQ(ierr);
  ierr = DMSetFromOptions(da); CHKERRQ(ierr);
  ierr = DMSetUp(da);CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(da, ctx.z_min, ctx.z_max, ctx.theta_min, ctx.theta_max, 0.0, 0.0); CHKERRQ(ierr);
  
  
  // --------------------------------------------
  // Create vectors
  //---------------------------------------------

  ierr = DMCreateGlobalVector(da,&Tau);CHKERRQ(ierr);
  ierr = VecDuplicate(Tau, &Cos); CHKERRQ(ierr);
  ierr = VecDuplicate(Tau, &R); CHKERRQ(ierr);
  ierr = VecDuplicate(Tau, &P); CHKERRQ(ierr);
  ierr = VecDuplicate(Tau, &Pn); CHKERRQ(ierr);
  ierr = VecDuplicate(Tau, &Pt); CHKERRQ(ierr);



  // -----------------------------------------------------------------
  // Loop over observer time and compute farfield pressure p' = p - po
  //------------------------------------------------------------------

  while(ctx.to <= ctx.tf){

    ierr = PetscPrintf(PETSC_COMM_WORLD, "to = %f\n", ctx.to); CHKERRQ(ierr);


    ierr = ComputeEmissionTime(da, Tau, R, Cos, &ctx); CHKERRQ(ierr);




    ctx.to+=ctx.dt;

  }


  // --------------------------------------------
  // Free all the memory, finalize MPI and exit   
  //---------------------------------------------

  DMDestroy(&da);
  VecDestroy(&Tau);
  VecDestroy(&Cos);
  VecDestroy(&R);
  VecDestroy(&P);
  VecDestroy(&Pn);
  VecDestroy(&Pt);


  return PetscFinalize();
}
