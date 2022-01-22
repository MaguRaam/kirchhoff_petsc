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


//----------------------------------------------------------------------------
// Structure to define 3D point or vector
//----------------------------------------------------------------------------
typedef struct{ 
  PetscReal x, y, z;
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
static inline PetscReal Magnitude(const Point* R){
  return PetscSqrtReal(R->x*R->x + R->y*R->y + R->z*R->z);
}


//----------------------------------------------------------------------------
// Distance between two points
//----------------------------------------------------------------------------
static inline PetscReal Distance(const Point* X, const Point* Y){

  Vector R;
  DistanceVector(X, Y, &R);
  return Magnitude(&R);
}



//----------------------------------------------------------------------------
// Normalize vector
//----------------------------------------------------------------------------
static inline void NormalizeVector(const Vector* R, Vector* r){

  PetscReal magR = Magnitude(R);

  r->x = (R->x)/magR;
  r->y = (R->y)/magR;
  r->z = (R->z)/magR;

}


//----------------------------------------------------------------------------
// Dot product <X, Y>
//----------------------------------------------------------------------------
static inline PetscReal Dot(const Vector* X, const Vector* Y){
  return X->x*Y->x + X->y*Y->y + X->z*Y->z;
}



// --------------------------------------------
// Emission time τ = t − r/c. where r = |x−y|
//---------------------------------------------
static inline PetscReal EmissionTime(PetscReal t, const Point *X, const Point *Y, PetscReal c){

  PetscReal r = Distance(X, Y);
  return t - r/c;
}



// --------------------------------------------
// Cos of angle between r = x−y and n vectors
//---------------------------------------------
static inline PetscReal CosOfNormalAndDistance(const Point *X, const Point *Y, const Point *n){

  Vector R, r;
  DistanceVector(X, Y, &R);
  NormalizeVector(&R, &r);
  return Dot(&r, n);
}


// --------------------------------------------
// Linear interpolation
//---------------------------------------------
static inline PetscReal LinearInterpolate(PetscReal tau, PetscReal t1, PetscReal t2, PetscReal p1, PetscReal p2){
  return p1 + ((p2 - p1)/(t2 - t1))*(tau - t1);
}


// --------------------------------------------
// Forward difference
//---------------------------------------------
static inline PetscReal ForwardDifference(PetscReal t1, PetscReal t2, PetscReal p1, PetscReal p2){
  return (p2 - p1)/(t2 - t1);
}



// -----------------------------------------------------------------------------------------------------------------------
// Get emission time τ = t − r/c, distance r = |x−y| and cos of angle between r = x−y and n vectors at each quad point Yq
//------------------------------------------------------------------------------------------------------------------------
PetscErrorCode ComputeEmissionDistanceAndCos(Vec Tauq, Vec Rq, Vec Cosq, PetscReal to, Point xo, PetscReal c){

  PetscErrorCode ierr;

  Vec                     Yq, Nq;
  const PetscScalar       *yq, *nq;
  PetscScalar             *rq, *tauq, *cosq;
  PetscInt                nqpts, q;
  PetscViewer             viewer;

  //create petsc vector to store quadrature points and  normal vectors in the current process:
  ierr = VecCreate(PETSC_COMM_WORLD, &Yq); CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD, &Nq); CHKERRQ(ierr);

  //set blocksize = 3, as the vectors have 3 components:
  ierr = VecSetBlockSize(Yq, 3); CHKERRQ(ierr);
  ierr = VecSetBlockSize(Nq, 3); CHKERRQ(ierr);

  //read petsc vector from a binary file
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"input/q_points.dat",FILE_MODE_READ,&viewer); CHKERRQ(ierr);

  //load quadrature points and normal vectors
  ierr = VecLoad(Yq, viewer); CHKERRQ(ierr);
  ierr = VecLoad(Nq, viewer); CHKERRQ(ierr);

  //Get no of quadpts in the current process:
  ierr = VecGetLocalSize(Tauq, &nqpts); CHKERRQ(ierr);

  //get array of quad points, normal vectors, emission time and cos at quadpoints:
  ierr = VecGetArrayRead(Yq, &yq); CHKERRQ(ierr);
  ierr = VecGetArrayRead(Nq, &nq); CHKERRQ(ierr);
  ierr = VecGetArray(Tauq, &tauq); CHKERRQ(ierr);
  ierr = VecGetArray(Cosq, &cosq); CHKERRQ(ierr);
  ierr = VecGetArray(Rq, &rq); CHKERRQ(ierr);

  //cast double pointer to Point pointer:
  Point* y = (Point*)yq;
  Point* n = (Point*)nq;

  //loop over quad points in the current process:
  for (q = 0; q < nqpts; ++q){

    //get distance r = x−y:
    rq[q] = Distance(&xo, &y[q]); CHKERRQ(ierr);

    //get emission time:
    tauq[q] = EmissionTime(to, &xo, &y[q], c); CHKERRQ(ierr);

    //get  cos of angle between r = x−y and n vectors
    cosq[q] = CosOfNormalAndDistance(&xo, &y[q], &n[q]); CHKERRQ(ierr);

  }


  //restore array
  ierr = VecRestoreArrayRead(Yq, &yq); CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Nq, &nq); CHKERRQ(ierr);
  ierr = VecRestoreArray(Tauq, &tauq); CHKERRQ(ierr);
  ierr = VecRestoreArray(Cosq, &cosq); CHKERRQ(ierr);
  ierr = VecRestoreArray(Rq, &rq); CHKERRQ(ierr);

  //Destroy objects
  VecDestroy(&Yq);
  VecDestroy(&Nq);
  PetscViewerDestroy(&viewer);

  return ierr;
}


// --------------------------------------------
// Compute P, Pn, Pt at emission time τ:
//---------------------------------------------
PetscErrorCode  ComputePressureAtEmissionTime(Vec P, Vec Pn, Vec Pt, Vec Tauq){

    PetscErrorCode ierr;
    PetscInt       it = 0, nqpts, q;
    PetscReal      t, dt, c;
    PetscReal      emission_begin, emission_end;
    PetscViewer    viewer1, viewer2;    

    Vec P1, Pn1;
    Vec P2, Pn2;

    PetscReal *p, *pn, *pt;
    const PetscReal *p1, *pn1;
    const PetscReal *p2, *pn2;
    const PetscReal *tauq;

    FILE *f1, *f2;
  
    //get emission begin and end time:
    ierr = VecMin(Tauq, NULL, &emission_begin); CHKERRQ(ierr);
    ierr = VecMax(Tauq, NULL, &emission_end); CHKERRQ(ierr);

    //print emission begin and emission end time:
    ierr = PetscPrintf(PETSC_COMM_WORLD, "\temission begin = %e   emission end = %e\n", emission_begin, emission_end); CHKERRQ(ierr);

    //Get no of quadpts in the current process:
    ierr = VecGetLocalSize(P, &nqpts); CHKERRQ(ierr);

    //create vectors:
    ierr = VecCreate(PETSC_COMM_WORLD, &P1); CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD, &Pn1); CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD, &P2); CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD, &Pn2); CHKERRQ(ierr);
    

    //get array from vector P, Pn, Pt and Tauq:
    ierr = VecGetArray(P, &p); CHKERRQ(ierr);
    ierr = VecGetArray(Pn, &pn); CHKERRQ(ierr);
    ierr = VecGetArray(Pt, &pt); CHKERRQ(ierr);
    ierr = VecGetArrayRead(Tauq, &tauq); CHKERRQ(ierr);

    //open t.dat  and dt.dat file in all process:
    f1 = fopen("input/t.dat", "r");
    f2 = fopen("input/dt.dat", "r");

    //throw file opening error:
    if (f1 == NULL){
        printf("could not open t.dat");
        exit(1);
    }

    if (f2 == NULL){
        printf("could not open dt.dat");
        exit(1);
    }


    //loop over time steps and interpolate p, pn, pt at emission time
    while ((c = fgetc(f1)) != EOF)
    {
        //get t and t + dt:
        fscanf(f1, "%lf", &t);
        fscanf(f2, "%lf", &dt);

        if (t > emission_begin && t < emission_end){

          
          //load P1, Pn1, Pt1 at current step and P2, Pn2, Pt2 next step: 
          char filename1[20], filename2[20];
          sprintf(filename1, "input/p-%05d.dat", it);
          sprintf(filename2, "input/p-%05d.dat", it+1);

          ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename1,FILE_MODE_READ,&viewer1); CHKERRQ(ierr);
          ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename2,FILE_MODE_READ,&viewer2); CHKERRQ(ierr);

          //PetscPrintf(PETSC_COMM_WORLD, "Reading p, pn, pt on Kirchoff surface from %s and %s\n", filename1, filename2);

          ierr = VecLoad(P1, viewer1); CHKERRQ(ierr);
          ierr = VecLoad(Pn1, viewer1); CHKERRQ(ierr);
          ierr = VecLoad(P2, viewer2); CHKERRQ(ierr);
          ierr = VecLoad(Pn2, viewer2); CHKERRQ(ierr);

          //get array from vector :
          ierr = VecGetArrayRead(P1, &p1); CHKERRQ(ierr);
          ierr = VecGetArrayRead(Pn1, &pn1); CHKERRQ(ierr);
          ierr = VecGetArrayRead(P2, &p2); CHKERRQ(ierr);
          ierr = VecGetArrayRead(Pn2, &pn2); CHKERRQ(ierr);


          //loop over quad points and interpolate p, pt, pn at emission time :
          for (q = 0; q < nqpts; ++q){

            //Check if the emission time is in the range (t, t + dt) !TODO is this right?
            if (tauq[q] > t && tauq[q] < t + dt){
              
              p[q] =   LinearInterpolate(tauq[q], t, t + dt, p1[q], p2[q]);
              pn[q] =  LinearInterpolate(tauq[q], t, t + dt, pn1[q], pn2[q]);
              pt[q] =  ForwardDifference(t, t + dt, p1[q], p2[q]);

            }

          }

          //restore array from vector:
          ierr = VecRestoreArrayRead(P1, &p1); CHKERRQ(ierr);
          ierr = VecRestoreArrayRead(Pn1, &pn1); CHKERRQ(ierr);
          ierr = VecRestoreArrayRead(P2, &p2); CHKERRQ(ierr);
          ierr = VecRestoreArrayRead(Pn2, &pn2); CHKERRQ(ierr);


          //Destroy veiwer:
          PetscViewerDestroy(&viewer1);
          PetscViewerDestroy(&viewer2);

        }
        
        //increment time step:
        it++;
    }

    //close file
    fclose(f1);
    fclose(f2);


    //restore array p, pn, pt, tauq:
    ierr = VecRestoreArray(P, &p); CHKERRQ(ierr);
    ierr = VecRestoreArray(Pn, &pn); CHKERRQ(ierr);
    ierr = VecRestoreArray(Pt, &pt); CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(Tauq, &tauq); CHKERRQ(ierr);


    //Destroy objects:
    VecDestroy(&P1);
    VecDestroy(&Pn1);
    VecDestroy(&P2);
    VecDestroy(&Pn2);
    
    return ierr;
}


// --------------------------------------------
// Compute Kirchhoff Integral:
//---------------------------------------------
PetscReal ComputeKirchhoffIntegral(Vec P, Vec Pn, Vec Pt, Vec Tauq, Vec Rq, Vec Cosq, PetscReal h, PetscReal to, Point xo, PetscReal c){

  PetscErrorCode    ierr;
  PetscInt          nqpts, q;
  const PetscReal   *p, *pn, *pt, *tauq, *rq, *cosq;
  PetscReal         JxW;                               
  PetscScalar       local_integral = 0.0, global_integral = 0.0;
  PetscScalar       cinv = 1.0/c;


  //Get array from vector
  ierr = VecGetArrayRead(P, &p); CHKERRQ(ierr);
  ierr = VecGetArrayRead(Pn, &pn); CHKERRQ(ierr);
  ierr = VecGetArrayRead(Pt, &pt); CHKERRQ(ierr);
  ierr = VecGetArrayRead(Tauq, &tauq); CHKERRQ(ierr);
  ierr = VecGetArrayRead(Rq, &rq); CHKERRQ(ierr);
  ierr = VecGetArrayRead(Cosq, &cosq); CHKERRQ(ierr);


  //Get no of quadpts in the current process:
  ierr = VecGetLocalSize(P, &nqpts); CHKERRQ(ierr);


  //Compute Jacobian weight:
  JxW = h*h*0.25;


  //loop over quad points in the current process and compute local kirchhoff integral:
  for (q = 0; q < nqpts; ++q)
    local_integral += (cinv*pt[q]*cosq[q] - pn[q])/rq[q] + (p[q]*cosq[q])/(rq[q]*rq[q]); 


  //sum kirchhoff integral from all the process:
  ierr = MPI_Allreduce(&local_integral, &global_integral,1,MPIU_REAL,MPIU_SUM, PETSC_COMM_WORLD);CHKERRQ(ierr);


  //Restore array from vector
  ierr = VecRestoreArrayRead(P, &p); CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Pn, &pn); CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Pt, &pt); CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Tauq, &tauq); CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Rq, &rq); CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Cosq, &cosq); CHKERRQ(ierr);

  return JxW*global_integral*(0.25/PETSC_PI);
}


// --------------------------------------------
// Main program
//---------------------------------------------
int main(int argc, char **argv)
{
    // --------------------------------------------
    // Initialize MPI 
    //---------------------------------------------
    PetscErrorCode ierr;
    ierr = PetscInitialize(&argc, &argv, NULL, help);CHKERRQ(ierr);


    // --------------------------------------------
    // Set important user defined parameters
    //---------------------------------------------
    //!The parameters should match the fluid solver:
    Point     xo              =   {0.34808, 0.34808, 0.34808};          //observer location
    PetscReal c               =   66.53119569044284;                    //acoustic wave speed
    PetscReal h               =   0.00304;                              //cell size
    PetscInt  Tnqpts          =   351384;                               //Total no of quad points

    
    PetscReal to              =   0.0;                                  // observer start time
    PetscReal tf              =   0.1;                                  // observer end time
    PetscReal dt              =   0.001;                                // observer time step
    

    // --------------------------------------------
    // Data members  
    //---------------------------------------------
    Vec                     Tauq, Rq, Cosq;                    
    Vec                     P, Pn, Pt;     
    PetscMPIInt             MyPID;                            // Rank of the current processor 
    PetscMPIInt             numProcs;                         // Size of the communicator
    

    // --------------------------------------------
    // Obtain the rank of the process and size of 
    // the communicator 
    //---------------------------------------------

    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&numProcs);CHKERRQ(ierr); 
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&MyPID);CHKERRQ(ierr); 
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Code running with %d processes\n", numProcs);CHKERRQ(ierr); 
    


    // --------------------------------------------
    // Create vectors
    //---------------------------------------------
    ierr = VecCreate(PETSC_COMM_WORLD, &Tauq); CHKERRQ(ierr);
    ierr = VecSetSizes(Tauq, PETSC_DECIDE, Tnqpts); CHKERRQ(ierr);
    ierr = VecSetUp(Tauq); CHKERRQ(ierr);
    ierr = VecDuplicate(Tauq, &Cosq); CHKERRQ(ierr);
    ierr = VecDuplicate(Tauq, &Rq); CHKERRQ(ierr);
    ierr = VecDuplicate(Tauq, &P); CHKERRQ(ierr);
    ierr = VecDuplicate(Tauq, &Pn); CHKERRQ(ierr);
    ierr = VecDuplicate(Tauq, &Pt); CHKERRQ(ierr);


    // -----------------------------------------------------------------
    // Loop over observer time and compute farfield pressure p' = p - po
    //------------------------------------------------------------------

    while(to <= tf){

	      ierr = PetscPrintf(PETSC_COMM_WORLD, "to = %f\n", to); CHKERRQ(ierr);
		 
        // -----------------------------------------------------------------------------------------------------------------------
        // Get emission time τ = t − r/c, distance r = |x−y| and cos of angle between r = x−y and n vectors at each quad point Yq
        //------------------------------------------------------------------------------------------------------------------------    
        ierr = PetscPrintf(PETSC_COMM_WORLD, "1 Compute tau, r, cos_theta at each quad point\n"); CHKERRQ(ierr);
        ierr = ComputeEmissionDistanceAndCos(Tauq, Rq, Cosq, to, xo, c); CHKERRQ(ierr);


        // -------------------------------------------------------------
        // Compute P', P'n, P't at emission time τ:
        //--------------------------------------------------------------

        ierr = PetscPrintf(PETSC_COMM_WORLD, "2 Compute p', p'n, p't at emission time tau\n"); CHKERRQ(ierr);
        ierr = ComputePressureAtEmissionTime(P, Pn, Pt, Tauq); CHKERRQ(ierr);


        // -------------------------------------------------------------
        // Compute Kirchhoff Integral
        //--------------------------------------------------------------
        PetscReal integral;
        ierr = PetscPrintf(PETSC_COMM_WORLD, "3 Compute Kirchhoff Integral\n"); CHKERRQ(ierr);
        integral = ComputeKirchhoffIntegral(P, Pn, Pt, Tauq, Rq, Cosq, h, to, xo, c);


        // -------------------------------------------------------------
        // Write perturbation pressure
        //--------------------------------------------------------------
        ierr = PetscPrintf(PETSC_COMM_WORLD, "4 Writing perturbation pressure at observer location xo and observer time to\n\n"); CHKERRQ(ierr);
        if (MyPID == 0) {
          FILE *f = fopen("p_kirchhoff.dat","a+");
          ierr = PetscFPrintf(PETSC_COMM_SELF, f, "%.16e\t%.16e\n", to, integral); CHKERRQ(ierr);
          fclose(f);
        }

        //increment time:
        to += dt;

    }
      


    //destroy petsc objects
    VecDestroy(&Tauq);
    VecDestroy(&Cosq);
    VecDestroy(&Rq);
    VecDestroy(&P);
    VecDestroy(&Pn);
    VecDestroy(&Pt);

    return PetscFinalize();
}
