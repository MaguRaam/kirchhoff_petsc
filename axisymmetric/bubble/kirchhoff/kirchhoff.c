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


// -------------------------------------------------
// Compute area of the cylindrical Kirchhoff surface
//--------------------------------------------------
PetscErrorCode  ComputeArea(DM da, PetscReal R, PetscReal hx, PetscReal hy){

  PetscErrorCode ierr;

  DM          coordDA;
  Vec         coordinates;
  DMDACoor2d  **coords;
  

  PetscReal  global_area = 0.0, local_area = 0.0;
  PetscInt    xs, ys, xm, ym, i, j;
  PetscReal   xc, yc;

  ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);
  ierr = DMGetCoordinateDM(da, &coordDA);CHKERRQ(ierr);
  ierr = DMGetCoordinates(da, &coordinates);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(coordDA, coordinates, &coords);CHKERRQ(ierr);

  for (j = ys; j < ys + ym; ++j)
    for (i = xs; i < xs + xm; ++i)
      local_area += 1.0;

  //accumulate area from all the process:
  ierr = MPI_Allreduce(&local_area, &global_area, 1, MPIU_REAL, MPIU_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);

  //scale area:
  global_area = R*global_area*hx*hy;

  ierr = PetscPrintf(PETSC_COMM_WORLD, "exact area = %.16g\n", 2.0*PETSC_PI*R*2.0*0.38); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "numerical area = %.16g\n", global_area); CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords);CHKERRQ(ierr);

  return ierr;
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
  // Set important user defined parameters
  //---------------------------------------------
  
  PetscReal x_min       =   -0.38;                                        //z direction begining                                        // z direction begining
  PetscReal x_max       =    0.38;                                        //z direction ending                                   // z direction ending
  PetscReal k_min       =    0.0;                                         // theta direction begining
  PetscReal k_max       =    2.0*PETSC_PI;                                // theta direction ending
  PetscInt  Nx          =    500;                                         // no of cells along z axis
  PetscInt  Ny          =    2000;                                        // no of cells along theta axis
  PetscReal hx          =    (x_max - x_min)/(PetscReal)(Nx);             // grid size along z axis
  PetscReal hy          =    (k_max - k_min)/(PetscReal)(Ny);             // grid size along theta axis
  PetscReal R           =    1.0;                                         // Location of Kirchhoff surface


  // --------------------------------------------
  // Data members  
  //---------------------------------------------  
  
  DM                      da;                                 // Grid object
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
  
  ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DMDA_STENCIL_BOX, Nx, Ny, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &da); CHKERRQ(ierr);
  ierr = DMSetUp(da);CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(da, x_min + 0.5*hx, x_max - 0.5*hx, k_min + 0.5*hy, k_max - 0.5*hy , 0.0, 0.0); CHKERRQ(ierr);




  // --------------------------------------------
  // Free all the memory, finalize MPI and exit   
  //---------------------------------------------

  ierr = DMDestroy(&da); CHKERRQ(ierr);

  return PetscFinalize();
}