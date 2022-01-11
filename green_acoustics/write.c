/*
 * monitor_function.c
 *      Author: sunder
 */ 

#include "hype.h"  

//----------------------------------------------------------------------------
// Write pressure data in vts format
//----------------------------------------------------------------------------


PetscErrorCode WriteVtk(DM da, Vec U, PetscInt step)
{
   PetscErrorCode ierr; 
   
   char filename[20];
   sprintf(filename, "sol-%05d.vtk", step); // 4 is the padding level, increase it for longer simulations 
   PetscViewer viewer;  
   PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
   PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_VTK);
   ierr = DMView(da, viewer);
   VecView(U, viewer);
   
   ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    
    return ierr; 
}
