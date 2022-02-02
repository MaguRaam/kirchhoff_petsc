/*
 * kirchhoff.c
 *      Author: Magu
 */ 
#include "hype.h"


//----------------------------------------------------------------------------------------------------
// Get i index of cells on Kirchhoff Line surface
//----------------------------------------------------------------------------------------------------

PetscErrorCode  GetCellIndexOnLineSurface(DM da, AppCtx* Ctx){

    PetscErrorCode  ierr;
    PetscInt        xs, ys, xm, ym, i, j;
    PetscScalar*    iq;
    PetscInt        index = 0;
    PetscViewer     viewer;

    ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);

    //get no of cells on Kirchhoff Line surface in the current process:
    Ctx->surface.ncells = 0;

    for (j = ys; j < ys+ym; ++j) {
        for (i = xs; i < xs+xm; ++i) {

            //check if the cell is located on the Kirchhoff line surface
            if (j == Ctx->surface.j)
                Ctx->surface.ncells++;  

        }
    }
    
    //create i index vector:
    ierr = VecCreate(PETSC_COMM_WORLD, &(Ctx->surface.i)); CHKERRQ(ierr);
    ierr = VecSetSizes(Ctx->surface.i, Ctx->surface.ncells, PETSC_DECIDE); CHKERRQ(ierr);
    ierr = VecSetUp(Ctx->surface.i); CHKERRQ(ierr);

    //get array from vector:
    ierr = VecGetArray(Ctx->surface.i, &iq); CHKERRQ(ierr);

    //collect i index on the Kirchhoff cells:
    for (j = ys; j < ys+ym; ++j) {
        for (i = xs; i < xs+xm; ++i) {

            //check if the cell is located on the Kirchhoff line surface
            if (j == Ctx->surface.j){
                iq[index] = i;
                index++;
            }
                

        }
    }

    //vec restore array:
    ierr = VecRestoreArray(Ctx->surface.i, &iq); CHKERRQ(ierr);


    //write i cell index of Kirchhoff line surface:
    ierr = PetscPrintf(PETSC_COMM_WORLD, "write i cell index of Kirchhoff line surface\n"); CHKERRQ(ierr);
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "kirchhoff/input/i.output", &viewer);
    ierr = VecView(Ctx->surface.i, viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    return ierr;
}


//----------------------------------------------------------------------------------------------------
// Initialize p', p'r vectors for Kirchhoff Line surface 
//----------------------------------------------------------------------------------------------------
PetscErrorCode InitializePressureOnLineSurface(DM da, AppCtx* Ctx){

    PetscErrorCode  ierr;

    //allocate memory for p' and p'r vectors stores on Kirchhoff Line surface:
    ierr = VecCreate(PETSC_COMM_WORLD, &(Ctx->surface.P)); CHKERRQ(ierr);
    ierr = VecSetSizes(Ctx->surface.P, Ctx->surface.ncells, PETSC_DECIDE); CHKERRQ(ierr);
    ierr = VecSetUp(Ctx->surface.P); CHKERRQ(ierr);
    ierr = VecDuplicate(Ctx->surface.P, &(Ctx->surface.Pr)); CHKERRQ(ierr);

    return  ierr;
}



//----------------------------------------------------------------------------------------------------
// Get perturbed pressure p' = p - p0 from conserved state variable
//----------------------------------------------------------------------------------------------------
PetscReal GetPerturbedPressure(const Field* q){

    PetscReal Q[nVar], V[nVar];

    for (int c = 0 ; c < nVar; ++c) Q[c] = q->comp[c];
    PDECons2Prim(Q, V);

    return V[3] - p0;
}



//----------------------------------------------------------------------------------------------------------------
// Interpolate p' and p'r at quadrature points on Kirchhoff line surface and Write to a file for Kirchhoff solver
//----------------------------------------------------------------------------------------------------------------

PetscErrorCode GetPressureOnLineSurface(Vec U, DM da, AppCtx* Ctx, PetscInt step){

    PetscErrorCode  ierr;
    PetscInt        xs, ys, xm, ym, i, j;
    Field           **u;
    PetscReal       pR, pL, pC;  
    PetscReal*      p;
    PetscReal*      pr;
    PetscInt        index = 0;
    PetscReal       h = Ctx->h;
    
    ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);
    

    // Scatter global->local to have access to the required ghost values 
    ierr = DMGlobalToLocalBegin(da, U, INSERT_VALUES, Ctx->localU);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da, U, INSERT_VALUES,   Ctx->localU);CHKERRQ(ierr);


    // Read the local solution to the array u  
    ierr = DMDAVecGetArrayRead(da, Ctx->localU, &u); CHKERRQ(ierr); 

     //get array from vector P and Pr:
    ierr = VecGetArray((Ctx->surface.P), &p); CHKERRQ(ierr);
    ierr = VecGetArray((Ctx->surface.Pr), &pr); CHKERRQ(ierr);

     //loop over cells in the current process:
    for (j = ys; j < ys+ym; ++j) {
        for (i = xs; i < xs+xm; ++i) {

            //check if the cell is located on the Kirchhoff line surface
            if (j == Ctx->surface.j){
                
                //collect perturbed pressure p' = p - p0 values from the neighbors in the stencil:                 
                pL = GetPerturbedPressure(&u[j-1][i]); 
                pC = GetPerturbedPressure(&u[j][i]); 
                pR = GetPerturbedPressure(&u[j+1][i]); 


                //compute p' and p'r at cell center:
                p[index]  = pC; 
                pr[index] = (pR - pL)/(2.0*h);
                

                //update cell index:
                index++;
            }

        }
    }

    //restore array
    ierr = VecRestoreArray((Ctx->surface.P), &p); CHKERRQ(ierr);
    ierr = VecRestoreArray((Ctx->surface.Pr), &pr); CHKERRQ(ierr);

    ierr = DMDAVecRestoreArrayRead(da, Ctx->localU, &u); CHKERRQ(ierr);
    
    return ierr;
}

//--------------------------------------------------------------------------------------
// Release all the memory related to Kirchhoff line surface
//--------------------------------------------------------------------------------------

PetscErrorCode  DestroyKirchoffLineSurface(AppCtx* Ctx){

    VecDestroy(&(Ctx->surface.i));
    VecDestroy(&(Ctx->surface.P));
    VecDestroy(&(Ctx->surface.Pr));

    return 0;
}
