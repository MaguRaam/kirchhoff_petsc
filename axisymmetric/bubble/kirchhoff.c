/*
 * kirchhoff.c
 *      Author: Magu
 */ 
#include "hype.h"


//----------------------------------------------------------------------------------------------------
// Get cell indices on Kirchhoff surface
//----------------------------------------------------------------------------------------------------

PetscErrorCode  GetCellIndexOnSurface(DM da, AppCtx* Ctx){

    PetscErrorCode  ierr;
    PetscInt        xs, ys, xm, ym, i, j;
    PetscScalar    *ia_curved, *ja_top, *ja_bottom;
    PetscInt        index_top = 0, index_bottom = 0, index_curved = 0;
    PetscViewer     viewer_top, viewer_bottom, viewer_curved;

    ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);

    //get no of cells on Kirchhoff top, bottom and curved surface in the current process:
    Ctx->surface.ncells_top = 0;
    Ctx->surface.ncells_bottom = 0;
    Ctx->surface.ncells_curved = 0;

    for (j = ys; j < ys+ym; ++j) {
        for (i = xs; i < xs+xm; ++i) {

            //check if the cell is located on the Kirchhoff top surface
            if (i == Ctx->surface.i_top && j <= Ctx->surface.j_curved)
                Ctx->surface.ncells_top++;

            //check if the cell is located on the Kirchhoff bottom surface
            if (i == Ctx->surface.i_bottom && j <= Ctx->surface.j_curved)
                Ctx->surface.ncells_bottom++;

            //check if the cell is located on the Kirchhoff curved surface
            if (i >= Ctx->surface.i_top && i <= Ctx->surface.i_bottom && j == Ctx->surface.j_curved)
                Ctx->surface.ncells_curved++;  


        }
    }


    //create j_top, j_bottom and i_curved index vector:
    ierr = VecCreate(PETSC_COMM_WORLD, &(Ctx->surface.i_curved)); CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD, &(Ctx->surface.j_top)); CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD, &(Ctx->surface.j_bottom)); CHKERRQ(ierr);

    ierr = VecSetSizes(Ctx->surface.i_curved, Ctx->surface.ncells_curved, PETSC_DECIDE); CHKERRQ(ierr);
    ierr = VecSetSizes(Ctx->surface.j_top, Ctx->surface.ncells_top, PETSC_DECIDE); CHKERRQ(ierr);
    ierr = VecSetSizes(Ctx->surface.j_bottom, Ctx->surface.ncells_bottom, PETSC_DECIDE); CHKERRQ(ierr);

    ierr = VecSetUp(Ctx->surface.i_curved); CHKERRQ(ierr);
    ierr = VecSetUp(Ctx->surface.j_top); CHKERRQ(ierr);
    ierr = VecSetUp(Ctx->surface.j_bottom); CHKERRQ(ierr);


    //get array from vector:
    ierr = VecGetArray(Ctx->surface.i_curved, &ia_curved); CHKERRQ(ierr);
    ierr = VecGetArray(Ctx->surface.j_top, &ja_top); CHKERRQ(ierr);
    ierr = VecGetArray(Ctx->surface.j_bottom, &ja_bottom); CHKERRQ(ierr);


    //collect indices on the Kirchhoff surface:
    for (j = ys; j < ys+ym; ++j) {
        for (i = xs; i < xs+xm; ++i) {

            //check if the cell is located on the Kirchhoff top surface
            if (i == Ctx->surface.i_top && j <= Ctx->surface.j_curved){ 
                ja_top[index_top] = j;
                index_top++;
            }

            //check if the cell is located on the Kirchhoff bottom surface
            if (i == Ctx->surface.i_bottom && j <= Ctx->surface.j_curved){
                ja_bottom[index_bottom] = j;
                index_bottom++;
            }
                

            //check if the cell is located on the Kirchhoff curved surface
            if (i >= Ctx->surface.i_top && i <= Ctx->surface.i_bottom && j == Ctx->surface.j_curved){
                ia_curved[index_curved] = i;
                index_curved++;
            }
                

        }
    }

    //vec restore array:
    ierr = VecRestoreArray(Ctx->surface.i_curved, &ia_curved); CHKERRQ(ierr);
    ierr = VecRestoreArray(Ctx->surface.j_top, &ja_top); CHKERRQ(ierr);
    ierr = VecRestoreArray(Ctx->surface.j_bottom, &ja_bottom); CHKERRQ(ierr);

    //write cell index of Kirchhoff curved surface:
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "kirchhoff/input/curved/i.output", &viewer_curved);
    ierr = VecView(Ctx->surface.i_curved, viewer_curved); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer_curved); CHKERRQ(ierr);

    //write cell index of Kirchhoff top surface:
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "kirchhoff/input/top/j.output", &viewer_top);
    ierr = VecView(Ctx->surface.j_top, viewer_top); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer_top); CHKERRQ(ierr);

    //write cell index of Kirchhoff bottom surface:
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "kirchhoff/input/bottom/j.output", &viewer_bottom);
    ierr = VecView(Ctx->surface.j_bottom, viewer_bottom); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer_bottom); CHKERRQ(ierr);

    //destroy vectors:
    VecDestroy(&(Ctx->surface.i_curved));
    VecDestroy(&(Ctx->surface.j_top));
    VecDestroy(&(Ctx->surface.j_bottom));


    return ierr;
}


//----------------------------------------------------------------------------------------------------
// Initialize p', p'r vectors for Kirchhoff  surface 
//----------------------------------------------------------------------------------------------------
PetscErrorCode InitializePressureOnSurface(AppCtx* Ctx){

    PetscErrorCode  ierr;

    //allocate memory for p' and p'r vectors stores on Kirchhoff curved surface:
    ierr = VecCreate(PETSC_COMM_WORLD, &(Ctx->surface.P_curved)); CHKERRQ(ierr);
    ierr = VecSetSizes(Ctx->surface.P_curved, 2*Ctx->surface.ncells_curved, PETSC_DECIDE); CHKERRQ(ierr);
    ierr = VecSetUp(Ctx->surface.P_curved); CHKERRQ(ierr);
    ierr = VecDuplicate(Ctx->surface.P_curved, &(Ctx->surface.Pr_curved)); CHKERRQ(ierr);

    //allocate memory for p' and p'z vectors stores on Kirchhoff top surface:
    ierr = VecCreate(PETSC_COMM_WORLD, &(Ctx->surface.P_top)); CHKERRQ(ierr);
    ierr = VecSetSizes(Ctx->surface.P_top, 2*Ctx->surface.ncells_top, PETSC_DECIDE); CHKERRQ(ierr);
    ierr = VecSetUp(Ctx->surface.P_top); CHKERRQ(ierr);
    ierr = VecDuplicate(Ctx->surface.P_top, &(Ctx->surface.Pz_top)); CHKERRQ(ierr);

    //allocate memory for p' and p'z vectors stores on Kirchhoff bottom surface:
    ierr = VecCreate(PETSC_COMM_WORLD, &(Ctx->surface.P_bottom)); CHKERRQ(ierr);
    ierr = VecSetSizes(Ctx->surface.P_bottom, 2*Ctx->surface.ncells_bottom, PETSC_DECIDE); CHKERRQ(ierr);
    ierr = VecSetUp(Ctx->surface.P_bottom); CHKERRQ(ierr);
    ierr = VecDuplicate(Ctx->surface.P_bottom, &(Ctx->surface.Pz_bottom)); CHKERRQ(ierr);


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

//----------------------------------------------------------------------------------------------------
// Reconstruct polynomial on a given cell
//----------------------------------------------------------------------------------------------------
void ReconstructOnGivenCell(Field** u, PetscInt i, PetscInt j, PetscReal coeffs[], PetscReal pcoeffs[]){

    PetscReal p_x_loc[s_width], p_y_loc[s_width], p_xy_loc[s_width];

    //get perturbed pressure values from the stencil:
    p_x_loc[0] = GetPerturbedPressure(&u[j][i-2]); 
    p_x_loc[1] = GetPerturbedPressure(&u[j][i-1]); 
    p_x_loc[2] = GetPerturbedPressure(&u[j][i]); 
    p_x_loc[3] = GetPerturbedPressure(&u[j][i+1]); 
    p_x_loc[4] = GetPerturbedPressure(&u[j][i+2]);

    p_y_loc[0] = GetPerturbedPressure(&u[j-2][i]); 
    p_y_loc[1] = GetPerturbedPressure(&u[j-1][i]); 
    p_y_loc[2] = GetPerturbedPressure(&u[j][i]); 
    p_y_loc[3] = GetPerturbedPressure(&u[j+1][i]); 
    p_y_loc[4] = GetPerturbedPressure(&u[j+2][i]);

    p_xy_loc[0] = GetPerturbedPressure(&u[j][i]);
    p_xy_loc[1] = GetPerturbedPressure(&u[j+1][i+1]);
    p_xy_loc[2] = GetPerturbedPressure(&u[j-1][i+1]);
    p_xy_loc[3] = GetPerturbedPressure(&u[j+1][i-1]);
    p_xy_loc[4] = GetPerturbedPressure(&u[j-1][i-1]);

    //reconstruct polynomial (get polynomial coefficients)
    weno(p_x_loc, p_y_loc, p_xy_loc, pcoeffs, coeffs);

}




//----------------------------------------------------------------------------------------------------------------
// Interpolate p' and p'r at quadrature points on Kirchhoff surface and Write to a file for Kirchhoff solver
//----------------------------------------------------------------------------------------------------------------

PetscErrorCode GetPressureOnSurface(Vec U, DM da, AppCtx* Ctx, PetscInt step){

    PetscErrorCode  ierr;
    PetscInt        xs, ys, xm, ym, i, j;
    Field           **u;
    PetscReal       coeffs[nDOF], pcoeffs[nDOF];
    PetscReal       *p_curved, *p_top, *p_bottom;
    PetscReal       *pr_curved, *pz_top, *pz_bottom;
    PetscInt        index_top = 0, index_bottom = 0, index_curved = 0;
    PetscReal       h = Ctx->h, p_r, p_z;
    PetscViewer     viewer_curved, viewer_top, viewer_bottom;
    
    ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);
    

    // Scatter global->local to have access to the required ghost values 
    ierr = DMGlobalToLocalBegin(da, U, INSERT_VALUES, Ctx->localU);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da, U, INSERT_VALUES,   Ctx->localU);CHKERRQ(ierr);


    // Read the local solution to the array u  
    ierr = DMDAVecGetArrayRead(da, Ctx->localU, &u); CHKERRQ(ierr); 

     //get array from vector P and Pr:
    ierr = VecGetArray((Ctx->surface.P_curved), &p_curved); CHKERRQ(ierr);
    ierr = VecGetArray((Ctx->surface.Pr_curved), &pr_curved); CHKERRQ(ierr);
    ierr = VecGetArray((Ctx->surface.P_top), &p_top); CHKERRQ(ierr);
    ierr = VecGetArray((Ctx->surface.Pz_top), &pz_top); CHKERRQ(ierr);
    ierr = VecGetArray((Ctx->surface.P_bottom), &p_bottom); CHKERRQ(ierr);
    ierr = VecGetArray((Ctx->surface.Pz_bottom), &pz_bottom); CHKERRQ(ierr);



     //loop over cells in the current process:
    for (j = ys; j < ys+ym; ++j) {
        for (i = xs; i < xs+xm; ++i) {

            //check if the cell is located on the Kirchhoff curved surface
            if (i >= Ctx->surface.i_top && i <= Ctx->surface.i_bottom && j == Ctx->surface.j_curved){

                //reconstruct polynomial
                ReconstructOnGivenCell(u, i, j, coeffs, pcoeffs);

                //loop over quadrature points on the line:
                for (PetscInt q = 0; q < N_gp2; ++q){

                    //compute p' and p'r:
                    p_curved[index_curved] = evaluate_polynomial(x_gp2[q], 0.0, coeffs);
                    
                   	evaluate_grad(coeffs, x_gp2[q], 0.0, h, &p_z, &p_r);
                    pr_curved[index_curved] = p_r;
                    index_curved++;

                }
            }

            
            //check if the cell is located on the Kirchhoff top surface
            if (i == Ctx->surface.i_top && j <= Ctx->surface.j_curved){

                //reconstruct polynomial
                ReconstructOnGivenCell(u, i, j, coeffs, pcoeffs);

                //loop over quadrature points on the line:
                for (PetscInt q = 0; q < N_gp2; ++q){
                
                    //compute p' and p'z:
                    p_top[index_top] = evaluate_polynomial(0.0, x_gp2[q], coeffs);

                    evaluate_grad(coeffs, 0.0, x_gp2[q], h, &p_z, &p_r);
                    pz_top[index_top] = -p_z;
                    index_top++;

                }
            }

            //check if the cell is located on the Kirchhoff bottom surface
            if (i == Ctx->surface.i_bottom && j <= Ctx->surface.j_curved){

                //reconstruct polynomial
                ReconstructOnGivenCell(u, i, j, coeffs, pcoeffs);

                //loop over quadrature points on the line:
                for (PetscInt q = 0; q < N_gp2; ++q){
                    
                    //compute p' and p'z:
                    p_bottom[index_bottom] = evaluate_polynomial(0.0, x_gp2[q], coeffs);

                    evaluate_grad(coeffs, 0.0, x_gp2[q], h, &p_z, &p_r);
                    pz_bottom[index_bottom] = +p_z;
                    index_bottom++;
                
                }         

            }


        }
    }

    //restore array
    ierr = VecRestoreArray((Ctx->surface.P_curved), &p_curved); CHKERRQ(ierr);
    ierr = VecRestoreArray((Ctx->surface.Pr_curved), &pr_curved); CHKERRQ(ierr);
    ierr = VecRestoreArray((Ctx->surface.P_top), &p_top); CHKERRQ(ierr);
    ierr = VecRestoreArray((Ctx->surface.Pz_top), &pz_top); CHKERRQ(ierr);
    ierr = VecRestoreArray((Ctx->surface.P_bottom), &p_bottom); CHKERRQ(ierr);
    ierr = VecRestoreArray((Ctx->surface.Pz_bottom), &pz_bottom); CHKERRQ(ierr);


    //write pressure and its normal data
    char filename_curved[50];
    sprintf(filename_curved, "kirchhoff/input/curved/p-%05d.dat", step); 
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename_curved, &viewer_curved);
    ierr = VecView(Ctx->surface.P_curved, viewer_curved); CHKERRQ(ierr);
    ierr = VecView(Ctx->surface.Pr_curved, viewer_curved); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer_curved); CHKERRQ(ierr);
    
    char filename_top[50];
    sprintf(filename_top, "kirchhoff/input/top/p-%05d.dat", step); 
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename_top, &viewer_top);
    ierr = VecView(Ctx->surface.P_top, viewer_top); CHKERRQ(ierr);
    ierr = VecView(Ctx->surface.Pz_top, viewer_top); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer_top); CHKERRQ(ierr);
    
    char filename_bottom[50];
    sprintf(filename_bottom, "kirchhoff/input/bottom/p-%05d.dat", step); 
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename_bottom, &viewer_bottom);
    ierr = VecView(Ctx->surface.P_bottom, viewer_bottom); CHKERRQ(ierr);
    ierr = VecView(Ctx->surface.Pz_bottom, viewer_bottom); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer_bottom); CHKERRQ(ierr);


    ierr = DMDAVecRestoreArrayRead(da, Ctx->localU, &u); CHKERRQ(ierr);
    
    return ierr;
}


//--------------------------------------------------------------------------------------
// Release all the memory related to Kirchhoff  surface
//--------------------------------------------------------------------------------------

PetscErrorCode  DestroyKirchoffSurface(AppCtx* Ctx){

    VecDestroy(&(Ctx->surface.P_curved));
    VecDestroy(&(Ctx->surface.Pr_curved));
    VecDestroy(&(Ctx->surface.P_top));
    VecDestroy(&(Ctx->surface.Pz_top));
    VecDestroy(&(Ctx->surface.P_bottom));
    VecDestroy(&(Ctx->surface.Pz_bottom));


    return 0;
}
