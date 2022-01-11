/*
 * kirchhoff.c
 *      Author: Magu
 */ 
#include "hype.h"

//--------------------------------------------------------------------------------------
// Compute no of cells on each face of the Kirchoff box surface in the current process
//--------------------------------------------------------------------------------------

PetscErrorCode  GetNumberOfCellsOnBoxFaces(DM da, AppCtx* Ctx){       
    
    PetscErrorCode ierr;
    PetscInt    xs, ys, zs, xm, ym, zm, i, j, k, f;
    PetscInt    index[3];
    PetscInt    f1[3] = {1, 0, 0};
    PetscInt    f2[3] = {2, 2, 1};

    //initialize ncells in each face in the current process to zero:
    for (f = 0; f < 6; ++f) Ctx->box.ncells[f] = 0;


    ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);

    //loop over cells in the current process and compute no of cells in each face of the box in the current process:
    
    for (k = zs; k < zs + zm; ++k){
        for (j = ys; j < ys + ym; ++j){
            for (i = xs; i < xs + xm; ++i){

                //get the index of the current cell:
                index[0] = i; index[1] = j; index[2] =  k;

                //loop over faces in each direction of the box
                for (f = 0; f < 3; ++f){

                    //check if the given cell index is on the face, if true increase the number of cells:
                    if ( index[f1[f]] >=  Ctx->box.lower[f1[f]] &&  index[f1[f]] <=  Ctx->box.upper[f1[f]] && index[f2[f]] >=  Ctx->box.lower[f2[f]] &&  index[f2[f]] <=  Ctx->box.upper[f2[f]]){
                        if (index[f] == Ctx->box.lower[f]) Ctx->box.ncells[2*f]++;
                        if (index[f] == Ctx->box.upper[f]) Ctx->box.ncells[2*f + 1]++;
                    }

                } 

            }
        }
    }

    return ierr;
}


//--------------------------------------------------------------------------------------
// Collect cell indices in each face of the Kirchoff box surface in the current process
//--------------------------------------------------------------------------------------

PetscErrorCode  GetCellIndicesOnBoxFaces(DM da, AppCtx* Ctx){       
    
    PetscErrorCode ierr;
    PetscInt    xs, ys, zs, xm, ym, zm, i, j, k, f;
    PetscInt    index[3];
    PetscInt    f1[3] = {1, 0, 0};
    PetscInt    f2[3] = {2, 2, 1};
    PetscInt    count[6] = {0};


    ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);

    //allocate memory for cell indices of each face of the Box surface:
    for (f = 0; f < 6; ++f) Ctx->box.indices[f] = (CellIndex*)malloc(Ctx->box.ncells[f]*sizeof(CellIndex));


    //loop over cells in the current process and collect the cell indices in each face of the box in the current process:
    
    for (k = zs; k < zs + zm; ++k){
        for (j = ys; j < ys + ym; ++j){
            for (i = xs; i < xs + xm; ++i){

                //get the index of the current cell:
                index[0] = i; index[1] = j; index[2] =  k;

                //loop over faces in each direction of the box
                for (f = 0; f < 3; ++f){

                    //check if the given cell index is on the face, if true collect the cell index:
                    if ( index[f1[f]] >=  Ctx->box.lower[f1[f]] &&  index[f1[f]] <=  Ctx->box.upper[f1[f]] && index[f2[f]] >=  Ctx->box.lower[f2[f]] &&  index[f2[f]] <=  Ctx->box.upper[f2[f]]){
                        
                        if (index[f] == Ctx->box.lower[f]){
                            Ctx->box.indices[2*f][count[2*f]].i = i;
                            Ctx->box.indices[2*f][count[2*f]].j = j;
                            Ctx->box.indices[2*f][count[2*f]].k = k;
                            count[2*f]++;
                        }

                        if (index[f] == Ctx->box.upper[f]){
                            Ctx->box.indices[2*f + 1][count[2*f + 1]].i = i;
                            Ctx->box.indices[2*f + 1][count[2*f + 1]].j = j;
                            Ctx->box.indices[2*f + 1][count[2*f + 1]].k = k;
                            count[2*f + 1]++;
                        }
                    }

                } 

            }
        }
    }


    return ierr;
}


//----------------------------------------------------------------------------------------------------
// Compute the quadrature points and normals on the box faces and Write to a file for Kirchhoff solver
//----------------------------------------------------------------------------------------------------

PetscErrorCode  GetQuadPointsOnBoxFaces(DM da, AppCtx* Ctx){

    PetscErrorCode  ierr;
    DM              coordDA;
    Vec             coordinates;
    DMDACoor3d      ***coords;
    PetscReal       xc, yc, zc, h = Ctx->h;
    PetscInt        nfaces = 0, i, j, k;
    PetscInt        f, c, q, iq = 0;
    PetscScalar*    xq;
    PetscScalar*    nq;
    PetscViewer     viewer;

    //get coordinates array:
    ierr = DMGetCoordinateDM(da, &coordDA);CHKERRQ(ierr);
    ierr = DMGetCoordinates(da, &coordinates);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(coordDA, coordinates, &coords);CHKERRQ(ierr);


    //compute total no of cell faces attached to the Kirchoff surface on the current process
    for (f = 0; f < 6; ++f) nfaces += Ctx->box.ncells[f];
    
    //compute total no of quad pts attached to the Kirchoff surface on the current process
    Ctx->box.nqpts = nfaces*N_gp2d; 

    //create petsc vector to store quadrature points in the current process
    ierr = VecCreate(PETSC_COMM_WORLD, &(Ctx->box.Xq)); CHKERRQ(ierr);
    ierr = VecSetBlockSize(Ctx->box.Xq, 3); CHKERRQ(ierr);
    ierr = VecSetSizes(Ctx->box.Xq, 3*(Ctx->box.nqpts), PETSC_DECIDE); CHKERRQ(ierr);
    ierr = VecSetUp(Ctx->box.Xq); CHKERRQ(ierr);
    
    //Create petsc vector to store normal vector at quadrature points in the current process
    ierr = VecDuplicate(Ctx->box.Xq, &(Ctx->box.Nq)); CHKERRQ(ierr);


    //get array of quad points
    ierr = VecGetArray(Ctx->box.Xq, &xq); CHKERRQ(ierr);
    Point* p = (Point*)xq;

    //get array of normal vector
    ierr = VecGetArray(Ctx->box.Nq, &nq); CHKERRQ(ierr);
    Point* n = (Point*)nq;


    //loop over faces on the Kirchhoff surface
    for (f = 0; f < 6; ++f){

        //loop over cells sharing the face
        for (c = 0; c < Ctx->box.ncells[f]; ++c){

            //get cell index:
            i = Ctx->box.indices[f][c].i;
            j = Ctx->box.indices[f][c].j;
            k = Ctx->box.indices[f][c].k;

            //get cell center:
            xc = coords[k][j][i].x; 
            yc = coords[k][j][i].y;
            zc = coords[k][j][i].z;

            //loop over quadrature points on the cell face attached to Kirchhoff surface and compute coordinates and normal

            for (q = 0; q < N_gp2d; ++q){

                switch (f)
                {
                case 0:
                    p[iq].x = xc - 0.5 * h;
                    p[iq].y = yc + x_gp2d[q] * h;
                    p[iq].z = zc + y_gp2d[q] * h;

                    n[iq].x = -1.0;
                    n[iq].y = 0.0;
                    n[iq].z = 0.0;
                    break;
                case 1:
                    p[iq].x = xc + 0.5*h;          
                    p[iq].y = yc + x_gp2d[q]*h;         
                    p[iq].z = zc + y_gp2d[q]*h;

                    n[iq].x     =  1.0;        
                    n[iq].y     =  0.0;         
                    n[iq].z     =  0.0;
                    break;
                case 2:
                    p[iq].x = xc + x_gp2d[q]*h;    
                    p[iq].y = yc - 0.5*h ;              
                    p[iq].z = zc + y_gp2d[q]*h;
                
                    n[iq].x     =  0.0;        
                    n[iq].y     = -1.0;         
                    n[iq].z     =  0.0;
                    break;
                case 3:
                    p[iq].x = xc + x_gp2d[q]*h;    
                    p[iq].y = yc + 0.5*h ;              
                    p[iq].z = zc + y_gp2d[q]*h;
                
                    n[iq].x     =  0.0;        
                    n[iq].y     =  1.0;         
                    n[iq].z     =  0.0;
                    break;
                case 4:
                    p[iq].x = xc + x_gp2d[q]*h;    
                    p[iq].y = yc + y_gp2d[q]*h;         
                    p[iq].z = zc - 0.5*h;
                    
                    n[iq].x     =  0.0;        
                    n[iq].y     =  0.0;         
                    n[iq].z     = -1.0;
                    break;
                case 5:
                    p[iq].x = xc + x_gp2d[q]*h;    
                    p[iq].y = yc + y_gp2d[q]*h;         
                    p[iq].z = zc + 0.5*h;
                    
                    n[iq].x     =  0.0;        
                    n[iq].y     =  0.0;         
                    n[iq].z     =  1.0;
                    break;
                }
                
                
                iq++;
            }
        }
    }

    //restore array
    ierr = VecRestoreArray(Ctx->box.Xq, &xq); CHKERRQ(ierr);
    ierr = VecRestoreArray(Ctx->box.Nq, &nq); CHKERRQ(ierr);

    //Compute total no of quad points:
    PetscInt tmp;
    ierr = VecGetSize(Ctx->box.Xq, &tmp); CHKERRQ(ierr);
    Ctx->box.Tnqpts = tmp/3;
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Total Nqpts = %d\n", Ctx->box.Tnqpts);CHKERRQ(ierr);

    //write quadrature points in a file:
    ierr = PetscPrintf(PETSC_COMM_WORLD, "writing quadrature points and normals on the Kirchoff surface\n"); CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"kirchhoff/q_points.dat",FILE_MODE_WRITE,&viewer); CHKERRQ(ierr);
    
    ierr = VecView(Ctx->box.Xq, viewer); CHKERRQ(ierr);
    ierr = VecView(Ctx->box.Nq, viewer); CHKERRQ(ierr);

    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords);CHKERRQ(ierr);

    return ierr;

}


//----------------------------------------------------------------------------------------------------
// Initialize p, pn vectors stored for Kirchhoff surface 
//----------------------------------------------------------------------------------------------------
PetscErrorCode InitializePressureOnBoxFaces(AppCtx* Ctx){

    PetscErrorCode  ierr;

    //create P, Pn vector defined on quad point on the box
    ierr = VecCreate(PETSC_COMM_WORLD, &(Ctx->box.P)); CHKERRQ(ierr);
    ierr = VecSetSizes(Ctx->box.P, Ctx->box.nqpts, PETSC_DECIDE); CHKERRQ(ierr);
    ierr = VecSetUp(Ctx->box.P); CHKERRQ(ierr);
    ierr = VecDuplicate(Ctx->box.P, &(Ctx->box.Pn)); CHKERRQ(ierr);

    return ierr;
}



//----------------------------------------------------------------------------------------------------
// Interpolate p and pn at quadrature points and Write to a file for Kirchhoff solver
//----------------------------------------------------------------------------------------------------
PetscErrorCode GetPressureOnBoxFaces(Vec U, DM da, AppCtx* Ctx, PetscInt step){

    PetscErrorCode  ierr;
    PetscReal*      p;
    PetscReal*      pn;
    PetscReal       u_x_loc[5], u_y_loc[5], u_z_loc[5], u_xy_loc[5], u_yz_loc[5], u_zx_loc[5], u_xyz_loc[8]; 
    PetscInt        i, j, k, f, c, q, iq = 0, z;
    Field           ***u;
    PetscReal       coeffs[2][dofs_per_cell];
    PetscViewer     viewer;
    PetscReal       h = Ctx->h, p_x, p_y, p_z;
    PetscInt		  oned_begin;
    
    // Scatter global->local to have access to the required ghost values 
    ierr = DMGlobalToLocalBegin(da, U, INSERT_VALUES, Ctx->localU);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da, U, INSERT_VALUES,   Ctx->localU);CHKERRQ(ierr);

    // Read the local solution to the array u  
    ierr = DMDAVecGetArrayRead(da, Ctx->localU, &u); CHKERRQ(ierr); 

    //get array from vector P, Pn:
    ierr = VecGetArray((Ctx->box.P), &p); CHKERRQ(ierr);
    ierr = VecGetArray((Ctx->box.Pn), &pn); CHKERRQ(ierr);
 

    //loop over faces on the Kirchhoff surface
    for (f = 0; f < 6; ++f){

        //loop over cells sharing the face
        for (c = 0; c < Ctx->box.ncells[f]; ++c){

            //get cell index:
            i = Ctx->box.indices[f][c].i;
            j = Ctx->box.indices[f][c].j;
            k = Ctx->box.indices[f][c].k;

            //get pressure values from the stencil:
            //!Note Remeber to change the component for pressure if you are using it in a fluid code:
            for (z = 0; z < nVar; ++z){
            
            	for (oned_begin = -2; oned_begin < 3; ++oned_begin) 
            	{
                  u_x_loc[oned_begin+2] = u[k][j][i+oned_begin].comp[z];
                  u_y_loc[oned_begin+2] = u[k][j+oned_begin][i].comp[z];
                  u_z_loc[oned_begin+2] = u[k+oned_begin][j][i].comp[z];
               }
            
                u_xy_loc[0] = u[k][j][i].comp[z];
                u_xy_loc[1] = u[k][j+1][i+1].comp[z];
                u_xy_loc[2] = u[k][j-1][i+1].comp[z];
                u_xy_loc[3] = u[k][j+1][i-1].comp[z];
                u_xy_loc[4] = u[k][j-1][i-1].comp[z];
                
                u_yz_loc[0] = u[k][j][i].comp[z];
                u_yz_loc[1] = u[k+1][j+1][i].comp[z];
                u_yz_loc[2] = u[k-1][j+1][i].comp[z];
                u_yz_loc[3] = u[k+1][j-1][i].comp[z];
                u_yz_loc[4] = u[k-1][j-1][i].comp[z];
                
                u_zx_loc[0] = u[k][j][i].comp[z];
                u_zx_loc[1] = u[k+1][j][i+1].comp[z];
                u_zx_loc[2] = u[k+1][j][i-1].comp[z];
                u_zx_loc[3] = u[k-1][j][i+1].comp[z];
                u_zx_loc[4] = u[k-1][j][i-1].comp[z];
                
                u_xyz_loc[0] = u[k+1][j+1][i+1].comp[z]; u_xyz_loc[1] = u[k+1][j+1][i-1].comp[z];
                u_xyz_loc[2] = u[k+1][j-1][i+1].comp[z]; u_xyz_loc[3] = u[k-1][j+1][i+1].comp[z];
                u_xyz_loc[4] = u[k+1][j-1][i-1].comp[z]; u_xyz_loc[5] = u[k-1][j+1][i-1].comp[z];
                u_xyz_loc[6] = u[k-1][j-1][i+1].comp[z]; u_xyz_loc[7] = u[k-1][j-1][i-1].comp[z];
                
                //reconstruct polynomial:
                weno(u_x_loc, u_y_loc, u_z_loc, u_xy_loc, u_yz_loc, u_zx_loc, u_xyz_loc, coeffs[z]);

            }

            //loop over quadrature points on the cell face attached to Kirchhoff surface and compute p, pn.
            for (q = 0; q < N_gp2d; ++q){
                
                switch (f)
                {
                case 0:
                    p[iq] = evaluate_polynomial(-0.5, x_gp2d[q], y_gp2d[q], coeffs[0]);
                    evaluate_grad(coeffs[0], -0.5, x_gp2d[q], y_gp2d[q], h, &p_x, &p_y, &p_z);
                    pn[iq] = -p_x;
                    break;
                case 1:
                    p[iq] = evaluate_polynomial(+0.5, x_gp2d[q], y_gp2d[q], coeffs[0]); 
                    evaluate_grad(coeffs[0], +0.5, x_gp2d[q], y_gp2d[q], h, &p_x, &p_y, &p_z);
                    pn[iq] = p_x;
                    break;
                case 2:
                    p[iq] = evaluate_polynomial(x_gp2d[q], -0.5, y_gp2d[q], coeffs[0]);
                    evaluate_grad(coeffs[0], x_gp2d[q], -0.5, y_gp2d[q], h, &p_x, &p_y, &p_z);
                    pn[iq] = -p_y;
                    break;
                case 3:
                    p[iq] = evaluate_polynomial(x_gp2d[q], +0.5, y_gp2d[q], coeffs[0]);
                    evaluate_grad(coeffs[0], x_gp2d[q], +0.5, y_gp2d[q], h, &p_x, &p_y, &p_z);
                    pn[iq] = p_y;
                    break;
                case 4:
                    p[iq] = evaluate_polynomial(x_gp2d[q], y_gp2d[q], -0.5, coeffs[0]);
                    evaluate_grad(coeffs[0], x_gp2d[q], y_gp2d[q], -0.5, h, &p_x, &p_y, &p_z);
                    pn[iq] = -p_z;
                    break;
                case 5:
                    p[iq] = evaluate_polynomial(x_gp2d[q], y_gp2d[q], +0.5, coeffs[0]);
                    evaluate_grad(coeffs[0], x_gp2d[q], y_gp2d[q], +0.5, h, &p_x, &p_y, &p_z);
                    pn[iq] = p_z;
                    break;
                }

                iq++;
            }


        }

    }

    //restore array
    ierr = VecRestoreArray((Ctx->box.P), &p); CHKERRQ(ierr);
    ierr = VecRestoreArray((Ctx->box.Pn), &pn); CHKERRQ(ierr);
 

    //write P, Pn in a file:
    char filename[30];
    sprintf(filename, "kirchhoff/p-%05d.dat", step); 
    ierr = PetscPrintf(PETSC_COMM_WORLD, "writing p, pn on Kirchoff surface %s\n", filename); CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer); CHKERRQ(ierr);

    ierr = VecView(Ctx->box.P, viewer); CHKERRQ(ierr);
    ierr = VecView(Ctx->box.Pn, viewer); CHKERRQ(ierr);
 

    //restore array:
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,Ctx->localU,&u);CHKERRQ(ierr);

    return ierr;
}





//--------------------------------------------------------------------------------------
// Release all the memory related to Kirchhoff box surface
//--------------------------------------------------------------------------------------

PetscErrorCode  DestroyKirchoffBoxSurface(AppCtx* Ctx){

    PetscInt f;

    for ( f = 0; f < 6; ++f) free(Ctx->box.indices[f]);
    VecDestroy(&(Ctx->box.Xq));
    VecDestroy(&(Ctx->box.Nq));
    VecDestroy(&(Ctx->box.P));
    VecDestroy(&(Ctx->box.Pn));
 

    return 0;
}
