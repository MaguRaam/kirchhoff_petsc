/*
 * rhs_function.c
 *      Author: sunder
 */ 
#include "hype.h" 

//----------------------------------------------------------------------------
// Compute the value of RHS for each cell in the domain using 
// conserved variables
//----------------------------------------------------------------------------

PetscErrorCode RHSFunction(TS ts, PetscReal t, Vec U, Vec RHS, void* ctx) {
    
    
    PetscErrorCode ierr;           
    AppCtx *Ctx = (AppCtx*)ctx; 
    DM da;                         
    PetscInt c,i,j,k,l,f,q,xs,ys,zs,xm,ym,zm,xsg,ysg,zsg,xmg,ymg,zmg,oned_begin,oned_end,irhs;         
    
    Field   ***u;                   
    Field   ***rhs;                 
    PetscReal s, s_max = 0.0;  
    PetscReal u_x_loc[5], u_y_loc[5], u_z_loc[5], u_xy_loc[5], u_yz_loc[5], u_zx_loc[5], u_xyz_loc[8];
    PetscReal dt, value, grad_x, grad_y, grad_z; 
    PetscReal coeffs[nDOF], sol[nDOF][nVar];
    PetscReal r1_h = 1./(Ctx->h); 
    PetscReal nx, ny, nz, xloc, yloc, zloc; 
    PetscInt local_i, local_j, local_k;
    PetscBool PAD;
    
    PetscReal Qgp[N_gp3d][nVar], grad_Qgp_x[N_gp3d][nVar], grad_Qgp_y[N_gp3d][nVar], grad_Qgp_z[N_gp3d][nVar];
    PetscReal QNode[N_node][nVar];
    PetscReal Q[nVar], gradQ_x[nVar], gradQ_y[nVar], gradQ_z[nVar], BgradQ[nVar];
    PetscReal QL[nVar], QR[nVar];
    PetscReal Flux[nVar], Flux_conv[nVar], Fluc[nVar], NCP[nVar]; 
    
    ierr = TSGetDM(ts,&da);CHKERRQ(ierr);
    
    // Scatter global->local to have access to the required ghost values 

    ierr = DMGlobalToLocalBegin(da, U, INSERT_VALUES, Ctx->localU);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da, U, INSERT_VALUES,   Ctx->localU);CHKERRQ(ierr);

    // Read the local solution to the array u  

    ierr = DMDAVecGetArrayRead(da, Ctx->localU, &u); CHKERRQ(ierr); 
    ierr = DMDAVecGetArray(da, RHS,&rhs);CHKERRQ(ierr);

    ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);
    ierr = DMDAGetGhostCorners(da, &xsg, &ysg, &zsg, &xmg, &ymg, &zmg);
    
    // 1) Apply boundary conditions (Periodic boundary conditions are taken care of by PETsc)

    oned_begin = 0; oned_end = Ctx->N_x-1;  // First in the x-direction 

    for (k = zsg; k < zsg+zmg; ++k) {
        for (j = ysg; j < ysg+ymg; ++j) {
            for (i = xsg; i < xsg+xmg; ++i) {
                
                if (i < 0)  { // Left boundary 
                    
                    if (Ctx->left_boundary == transmissive) { // Transmissive/Outflow Boundary
                        irhs = oned_begin; 
                        for (c = 0; c < nVar; ++c)
                            u[k][j][i].comp[c] = u[k][j][irhs].comp[c];
                    }
                    
                    if (Ctx->left_boundary == reflective) { // Reflective Boundary 
                        irhs = oned_begin - 1 - i; 
                        for (c = 0; c < nVar; ++c)
                            u[k][j][i].comp[c] = u[k][j][irhs].comp[c];
                        
                        u[k][j][i].comp[1] = -u[k][j][i].comp[1]; // Reflect x-momentum component 
                    }
                }
                
                if (i >= Ctx->N_x) { // Right Boundary 
                    
                    if (Ctx->right_boundary == transmissive) { // Transmissive/Outflow Boundary 
                        irhs = oned_end; 
                        for (c = 0; c < nVar; ++c)
                            u[k][j][i].comp[c] = u[k][j][irhs].comp[c];
                    }
                    
                    if (Ctx->right_boundary == reflective) { // Reflective Boundary
                        irhs = 2*oned_end - i + 1; 
                        for (c = 0; c < nVar; ++c)
                            u[k][j][i].comp[c] = u[k][j][irhs].comp[c];
                        
                        u[k][j][i].comp[1] = -u[k][j][i].comp[1]; // Reflect x-momentum component 
                    }
                }
            }
        }
    }

    oned_begin = 0; oned_end = Ctx->N_y-1; // Next in the y-direction 
	
    for (k = zsg; k < zsg+zmg; ++k) {
        for (j = ysg; j < ysg+ymg; ++j) {
            for (i = xsg; i < xsg+xmg; ++i) {
                
                if (j < 0) { // Bottom boundary 
                    
                    if (Ctx->bottom_boundary == transmissive) { // Transmissive/Outflow Boundary
                        irhs = oned_begin;
                        for (c = 0; c < nVar; ++c)
                            u[k][j][i].comp[c] = u[k][irhs][i].comp[c]; 
                    }
                    
                    if (Ctx->bottom_boundary == reflective) { // Reflective Boundary
                        irhs = oned_begin - 1 - j;
                        for (c = 0; c < nVar; ++c)
                            u[k][j][i].comp[c] = u[k][irhs][i].comp[c]; 
                        
                        u[k][j][i].comp[2] = -u[k][j][i].comp[2]; // Reflect y-momentum component
                    }
                }
                
                if (j >= Ctx->N_y) { // Top boundary 
                    
                    if (Ctx->top_boundary == transmissive) { // Transmissive/Outflow Boundary
                        irhs = oned_end;
                        for (c = 0; c < nVar; ++c)
                            u[k][j][i].comp[c] = u[k][irhs][i].comp[c];
                    }
                    
                    if (Ctx->top_boundary == reflective) { // Reflective Boundary
                        irhs = 2*oned_end - j + 1;
                        for (c = 0; c < nVar; ++c)
                            u[k][j][i].comp[c] = u[k][irhs][i].comp[c];
                        
                        u[k][j][i].comp[2] = -u[k][j][i].comp[2]; // Reflect y-momentum component
                    }
                }
            }
        }
    }
	
    oned_begin = 0; oned_end = Ctx->N_z-1; // Finally in the z-direction 
	
    for (k = zsg; k < zsg+zmg; ++k) {
        for (j = ysg; j < ysg+ymg; ++j) {
            for (i = xsg; i < xsg+xmg; ++i) {
                
                if (k < 0) { // Back boundary 
                    
                    if (Ctx->back_boundary == transmissive) { // Transmissive/Outflow Boundary
                        irhs = oned_begin;
                        for (c = 0; c < nVar; ++c)
                            u[k][j][i].comp[c] = u[irhs][j][i].comp[c]; 
                    }
                    
                    if (Ctx->back_boundary == reflective) { // Reflective Boundary
                        irhs = oned_begin - 1 - k;
                        for (c = 0; c < nVar; ++c)
                            u[k][j][i].comp[c] = u[irhs][j][i].comp[c]; 
                        
                        u[k][j][i].comp[3] = -u[k][j][i].comp[3]; // Reflect z-momentum component
                    }
                }
                
                if (k >= Ctx->N_z) { // Front boundary 
                    
                    if (Ctx->front_boundary == transmissive) { // Transmissive/Outflow Boundary
                        irhs = oned_end;
                        for (c = 0; c < nVar; ++c)
                            u[k][j][i].comp[c] = u[irhs][j][i].comp[c];
                    }
                    
                    if (Ctx->front_boundary == reflective) { // Reflective Boundary
                        irhs = 2*oned_end - k + 1;
                        for (c = 0; c < nVar; ++c)
                            u[k][j][i].comp[c] = u[irhs][j][i].comp[c];
                        
                        u[k][j][i].comp[3] = -u[k][j][i].comp[3]; // Reflect z-momentum component
                    }
                }
            }
        }
    }
    
    // 2) Find the boundary extrapolated values of the conserved variables in all the cells 
    
    for ( k = zs-1; k < zs+zm+1; ++k ) {
        for ( j = ys-1; j < ys+ym+1; ++j ) { 
            for ( i = xs-1; i < xs+xm+1; ++i ) {
            
                
                local_k = k - (zs-1);
                local_j = j - (ys-1); 
                local_i = i - (xs-1); 
                PAD = PETSC_TRUE;
                
                for (c = 0; c < nVar; ++c) {
                    
                    for (oned_begin = -2; oned_begin < 3; ++oned_begin) {
                        u_x_loc[oned_begin+2] = u[k][j][i+oned_begin].comp[c];
                        u_y_loc[oned_begin+2] = u[k][j+oned_begin][i].comp[c];
                        u_z_loc[oned_begin+2] = u[k+oned_begin][j][i].comp[c];
                    }
                    
                    u_xy_loc[0] = u[k][j][i].comp[c];
                    u_xy_loc[1] = u[k][j+1][i+1].comp[c];
                    u_xy_loc[2] = u[k][j-1][i+1].comp[c];
                    u_xy_loc[3] = u[k][j+1][i-1].comp[c];
                    u_xy_loc[4] = u[k][j-1][i-1].comp[c];
                    
                    u_yz_loc[0] = u[k][j][i].comp[c];
                    u_yz_loc[1] = u[k+1][j+1][i].comp[c];
                    u_yz_loc[2] = u[k-1][j+1][i].comp[c];
                    u_yz_loc[3] = u[k+1][j-1][i].comp[c];
                    u_yz_loc[4] = u[k-1][j-1][i].comp[c];
                    
                    u_zx_loc[0] = u[k][j][i].comp[c];
                    u_zx_loc[1] = u[k+1][j][i+1].comp[c];
                    u_zx_loc[2] = u[k+1][j][i-1].comp[c];
                    u_zx_loc[3] = u[k-1][j][i+1].comp[c];
                    u_zx_loc[4] = u[k-1][j][i-1].comp[c];
                    
                    u_xyz_loc[0] = u[k+1][j+1][i+1].comp[c]; u_xyz_loc[1] = u[k+1][j+1][i-1].comp[c];
                    u_xyz_loc[2] = u[k+1][j-1][i+1].comp[c]; u_xyz_loc[3] = u[k-1][j+1][i+1].comp[c];
                    u_xyz_loc[4] = u[k+1][j-1][i-1].comp[c]; u_xyz_loc[5] = u[k-1][j+1][i-1].comp[c];
                    u_xyz_loc[6] = u[k-1][j-1][i+1].comp[c]; u_xyz_loc[7] = u[k-1][j-1][i-1].comp[c];
                    
                    weno(u_x_loc, u_y_loc, u_z_loc, u_xy_loc, u_yz_loc, u_zx_loc, u_xyz_loc, coeffs);
                    
                    // Store the coefficients 
                
                    for (l = 0; l < nDOF; ++l) {
                        sol[l][c] = coeffs[l]; 
                    }
                    
                    // Evaluate solution at various nodes to check physical admissibility of the solution  
                
                    for (q = 0; q < N_node; ++q) {
                        QNode[q][c] = 0.0;
                        
                        for (l = 0; l < nDOF; ++l)
                            QNode[q][c] += coeffs[l]*get_element_2d(Ctx->phiNode,q,l);
                    }
                }
                
                // Check physical admissibility of the solution  and reduce to TVD if necessary 
                
                for (q = 0; q < N_node; ++q) {
                    
                    for (c = 0; c < nVar; ++c) {
                        Q[c] = QNode[q][c];
                    }
                    
                    PAD = PDECheckPAD(Q);
                    
                    if (PAD == PETSC_FALSE)
                        break; 
                }
                
                if (PAD == PETSC_FALSE) {
                
                    for (c = 0 ; c < nVar; ++c) {
                        sol[0][c] = u[k][j][i].comp[c];
                        sol[1][c] = minmod(u[k][j][i+1].comp[c]-u[k][j][i].comp[c],u[k][j][i].comp[c]-u[k][j][i-1].comp[c]);
                        sol[2][c] = minmod(u[k][j+1][i].comp[c]-u[k][j][i].comp[c],u[k][j][i].comp[c]-u[k][j-1][i].comp[c]);
                        sol[3][c] = minmod(u[k+1][j][i].comp[c]-u[k][j][i].comp[c],u[k][j][i].comp[c]-u[k-1][j][i].comp[c]);
                        
                        for (l = 4; l < nDOF; ++l) {
                            sol[l][c] = 0.0; 
                        }
                    }
                }
                
                // Find the values of conserved variables at face quadrature points 
                
                for (c = 0 ; c < nVar; ++c) {
                    
                    // Get coefficients 
                    
                    for (f = 0; f < 6; ++f) {
                        for (q = 0; q < N_gp2d; ++q) {
                            
                            value = 0.0; 
                            
                            for (l = 0; l < nDOF; ++l)
                                value += sol[l][c]*get_element_3d(Ctx->phiFace,f,q,l);
                            
                            set_element_6d(Ctx->u_bnd, local_k, local_j, local_i, c, f, q, value);
                        }
                    }
                    
                    
                    
                    for (q = 0; q < N_gp3d; ++q) {
                    
                        value = 0.0; grad_x = 0.0; grad_y = 0.0; grad_z = 0.0; 
                        
                        for (l = 0; l < nDOF; ++l) {
                            value  += sol[l][c]*get_element_2d(Ctx->phiVol,q,l);
                            grad_x += sol[l][c]*get_element_2d(Ctx->gradphiVol_x,q,l);
                            grad_y += sol[l][c]*get_element_2d(Ctx->gradphiVol_y,q,l);
                            grad_z += sol[l][c]*get_element_2d(Ctx->gradphiVol_z,q,l);
                        }
                        
                        Qgp[q][c] = value;
                        grad_Qgp_x[q][c] = r1_h*grad_x; 
                        grad_Qgp_y[q][c] = r1_h*grad_y;
                        grad_Qgp_z[q][c] = r1_h*grad_z;
                    }
                    
                
                }
                
                // Add smooth part of the non-conservative product 
            
                if (i >= xs && i < xs+xm && j >= ys && j < ys+ym && k >= zs && k < zs+zm ) {
                    
                    for (c = 0 ; c < nVar; ++c) 
                        rhs[k][j][i].comp[c] = 0.0; 
                    
                    for (q = 0; q < N_gp3d; ++q) {
                    
                        for (c = 0 ; c < nVar; ++c) {
                            Q[c] = Qgp[q][c];
                            gradQ_x[c] = grad_Qgp_x[q][c];
                            gradQ_y[c] = grad_Qgp_y[q][c];
                            gradQ_z[c] = grad_Qgp_z[q][c];
                        }
                        
                        PDENCP(Q, gradQ_x, gradQ_y, gradQ_z, BgradQ);
                        
                        for (c = 0 ; c < nVar; ++c) 
                            rhs[k][j][i].comp[c] += -w_gp3d[q]*BgradQ[c];
                    }
                }
                
            }
        }
    }

    // 3) Find the upwind fluxes 
    
    // in x-direction 
    
    nx = 1.0; ny = 0.0; nz = 0.0; 
    
    for (k=zs; k<zs+zm; ++k) {
        for (j=ys; j<ys+ym; ++j) {
            for (i=xs; i <xs+xm+1; ++i) {
                
                xloc = Ctx->x_min+(PetscReal)i*(Ctx->h); yloc = Ctx->y_min+(PetscReal)j*(Ctx->h); zloc = Ctx->z_min+(PetscReal)k*(Ctx->h);
                        
                local_k = k - (zs-1);
                local_j = j - (ys-1); 
                local_i = i - (xs-1); 
            
                for (c = 0; c < nVar; ++c) {
                    Flux[c] = 0.0; 
                    NCP[c] = 0.0; 
                }
 
                for (q = 0; q < N_gp2d; ++q) {
                
                    for (c = 0; c < nVar; ++c) {
                        
                        QL[c] = get_element_6d(Ctx->u_bnd, local_k, local_j, local_i-1, c, 1, q);
                        QR[c] = get_element_6d(Ctx->u_bnd, local_k, local_j, local_i,   c, 0, q);
                    }
                
                    s = LLFRiemannSolver(QL, QR, nx, ny, nz, xloc, yloc, zloc, Flux_conv, Fluc); if (s>s_max) s_max = s; 
                 
                    for (c = 0; c < nVar; ++c) {
                        Flux[c] += w_gp2d[q]*(Flux_conv[c]);
                        NCP[c] += w_gp2d[q]*(Fluc[c]);
                    }
                }
        
                for (c = 0; c < nVar; ++c) {
                    set_element_4d(Ctx->Fx, k-zs, j-ys, i-xs, c, Flux[c]); 
                    set_element_4d(Ctx->Dx, k-zs, j-ys, i-xs, c, NCP[c]);
                }
            }
        }
    }
    
    // in y-direction 
    
    nx = 0.0; ny = 1.0; nz = 0.0; 
    
    for (k=zs; k<zs+zm; ++k) {
        for (j=ys; j<ys+ym+1; ++j) {
            for (i=xs; i <xs+xm; ++i) {
                
                xloc = Ctx->x_min+(PetscReal)i*(Ctx->h); yloc = Ctx->y_min+(PetscReal)j*(Ctx->h); zloc = Ctx->z_min+(PetscReal)k*(Ctx->h);
                        
                local_k = k - (zs-1);
                local_j = j - (ys-1); 
                local_i = i - (xs-1); 
            
                for (c = 0; c < nVar; ++c) {
                    Flux[c] = 0.0; 
                    NCP[c] = 0.0; 
                }
 
                for (q = 0; q < N_gp2d; ++q) {
                
                    for (c = 0; c < nVar; ++c) {
                        
                        QL[c] = get_element_6d(Ctx->u_bnd, local_k, local_j-1, local_i, c, 3, q);
                        QR[c] = get_element_6d(Ctx->u_bnd, local_k, local_j, local_i,   c, 2, q);
                    }
                
                    s = LLFRiemannSolver(QL, QR, nx, ny, nz, xloc, yloc, zloc, Flux_conv, Fluc); if (s>s_max) s_max = s; 
                 
                
                    for (c = 0; c < nVar; ++c) {
                        Flux[c] += w_gp2d[q]*(Flux_conv[c]);
                        NCP[c] += w_gp2d[q]*(Fluc[c]);
                    }
                }
        
                for (c = 0; c < nVar; ++c) {
                    set_element_4d(Ctx->Fy, k-zs, j-ys, i-xs, c, Flux[c]);
                    set_element_4d(Ctx->Dy, k-zs, j-ys, i-xs, c, NCP[c]);
                }
            }
        }
    }
    
    // in z-direction 
    
    nx = 0.0; ny = 0.0; nz = 1.0; 
    
    for (k=zs; k<zs+zm+1; ++k) {
        for (j=ys; j<ys+ym; ++j) {
            for (i=xs; i <xs+xm; ++i) {
                
                xloc = Ctx->x_min+(PetscReal)i*(Ctx->h); yloc = Ctx->y_min+(PetscReal)j*(Ctx->h); zloc = Ctx->z_min+(PetscReal)k*(Ctx->h);
                        
                local_k = k - (zs-1);
                local_j = j - (ys-1); 
                local_i = i - (xs-1); 
            
                for (c = 0; c < nVar; ++c) {
                    Flux[c] = 0.0; 
                    NCP[c] = 0.0; 
                }
 
                for (q = 0; q < N_gp2d; ++q) {
                
                    for (c = 0; c < nVar; ++c) {
                        
                        QL[c] = get_element_6d(Ctx->u_bnd, local_k-1, local_j, local_i, c, 5, q);
                        QR[c] = get_element_6d(Ctx->u_bnd, local_k, local_j, local_i,   c, 4, q);
                    }
                
                    s = LLFRiemannSolver(QL, QR, nx, ny, nz, xloc, yloc, zloc, Flux_conv, Fluc); if (s>s_max) s_max = s; 
                 
                    for (c = 0; c < nVar; ++c) {
                        Flux[c] += w_gp2d[q]*(Flux_conv[c]);
                        NCP[c] += w_gp2d[q]*(Fluc[c]);
                    }
                }
        
                for (c = 0; c < nVar; ++c) {
                    set_element_4d(Ctx->Fz, k-zs, j-ys, i-xs, c, Flux[c]);
                    set_element_4d(Ctx->Dz, k-zs, j-ys, i-xs, c, NCP[c]);
                }
            }
        }
    }
    
    // 4) Find the rhs in each cell 


    for(k=zs; k<zs+zm; ++k) {
        for (j=ys; j<ys+ym; ++j) {
            for (i=xs; i<xs+xm; ++i) {

                for (c = 0 ; c < nVar; ++c) {
        
                    rhs[k][j][i].comp[c] += -r1_h*( get_element_4d(Ctx->Fx, k-zs, j-ys, i+1-xs, c) - get_element_4d(Ctx->Fx, k-zs, j-ys, i-xs, c) + 
                                              0.5*( get_element_4d(Ctx->Dx, k-zs, j-ys, i+1-xs, c) + get_element_4d(Ctx->Dx, k-zs, j-ys, i-xs, c) ) )
                                            -r1_h*( get_element_4d(Ctx->Fy, k-zs, j+1-ys, i-xs, c) - get_element_4d(Ctx->Fy, k-zs, j-ys, i-xs, c) + 
                                              0.5*( get_element_4d(Ctx->Dy, k-zs, j+1-ys, i-xs, c) + get_element_4d(Ctx->Dy, k-zs, j-ys, i-xs, c) ) )                          
                                            -r1_h*( get_element_4d(Ctx->Fz, k+1-zs, j-ys, i-xs, c) - get_element_4d(Ctx->Fz, k-zs, j-ys, i-xs, c) + 
                                              0.5*( get_element_4d(Ctx->Dz, k+1-zs, j-ys, i-xs, c) + get_element_4d(Ctx->Dz, k-zs, j-ys, i-xs, c) ) );
                }


            }
        }
    }
    
    ierr = DMDAVecRestoreArray(da,Ctx->localU,&u);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,RHS,&rhs);CHKERRQ(ierr);

    dt = (Ctx->CFL*Ctx->h)/(3.0*s_max); // 3.0 is for three dimensions

    ierr = MPI_Allreduce(&dt, &Ctx->dt, 1, MPIU_REAL,MPIU_MIN, PetscObjectComm((PetscObject)da)); CHKERRQ(ierr);

    return ierr; 
}
