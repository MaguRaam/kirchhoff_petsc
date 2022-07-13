/*
 * rhs_function_prim.c
 *      Author: sunder
 */ 
#include "hype.h" 

//----------------------------------------------------------------------------
// Compute the value of RHS for each cell in the domain using 
// primitive variables
//----------------------------------------------------------------------------

PetscErrorCode RHSFunctionPrim(TS ts, PetscReal t, Vec U, Vec RHS, void* ctx) {
    
    PetscErrorCode ierr;           
    AppCtx *Ctx = (AppCtx*)ctx; 
    DM da;                         
    PetscInt c,i,j,k,l,f,q,xs,ys,zs,xm,ym,zm,xsg,ysg,zsg,xmg,ymg,zmg,oned_begin,oned_end,irhs;         
    
    Field   ***w;                   
    Field   ***rhs;                 
    PetscReal s, s_max = 0.0;  
    PetscReal w_x_loc[5], w_y_loc[5], w_z_loc[5], w_xy_loc[5], w_yz_loc[5], w_zx_loc[5], w_xyz_loc[8];
    PetscReal dt, value, grad_x, grad_y, grad_z; 
    PetscReal coeffs[nDOF], sol[nDOF][nVar];
    PetscReal r1_h = 1./(Ctx->h); 
    PetscReal nx, ny, nz, xloc, yloc, zloc; 
    PetscInt local_i, local_j, local_k;
    PetscBool PAD;
    
    PetscReal Vgp[N_gp3d][nVar], grad_Vgp_x[N_gp3d][nVar], grad_Vgp_y[N_gp3d][nVar], grad_Vgp_z[N_gp3d][nVar];
    PetscReal VNode[N_node][nVar];
    PetscReal V[nVar], gradV_x[nVar], gradV_y[nVar], gradV_z[nVar], BgradQ[nVar];
    PetscReal VL[nVar], VR[nVar];
    PetscReal Flux[nVar], Flux_conv[nVar], Fluc[nVar], NCP[nVar]; 
    
    ierr = TSGetDM(ts,&da);CHKERRQ(ierr);
    
    // Calculate primitive variables in each cell 
    
    ierr = ComputePrimitiveVariables(U, Ctx->W, da);CHKERRQ(ierr);
    
    // Scatter global->local to have access to the required ghost values 

    ierr = DMGlobalToLocalBegin(da, Ctx->W, INSERT_VALUES, Ctx->localU);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da, Ctx->W, INSERT_VALUES,   Ctx->localU);CHKERRQ(ierr);

    // Read the local solution to the array u  

    ierr = DMDAVecGetArrayRead(da, Ctx->localU, &w); CHKERRQ(ierr); 
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
                            w[k][j][i].comp[c] = w[k][j][irhs].comp[c];
                    }
                    
                    if (Ctx->left_boundary == reflective) { // Reflective Boundary 
                        irhs = oned_begin - 1 - i; 
                        for (c = 0; c < nVar; ++c)
                            w[k][j][i].comp[c] = w[k][j][irhs].comp[c];
                        
                        w[k][j][i].comp[1] = -w[k][j][i].comp[1]; // Reflect x-momentum component 
                    }
                }
                
                if (i >= Ctx->N_x) { // Right Boundary 
                    
                    if (Ctx->right_boundary == transmissive) { // Transmissive/Outflow Boundary 
                        irhs = oned_end; 
                        for (c = 0; c < nVar; ++c)
                            w[k][j][i].comp[c] = w[k][j][irhs].comp[c];
                    }
                    
                    if (Ctx->right_boundary == reflective) { // Reflective Boundary
                        irhs = 2*oned_end - i + 1; 
                        for (c = 0; c < nVar; ++c)
                            w[k][j][i].comp[c] = w[k][j][irhs].comp[c];
                        
                        w[k][j][i].comp[1] = -w[k][j][i].comp[1]; // Reflect x-momentum component 
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
                            w[k][j][i].comp[c] = w[k][irhs][i].comp[c]; 
                    }
                    
                    if (Ctx->bottom_boundary == reflective) { // Reflective Boundary
                        irhs = oned_begin - 1 - j;
                        for (c = 0; c < nVar; ++c)
                            w[k][j][i].comp[c] = w[k][irhs][i].comp[c]; 
                        
                        w[k][j][i].comp[2] = -w[k][j][i].comp[2]; // Reflect y-momentum component
                    }
                }
                
                if (j >= Ctx->N_y) { // Top boundary 
                    
                    if (Ctx->top_boundary == transmissive) { // Transmissive/Outflow Boundary
                        irhs = oned_end;
                        for (c = 0; c < nVar; ++c)
                            w[k][j][i].comp[c] = w[k][irhs][i].comp[c];
                    }
                    
                    if (Ctx->top_boundary == reflective) { // Reflective Boundary
                        irhs = 2*oned_end - j + 1;
                        for (c = 0; c < nVar; ++c)
                            w[k][j][i].comp[c] = w[k][irhs][i].comp[c];
                        
                        w[k][j][i].comp[2] = -w[k][j][i].comp[2]; // Reflect y-momentum component
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
                            w[k][j][i].comp[c] = w[irhs][j][i].comp[c]; 
                    }
                    
                    if (Ctx->back_boundary == reflective) { // Reflective Boundary
                        irhs = oned_begin - 1 - k;
                        for (c = 0; c < nVar; ++c)
                            w[k][j][i].comp[c] = w[irhs][j][i].comp[c]; 
                        
                        w[k][j][i].comp[3] = -w[k][j][i].comp[3]; // Reflect z-momentum component
                    }
                }
                
                if (k >= Ctx->N_z) { // Front boundary 
                    
                    if (Ctx->front_boundary == transmissive) { // Transmissive/Outflow Boundary
                        irhs = oned_end;
                        for (c = 0; c < nVar; ++c)
                            w[k][j][i].comp[c] = w[irhs][j][i].comp[c];
                    }
                    
                    if (Ctx->front_boundary == reflective) { // Reflective Boundary
                        irhs = 2*oned_end - k + 1;
                        for (c = 0; c < nVar; ++c)
                            w[k][j][i].comp[c] = w[irhs][j][i].comp[c];
                        
                        w[k][j][i].comp[3] = -w[k][j][i].comp[3]; // Reflect z-momentum component
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
                        w_x_loc[oned_begin+2] = w[k][j][i+oned_begin].comp[c];
                        w_y_loc[oned_begin+2] = w[k][j+oned_begin][i].comp[c];
                        w_z_loc[oned_begin+2] = w[k+oned_begin][j][i].comp[c];
                    }
                    
                    w_xy_loc[0] = w[k][j][i].comp[c];
                    w_xy_loc[1] = w[k][j+1][i+1].comp[c];
                    w_xy_loc[2] = w[k][j-1][i+1].comp[c];
                    w_xy_loc[3] = w[k][j+1][i-1].comp[c];
                    w_xy_loc[4] = w[k][j-1][i-1].comp[c];
                    
                    w_yz_loc[0] = w[k][j][i].comp[c];
                    w_yz_loc[1] = w[k+1][j+1][i].comp[c];
                    w_yz_loc[2] = w[k-1][j+1][i].comp[c];
                    w_yz_loc[3] = w[k+1][j-1][i].comp[c];
                    w_yz_loc[4] = w[k-1][j-1][i].comp[c];
                    
                    w_zx_loc[0] = w[k][j][i].comp[c];
                    w_zx_loc[1] = w[k+1][j][i+1].comp[c];
                    w_zx_loc[2] = w[k+1][j][i-1].comp[c];
                    w_zx_loc[3] = w[k-1][j][i+1].comp[c];
                    w_zx_loc[4] = w[k-1][j][i-1].comp[c];
                    
                    w_xyz_loc[0] = w[k+1][j+1][i+1].comp[c]; w_xyz_loc[1] = w[k+1][j+1][i-1].comp[c];
                    w_xyz_loc[2] = w[k+1][j-1][i+1].comp[c]; w_xyz_loc[3] = w[k-1][j+1][i+1].comp[c];
                    w_xyz_loc[4] = w[k+1][j-1][i-1].comp[c]; w_xyz_loc[5] = w[k-1][j+1][i-1].comp[c];
                    w_xyz_loc[6] = w[k-1][j-1][i+1].comp[c]; w_xyz_loc[7] = w[k-1][j-1][i-1].comp[c];
                    
                    weno(w_x_loc, w_y_loc, w_z_loc, w_xy_loc, w_yz_loc, w_zx_loc, w_xyz_loc, coeffs);
                    
                    // Store the coefficients 
                
                    for (l = 0; l < nDOF; ++l) {
                        sol[l][c] = coeffs[l]; 
                    }
                    
                    // Evaluate solution at various nodes to check physical admissibility of the solution  
                
                    for (q = 0; q < N_node; ++q) {
                        VNode[q][c] = 0.0;
                        
                        for (l = 0; l < nDOF; ++l)
                            VNode[q][c] += coeffs[l]*get_element_2d(Ctx->phiNode,q,l);
                    }
                }
                
                // Check physical admissibility of the solution  and reduce to TVD if necessary 
                
                for (q = 0; q < N_node; ++q) {
                    
                    for (c = 0; c < nVar; ++c) {
                        V[c] = VNode[q][c];
                    }
                    
                    PAD = PDECheckPADPrim(V);
                    
                    if (PAD == PETSC_FALSE)
                        break; 
                }
                
                if (PAD == PETSC_FALSE) {
                
                    for (c = 0 ; c < nVar; ++c) {
                        sol[0][c] = w[k][j][i].comp[c];
                        sol[1][c] = minmod(w[k][j][i+1].comp[c]-w[k][j][i].comp[c],w[k][j][i].comp[c]-w[k][j][i-1].comp[c]);
                        sol[2][c] = minmod(w[k][j+1][i].comp[c]-w[k][j][i].comp[c],w[k][j][i].comp[c]-w[k][j-1][i].comp[c]);
                        sol[3][c] = minmod(w[k+1][j][i].comp[c]-w[k][j][i].comp[c],w[k][j][i].comp[c]-w[k-1][j][i].comp[c]);
                        
                        for (l = 4; l < nDOF; ++l) {
                            sol[l][c] = 0.0; 
                        }
                    }
                    
                }
                
                // Find the values of primitive variables at face quadrature points 
                
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
                        
                        Vgp[q][c] = value;
                        grad_Vgp_x[q][c] = r1_h*grad_x; 
                        grad_Vgp_y[q][c] = r1_h*grad_y;
                        grad_Vgp_z[q][c] = r1_h*grad_z;
                    }
                
                }
                
                // Add smooth part of the non-conservative product 
            
                if (i >= xs && i < xs+xm && j >= ys && j < ys+ym && k >= zs && k < zs+zm ) {
                    
                    for (c = 0 ; c < nVar; ++c) 
                        rhs[k][j][i].comp[c] = 0.0; 
                    
                    for (q = 0; q < N_gp3d; ++q) {
                    
                        for (c = 0 ; c < nVar; ++c) {
                            V[c] = Vgp[q][c];
                            gradV_x[c] = grad_Vgp_x[q][c];
                            gradV_y[c] = grad_Vgp_y[q][c];
                            gradV_z[c] = grad_Vgp_z[q][c];
                        }
                        
                        PDENCPPrim(V, gradV_x, gradV_y, gradV_z, BgradQ);
                        
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
                        
                        VL[c] = get_element_6d(Ctx->u_bnd, local_k, local_j, local_i-1, c, 1, q);
                        VR[c] = get_element_6d(Ctx->u_bnd, local_k, local_j, local_i,   c, 0, q);
                    }
                
                    s = LLFRiemannSolverPrim(VL, VR, nx, ny, nz, xloc, yloc, zloc, Flux_conv, Fluc); if (s>s_max) s_max = s; 
                 
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
                        
                        VL[c] = get_element_6d(Ctx->u_bnd, local_k, local_j-1, local_i, c, 3, q);
                        VR[c] = get_element_6d(Ctx->u_bnd, local_k, local_j, local_i,   c, 2, q);
                    }
                
                    s = LLFRiemannSolverPrim(VL, VR, nx, ny, nz, xloc, yloc, zloc, Flux_conv, Fluc); if (s>s_max) s_max = s; 
                 
                
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
                        
                        VL[c] = get_element_6d(Ctx->u_bnd, local_k-1, local_j, local_i, c, 5, q);
                        VR[c] = get_element_6d(Ctx->u_bnd, local_k, local_j, local_i,   c, 4, q);
                    }
                
                    s = LLFRiemannSolverPrim(VL, VR, nx, ny, nz, xloc, yloc, zloc, Flux_conv, Fluc); if (s>s_max) s_max = s; 
                 
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



    
    ierr = DMDAVecRestoreArray(da,Ctx->localU,&w);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,RHS,&rhs);CHKERRQ(ierr);

    dt = (Ctx->CFL*Ctx->h)/(3.0*s_max); // 3.0 is for three dimensions

    ierr = MPI_Allreduce(&dt, &Ctx->dt, 1, MPIU_REAL,MPIU_MIN, PetscObjectComm((PetscObject)da)); CHKERRQ(ierr);

    return ierr; 
}
