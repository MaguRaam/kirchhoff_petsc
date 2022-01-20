/*
 * pde_prim.c
 *      Author: sunder
 */ 

#include "hype.h"

/* Definitions of various PDE functions */

//----------------------------------------------------------------------------
// Find the conservative flux components F in the given normal direction
//----------------------------------------------------------------------------

PetscReal PDEFluxPrim (const PetscReal* V, 
                       const PetscReal nx, const PetscReal ny, const PetscReal nz, 
                       const PetscReal x,  const PetscReal y,  const PetscReal z, 
                       PetscReal* F) {


    PetscReal phi = V[5];
    PetscReal g = 1.0 + (g1-1.0)*(g2-1.0)/((1.0-phi)*(g1 -1.0) + phi*(g2-1.0));
    PetscReal P_inf = ((g -1.0)/g)*( g1*p1*phi/(g1 - 1.0) + g2*p2*(1.0 - phi)/(g2 - 1.0) );
    PetscReal rho = V[0];
    PetscReal u = V[1];
    PetscReal v = V[2];
    PetscReal w = V[3];
    PetscReal un = u*nx + v*ny + w*nz;
    PetscReal rhoun = rho*un;  
    PetscReal p = V[4]; 
    PetscReal E = (p + g*P_inf)/(g - 1.0) + 0.5*rho*(u*u + v*v + w*w);
    
    // Check if the input state is physically admissible 

    if (rho < rho_floor) {
        printf("Negative density = %f\n", rho);
        printf("At x = %f, y = %f, z = %f\n", x, y, z);  
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    if ((p + P_inf)  < prs_floor) {
        printf("Negative pressure, p = %f\n", p + P_inf);
        printf("At x = %f, y = %f, z = %f\n", x, y, z);  
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    // Now find the fluxes 
    
    F[0] = rhoun;
    F[1] = rhoun*u + nx*p;
    F[2] = rhoun*v + ny*p;
    F[3] = rhoun*w + nz*p;
    F[4] = (E + p)*un;
    F[5] = 0.0; 
    
    // Also obtain the maximum eigen value

    PetscReal s_max = PetscAbsReal(un) + PetscSqrtReal(g*(p+P_inf)/rho);

    return s_max;
}

//----------------------------------------------------------------------------
// Non-conservative product BgradQ
//----------------------------------------------------------------------------

void PDENCPPrim(const PetscReal *V, const PetscReal *grad_V_x, const PetscReal *grad_V_y, const PetscReal *grad_V_z, PetscReal *BgradQ) {
    
    PetscReal phi_x = grad_V_x[5]; 
    PetscReal phi_y = grad_V_y[5];
    PetscReal phi_z = grad_V_z[5];

    BgradQ[0] = 0.0;
    BgradQ[1] = 0.0;
    BgradQ[2] = 0.0; 
    BgradQ[3] = 0.0; 
    BgradQ[4] = 0.0; 
    BgradQ[5] = (V[1]*phi_x + V[2]*phi_y + V[3]*phi_z);
}

//----------------------------------------------------------------------------
// Check physical admissibility of the input state 
//----------------------------------------------------------------------------

PetscBool PDECheckPADPrim(const PetscReal *V) {
    
    PetscBool PAD = PETSC_TRUE;
    
    PetscReal rho = V[0];
    PetscReal phi = V[5]; 
    PetscReal g = 1.0 + (g1-1.0)*(g2-1.0)/((1.0-phi)*(g1 -1.0) + phi*(g2-1.0));
    PetscReal P_inf = ((g -1.0)/g)*( g1*p1*phi/(g1 - 1.0) + g2*p2*(1.0 - phi)/(g2 - 1.0) );
    PetscReal p = V[4];
    
    // Check if the input state is physically admissible 

    if (rho < rho_floor) {    
        PAD = PETSC_FALSE;
    }

    if (p + P_inf  < prs_floor) {
        PAD = PETSC_FALSE;
    }
    
    if (phi < 0.0 || phi > 1.0) {
        PAD = PETSC_FALSE;
    }
    
    return PAD; 
}

//----------------------------------------------------------------------------
// Rusanov/LLF Riemann Solver 
//----------------------------------------------------------------------------

PetscReal LLFRiemannSolverPrim(const PetscReal *VL, const PetscReal *VR,
                           PetscReal nx, PetscReal ny, PetscReal nz,
                           PetscReal x,  PetscReal y, PetscReal z, 
                           PetscReal *F, PetscReal *D) {
        
    PetscInt c;
    PetscReal FL[nVar], FR[nVar], grad_V_x[nVar], grad_V_y[nVar], grad_V_z[nVar];
    PetscReal QL[nVar], QR[nVar];
    PetscReal dv; 
    PetscReal V_av[nVar];
    
    PDEPrim2Cons(VL,QL); PDEPrim2Cons(VR,QR);

    PetscReal s_max_l = PDEFluxPrim(VL, nx, ny, nz, x, y, z, FL); 
    PetscReal s_max_r = PDEFluxPrim(VR, nx, ny, nz, x, y, z, FR);

    PetscReal s_max = PetscMax(s_max_l, s_max_r);

    for (c = 0; c < nVar; ++c) {
        V_av[c] = 0.5*(VR[c]+VL[c]); 
        dv = VR[c]-VL[c];
        grad_V_x[c] = nx*dv;
        grad_V_y[c] = ny*dv;
        grad_V_z[c] = nz*dv;
        F[c] = 0.5*( FR[c] + FL[c] - s_max*(QR[c]-QL[c]) );
    }
    
    PDENCPPrim(V_av, grad_V_x, grad_V_y, grad_V_z, D);
    
    return s_max; 
}
