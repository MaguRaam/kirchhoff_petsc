/*
 * pde.c
 *      Author: sunder
 */ 

#include "hype.h"

/* Definitions of various PDE functions */

//----------------------------------------------------------------------------
// Convert a conserved variable to primitive variable
//----------------------------------------------------------------------------

void PDECons2Prim(const PetscReal* Q, PetscReal* V) {
    
    PetscReal g = 1.0 + (g1-1.0)*(g2-1.0)/((1.0-Q[5])*(g1 -1.0) + Q[5]*(g2-1.0));
    PetscReal P_inf = ((g -1.0)/g)*( g1*p1*Q[5]/(g1 - 1.0) + g2*p2*(1.0 - Q[5])/(g2 - 1.0) );
    PetscReal irho = 1.0/Q[0]; 
    
    V[0] = Q[0];
    V[1] = Q[1]*irho;
    V[2] = Q[2]*irho;
    V[3] = Q[3]*irho; 
    V[4] = (g -1.0)*( Q[4] - 0.5*irho*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]) )  - g*P_inf;
    V[5] = Q[5];
}

//----------------------------------------------------------------------------
// Convert a primitive variable to conserved variable
//----------------------------------------------------------------------------

void PDEPrim2Cons(const PetscReal* V, PetscReal* Q) {
    
    PetscReal g = 1.0 + (g1-1.0)*(g2-1.0)/((1.0-V[5])*(g1 -1.0) + V[5]*(g2-1.0));
    PetscReal P_inf = ((g -1.0)/g)*( g1*p1*V[5]/(g1 - 1.0) + g2*p2*(1.0 - V[5])/(g2 - 1.0) );

    PetscReal e = (V[4] + g*P_inf)/(g - 1.0);
    PetscReal k = 0.5*V[0]*(V[1]*V[1] + V[2]*V[2] + V[3]*V[3]);

    Q[0] = V[0];
    Q[1] = V[0]*V[1];
    Q[2] = V[0]*V[2];
    Q[3] = V[0]*V[3]; 
    Q[4] = k + e;
    Q[5] = V[5];
}

//----------------------------------------------------------------------------
// Find the conservative flux components F in the given normal direction
//----------------------------------------------------------------------------

PetscReal PDEFlux (const PetscReal* Q, 
                   const PetscReal nx, const PetscReal ny, const PetscReal nz, 
                   const PetscReal x,  const PetscReal y,  const PetscReal z, 
                   PetscReal* F) {

    
    PetscReal phi = Q[5];
    PetscReal g = 1.0 + (g1-1.0)*(g2-1.0)/((1.0-phi)*(g1 -1.0) + phi*(g2-1.0));
    PetscReal P_inf = ((g -1.0)/g)*( g1*p1*phi/(g1 - 1.0) + g2*p2*(1.0 - phi)/(g2 - 1.0) );
    PetscReal rho = Q[0];
    PetscReal irho = 1.0/rho;  
    PetscReal u = Q[1]*irho;
    PetscReal v = Q[2]*irho;
    PetscReal w = Q[3]*irho;
    PetscReal un = u*nx + v*ny + w*nz;
    PetscReal rhoun = rho*un;  
    PetscReal E = Q[4]; 
    PetscReal p = (g -1.0)*( E - 0.5*rho*(u*u + v*v + w*w) )  - g*P_inf;
    
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

    PetscReal s_max = PetscAbsReal(un) + PetscSqrtReal(g*(p+P_inf)*irho);

    return s_max;
}

//----------------------------------------------------------------------------
// Non-conservative product BgradQ
//----------------------------------------------------------------------------

void PDENCP(const PetscReal *Q, const PetscReal *grad_Q_x, const PetscReal *grad_Q_y, const PetscReal *grad_Q_z, PetscReal *BgradQ) {
    
    PetscReal phi_x = grad_Q_x[5]; 
    PetscReal phi_y = grad_Q_y[5];
    PetscReal phi_z = grad_Q_z[5];

    BgradQ[0] = 0.0;
    BgradQ[1] = 0.0;
    BgradQ[2] = 0.0; 
    BgradQ[3] = 0.0; 
    BgradQ[4] = 0.0; 
    BgradQ[5] = (Q[1]*phi_x + Q[2]*phi_y + Q[3]*phi_z)/Q[0];
}

//----------------------------------------------------------------------------
// Check physical admissibility of the input state 
//----------------------------------------------------------------------------

PetscBool PDECheckPAD(const PetscReal *Q) {
    
    PetscBool PAD = PETSC_TRUE;
    
    PetscReal phi = Q[5]; 
    PetscReal g = 1.0 + (g1-1.0)*(g2-1.0)/((1.0-phi)*(g1 -1.0) + phi*(g2-1.0));
    PetscReal P_inf = ((g -1.0)/g)*( g1*p1*phi/(g1 - 1.0) + g2*p2*(1.0 - phi)/(g2 - 1.0) );
    PetscReal rho = Q[0];
    PetscReal irho = 1.0/rho;  
    PetscReal u = Q[1]*irho;
    PetscReal v = Q[2]*irho;
    PetscReal w = Q[3]*irho;
    PetscReal p = (g -1.0)*( Q[4] - 0.5*rho*(u*u + v*v + w*w) )  - g*P_inf;
    
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

PetscReal LLFRiemannSolver(const PetscReal *QL, const PetscReal *QR,
                           PetscReal nx, PetscReal ny, PetscReal nz,
                           PetscReal x,  PetscReal y, PetscReal z, 
                           PetscReal *F, PetscReal *D) {
        
    PetscInt c;
    PetscReal FL[nVar], FR[nVar], grad_Q_x[nVar], grad_Q_y[nVar], grad_Q_z[nVar];
    PetscReal Q_av[nVar], dq;

    PetscReal s_max_l = PDEFlux(QL, nx, ny, nz, x, y, z, FL); 
    PetscReal s_max_r = PDEFlux(QR, nx, ny, nz, x, y, z, FR);

    PetscReal s_max = PetscMax(s_max_l, s_max_r);

    for (c = 0; c < nVar; ++c) {
        Q_av[c] = 0.5*(QR[c]+QL[c]);
        dq = QR[c]-QL[c]; 
        grad_Q_x[c] = nx*dq;
        grad_Q_y[c] = ny*dq;
        grad_Q_z[c] = nz*dq;
        F[c] = 0.5*( FR[c] + FL[c] - s_max*dq );
    }
    
    PDENCP(Q_av, grad_Q_x, grad_Q_y, grad_Q_z, D);
    
    return s_max; 
}
