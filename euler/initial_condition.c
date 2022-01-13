/*
 * initial_condition.c
 *      Author: sunder
 */
#include "hype.h"



//----------------------------------------------------------------------------
// Initial condition function
//----------------------------------------------------------------------------

//TODO Initial condition works for only unit domain:

Field InitialCondition(PetscReal x, PetscReal y, PetscReal z)
{

    Field Q0;
    Field V0;

    //-----------------------------------------------------------------------

    //acoustic perturbation
    PetscReal rho_, u_, v_, w_, p_;                         
    PetscReal r = PetscSqrtReal((x-1.0)*(x-1.0) + (y-1.0)*(y-1.0) + (z-1.0)*(z-1.0));
    
    //amplitude of density perturbation:
    PetscReal A = 0.01;

    //density and pressure perturbation:
    if (r < .125){
        rho_ = A*PetscExpReal(-30.0*r*r*r);
        p_   = c0*c0*rho_;
    }
    else {
        rho_ = 0.0;
        p_   = 0.0;
    }

    //velocity perturbation:
    u_ = 0.0;
    v_ = 0.0;
    w_ = 0.0;


    //inital state = mean_state + perturbed_state
    V0.comp[0] = rho0 + rho_;
    V0.comp[1] = u0 + u_;
    V0.comp[2] = v0 + v_;
    V0.comp[3] = w0 + w_;
    V0.comp[4] = p0 + p_;

    //-----------------------------------------------------------------------

    PDEPrim2Cons(&V0, &Q0);

    return Q0;
}
