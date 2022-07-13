/*
 * initial_condition.c
 *      Author: sunder
 */ 
#include "hype.h" 

//----------------------------------------------------------------------------
// Initial condition function 
//----------------------------------------------------------------------------

void InitialCondition(PetscReal x, PetscReal y, PetscReal z, PetscReal* Q0) {

    PetscReal V0[nVar];

    //----------------------------------------------------------------------- 
    PetscReal R0 = 0.038;
    const PetscReal x0 = 0.0; const PetscReal y0 = 0.0; const PetscReal z0 = 0.0;
    PetscReal R = PetscSqrtReal((x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0));

    // Outside the bubble
    PetscReal rhoo   = 1000.0;  
    PetscReal po     = 1.0e6;
    PetscReal alphao = 1.0;


    // Inside the bubble 
    PetscReal rhoi   = 1.0;
    PetscReal pi     = 1.0e5;
    PetscReal alphai = 0.0;

    PetscReal h = 20.0*R0/250.0;     
    PetscReal smear = 2.0;


    V0[0] = 0.5*((rhoi+rhoo) + (rhoo-rhoi)*PetscTanhReal((R-R0)/(smear*h)));
    V0[1] = 0.0;
    V0[2] = 0.0; 
    V0[3] = 0.0; 
    V0[4] = 0.5*((pi+po) + (po-pi)*PetscTanhReal((R-R0)/(smear*h)));
    V0[5] = 0.5*((alphai+alphao) + (alphao-alphai)*PetscTanhReal((R-R0)/(smear*h)));


    //----------------------------------------------------------------------- 
    PDEPrim2Cons(V0, Q0);
}


