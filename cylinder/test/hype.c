static char help[] = "Write pressure data from acoustic monopole source to validate Kirchhoff solver.\n\n";

#include "hype.h" 

//----------------------------------------------------------------------------
// Main function of the code 
//----------------------------------------------------------------------------

int main(int argc,char **argv){

    // --------------------------------------------
    // Initialize MPI 
    //---------------------------------------------

    PetscErrorCode ierr;                    /* For catching PETSc errors */ 
    PetscLogDouble start_time, end_time;    /* For logging the time values */

    ierr = PetscInitialize(&argc, &argv, (char*)0, help);CHKERRQ(ierr);
    ierr =  PetscTime(&start_time);CHKERRQ(ierr); 

    // --------------------------------------------
    // Set important user defined parameters  
    //---------------------------------------------

    AppCtx Ctx; 

    Ctx.x_min           = -20.0;
    Ctx.x_max           =  20.0;
    Ctx.y_min           =  0.0;
    Ctx.y_max           =  20.0;
    Ctx.N_x             =  400;
    Ctx.N_y             =  200;
    Ctx.dt              =  0.001;
    Ctx.InitialStep     =  0;
    Ctx.InitialTime     =  0.0;
    Ctx.FinalTime       =  0.1;
    Ctx.left_boundary   =  transmissive;
    Ctx.right_boundary  =  transmissive;
    Ctx.bottom_boundary =  transmissive;
    Ctx.top_boundary    =  transmissive;
    Ctx.Restart         =  PETSC_FALSE;
    Ctx.h               =  (Ctx.x_max - Ctx.x_min)/(PetscReal)(Ctx.N_x);

    //index of Kirchhoff surface: 
    Ctx.surface.i_top    =   50;
    Ctx.surface.i_bottom =   349;
    Ctx.surface.j_curved =   9; 

    //observer grid point:
    Ctx.xo               =   3.0;
    Ctx.yo               =   3.0;
    Ctx.zo               =   3.0;

    // --------------------------------------------
    // Data members  
    //---------------------------------------------

    Vec U;                           // Solution Vector (Conserved variables)
    DM da;                           // Grid object 
    PetscMPIInt MyPID;               // Rank of the current processor 
    PetscMPIInt numProcs;            // Size of the communicator

    // --------------------------------------------
    // Obtain the rank of the process and size of 
    // the communicator 
    //---------------------------------------------

    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&numProcs);CHKERRQ(ierr); 
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&MyPID);CHKERRQ(ierr); 
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Code running with %d processes\n", numProcs);CHKERRQ(ierr); 

    // --------------------------------------------
    // Initialize the grid and set field names
    //---------------------------------------------

    DMBoundaryType x_boundary;
    DMBoundaryType y_boundary;

    if (Ctx.left_boundary == periodic || Ctx.right_boundary == periodic)
        x_boundary = DM_BOUNDARY_PERIODIC;
    else
        x_boundary = DM_BOUNDARY_GHOSTED; 

    if (Ctx.bottom_boundary == periodic || Ctx.top_boundary == periodic)
        y_boundary = DM_BOUNDARY_PERIODIC;
    else
        y_boundary = DM_BOUNDARY_GHOSTED; 

    ierr = DMDACreate2d(PETSC_COMM_WORLD, // Global communicator      
                        x_boundary,       // Boundary conditions in x-direction 
                        y_boundary,       // Boundary conditions in y-direction
                        DMDA_STENCIL_BOX, // Stencil type (other is star type)
                        Ctx.N_x,          // No. of cells in x-direction 
                        Ctx.N_y,          // No. of cells in y-direction
                        PETSC_DECIDE,     // Domain decomposition in x-direction 
                        PETSC_DECIDE,     // Domain decomposition in y-direction
                        nVar,             // No. of dofs per cell 
                        3,                // Width of the stencil
                        NULL,
                        NULL,
                        &da);CHKERRQ(ierr); // da object 

    ierr = DMSetUp(da);CHKERRQ(ierr);

    // Now create various global vectors 
    ierr = DMCreateGlobalVector(da, &U);CHKERRQ(ierr);
    ierr = DMCreateLocalVector(da,&Ctx.localU);CHKERRQ(ierr);

    // Set coordinates of cell centers 
    ierr = DMDASetUniformCoordinates(da,
                                         Ctx.x_min + 0.5*Ctx.h, Ctx.x_max - 0.5*Ctx.h,
                                        Ctx.y_min + 0.5*Ctx.h, Ctx.y_max - 0.5*Ctx.h,
                                        0.0,0.0);CHKERRQ(ierr);

    // Allocate Memory for Kirchhoff line surface:
    ierr = GetCellIndexOnSurface(da, &Ctx); CHKERRQ(ierr);
    ierr = InitializePressureOnSurface(&Ctx); CHKERRQ(ierr);

    // Set names of the fields
    ierr = PetscObjectSetName((PetscObject)U,"sol");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,0,"p");CHKERRQ(ierr);
    
    // Initialize solution vectors with zero    
    ierr = VecSet(U, 0.0); CHKERRQ(ierr);

    // --------------------------------------------
    // Advance solution in time   
    //---------------------------------------------
    PetscInt it = 0;
    while (it*Ctx.dt < Ctx.FinalTime){

        ierr = ComputePressureExact(U, da, it*Ctx.dt, Ctx); CHKERRQ(ierr);
        ierr = GetPressureOnSurface(U, da, &Ctx, it); CHKERRQ(ierr);

         //Write pressure at observer location:
        if (MyPID == 0) {
            FILE *f = fopen("kirchhoff/p_exact.dat","a+");
            ierr = PetscFPrintf(PETSC_COMM_SELF, f, "%.8e\t%.8e\n", it*Ctx.dt, Pressure(Ctx.xo, Ctx.yo, Ctx.zo, it*Ctx.dt).comp[0]); CHKERRQ(ierr);
            fclose(f);
        }


        //write time.dat and time_step.dat
        if (MyPID == 0) {
        
          FILE *f = fopen("kirchhoff/input/t.dat","a+");
          ierr = PetscFPrintf(PETSC_COMM_SELF, f, "%.16e\n", it*Ctx.dt); CHKERRQ(ierr);
          fclose(f);
          
        }
				
		    //write data vtu data:
		    ierr = WriteVtk(da, U, it); CHKERRQ(ierr);

        it++;
    }


    // --------------------------------------------
    // Free all the memory, finalize MPI and exit   
    //---------------------------------------------
    
    ierr = VecDestroy(&U);CHKERRQ(ierr);
    ierr = VecDestroy(&Ctx.localU);CHKERRQ(ierr);
    ierr = DMDestroy(&da);CHKERRQ(ierr);
    ierr = DestroyKirchoffSurface(&Ctx); CHKERRQ(ierr);

    // --------------------------------------------
    // Print the time taken for simulation       
    //---------------------------------------------

    ierr =  PetscTime(&end_time);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Time taken =  %g\n",(double)(end_time - start_time));CHKERRQ(ierr);

    ierr = PetscFinalize();CHKERRQ(ierr);
    return ierr;
}
