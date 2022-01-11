/*
 * hype.c
 *      Author: sunder
 */ 

static char help[] = "Compute pressure radiated by a monopole source place at (0,0,0) using Green's function and farfield pressure is obtained by Kirchhoff solver.\n\n";

#include "hype.h" 

//----------------------------------------------------------------------------
// Main function of the code 
//----------------------------------------------------------------------------

int main(int argc,char **argv) {

    // --------------------------------------------
    // Initialize MPI 
    //---------------------------------------------
    
    PetscErrorCode ierr;                    // For catching PETSc errors  
    PetscLogDouble start_time, end_time;    // For logging the time values 
    
    ierr = PetscInitialize(&argc, &argv, (char*)0, help);CHKERRQ(ierr);
    
    ierr =  PetscTime(&start_time);CHKERRQ(ierr);  
    
    // --------------------------------------------
    // Set important user defined parameters  
    //---------------------------------------------

    AppCtx Ctx; 

    //!Note: Ensure that point (0,0,0) is inside the domain:
    Ctx.x_min           = -5.0;                            
    Ctx.x_max           =  5.0;                            
    Ctx.y_min           = -5.0;                           
    Ctx.y_max           =  5.0;
    Ctx.z_min           = -5.0; 
    Ctx.z_max           =  5.0;
    Ctx.N_x             =  200;
    Ctx.N_y             =  200;
    Ctx.N_z             =  200;
    Ctx.h               =  0.05; 
    Ctx.InitialStep     =  0; 
    Ctx.InitialTime     =  0.0;
    Ctx.dt              =  0.0001;
    Ctx.Nt              =  1000;
    Ctx.FinalTime       =  Ctx.dt*Ctx.Nt;                                                        
    Ctx.WriteInterval   =  50;      
    Ctx.RestartInterval =  1000;
    Ctx.left_boundary   =  periodic;                   
    Ctx.right_boundary  =  periodic;                   
    Ctx.bottom_boundary =  periodic;                     
    Ctx.top_boundary    =  periodic;    
    Ctx.front_boundary  =  periodic;                     
    Ctx.back_boundary   =  periodic; 
    Ctx.Restart         =  PETSC_FALSE;   
    
    Ctx.box.lower[0]    =   80;
    Ctx.box.lower[1]    =   80;
    Ctx.box.lower[2]    =   80;
    Ctx.box.upper[0]    =   120;
    Ctx.box.upper[1]    =   120;
    Ctx.box.upper[2]    =   120;

    Ctx.xo              =   3.0;
    Ctx.yo              =   3.0;
    Ctx.zo              =   3.0; 

    // --------------------------------------------
    // Data members  
    //---------------------------------------------

    Vec U;                           // Solution Vector
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
    DMBoundaryType z_boundary; 
    
    if (Ctx.left_boundary == periodic || Ctx.right_boundary == periodic)
        x_boundary = DM_BOUNDARY_PERIODIC;
    else
        x_boundary = DM_BOUNDARY_GHOSTED; 

    if (Ctx.bottom_boundary == periodic || Ctx.top_boundary == periodic)
        y_boundary = DM_BOUNDARY_PERIODIC;
    else
        y_boundary = DM_BOUNDARY_GHOSTED;
    
    if (Ctx.front_boundary == periodic || Ctx.back_boundary == periodic)
        z_boundary = DM_BOUNDARY_PERIODIC;
    else
        z_boundary = DM_BOUNDARY_GHOSTED;

    ierr = DMDACreate3d(PETSC_COMM_WORLD,x_boundary,y_boundary,z_boundary,DMDA_STENCIL_BOX,Ctx.N_x,Ctx.N_y,Ctx.N_z,
                        PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,                             
                        nVar,3,NULL,NULL,NULL,&da);CHKERRQ(ierr);  
                        
    ierr = DMSetUp(da);CHKERRQ(ierr);
    
    // Now create various global vectors 

    ierr = DMCreateGlobalVector(da, &U); CHKERRQ(ierr);

    // Set coordinates of cell centers 
    
    ierr = DMDASetUniformCoordinates(da,Ctx.x_min + 0.5*Ctx.h, Ctx.x_max + 0.5*Ctx.h,
                                        Ctx.y_min + 0.5*Ctx.h, Ctx.y_max + 0.5*Ctx.h,
                                        Ctx.z_min + 0.5*Ctx.h, Ctx.z_max + 0.5*Ctx.h);CHKERRQ(ierr);
    
    
    //Kirchhoff Box Surface:
    ierr = GetNumberOfCellsOnBoxFaces(da, &Ctx); CHKERRQ(ierr);
    ierr = GetCellIndicesOnBoxFaces(da, &Ctx); CHKERRQ(ierr);
    ierr = GetQuadPointsOnBoxFaces(da, &Ctx); CHKERRQ(ierr);
    ierr = InitializePressureOnBoxFaces(&Ctx); CHKERRQ(ierr);


    // Set names of the fields
    ierr = DMDASetFieldName(da,0,"p");CHKERRQ(ierr);
    ierr = DMCreateLocalVector(da,&Ctx.localU);CHKERRQ(ierr);    

    // Initialize solution vectors with zero    
    ierr = VecSet(U, 0.0); CHKERRQ(ierr);


    // --------------------------------------------
    // Advance solution in time   
    //---------------------------------------------
    for (int it = 0; it < Ctx.Nt; ++it){

        ierr = ComputePressureExact(U, da, it*Ctx.dt, &Ctx); CHKERRQ(ierr);
        ierr = GetPressureOnBoxFaces(U, da, &Ctx, it); CHKERRQ(ierr);
    


        //Write pressure at observer location:
        if (MyPID == 0) {
            FILE *f = fopen("kirchhoff/p_exact.dat","a+");
            ierr = PetscFPrintf(PETSC_COMM_SELF, f, "%.8e\t%.8e\n", it*Ctx.dt, Pressure(Ctx.xo, Ctx.yo, Ctx.zo, it*Ctx.dt).comp[0]); CHKERRQ(ierr);
            fclose(f);
        }


        //write time.dat and time_step.dat
        if (MyPID == 0) {
        
          FILE *f1 = fopen("kirchhoff/t.dat","a+");
          FILE *f2 = fopen("kirchhoff/dt.dat","a+");

          ierr = PetscFPrintf(PETSC_COMM_SELF, f1, "%.16e\n", it*Ctx.dt); CHKERRQ(ierr);
          ierr = PetscFPrintf(PETSC_COMM_SELF, f2, "%.16e\n", Ctx.dt); CHKERRQ(ierr);
          
          fclose(f1);
          fclose(f2);
        }



        //if (it % 100 == 0) ierr = WriteVtk(da, U, it); CHKERRQ(ierr);
    }

    // --------------------------------------------
    // Free all the memory, finalize MPI and exit   
    //---------------------------------------------
    
    ierr = VecDestroy(&U);CHKERRQ(ierr);
    ierr = DMDestroy(&da);CHKERRQ(ierr);
    ierr = VecDestroy(&Ctx.localU);CHKERRQ(ierr);
    ierr = DestroyKirchoffBoxSurface(&Ctx); CHKERRQ(ierr);


    // --------------------------------------------
    // Print the time taken for simulation       
    //---------------------------------------------
    
    ierr =  PetscTime(&end_time);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Time taken =  %g\n",(double)(end_time - start_time));CHKERRQ(ierr);
    
    ierr = PetscFinalize();CHKERRQ(ierr);
}
