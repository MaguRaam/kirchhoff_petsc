/*
 * monitor_function.c
 *      Author: sunder
 */ 

#include "hype.h"  

//----------------------------------------------------------------------------
// Monitor function for additional processing in the intermediate time steps 
//----------------------------------------------------------------------------

PetscErrorCode MonitorFunction (TS ts,PetscInt step, PetscReal time, Vec U, void *ctx) {

    PetscErrorCode ierr; 
    AppCtx *Ctx = (AppCtx*)ctx;
    DM da;
    PetscMPIInt MyPID;

    ierr = TSGetDM(ts,&da); CHKERRQ(ierr);

    // Set the time step based on CFL condition 

    ierr = TSSetTimeStep(ts, Ctx->dt); CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,"%d t = %f\n", step, time); CHKERRQ(ierr);

    // Plot the solution at the required time interval 

    if (Ctx->WriteInterval != 0) {

        if(step%Ctx->WriteInterval == 0) {
            
            char filename[30];
            sprintf(filename, "plot/sol-%08d.vts", step); // 4 is the padding level, increase it for longer simulations 
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing data in vts format to %s at t = %f, step = %d\n", filename, time, step); CHKERRQ(ierr);
            PetscViewer viewer;  
            ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);CHKERRQ(ierr);
            ierr = DMView(da, viewer);
            VecView(U, viewer);
            
            ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
        
        }
    }

    if (Ctx->RestartInterval != 0) {
        
        PetscViewer    viewer_binary;

        if(step%Ctx->RestartInterval == 0) {
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing data in binary to restart1.bin at t = %f\n", time); CHKERRQ(ierr);
            ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"restart1.bin",FILE_MODE_WRITE, &viewer_binary); CHKERRQ(ierr);
            ierr = VecView(U,viewer_binary); CHKERRQ(ierr);
            ierr = PetscViewerDestroy(&viewer_binary); CHKERRQ(ierr);
        }
        
        if(step%(Ctx->RestartInterval + 7) == 0) {
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing data in binary to restart2.bin at t = %f\n", time); CHKERRQ(ierr);
            ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"restart2.bin",FILE_MODE_WRITE, &viewer_binary); CHKERRQ(ierr);
            ierr = VecView(U,viewer_binary); CHKERRQ(ierr);
            ierr = PetscViewerDestroy(&viewer_binary); CHKERRQ(ierr);
        }
    }




    //Write the perturbation pressure p' = (p - p0) at the observer grid point:
    Field   ***u; 
    Field   V;
    PetscInt i, j, k, xs, ys, zs, xm, ym, zm;
    ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da, U, &u); CHKERRQ(ierr); 

    for (k = zs; k < zs + zm; ++k){
        for (j = ys; j < ys + ym; ++j){
            for (i = xs; i < xs + xm; ++i){

                if (i == Ctx->iox && j == Ctx->ioy && k == Ctx->ioz)
                {
                    //get primitive variables from conserved variables:
                    PDECons2Prim(&u[k][j][i], &V);

                    //write p' = (p - p0):
                    FILE *f = fopen("kirchhoff/p_weno.dat", "a+");
                    ierr = PetscFPrintf(PETSC_COMM_SELF, f, "%e\t%.12e\n", time, V.comp[4] - p0); CHKERRQ(ierr);
                    fclose(f);
                }
            }
        }
    }

    ierr = DMDAVecRestoreArrayRead(da, U, &u); CHKERRQ(ierr);

    //write time.dat and time_step.dat
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&MyPID);CHKERRQ(ierr);

    if (MyPID == 0) {
    
      FILE *f1 = fopen("kirchhoff/input/t.dat","a+");
      FILE *f2 = fopen("kirchhoff/input/dt.dat","a+");

      ierr = PetscFPrintf(PETSC_COMM_SELF, f1, "%.16e\n", time); CHKERRQ(ierr);
      ierr = PetscFPrintf(PETSC_COMM_SELF, f2, "%.16e\n", Ctx->dt); CHKERRQ(ierr);
      
      fclose(f1);
      fclose(f2);
    }


    //Write p', p'n on Kirchhoff surface:
    ierr = GetPressureOnBoxFaces(U, da, Ctx, step); CHKERRQ(ierr);


    return ierr; 
} 

