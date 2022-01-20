import numpy as np 
import os
from datetime import datetime

dim = 3
nVar = 6 

def plot_tecplot(in_file, out_file):
    
    infile = open(in_file, "r")
    
    # Read header lines 

    line = infile.readline()
    line = infile.readline()
    line = infile.readline()
    line = infile.readline()
    line = infile.readline()
    split_line = line.split()

    N_x = int(split_line[1])
    N_y = int(split_line[2])
    N_z = int(split_line[3])


    line = infile.readline()

    # Get number of cells in the mesh 

    split_line = line.split()
    n_cells = int(split_line[1])

    W = np.zeros((n_cells,nVar))

    coords = np.zeros((n_cells,dim))

    # Read the coordinates 

    for i in range(0, n_cells):
        line = infile.readline()
        split_line = line.split()
        coords[i, 0] = np.float64(split_line[0])
        coords[i, 1] = np.float64(split_line[1])
        coords[i, 2] = np.float64(split_line[2])

    # Read some more headers 
        
    line = infile.readline()
    line = infile.readline()
    line = infile.readline()

    for i in range(0, n_cells):
        line = infile.readline()
        split_line = line.split()
        for c in range (0, nVar):
            W[i,c] = np.float64(split_line[c])
        
    infile.close()

    # Reshape arrays to fit rectangular coordinates 

    coords = coords.reshape((N_z, N_y, N_x, dim))
    W = W.reshape((N_z, N_y, N_x, nVar))

    # Start plotting 
            
    outfile = open(out_file, "w+")
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    
    print("TITLE = \"Euler Equations. File created on: ", dt_string, "\"",  file = outfile)
    outfile.write("\nVARIABLES = \"x\", \"y\", \"z\", \"Density\", \"V_x\", \"V_y\", \"V_z\", \"Pressure\", \"Phi\" ")
    print('\nZone I = ', N_x, "J = ", N_y, "K = ", N_z, file = outfile) 
    
    for k in range(0, N_z):
        for j in range(0, N_y):
            for i in range (0, N_x):
                print(coords[k,j,i,0], coords[k,j,i,1], coords[k,j,i,2], W[k,j,i,0], W[k,j,i,1], W[k,j,i,2], W[k,j,i,3], W[k,j,i,4], W[k,j,i,5], file = outfile)

            
    outfile.close()
###############################################################################################################################

for in_file in os.listdir('.'):
    if in_file.endswith('.vtk'):
        print(in_file)
        out_file = in_file.replace("vtk", "dat")
        plot_tecplot(in_file, out_file)

