
from numpy import double
import math


#wave speed
c = 250

#observer point
xo = 3.0
yo = 3.0
zo = 3.0

# axisymmetric flow domain
x_min = -20.0
x_max = 20.0
y_min = 0.0
y_max = 20.0
N_x = 400
N_y = 200

# cell size
h = (y_max - y_min)/double(N_y)

# simulation time
InitialTime = 0.0
FinalTime = 0.1
dt = 0.001

#Kirchhoff grid
i_top    =   50
i_bottom =   349
j_curved =   9 

ncells_r = j_curved + 1
ncells_z = (i_bottom - i_top) + 1
ncells_theta = int((2.0*math.pi)/double(h))

z_min = x_min + i_top*h
z_max = x_min + (i_bottom + 1)*h
r_min = 0.0
r_max = y_min + (j_curved + 1)*h
theta_min = 0.0
theta_max = 2.*math.pi


print("Cylinder geometry:\n")
print("cell size h = ", h)
print("z_min = ", z_min)
print("z_max = ", z_max)
print("r_max = ", r_max)
print("ncells_r = ", ncells_r)
print("ncells_z = ", ncells_z)
print("ncells_theta = ", ncells_theta)



