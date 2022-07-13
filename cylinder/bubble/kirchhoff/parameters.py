#simulation parameters!
from numpy import double
import math

#----------------------Scaled Values----------------------------------------#

#properties of multiphase system (air bubble in water)

#water medium  (phi = 1):
g1                = 4.4          #Specific heat ratio of water
p_inf1            = 6000.0       #Stiffness constant of water
rho1              = 1000.0       #initial density of water
p1                = 1.0e6        #initial pressure of water


#air bubble (phi = 0):
g2                = 1.4          #Specific heat ratio of air
p_inf2            = 0.0          #Stiffness constant of air
rho2              = 1.0          #initial density of air
p2                = 1.0e5        #initial pressure of air


#------------------------------------------------------#

#bubble radius:
R  =  0.038



#Grid: (z, r)
xmin    =   -10.0*R
xmax    =   10.0*R
ymin    =   0
ymax    =   10.0*R
Nx      =   500
Ny      =   250
h       =   (xmax - xmin)/Nx
Tf      =   0.2

#ncells/radius
print("ncells/radius = ", R/h)


#cylinder surface:
j_curved        =  150
i_top           =  25
i_bottom        =  475

Ncells_r        =   j_curved + 1
Ncells_z        =   (i_bottom - i_top) + 1
Ncells_theta    =   int((2.0*math.pi)/double(h))

z_min = xmin + i_top*h              #-9R
z_max = xmin + (i_bottom + 1)*h     #+9R
r_min = 0.0                         
r_max = ymin + (j_curved + 1)*h     #+6R
theta_min = 0.0
theta_max = 2.*math.pi

#observer point
io = 250
jo = 225
z_o = xmin + (io + 0.5)*h           #0.R
r_o = ymin + (jo + 0.5)*h           #9.0R

print("zo =", z_o)
print("ro =", r_o)

print("r_min = ", r_min)
print("r_max = ", r_max)
print("theta_min = ", theta_min)
print("theta_max = ", theta_max)
print("z_min = ", z_min)
print("z_max = ", z_max)
print("ncells_r = ", Ncells_r)
print("ncells_z = ", Ncells_z)
print("ncells_theta = ", Ncells_theta)



print("speed of sound in water = ", math.sqrt(g1*(p1 + p_inf1)/rho1) );
print("Grid size h = ", h)
print("Final time Tf = ", Tf)


#---------------------------SI units-------------------------------------------#













