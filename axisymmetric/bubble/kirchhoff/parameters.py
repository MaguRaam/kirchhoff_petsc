#simulation parameters!

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
Tf      =   8.2410553987786652e-03

#cylinder surface:
j_curved        =  230
i_top           =  10
i_bottom        =  490

Ncells_r        =   j_curved + 1
Ncells_z        =   (i_bottom - i_top) + 1
Ncells_theta    =   200

r_min           =   0.0
r_max           =   (j_curved + 0.5)*h
theta_min       =   0.0
theta_max       =   2*math.pi
z_min           =   xmin + (i_top + 0.5)*h 
z_max           =   xmin + (i_bottom + 0.5)*h

#observer point
io = 250
jo = 240
z_o = xmin + (io + 0.5)*h
r_o = (jo + 0.5)*h

print("zo =", z_o)
print("ro =", r_o)

print("r_min = ", r_min)
print("r_max = ", r_max)
print("theta_min = ", theta_min)
print("theta_max = ", theta_max)
print("z_min = ", z_min)
print("z_max = ", z_max)



print("speed of sound in water = ", math.sqrt(g1*(p1 + p_inf1)/rho1) );
print("Grid size h = ", h)
print("Final time Tf = ", Tf)


#---------------------------SI units-------------------------------------------#













