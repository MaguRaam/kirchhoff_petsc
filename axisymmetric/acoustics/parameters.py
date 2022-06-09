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
Tf      =   3.0*3.5e-4



print("speed of sound in water = ", math.sqrt(g1*(p1 + p_inf1)/rho1) );
print("Grid size h = ", h)
print("Final time Tf = ", Tf)


#---------------------------SI units-------------------------------------------#













