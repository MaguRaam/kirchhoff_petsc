#simulation parameters!

import math


#properties of multiphase system (air bubble in water)

#water medium  (phi = 1):
g1                = 4.4      	#Specific heat ratio of water
p_inf1            = 6000.0   	#Stiffness constant of water
rho1              = 1000.0      #initial density of water
p1                = 1.0e6      	#initial pressure of water


#air bubble (phi = 0):
g2                = 1.4      	#Specific heat ratio of air
p_inf2            = 0.0      	#Stiffness constant of air
rho2              = 1.0   	 	#initial density of air
p2                = 1.0e5    	#initial pressure of air


#------------------------------------------------------#

#bubble radius:
R  =  0.038


#Grid:
N = 250.0
xmin = -10.0*R
xmax = +10.0*R
h = (xmax - xmin)/N

#Kirchhoff box
lower    =   30
upper    =   220

#Reciever location 
iox  =   240

#-----------------------------------------------------------#


#Reciever location (coordinate of cell center)
ox   =  xmin + (iox)*h - 0.5*h

#Kirchhoff box coordinate:
lowerx = xmin + (lower)*h - h
upperx = xmin + (upper)*h


print("speed of sound in water = ", math.sqrt(g1*(p1 + p_inf1)/rho1) );
print('Kirchhoff box ranges from [',lowerx, upperx,']')
print("observer location = ",ox)
print("Grid size h = ", h)
print("no of grid points in bubble = ", 2.0*R/h)



