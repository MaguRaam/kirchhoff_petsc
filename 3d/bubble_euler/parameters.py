#simulation parameters!

import math


#properties of multiphase system (air bubble in water)

#water medium  (phi = 1):
g1                = 4.4      #Specific heat ratio of water
p_inf1            = 6000.0   #Stiffness constant of water
rho1              = 1.0      #initial density of water
p1                = 1.0      #initial pressure of water


#air bubble (phi = 0):
g2                = 1.4      #Specific heat ratio of air
p_inf2            = 0.0      #Stiffness constant of air
rho2              = 1.0e-3   #initial density of air
p2                = 1.0e-1   #initial pressure of air


#------------------------------------------------------#

#Grid:
N = 200.0
xmin = -5.0
xmax = 5.0
h = (xmax - xmin)/N

#Kirchhoff box
lower    =   40
upper    =   160

#Reciever location 
iox  =   180

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




