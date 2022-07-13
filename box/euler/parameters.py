#simulation parameters!

import math


#speed of sound using Newton–Laplace equation
def sound_speed(K, rho):
    return math.sqrt(K/rho)



#Acoustic medium (ideal gas)
gamma = 1.4                     #specific heat ratio
rho0 = 1.0                      #mean density
p0 = 1.0                        #mean pressure
K0 = gamma*p0                   #bulk modulus
c0 = sound_speed(K0, rho0)      #speed of sound

print("speed of sound = ", end="")
print('%.17f'%c0)

#------------------------------------------------------#

#Grid:
N = 200.0
xmin = 0.0
xmax = 2.0
h = (xmax - xmin)/N

#Kirchhoff box
lower    =   40
upper    =   160

#Reciever location (cell index)
iox  =   170

#-----------------------------------------------------------#


#Reciever location (coordinate of cell center)
ox   =  xmin + (iox)*h - 0.5*h;

#Kirchhoff box coordinate:
lowerx = xmin + (lower)*h - h
upperx = xmin + (upper)*h

print('Kirchhoff box ranges from [',lowerx, upperx,']')
print("observer location = ",ox)
print("Grid size h = ", h)











