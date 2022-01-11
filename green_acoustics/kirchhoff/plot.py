#Plot Green's function solution and Kirchoff solution at observer location 

import numpy as np
import matplotlib.pyplot as plt

#source frequency
f0 = 100.0
t0 = 4.0/ f0

#time step
dt = 0.0001
nt = 1000

#time axis
time = np.linspace(0 * dt, nt * dt, nt)

#source function
src  = -2. * (time - t0) * (f0 ** 2) * (np.exp(-1.0 * (f0 ** 2) * (time - t0) ** 2))


#plot source time function
plt.figure(figsize=(9,9))
plt.plot(time, src)
plt.title('Source time function')
plt.xlabel('Time')
plt.ylabel('Amplitude')
plt.grid()
plt.savefig("Source")


#Plot Green's function solution and Kirchoff solution at observer location 
p_exact = np.loadtxt("p_exact.dat")
p_kirchhoff = np.loadtxt("p_kirchhoff.dat")


fig,ax = plt.subplots()
plt.figure(figsize=(9,9))
plt.rcParams.update({'font.size': 17})
plt.xlabel('observer time, $t_{0}$')
plt.ylabel('acoustic pressure, $p - p_{0}$')
leg1,= plt.plot(p_kirchhoff[:,0], p_kirchhoff[:,1],"r",marker = 'o')
leg2,= plt.plot(p_exact[:,0],p_exact[:,1],"b-")
plt.legend((leg1, leg2), ('Kirchhoff sol.', 'Analytical sol.'), loc='upper right', fontsize=15, numpoints=1)
plt.savefig("Pressure")

