#Plot exact solution and Kirchoff solution at observer location 

import numpy as np
import matplotlib.pyplot as plt


p_exact = np.loadtxt("p_exact.dat")
p_kirchhoff = np.loadtxt("p_kirchhoff.dat")


fig,ax = plt.subplots()
plt.figure(figsize=(14,14))
plt.rcParams.update({'font.size': 21})
plt.xlabel('observer time, $t_{0}$')
plt.ylabel('acoustic pressure, $p$')
leg1,= plt.plot(p_kirchhoff[:,0], p_kirchhoff[:,1],"b",marker = 'o', linewidth=2)
leg2,= plt.plot(p_exact[:,0],p_exact[:,1],"r-", linewidth=2)
plt.legend((leg1, leg2), ('Kirchhoff sol.', 'Analytical sol.'), loc='upper right', fontsize=22, numpoints=1)
plt.savefig("Pressure")
