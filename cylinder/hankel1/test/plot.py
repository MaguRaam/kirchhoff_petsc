#plot cylindrical wave 
import matplotlib.pyplot as plt
import numpy as np



#load data
data = np.loadtxt("plot.dat")

#plot
fig,ax = plt.subplots()
plt.figure(figsize=(14,14))
plt.rcParams.update({'font.size': 30})
plt.plot(data[:,0], data[:,1],linestyle="-", color='red', linewidth=3)

plt.xlabel('radius, $r$')
plt.ylabel('acoustic pressure, $p - p_{0}$')
plt.savefig("cylindrical.png")
