#plot cylindrical wave 
import matplotlib.pyplot as plt
import numpy as np



#load data
data = np.loadtxt("plot.dat")

#plot
plt.plot(data[:,0], data[:,1],linestyle="-", color='red')


plt.title("cylindrical wave at t = 1")
plt.xlabel("r")
plt.ylabel("p")
plt.savefig("t_1.png")
