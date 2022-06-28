#plot cylinder
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np


#load data
data = np.loadtxt("cylinder.dat")


# Creating figure
fig = plt.figure()
ax = plt.axes(projection ="3d")
ax.set_xlabel('X', fontweight ='bold')
ax.set_ylabel('Y', fontweight ='bold')
ax.set_zlabel('Z', fontweight ='bold')

# Creating plot
ax.scatter3D(data[:,0], data[:,1], data[:,2], color = "blue")
plt.savefig("quadpoints.png")
plt.show()
