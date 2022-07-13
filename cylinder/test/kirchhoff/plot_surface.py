import numpy as np
import matplotlib.pyplot as plt

# create empty pressure array:
p = []

# read kirchhoff surface data from input/top or input/bottom or input/curved
filename = "input/curved/p-00100.dat"

file = open(filename, "r")

# skip first two lines
line = file.readline()
line = file.readline()

for line in file:

    # break when you reach end of pressure vector:
    if "Vec" in line:
        break

    # skip process
    if "Process" in line:
        continue
    
    #append data to array
    p.append(float(line))

file.close()

plt.plot(p)
fig = filename.replace(".dat",".png")
plt.savefig(fig)