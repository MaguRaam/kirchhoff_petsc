import numpy
from matplotlib import pyplot, cm
from math import log



def convergence_plot(h,error,filename):
  fig,ax = pyplot.subplots()
  pyplot.figure(figsize=(35,35))
  pyplot.rcParams.update({'font.size': 30})
  pyplot.loglog(h, error,'-b',label = 'error',marker = 'o',linewidth=10)
  leg = ax.legend()
  pyplot.title("Convergence of Kirchhoff Integral")
  pyplot.xlabel("h")
  pyplot.ylabel("Error")
  pyplot.savefig(filename)

def convergence_rate(error):
  for i in range(1,numpy.size(error)):
    covergence_rate = (log(error[i-1]) - log(error[i]))/log(2)
    print(covergence_rate)



#load h and error abs(p - pexact)
data = numpy.loadtxt("error.dat")
convergence_plot(data[:,0],data[:,1],"convergence")


#linfty rate of convergence:
print("Rate of convergence")
convergence_rate(data[:,1])
