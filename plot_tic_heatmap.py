import numpy as np
import pyemma.coordinates
import pyemma
import pyemma.plots as mplt
from matplotlib import pylab as plt
Y=np.loadtxt('tica.dat')
# histogram data
z,x,y = np.histogram2d(Y[:,0],Y[:,1], bins=50)
# compute free energies
F = -np.log(z)
# contour plot
extent = [x[0], x[-1], y[0], y[-1]]
plt.contourf(F.T, 50, cmap=plt.cm.hot, extent=extent)
plt.show()
