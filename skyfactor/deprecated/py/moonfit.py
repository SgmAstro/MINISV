import numpy as np
import pylab as pl

from   diff2          import diff2
from   scipy.optimize import minimize


def X0(airmass, params):
      ##  Dependence on moonfrac for a given airmass.                                                                                                                                           
      return  params[0] + params[1] * airmass ** 0.75

def Y0(moonfrac, params):
      ##  Dependence on moonfrac for a given airmass. 
      Y = moonfrac - 0.1
      
      return  params[0] * np.exp(Y / 0.29) - params[0]


_airmass  = [1.0, 1.2, 1.6, 2.0, 2.4, 2.8, 3.0]
fits      = []

index     =  2

for airmass in _airmass:
  dat = np.loadtxt('dat/moonfit_{:.1f}_{:d}.txt'.format(airmass, index))
  pl.plot(dat[:,0], dat[:,1], label='{:.2f}'.format(airmass))

  res = minimize(diff2, np.array([1.0]), args=([dat[:,0], dat[:,1], Y0]), method='Nelder-Mead', tol=1e-6)

  # print(airmass, res.x)

  pl.plot(dat[:,0], Y0(dat[:,0], res.x), 'k-')

  fits.append([airmass, res.x])
  
pl.legend(frameon=False)
pl.show()
pl.clf()

##
fits = np.array(fits)

pl.plot(fits[:,0], fits[:,1], 'k-')

_abs = np.arange(0.0, 3.0, 0.01)

res  = minimize(diff2, np.array([0.0, 1.0]), args=([fits[:,0], fits[:,1], X0]), method='Nelder-Mead', tol=1e-6)

print(res.x)

pl.plot(fits[:,0], X0(fits[:,0], res.x), 'c--')

pl.show()
