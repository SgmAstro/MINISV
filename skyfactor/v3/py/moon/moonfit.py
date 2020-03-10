import sys
import numpy           as      np
import pylab           as      pl

from   diff2           import  diff2
from   scipy.optimize  import  minimize


def Y0(moonfrac, params):
  ##  Dependence of Z0 & Z1 on moonfrac for a given airmass. 
  Y = moonfrac - 0.1
      
  return  params[0] * np.exp(Y / 0.29) - params[0]

##  
_airmass  = [1.0, 1.2, 1.6, 2.0, 2.4, 2.8, 3.0]
fits      = []

##  Which parameter to model, [2, 3]. 
index     =  np.int(sys.argv[1])

for airmass in _airmass:
  ##  Fit for dependence of common parameter to Z0 and Z1 on moonfrac., for given airmass, for both 
  ##  parameters to moonmodel, indexed by 2 or 3.
  dat = np.loadtxt('../../dat/moonfit_{:.1f}_{:d}.txt'.format(airmass, index))
  pl.plot(dat[:,0], dat[:,1], label='{:.2f}'.format(airmass))

  res = minimize(diff2, np.array([1.0]), args=([dat[:,0], dat[:,1], Y0]), method='Nelder-Mead', tol=1e-6)

  # print(airmass, res.x)

  pl.plot(dat[:,0], Y0(dat[:,0], res.x), 'k-')

  fits.append([airmass, res.x])

pl.title('INDEX {:d}'.format(index))  
pl.legend(frameon=False)
pl.show()
pl.clf()

##  
##  Remaining dependence of the parameter to Y0 on airmass.  Later to be X0. 
def X0(airmass, params):
  ##  Dependence on moonfrac for a given airmass.                                                                                                                                                                                   
  return  params[0] + params[1] * airmass ** 0.75

##
fits = np.array(fits)

pl.plot(fits[:,0], fits[:,1], 'k-')

_abs = np.arange(0.0, 3.0, 0.01)
res  = minimize(diff2, np.array([0.0, 1.0]), args=([fits[:,0], fits[:,1], X0]), method='Nelder-Mead', tol=1e-6)

print(index, res.x)

pl.plot(fits[:,0], X0(fits[:,0], res.x), 'c--')
pl.xlabel('AIRMASS')
pl.ylabel('X0')
pl.title('INDEX: {}'.format(index))
pl.show()
