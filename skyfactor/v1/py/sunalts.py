import numpy as np
import pylab as pl

from   diff2          import diff2
from   scipy.optimize import minimize


def Y0(airmass, params):
  ##  [7.05246892e+02   1.82683325e-01]                                                                                                                                                                                       
  return  params[0] * np.exp(-airmass / params[1])

def Y1(airmass, params):
  ##  [1.84684173  1.92248671]                                                                                                                                                                                                
  return  params[0] + params[1] / airmass ** 2.5

def Z(sunalt, params):
  ##  Dependence on sunalt. for a given airmass.                                                                                                                                                                            
  return  np.exp(5. + sunalt / params[1]) - params[0]

## 
fits = []

for airmass in np.arange(1.0, 2.8, 0.2):
    dat = np.loadtxt('dat/sunalt_{:.1f}.txt'.format(airmass))

    pl.plot(dat[:,0], dat[:,1], label=airmass)

    res = minimize(diff2, np.array([np.exp(2.), 1.0]), args=([dat[:,0], dat[:,1], Z]), method='Nelder-Mead', tol=1e-6)
    
    print(airmass, res.x)
    
    pl.plot(dat[:,0], Z(dat[:,0], res.x), 'r-')

    fits.append([airmass, res.x[0], res.x[1]])
    
##
pl.legend()    
pl.show()
'''
pl.clf()

fits = np.array(fits)

pl.plot(fits[:,0], fits[:,1])
pl.plot(fits[:,0], fits[:,2])

res = minimize(diff2, np.array([1., 1.]), args=([fits[:,0], fits[:,1], Y0]), method='Nelder-Mead', tol=1e-6)                                                                                                
pl.plot(fits[:,0], Y0(fits[:,0], res.x), 'c--')

print(res.x)

res = minimize(diff2, np.array([1., 1.]), args=([fits[:,0], fits[:,2], Y1]), method='Nelder-Mead', tol=1e-6)
pl.plot(fits[:,0], Y1(fits[:,0], res.x), 'g--')

print(res.x)

pl.show()

##  Test.
sunalt  = -18.

for airmass in np.arange(1.0, 2.8, 0.2):
    _Y0 = Y0(airmass, np.array([7.05246892e+02, 1.82683325e-01]))
    _Y1 = Y1(airmass, np.array([1.84684173,     1.92248671]))

    _Z  =  Z(sunalt, np.array([_Y0, _Y1]))

    print(airmass, sunalt, _Y0, _Y1, _Z)
'''
