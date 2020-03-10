import  numpy           as      np
import  pylab           as      pl

from    diff2           import  diff2
from    scipy.optimize  import  minimize


def Z(sunalt, params):
  ##  Dependence on sunalt. for a given airmass.                                                                                                                                                                            
  return  np.exp(5. + sunalt / params[1]) - params[0]

## 
fits = []

for airmass in [1.0, 1.2, 1.6, 2.0, 2.4, 2.8, 3.0]:
    dat = np.loadtxt('../../dat/sunalt_{:.1f}.txt'.format(airmass))

    pl.plot(dat[:,0], dat[:,1], label=airmass)

    res = minimize(diff2, np.array([np.exp(2.), 1.0]), args=([dat[:,0], dat[:,1], Z]), method='Nelder-Mead', tol=1e-6)
    
    print(airmass, res.x)
    
    pl.plot(dat[:,0], Z(dat[:,0], res.x), 'k-', alpha=0.6)

    fits.append([airmass, res.x[0], res.x[1]])
    
##
pl.xlabel('SUN ALT.')

pl.legend(frameon=False, ncol=2, loc=2)    

pl.show()

exit(1)

pl.clf()

fits = np.array(fits)

pl.plot(fits[:,0], fits[:,1])
pl.plot(fits[:,0], fits[:,2])

def Y0(airmass, params):
  ##  [6.9489374,  -3.96722344]                                                                                                                                                                                                       
  return  params[0] + params[1] * airmass

def Y1(airmass, params):
  ##  [4.07622358, -0.09630329]                                                                                                                                                                                                       
  return  params[0] + params[1] * airmass

res = minimize(diff2, np.array([1., 1.]), args=([fits[:,0], fits[:,1], Y0]), method='Nelder-Mead', tol=1e-6)                                                                                                
pl.plot(fits[:,0], Y0(fits[:,0], res.x), 'c--')

print(res.x)

res = minimize(diff2, np.array([1., 1.]), args=([fits[:,0], fits[:,2], Y1]), method='Nelder-Mead', tol=1e-6)
pl.plot(fits[:,0], Y1(fits[:,0], res.x), 'g--')

print(res.x)

pl.show()
