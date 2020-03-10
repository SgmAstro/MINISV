import numpy             as np
import pylab             as pl

import matplotlib.pyplot as plt

from   astropy.table     import Table
from   scipy.optimize    import minimize
from   diff2             import diff2


def model(moonsep, params):
    ##  Dependence on moonsep for given moonalt, moonfrac and airmass. 
    return  params[0] + params[1] * np.abs(moonsep - 96.) ** 1.86328519

## 
airmass   = 2.0

##  airmass, moonfrac, moonalt, moonsep, expfac
dat       = np.loadtxt('dat/moons_{:.1f}.txt'.format(airmass))
dat       = Table(dat, names=['AIRMASS', 'MOONFRAC', 'MOONALT', 'MOONSEP', 'EXPFAC'])

moonalts  = np.unique(dat['MOONALT'].quantity)
moonfracs = np.unique(dat['MOONFRAC'].quantity)

fits      = []

for moonalt in moonalts[::-1]:
    for moonfrac in moonfracs[::-1]:
        ##  Moon separation dependence for fixed moonalt and moonfrac. 
        sample    = dat[(dat['MOONALT'] == moonalt) & (dat['MOONFRAC'] == moonfrac)]

        ##  pl.plot(sample['MOONSEP'], sample['EXPFAC'], 'k-')

        res = minimize(diff2, np.array([10., 1.e-3]), args=([sample['MOONSEP'], sample['EXPFAC'], model]), method='Nelder-Mead', tol=1e-6)

        ##  print(moonalt, moonfrac, model(sample['MOONSEP'], res.x))

        ##  pl.plot(sample['MOONSEP'], np.clip(10. + (sample['MOONSEP'] - 90.)**2. / 1.5e2, 1., None), 'k-')   

        print(res.x)
        
        ##  pl.plot(sample['MOONSEP'], np.clip(model(sample['MOONSEP'], res.x), 1., None), 'c-')

        fits.append([moonalt, moonfrac, res.x[0], res.x[1]])
        
##  pl.show()

##
pl.clf

fits = np.array(fits)
_abs = np.arange(0.1, 90., 1.)

def Y0(moonalt, params):
    ##  To be clipped below at unity.
    return  1. + params[0] * np.log(moonalt / 0.54)

def Y1(moonalt, params):
    ##  Tends to zero at moonalt of 0.0    
    return       params[0] * np.log(moonalt / 0.54)

index  =  3
result = []

for moonfrac in moonfracs[::-1]:
    ##  Dependence on moonalt for fixed moonfrac. 
    sample = fits[fits[:,1] == moonfrac]

    if index == 2:
      res    = minimize(diff2, np.array([2.6]), args=([sample[sample[:,0] > 0.0][:,0], sample[sample[:,0] > 0.0][:,index], Y0]), method='Nelder-Mead', tol=1e-6)
      pl.plot(_abs, np.clip(Y0(_abs, res.x), 1., None), 'k-')   

    if index == 3:
      res   = minimize(diff2, np.array([0.002]), args=([sample[sample[:,0] > 0.0][:,0], sample[sample[:,0] > 0.0][:,index], Y1]), method='Nelder-Mead', tol=1e-6)
      pl.plot(_abs, np.clip(Y1(_abs, res.x), 0.0, None), 'k-')   
    
    pl.plot(sample[:,0], sample[:,index])

    print(res.x)

    result.append([moonfrac, res.x])
    
pl.show()

# pl.clf()

result = np.array(result)

# pl.plot(result[:,0], result[:,1])
# pl.show()

##  np.savetxt('dat/moonfit_{:.1f}_{:d}.txt'.format(airmass, index), result)
