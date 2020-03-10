import sys
import numpy             as     np
import pylab             as     pl

import matplotlib.pyplot as     plt

from   astropy.table     import Table
from   scipy.optimize    import minimize
from   diff2             import diff2


def moonmodel(moonsep, params):
    ##  2-parameter model for dependence on moonsep for given moonalt, 
    ##  moonfrac and airmass. 
    ##
    ##  Note:  assumed fixed center for moonsep, and powerlaw index following 
    ##         trial and error with these parameters. 
    ##
    ##         params[0] ->  dependence on moon alt. -> Z0.
    ##         params[1] ->  dependence on moon alt. -> Z1.      

    return  params[0] + params[1] * np.abs(moonsep - 96.) ** 1.863

##  Airmass to run provided by command line. 
all_airmass = [1.0, 1.2, 1.6, 2.0, 2.4, 2.8, 3.0]
airmass     = all_airmass[np.int(sys.argv[1])]

##  CHH Moon model called for a grid of conditions. 
dat         = np.loadtxt('../../dat/moons_{:.1f}.txt'.format(airmass))
dat         = Table(dat, names=['AIRMASS', 'MOONFRAC', 'MOONALT', 'MOONSEP', '4000', '6700', '7160', '8070', '9060'])

##  Get normalisation.
norm        = np.loadtxt('../../dat/moon_norm.txt')
norm        = Table(norm, names=['AIRMASS', 'MOONFRAC', 'MOONALT', 'MOONSEP', '4000', '6700', '7160', '8070', '9060'])

##  Normalisation dat. 
for band in ['4000', '6700', '7160', '8070', '9060']:
    dat[band] /= norm[band]

    dat[band]  = np.clip(dat[band], 1.0, None)

print(dat)

##  Choose a band to define the expfac. 
dat['EXPFAC'] = dat['4000']

assert np.all(dat['EXPFAC'] >= 1.0)

##  Assumed by moon model. 
assert np.all(dat[dat['MOONALT'] < 0.5]['EXPFAC'] == 1.0)

##  Now start the fitting. 
moonalts      = np.unique(dat['MOONALT'].quantity)
moonfracs     = np.unique(dat['MOONFRAC'].quantity)

fits = []

for moonalt in moonalts[::-1]:
  for moonfrac in moonfracs[::-1]:
    ##  Fit for the moon separation dependence using 2-paramters, which will be dependent on moonalt, moonfrac and airmass. 
    sample = dat[(dat['MOONALT'] == moonalt) & (dat['MOONFRAC'] == moonfrac)]

    res    = minimize(diff2, np.array([10., 1.e-3]), args=([sample['MOONSEP'], sample['EXPFAC'], moonmodel]), method='Nelder-Mead', tol=1e-6)

    print(moonalt, moonfrac, res.x)

    ##  Sanity check plots on good fit. 
    ##  
    pl.plot(sample['MOONSEP'], sample['EXPFAC'], 'k-')   
    pl.plot(sample['MOONSEP'], np.clip(moonmodel(sample['MOONSEP'], res.x), 1., None), 'c--')

    fits.append([moonalt, moonfrac, res.x[0], res.x[1]])
        
pl.xlabel('MOON SEP')
pl.ylabel('EXPFAC')

pl.show()

##
pl.clf

fits = np.array(fits)
_abs = np.arange(0.1, 90., 1.)

##  For the two parameter model for dependence on moon separation above, the dependence of these 
##  parameters on moonalt. is given by Z0 & Z1 for fixed moonfrac and airmass.  
def Z0(moonalt, params):
    ##  To be clipped below at unity.
    return  1. + params[0] * np.log(moonalt / 0.6)

def Z1(moonalt, params):
    ##  Tends to zero at moonalt of 0.0    
    return       params[0] * np.log(moonalt / 0.6)

## 
for index in [2, 3]:
  pl.clf()

  ##  Which parameter to model now, [2, 3], index to fits. 
  result = []

  ##  Find the best fitting paramter common to both Z0 and Z1 that models the dependence on moon frac., for given airmass.  
  for moonfrac in moonfracs[::-1]:
    ##  Dependence on moonalt for fixed moonfrac. 
    sample = fits[fits[:,1] == moonfrac]

    if index == 2:
      res    = minimize(diff2, np.array([2.6]), args=([sample[sample[:,0] > 0.0][:,0], sample[sample[:,0] > 0.0][:,index],  Z0]), method='Nelder-Mead', tol=1e-6)
      pl.plot(_abs, np.clip(Z0(_abs, res.x), 1., None), 'k-')   

    if index == 3:
      res   = minimize(diff2, np.array([0.002]), args=([sample[sample[:,0] > 0.0][:,0], sample[sample[:,0] > 0.0][:,index], Z1]), method='Nelder-Mead', tol=1e-6)
      pl.plot(_abs, np.clip(Z1(_abs, res.x), 0.0, None), 'k-')   

    print(res.x)

    ##  Sanity check on good fit. 
    pl.plot(sample[:,0], sample[:,index], '^', markersize=4)

    result.append([moonfrac, res.x[0]])

  pl.xlabel('MOON ALT.')
  pl.ylabel('Z')

  pl.title('AIRMASS {:.1f} and INDEX {:d}'.format(airmass, index))
    
  pl.show()

  result = np.array(result)

  ##  pl.plot(result[:,0], result[:,1])
  ##  pl.title('AIRMASS {:.1f} and INDEX {:d}'.format(airmass, index))
  ##  pl.show()

  ##  Save data later to become Y0 and Y1, i.e. common parameter to Z0 & Z1, which depends on moon frac. and airmass. 
  np.savetxt('../../dat/moonfit_{:.1f}_{:d}.txt'.format(airmass, index), result)
