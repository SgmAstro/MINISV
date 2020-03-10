import time
import specsim
import pylab                    as     pl
import numpy                    as     np
import matplotlib.pyplot        as     plt

from   feasibgs.skymodel        import Isky_newKS_twi, _cI_twi
from   scipy.interpolate        import interp1d
from   desisurvey.etc           import schlegel_twi
from   scipy.optimize           import minimize
from   diff2                    import diff2


dwave       = 10.

##  Run for each airmass.                                                                                                                                                                                                              
all_airmass = [1.0, 1.2, 1.6, 2.0, 2.4, 2.8, 3.0]

sunalts     = np.arange(-18., -11., 1.0)
sunseps     = np.arange( 10., 180., 5.0)

bandpasses  = {'4000': [4000., 5350.], '6700': [6700., 6820.], '7160': [7160., 7220.], '8070': [8070., 8260.], '9060': [9060., 9280.]}

##  
nruns       = len(sunalts) * len(sunseps)
t0          = time.time()

for airmass in all_airmass:
  result    = []
  count     =  0

  for sunalt in sunalts:
    for sunsep in sunseps:
      row   = []

      ##  Note:  not valid below -20, or above -14, in alt. 
      wave, twi = _cI_twi(sunalt, sunsep, airmass)

      ##  Angstroms sampled at every 10.
      wave      = 10. * wave
      twi       = twi / np.pi
      
      for band in bandpasses.keys():
        lolam, hilam = bandpasses[band][0], bandpasses[band][1]

        sample       = twi[(lolam <= wave) & (wave <= hilam)]
        flux         = np.sum(sample) * dwave
        
        row.append(flux)

      toappend  = [sunalt, sunsep, airmass] + row

      if np.all(np.array(row) > 0.0):            
        result.append(toappend)

      count  += 1

      print('AIRMASS:  {};  SUNALT:  {};  SUNSEP:  {}'.format(airmass, sunalt, sunsep))


  ##
  diff          = time.time() - t0
  result        = np.array(result)

  np.savetxt('../dat/suns_{:.1f}.txt'.format(airmass), result, fmt='%.6lf')
