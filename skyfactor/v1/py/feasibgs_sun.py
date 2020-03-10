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


dwave      = 10.

##  Isky_newKS_twi(airmass, moonill, moonalt, moonsep, sunalt, sunsep)
## 
##  returns  AA, $10^{-17} erg/cm2/s/arcsec2
wave, dark = _cI_twi(-19., 50., 1.0)

assert np.all(dark > 0.0)

##  Angstroms sampled at every 10.                                                                                                                                                                         
wave       = 10. * wave
dark       = dark / np.pi

dark_g     = np.sum(dark[(4000. <= wave) & (wave <= 5500.)]) * dwave

airmass    = 2.8

sunalts    = np.arange(-18., -11., 1.0)
sunseps    = np.arange( 10., 180., 5.0)

nruns      = len(sunalts) * len(sunseps)

t0         = time.time()

result     = []

for sunalt in sunalts:
  for sunsep in sunseps:
    ##  Note:  not valid below -20 in alt. 
    _w_twi, I_twi = _cI_twi(sunalt, sunsep, airmass)

    ##  Angstroms sampled at every 10.
    w_twi         = 10. * _w_twi
    sky           = I_twi / np.pi
        
    ##  print(w_twi)
    ##  print(sky)
            
    if not np.all(I_twi > 0.0):
      pass
            ##  print(sunalt, sunsep, airmass)
            
    else:
      flux_g    = np.sum(sky[(4000. <= w_twi) & (w_twi <= 5500.)]) * 10.0
      expfac    = flux_g / dark_g
                
      result.append([sunalt, sunsep, airmass, expfac])

##
diff          = time.time() - t0

result        = np.array(result)

norm          = np.median(result[result[:,0] == -18.][:,3])
result[:,3]  /= norm

##  Enforce clip on dark time. 
result[result[:,3] < 1.0][:,3] = 1.0

np.savetxt('dat/suns_{:.1f}.txt'.format(airmass), result, fmt='%.6lf')
