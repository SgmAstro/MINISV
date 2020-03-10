import time
import specsim
import pylab                    as     pl
import numpy                    as     np
import matplotlib.pyplot        as     plt

from   feasibgs.skymodel        import Isky_newKS_twi, _cI_twi
from   scipy.interpolate        import interp1d
from   desisurvey.etc           import schlegel_twi
from   moonmodel                import moonmodel


dwave      = 0.1

##  Isky_newKS_twi(airmass, moonill, moonalt, moonsep, sunalt, sunsep)
## 
##  returns  AA, $10^{-17} erg/cm2/s/arcsec2
wave, dark = Isky_newKS_twi(1.0, 0.0, -90., 120., -90., 120.)

##  Defines expfac = 1.0
dark_g     = np.sum(dark[(4000. <= wave.value) & (wave.value <= 5500.)]) * dwave

##  Run for each airmass. 
airmass    = 3.0

moonfracs  = np.arange( 0.1,  1.0, 0.1)
moonalts   = np.arange(  0.,  90., 5.0)
moonseps   = np.arange( 50., 140., 5.0)

nruns      = len(moonfracs) * len(moonalts) * len(moonseps)

t0         = time.time()

result     = []
count      =  0.

for moonfrac in moonfracs:
  for moonalt in moonalts:
    for moonsep in moonseps:
      wave, sky = Isky_newKS_twi(airmass, moonfrac, moonalt, moonsep, -90., 120.)
      flux_g    = np.sum(sky[(4000. <= wave.value) & (wave.value <= 5500.)]) * dwave
      expfac    = flux_g / dark_g

      complete  = 100. * count / nruns

      print(complete, airmass, moonfrac, moonalt, moonsep, expfac)
        
      result.append([airmass, moonfrac, moonalt, moonsep, expfac])

      count     += 1
      
    output = np.array(result)
    np.savetxt('dat/moons_{:.1f}.txt'.format(airmass), output, fmt='%.6lf')

##
result = np.array(result)
        
## 
diff   = time.time() - t0

print(result)

np.savetxt('dat/moons_{:.1f}.txt'.format(airmass), result, fmt='%.6lf')
