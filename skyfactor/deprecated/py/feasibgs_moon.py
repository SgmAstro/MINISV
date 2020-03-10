import time
import specsim
import pylab                    as     pl
import numpy                    as     np
import matplotlib.pyplot        as     plt

from   feasibgs.skymodel        import Isky_newKS_twi, _cI_twi
from   scipy.interpolate        import interp1d
from   desisurvey.etc           import schlegel_twi
from   moonmodel                import moonmodel


## 
##  https://desi.lbl.gov/DocDB/cgi-bin/private/RetrieveFile?docid=3985;filename=Sky%20Camera%20FDR%20v1.pdf;version=1
##  https://desi.lbl.gov/DocDB/cgi-bin/private/RetrieveFile?docid=3422;filename=SkyCameraPDR-v2.pdf;version=2
##

dwave       = 0.1

##  Run for each airmass. 
all_airmass = [1.0, 1.2, 1.6, 2.0, 2.4, 2.8, 3.0]

moonfracs   = np.arange( 0.1,  1.0, 0.1)
moonalts    = np.arange(  0.,  90., 5.0)
moonseps    = np.arange( 50., 140., 5.0)

bandpasses  = {'4000': [4000., 5350.], '6700': [6700., 6820.], '7160': [7160., 7220.], '8070': [8070., 8260.], '9060': [9060., 9280.]}

nruns       = len(moonfracs) * len(moonalts) * len(moonseps)

t0          = time.time()

result      = []
count       =  0.
 
for airmass in all_airmass:
  for moonfrac in moonfracs:
    for moonalt in moonalts:
      for moonsep in moonseps:
        row  = []

        complete       = 100. * count / nruns

        for band in bandpasses.keys():
          lolam, hilam = bandpasses[band][0], bandpasses[band][1]

          wave, sky    = Isky_newKS_twi(airmass, moonfrac, moonalt, moonsep, -90., 120.)
          flux         = np.sum(sky[(lolam <= wave.value) & (wave.value <= hilam)]) * dwave
         
          row.append(flux)
      
        print(complete, airmass, moonfrac, moonalt, moonsep, expfac)

        count     += 1

        result.append([airmass, moonfrac, moonalt, moonsep] + row)


  output = np.array(result)
  np.savetxt('dat/moons_{:.1f}.txt'.format(airmass), output, fmt='%.6lf')

##
result = np.array(result)
        
## 
diff   = time.time() - t0

print(result)

np.savetxt('dat/moons_{:.1f}.txt'.format(airmass), result, fmt='%.6lf')
