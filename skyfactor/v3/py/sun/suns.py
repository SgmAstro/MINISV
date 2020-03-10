import  sys
import  time
import  specsim
import  pylab                         as      pl
import  numpy                         as      np
import  matplotlib.pyplot             as      plt

from    feasibgs.skymodel             import  Isky_newKS_twi, _cI_twi
from    scipy.interpolate             import  interp1d
from    scipy.optimize                import  minimize
from    diff2                         import  diff2
from    astropy.table                 import  Table
from    indices                       import  index_model


##  Airmass to run provided by command line.                                                                                                                                               
all_airmass = [1.0, 1.2, 1.4, 1.6, 2.0, 2.4]
airmass     = all_airmass[np.int(sys.argv[1])]

def sunmodel(sunsep, params):
  ##  One parameter fit for the dependence on sunsep for fixed sunalt and airmass.                                                                                                                            
  _params    = np.zeros(4)
  
  _params[0] = params[0]
  _params[1] = 2.25415081e-03
  _params[2] = 1.2356e+02
  _params[3] = index_model(airmass, np.array([2.04697753e+00, -3.77255962e-02, -2.08474058e-04]))
  
  return  np.ones_like(sunsep) + _params[0] + _params[1] * np.abs(sunsep - _params[2]) ** _params[3]

##  Sunalt, sunsep, airmass
dat  = np.loadtxt('../../dat/suns_{:.1f}.txt'.format(airmass))
dat  = Table(dat, names=['SUNALT', 'SUNSEP', 'AIRMASS', '4000', '6700', '7160', '8070', '9060'])
dat  = dat[dat['SUNALT'] <= -14.]

##  Get normalisation.                                                                                                                                                                                          
norm = np.loadtxt('../../dat/sun_norm.txt')
norm = Table(norm, names=['SUNALT', 'SUNSEP', 'AIRMASS', '4000', '6700', '7160', '8070', '9060'])

##  Normalisation.                                                                                                                                                                                          
for band in ['4000', '6700', '7160', '8070', '9060']:
  dat[band] /= norm[band]
  dat[band]  = np.clip(dat[band], 1.0, None)

##  Choose a band to define the expfac.                                                                                                                                                                         
dat['EXPFAC'] = dat['4000']

assert np.all(dat['EXPFAC'] >= 1.0)

## 
print(dat)

##                                                                                                                                                                                                                                     
pl.clf()

abscissae = np.arange(0.0, 180., 0.1)
indices   = []
fits      = []

for sunalt in np.unique(dat['SUNALT'].quantity):                                                                                                                                                                                     
  ##  np.array([1.2, 2.25e-3, 124, 1.5])
  res     = minimize(diff2, np.array([1.2]), args=([dat[dat['SUNALT'] == sunalt]['SUNSEP'].quantity, dat[dat['SUNALT'] == sunalt]['EXPFAC'].quantity, sunmodel]), method='Nelder-Mead', tol=1e-6)                       
  pl.plot(dat[dat['SUNALT'] == sunalt]['SUNSEP'], dat[dat['SUNALT'] == sunalt]['EXPFAC'], 'k.')
  pl.plot(abscissae, np.clip(sunmodel(abscissae, res.x), 1.0, None), 'k-')                                                                                                     
    
  ##  For determing the power-law index dependence on airmass. 
  fits.append([sunalt.value, res.x[0]])                                                                                                                                                                                                 
  
  print(sunalt, res.x)

  ##  indices.append([airmass, sunalt, res.x[1]])

pl.xlabel('SUN SEP.')
pl.ylabel('EXP FAC.')
pl.title('AIRMASS {}'.format(airmass))

pl.show()

##  np.savetxt('../../dat/indices_{:.2f}.txt'.format(airmass), np.array(indices))

##  
pl.clf()

fits = np.array(fits)
                                                                                                                                                                                                                                  
pl.plot(fits[:,0], fits[:,1])                                                                                                                                                                                                                                                                                                                                                                                                                                               
pl.show()                                                                                                                                                                                                                               

np.savetxt('../../dat/sunalt_{:.1f}.txt'.format(airmass), fits, fmt='%.6lf')
