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
all_airmass = [1.0, 1.2, 1.6, 2.0, 2.4, 2.8, 3.0]
airmass     = all_airmass[np.int(sys.argv[1])]

def sunmodel(sunsep, param):
  ##  One parameter fit for the dependence on sunsep for fixed sunalt and airmass.                                                                                                                             
  _params    = np.zeros(4)

  _params[0] = param
  _params[1] = 3.79912194e-04
  _params[2] = 1.2356e+02
  _params[3] = index_model(airmass, np.array([2.00688895, -0.01108773, -0.00291262]))

  return  np.ones_like(sunsep) + _params[0] + _params[1] * np.abs(sunsep - _params[2]) ** _params[3]

##  Sunalt, sunsep, airmass
dat  = np.loadtxt('../../dat/suns_{:.1f}.txt'.format(airmass))
dat  = Table(dat, names=['SUNALT', 'SUNSEP', 'AIRMASS', '4000', '6700', '7160', '8070', '9060'])

##  Get normalisation.                                                                                                                                                                                          
norm = np.loadtxt('../../dat/sun_norm.txt')
norm = Table(norm, names=['SUNALT', 'SUNSEP', 'AIRMASS', '4000', '6700', '7160', '8070', '9060'])

##  Normalisation.                                                                                                                                                                                          
for band in ['4000', '6700', '7160', '8070', '9060']:
  dat[band] /= norm[band]
  dat[band]  = np.clip(dat[band], 1.0, None)

##  Choose a band to define the expfac.                                                                                                                                                                         
dat['EXPFAC'] = dat['7160']

assert np.all(dat['EXPFAC'] >= 1.0)

## 
print(dat)

##                                                                                                                                                                                                                                     
pl.clf()

abscissae = np.arange(0.0, 180., 0.1)
indices   = []
fits      = []

for sunalt in np.unique(dat['SUNALT'].quantity):                                                                                                                                                                                     
  ##  np.array([1.2, 2.0])
  res     = minimize(diff2, np.array([1.2]), args=([dat[dat['SUNALT'] == sunalt]['SUNSEP'].quantity, dat[dat['SUNALT'] == sunalt]['EXPFAC'].quantity, sunmodel]), method='Nelder-Mead', tol=1e-6)                                      

  pl.plot(dat[dat['SUNALT'] == sunalt]['SUNSEP'], dat[dat['SUNALT'] == sunalt]['EXPFAC'], 'k.')
  pl.plot(abscissae, np.clip(sunmodel(abscissae, res.x), 1.0, None), 'k-')                                                                                                     
    
  ##  For determing the power-law index dependence on airmass. 
  fits.append([sunalt.value, res.x[0]])                                                                                                                                                                                                 

  print(sunalt, res.x)

  ##  indices.append([airmass, sunalt, res.x[1]])

pl.xlabel('SUN SEP.')
pl.ylabel('EXP FAC.')

##  pl.show()

##  np.savetxt('../../dat/indices_{:.2f}.txt'.format(airmass), np.array(indices))

##  
pl.clf()

fits = np.array(fits)
                                                                                                                                                                                                                                  
pl.plot(fits[:,0], fits[:,1])                                                                                                                                                                                                                                                                                                                                                                                                                                                   
def model(sunalt, params):                                                                                                                                                                                                                return params[0] + params[1] * np.exp(sunalt / params[2])                                                                                                                                                                           

res = minimize(diff2, np.array([-2.3, 278.5, 3.2]), args=([fits[:,0], fits[:,1], model]), method='Nelder-Mead', tol=1e-6)                                                                                                            
pl.plot(np.arange(-18., -12., 0.1), model(np.arange(-18., -12., 0.1), res.x), 'c-')                                                                                                                                                                                
print(res.x)                                                                                                                                                                                                                             
pl.show()                                                                                                                                                                                                                               

np.savetxt('../../dat/sunalt_{:.1f}.txt'.format(airmass), fits, fmt='%.6lf')
