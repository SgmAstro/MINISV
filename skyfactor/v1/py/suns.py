import time
import specsim
import pylab                        as     pl
import numpy                        as     np
import matplotlib.pyplot            as     plt

from   feasibgs.skymodel            import Isky_newKS_twi, _cI_twi
from   scipy.interpolate            import interp1d
from   desisurvey.etc               import schlegel_twi
from   scipy.optimize               import minimize
from   diff2                        import diff2
from   desisurvey.skymodel.sunmodel import sunmodel


def model(sunsep, param, airmass):
    ##  Dependence on sunsep for fixed sunalt and airmass.                                                                                                                                                                             
    '''                                                                                                                                                                                                                                
    -18.0 -1.55351352692                                                                                                                                                                                                                   -17.0 -1.35187339783                                                                                                                                                                                                               
    -16.0 -0.907573699951                                                                                                                                                                                                              
    -15.0 -0.220613479614                                                                                                                                                                                                                  -14.0  0.709006309509                                                                                                                                                                                                              
    -13.0  1.88128757477                                                                                                                                                                                                               
    -12.0  3.29622745514                                                                                                                                                                                                               
    '''

    _params    = np.zeros(4)

    _params[0] = param
    _params[1] = 1.91579937e-04
    _params[2] = 1.22515026e+02
    _params[3] = 1.00131852 + 1.11553212 / airmass -0.11919901 / airmass**2.

    return np.ones_like(sunsep) + _params[0] + _params[1] * np.abs(sunsep - _params[2])**_params[3]

airmass = 1.6

# sunalt, sunsep, airmass
result  = np.loadtxt('dat/suns_{:.1f}.txt'.format(airmass))
sunalts = result[:,0]
                                                                                                                                                                                                                               
pl.clf()

_xs = np.arange(0., 10., 0.01)                                                                                                                                                                                                       
pl.plot(_xs, _xs, 'k-')

for i, x in enumerate(result):
  # sunsep, airmass, sunalt
  expfac    = sunmodel(x[1], x[2], x[0])
  expfac    = np.clip(expfac, 1., None)                                                                                                                                                                                                
                                                                                                                                                                                                                                       
  pl.plot(x[3], expfac, '^', c='k', markersize=5)                                                                                                                                                                                      
                                                                                                                                                                                                                                       
  #  print(x[1], x[2], x[0], x[3], expfac)
  
pl.xlim(0., 5.)
pl.ylim(0., 5.)

pl.show()
pl.clf()

                                                                                                                                                                                                                                    
##                                                                                                                                                                                                                                    
pl.axhline(y=1.,  xmin=0., xmax=1., alpha=0.5, c='k')                                                                                                                                                                                  
pl.axhline(y=2.5, xmin=0., xmax=1., alpha=1.0, c='k')

pl.plot(result[:,0], schlegel_twi(result[:,0]), 'r-', label='DS')
plt.scatter(result[:,0], result[:,3], c=result[:,1], marker='.', label='CHH')

pl.colorbar(label='Sun separation [deg.]')

pl.yscale('log')

pl.xlim(-18.5, -12.5)                                                                                                                                                                                                                  
pl.ylim(  0.1,   10.)

pl.xlabel('Sun alt. [deg.]')                                                                                                                                                                                                           
pl.ylabel('Twilight factor')                                                                                                                                                                                                          
pl.legend(frameon=False)
pl.show()

##                                                                                                                                                                                                                                     
pl.clf()

pl.axhline(y=1.,  xmin=0., xmax=1., alpha=1.0, c='k')
pl.axhline(y=2.5, xmin=0., xmax=1., alpha=0.2, c='k')

plt.scatter(result[:,1], result[:,3], c=result[:,0], marker='.', label='CHH')

abscissae = np.arange(0.0, 180., 0.1)

fits      = []

for sunalt in sunalts:
    ##                                                                                                                                                                                                                                  
    ##  res = minimize(diff2, np.array([5.]), args=([result[result[:,0] == sunalt][:,1], result[result[:,0] == sunalt][:,3], model, airmass]), method='Nelder-Mead', tol=1e-6)                                                          
    ##  pl.plot(abscissae, np.clip(model(abscissae, res.x, airmass), 1.0, None), 'k-')                                                                                                                                                   
    pl.plot(abscissae, np.clip(sunmodel(abscissae, airmass, sunalt), 1.0, None), 'k-')

    ##  print(sunalt, res.x)                                                                                                                                                                                                            
    ##  fits.append([sunalt, res.x[0]])                                                                                                                                                                                                  
pl.colorbar(label='Sun altitude [deg.]')

pl.yscale('log')

pl.xlim(0.0, 180.)
pl.ylim(0.1,  10.)

pl.xlabel('Sun sep. [deg.]')
pl.ylabel('Twilight factor')

pl.legend(frameon=False)

pl.show()

exit(1)

pl.clf()

fits = np.array(fits)

'''                                                                                                                                                                                                                                  
pl.plot(fits[:,0], fits[:,1])                                                                                                                                                                                                                                                                                                                                                                                                                                                   
def model(sunalt, params, dummy):                                                                                                                                                                                                       
    ##                                                                                                                                                                                                                                  
    return params[0] + params[1] * np.exp(-sunalt / params[2])                                                                                                                                                                           
res = minimize(diff2, np.array([1., 0.1, 10.]), args=([fits[:,0], fits[:,1], model, None]), method='Nelder-Mead', tol=1e-6)                                                                                                            
pl.plot(fits[:,0], model(fits[:,0], res.x, None), 'c--')                                                                                                                                                                                

print(res.x)                                                                                                                                                                                                                             
pl.show()                                                                                                                                                                                                                               
'''
##  print(airmass)

np.savetxt('sunalt_{:.1f}.txt'.format(airmass), fits, fmt='%.6lf')
