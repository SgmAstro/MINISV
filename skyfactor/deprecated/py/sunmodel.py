import numpy          as     np
import pylab          as     pl

from   astropy.table  import Table
from   desisurvey.etc import schlegel_twi


def Y0(airmass, params):
    ##  [7.05246892e+02   1.82683325e-01]                                                                                                                                                                                              
    return  params[0] * np.exp(-airmass / params[1])

def Y1(airmass, params):
    ##  [1.84684173  1.92248671]                                                                                                                                                                                                       
    return  params[0] + params[1] / airmass ** 2.5

def Z(sunalt, params):
    ##  Dependence on sunalt. for a given airmass.                                                                                                                                                                                     
    return  np.exp(5. + sunalt / params[1]) - params[0]

def sunmodel(sunsep, airmass, sunalt):
    #  Airmass dependence of sunalt shape;  See sunalts.py;                                                                                                                                                                            
    if sunalt < -18.:
        return  np.ones_like(sunsep)

    if sunalt > -14:
        ##  Invalid above -14 deg. 
        return  1.e99

    _Y0    =  Y0(airmass, np.array([7.05246892e+02, 1.82683325e-01]))
    _Y1    =  Y1(airmass, np.array([1.84684173, 1.92248671]))

    # One parameter model for dependence on sun alt. for given airmass, Z.                                                                                                                                                             
    _Z     =  Z(sunalt,  np.array([_Y0, _Y1]))

    index  =  1.00131852 + 1.11553212 / airmass -0.11919901 / airmass**2.

    # Dependence on sunsep for fixed sunalt and airmass.                                                                                                                                                                               
    result =  np.ones_like(sunsep) + _Z + 1.91579937e-04 * np.abs(sunsep - 1.22515026e+02) ** index
    result =  np.clip(result, 1., None) 

    return  result
    
  
if __name__ == '__main__':
  cols      = ['SUNALT', 'SUNSEP', 'AIRMASS', 'EXPFAC']

  _airmass  = [1.0, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8]

  _abs      =  np.arange(0.0, 10., 0.01)
  pl.plot(_abs, _abs, 'k-', alpha=0.2)

  for airm in _airmass:
    dat  =  np.loadtxt('dat/suns_{:.1f}.txt'.format(airm))
    dat  =  Table(dat, names=cols)

    for x in dat:
      result = sunmodel(x['SUNSEP'], x['AIRMASS'], x['SUNALT'])

      ##  print(x['MOONSEP'], x['MOONALT'], x['MOONFRAC'], x['AIRMASS'], x['EXPFAC'], result)                                                                                                                                          

      pl.plot(x['EXPFAC'], result, 'k.', markersize=1, alpha=0.3)

  pl.xlim(0., 5.)
  pl.ylim(0., 5.)
  pl.show()
  pl.clf()
  
  for airm in _airmass:  
    dat  =  np.loadtxt('dat/suns_{:.1f}.txt'.format(airm))
    dat  =  Table(dat, names=cols)
    dat  =  dat[dat['SUNALT'] <= -14.]

    for x in dat:
      result = sunmodel(x['SUNSEP'], x['AIRMASS'], x['SUNALT'])
      
      ##  print(x['MOONSEP'], x['MOONALT'], x['MOONFRAC'], x['AIRMASS'], x['EXPFAC'], result)

      pl.plot(x['SUNALT'],  x['EXPFAC'], 'k.', markersize=1, alpha=0.3)
      #pl.plot(x['SUNALT'] + 0.2, result, 'b.', markersize=1, alpha=0.3)

  pl.plot(dat['SUNALT'].quantity, schlegel_twi(dat['SUNALT'].quantity), 'r-', markersize=1)
      
  pl.axhline(y=1.0, xmin=0, xmax=1, c='k', alpha=0.5)
  pl.axhline(y=2.5, xmin=0, xmax=1, c='k', alpha=1.0)
  
  pl.xlim(-18., -13.5)
  pl.ylim(  0.,  20.0)
  
  pl.yscale('log')
  pl.show()
