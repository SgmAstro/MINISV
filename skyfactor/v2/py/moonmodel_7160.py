import  numpy              as      np
import  pylab              as      pl
import  matplotlib.pyplot  as      plt

from    astropy.table      import  Table


def X0(airmass):  
  ##  v1:         return  -0.01363498  + 0.15419213 * airmass ** 0.75 
  ##  v2 (7160):  
  return  -0.01474017 + 0.07491027 * airmass ** 0.75

def X1(airmass):
  ##  v1:         return  -7.68546710e-06 + 5.75539255e-05 * airmass ** 0.75
  ##  v2 (7160):  
  return  -6.84340500e-06 + 2.79721586e-05 * airmass ** 0.75

def Y0(moonfrac, airmass):
  ##  Dependence on moonfrac for a given airmass.                                                                                                                                                
  _X0 = X0(airmass)

  Y   = moonfrac - 0.1

  return  _X0 * np.exp(Y / 0.29) - _X0

def Y1(moonfrac, airmass):
  ##  Dependence on moonfrac for a given airmass.                                                                                                                                               
  _X1 = X1(airmass)
  Y   = moonfrac - 0.1

  return  _X1 * np.exp(Y / 0.29) - _X1

def Z0(moonalt, moonfrac, airmass):
  ##  Dependence on moonalt for fixed moonfrac and airmass.                                                                                                                                      
  ##  Note:  np.clip(Z0(_abs, res.x), 1., None)                                                                                                                                                  
  _Y0    = Y0(moonfrac, airmass)

  result =  1. + _Y0 * np.log(moonalt / 0.07)
  
  return  result
  
def Z1(moonalt, moonfrac, airmass):
  ##  Dependence on moonalt for fixed moonfrac and airmass.                                                                                                                                      
  ##  Note:  np.clip(Y1(_abs, res.x), 0.0, None)                                                                                                                                                 
  _Y1 = Y1(moonfrac, airmass)

  result = _Y1 * np.log(moonalt / 0.07)
  
  return  result
  
def moonmodel(moonsep, moonalt, moonfrac, airmass):
  ##  Dependence on separation for fixed moon alt., moon frac. and airmass.                                                                                                                      
  ##  Note:  np.clip(model(sample['MOONSEP'], res.x), 1., None)                                                                                                                                  
  if moonalt < 0.5:
    return  np.ones_like(moonsep)

  _Z0    = np.clip(Z0(moonalt, moonfrac, airmass), 1., None)
  _Z1    = np.clip(Z1(moonalt, moonfrac, airmass), 0., None)

  result = _Z0 + _Z1 * np.abs(moonsep - 96.) ** 1.86328519
  result = np.clip(result, 1., None)
  
  return  result


if __name__ == '__main__':
  cols       = ['AIRMASS', 'MOONFRAC', 'MOONALT', 'MOONSEP', '4000', '6700', '7160', '8070', '9060']
  _airmass   = [1.0, 1.6, 2.0, 2.4, 2.8]

  ##  Get normalisation.                                                                                                                                                                                                            
  norm       = np.loadtxt('../dat/moon_norm.txt')
  norm       = Table(norm, names=cols)

  ##  
  colors     = plt.rcParams['axes.prop_cycle'].by_key()['color']

  ##  
  fig, axarr = plt.subplots(nrows=1, ncols=5, figsize=(20, 4))

  for i, airm in enumerate(_airmass):  
    dat  =  np.loadtxt('../dat/moons_{:.1f}.txt'.format(airm))
    dat  =  Table(dat, names=cols)

    for band in ['4000', '6700', '7160', '8070', '9060']:
      dat[band] /= norm[band]
      dat[band]  = np.clip(dat[band], 1.0, None)

    print(dat)

    ##  Choose a band to define the expfac.                                                                                                                                                                              
    dat['EXPFAC'] = dat['7160']

    assert np.all(dat['EXPFAC'] >= 1.0)

    dat  =  dat[dat['MOONFRAC'] > 0.6]

    _abs =  np.arange(0.0, 40., 0.01)
  
    axarr[i].plot(_abs, _abs, 'k-')

    analytic = []
    
    for x in dat:
      analytic.append(moonmodel(x['MOONSEP'], x['MOONALT'], x['MOONFRAC'], x['AIRMASS']))
    
    ##  print(x['MOONSEP'], x['MOONALT'], x['MOONFRAC'], x['AIRMASS'], x['EXPFAC'], result)

    im = axarr[i].scatter(dat['EXPFAC'], np.array(analytic), s=1., c=dat['MOONALT'], vmin=0., vmax=90.)

    ##  axarr[i].axis('equal')
    axarr[i].set_xlabel('CHH EXPFAC')

    if i == 0:
      axarr[i].set_ylabel('ANALYTIC')

    axarr[i].set_title('AIRMASS: {:.1f}'.format(airm))

    axarr[i].set_xlim((1., 25.))
    axarr[i].set_ylim((1., 25.))

  fig.colorbar(im, orientation='vertical', label=r'MOON ALT.')

  ## 
  pl.show()
    
