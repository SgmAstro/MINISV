import  numpy              as      np
import  pylab              as      pl
import  matplotlib.pyplot  as      plt

from    astropy.table      import  Table


def index_model(airmass, params):
  ##  [ 2.00688895 -0.01108773 -0.00291262]                                                                                                                                              
  return  params[0] + params[1] * airmass + params[2] * (airmass - 1.5)**2.

def Y0(airmass, params):
  ##  [35.38645069, -100.59010158, 10.14106723].                                                                                                                                           
  return  params[0] + params[1] * (airmass - 1.) + params[2] * (airmass - 1.7)**2.

def Y1(airmass, params):
  ##  [5.82270117, -0.52464915].                                                                                                                                                            
  return  params[0] + params[1] * (airmass - 1.)

def Z(sunalt, params):
  ##  Dependence on sunalt. for a given airmass.                                                                                                                                                                              
  return  np.exp(9.45 + sunalt / params[1]) - params[0]

def sunmodel(sunsep, airmass, sunalt):
  #  Airmass dependence of sunalt shape;  See sunalts.py;                                                                                                                                                                            
  if sunalt < -18.:
    return  np.ones_like(sunsep)

  if (sunalt > -14) | (airmass > 2.0):
    ##  Invalid above -14 deg. 
    return  1.e99

  _Y0    =  Y0(airmass, np.array([10.97556349, -94.02393911, 8.89981454]))
  _Y1    =  Y1(airmass, np.array([2.40058904, -0.07732892]))
  
  # One parameter model for dependence on sun alt. for given airmass, Z.                                                                                                                                                             
  _Z     =  Z(sunalt,   np.array([_Y0, _Y1]))

  index  =  index_model(airmass, np.array([2.04697753e+00, -3.77255962e-02, -2.08474058e-04]))

  # Dependence on sunsep for fixed sunalt and airmass.                                                                                                                                                                               
  result =  np.ones_like(sunsep) + _Z + 2.25415081e-03 * np.abs(sunsep - 1.2356e+02) ** index
  result =  np.clip(result, 1., None) 

  return  result
    
  
if __name__ == '__main__':
  _airmass  = [1.0, 1.2, 1.4, 1.6, 2.0]

  _abs      = np.arange(0.0, 140., 0.01)

  ##  Get normalisation.                                                                                                                                                                                                               
  norm      = np.loadtxt('../dat/sun_norm.txt')
  norm      = Table(norm, names=['SUNALT', 'SUNSEP', 'AIRMASS', '4000', '6700', '7160', '8070', '9060'])

  ##                                                                                                                                                                                                                                  
  colors     = plt.rcParams['axes.prop_cycle'].by_key()['color']

  ##                                                                                                                                                                                                                                  
  fig, axarr = plt.subplots(nrows=1, ncols=5, figsize=(20, 4))

  for i, airm in enumerate(_airmass):
    dat  = np.loadtxt('../dat/suns_{:.1f}.txt'.format(airm))
    dat  = Table(dat, names=['SUNALT', 'SUNSEP', 'AIRMASS', '4000', '6700', '7160', '8070', '9060'])
    dat  = dat[dat['SUNALT'] <= -14.]

    for band in ['4000', '6700', '7160', '8070', '9060']:
      dat[band] /= norm[band]
      dat[band]  = np.clip(dat[band], 1.0, None)

    ##  Choose a band to define the expfac.                                                                                                                                                                                  
    dat['EXPFAC'] = dat['4000']

    assert  np.all(dat['EXPFAC'] >= 1.0)
 
    ##                                                                                                                                                                                                                    
    analytic = []

    for x in dat:
      analytic.append(sunmodel(x['SUNSEP'], x['AIRMASS'], x['SUNALT']))

    analytic        = np.array(analytic)
    dat['ANALYTIC'] = analytic 
    
    ##
    axarr[i].plot(_abs, _abs, 'k-', alpha=.2)

    ##  im = axarr[i].scatter(dat['EXPFAC'], dat['ANALYTIC'], s=1., vmin=-20., vmax=-12., c=dat['SUNALT'])
    im = axarr[i].scatter(dat['EXPFAC'], dat['ANALYTIC'], s=1., vmin=10., vmax=165., c=dat['SUNSEP'])

    axarr[i].set_xlabel('CHH EXPFAC')

    if i == 0:
      axarr[i].set_ylabel('ANALYTIC')

    axarr[i].set_title('AIRMASS: {:.1f}'.format(airm))

    axarr[i].set_xlim((1., 35.))
    axarr[i].set_ylim((1., 35.))

    dat.sort('EXPFAC')
    
    print(dat)
  
  fig.colorbar(im, orientation='vertical', label=r'SUN ALT. @ $4000\AA$')
  ##  fig.colorbar(im, orientation='vertical', label=r'SUN SEP.')

  pl.show()
  
