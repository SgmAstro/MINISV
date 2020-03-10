import numpy           as np
import pylab           as pl

from   scipy.optimize  import minimize
from   astropy.table   import Table


def index_model(airmass, params):
    ##  [ 2.00688895 -0.01108773 -0.00291262]
    return  params[0] + params[1] * airmass + params[2] * (airmass - 1.5)**2.

if __name__ == '__main__':
    from   diff2           import diff2

    all_airmass = [1.0, 1.2, 1.6, 2.0, 2.4]

    dat = []

    for airmass in all_airmass:
      ##  Fit for the dependence of the powerlaw index in sun separation with airmass.  See sunmodel. 
      sunalt, __, indices = np.loadtxt('../../dat/indices_{:.2f}.txt'.format(airmass), unpack=True)
      
      mean = np.mean(indices)
      std  =  np.std(indices)

      dat.append([airmass, mean, std])
  
    dat = Table(np.array(dat), names=['AIRMASS', 'MEAN', 'STD'])
    
    print(dat)
    
    ##  Effectively no standard deviation across altitude. 
    pl.plot(dat['AIRMASS'], dat['MEAN'])

    res = minimize(diff2, np.array([1., 0.1, 0.1]), args=([dat['AIRMASS'], dat['MEAN'], index_model]), method='Nelder-Mead', tol=1e-6)                                                                                                    
    print(res.x)

    abs = np.arange(1., 3.0, 0.01)
    pl.plot(abs, index_model(abs, res.x), 'c--')                                                                                                                                                                           

    pl.show()
