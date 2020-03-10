import numpy          as np
import pylab          as pl

from   scipy.optimize import minimize
from   diff2          import diff2


def model(airmass, params):
    ##  [ 1.00131852  1.11553212 -0.11919901]    
    return params[0] + params[1] / airmass + params[2] / airmass**2.


##  Fit for the dependence of the powerlaw in sun separation.  See sunmodel. 
dat = np.loadtxt('indices.txt')

pl.plot(dat[:,0], dat[:,1])

res = minimize(diff2, np.array([1., 0.1, 0.01]), args=([dat[:,0], dat[:,1], model]), method='Nelder-Mead', tol=1e-6)                                                                                                       

pl.plot(dat[:,0], model(dat[:,0], res.x), 'c--')                                                                                                                                                                           
                                                                                                                                                                                                                            
print(res.x)

x = np.arange(1.0, 6.0, 0.1)
p = np.array([1.00131852, 1.11553212, -0.11919901])

pl.plot(x, model(x, p))
pl.show()
