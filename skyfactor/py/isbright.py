import numpy             as      np
import pylab             as      pl
import matplotlib.pyplot as      plt

from   moonmodel_7160    import  moonmodel
from   sunmodel          import  sunmodel


moonalts       = np.arange(0.0, 90.,   1.)
moonfracs      = np.arange( 0.,  1., 0.01)

all_airmass    = np.array([1.0, 1.6, 2.0, 2.4, 2.8])
moonseps       = np.array([50., 96.])

X, Y           = np.meshgrid(moonalts, moonfracs)

##                                                                                                                                                                                                                                   
fig, axarr     = plt.subplots(nrows=1, ncols=5, figsize=(20, 4))

for k, airmass in enumerate(all_airmass):
  for moonsep, alpha in zip(moonseps, [1., 0.5]):
    Z          = np.ones_like(X)

    print('Solving for airmass {} and moon sep. {}'.format(airmass, moonsep))

    for i, moonfrac in enumerate(moonfracs):
      for j, moonalt in enumerate(moonalts):
        Z[i,j] = moonmodel(moonsep, moonalt, moonfrac, airmass)

    CS = axarr[k].contour(X, Y, Z, levels=[2.5, 5.0, 7.5, 10.0], alpha=alpha)
    
    axarr[k].clabel(CS, [2.5, 5.0, 7.5, 10.0], inline=1, fontsize=5)

    axarr[k].set_ylim(0.,  1.)
    axarr[k].set_xlim(0., 90.)

    axarr[k].set_title('Airmass:  {}'.format(airmass))

##
axarr[0].set_ylabel('Moon frac.')

for i in range(5):
  axarr[i].set_xlabel('Moon alt.')

pl.show()                                                                                                                                                                                                                        

