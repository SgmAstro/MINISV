import numpy as np


def diff2(params, args):
    _abs    = args[0]
    data    = args[1]
    model   = args[2]

    return  np.sum((data - model(_abs, params))**2.)
