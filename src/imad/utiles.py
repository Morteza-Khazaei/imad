import platform
import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer


import warnings
warnings.filterwarnings("ignore")

if platform.system() == 'Windows':
    lib = ctypes.cdll.LoadLibrary('prov_means.dll')
elif platform.system() == 'Linux':
    lib = ctypes.cdll.LoadLibrary('libprov_means.so')
elif platform.system() == 'Darwin':
    lib = ctypes.cdll.LoadLibrary('libprov_means.dylib')
provmeans = lib.provmeans
provmeans.restype = None
c_double_p = ctypes.POINTER(ctypes.c_double)
provmeans.argtypes = [ndpointer(np.float64),
                      ndpointer(np.float64),
                      ctypes.c_int,
                      ctypes.c_int,
                      c_double_p,
                      ndpointer(np.float64),
                      ndpointer(np.float64)]

# -----------------
# provisional means
# -----------------

class CPM(object):
    '''Provisional means algorithm'''
    def __init__(self, N):
        self.mn = np.zeros(N)
        self.cov = np.zeros((N, N))
        self.sw = 0.0000001

    def update(self, Xs, Ws=None):
        n, N = np.shape(Xs)
        if Ws is None:
            Ws = np.ones(n)
        sw = ctypes.c_double(self.sw)
        mn = self.mn
        cov = self.cov
        provmeans(Xs, Ws, N, n, ctypes.byref(sw), mn, cov)
        self.sw = sw.value
        self.mn = mn
        self.cov = cov

    def covariance(self):
        c = np.mat(self.cov/(self.sw-1.0))
        d = np.diag(np.diag(c))
        return c + c.T - d

    def means(self):
        return self.mn