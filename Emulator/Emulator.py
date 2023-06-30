import numpy as np
import scipy as sp
from scipy.interpolate import CubicSpline

class Emulator:
    def __init__(self,axis,matrix,cutoff):
        self.axis   = axis
        self.matrix = matrix
        self.cutoff = cutoff
        self.height,self.width = self.matrix.shape
        self.build()
    
    def normalize(self):
        ''' Normalize the data '''
        self.means = [ np.average(self.matrix[:,i]) for i in range(self.width) ]
        self.sdevs = [ np.std(self.matrix[:,i])     for i in range(self.width) ]
        self.renor = [ [
                                ( self.matrix[i][j]-self.means[j] )/self.sdevs[j] 
                            if self.sdevs[j] != 0 else 0 
                        for j in range(self.width) ] 
                      for i in range(self.height) ]
        self.renor = np.array(self.renor)
    
    def build(self):
        self.normalize()
        self.u,self.s,self.v = np.linalg.svd(self.renor)
        self.splines = [ sp.interpolate.CubicSpline(self.axis,self.u[:,i]) for i in range(self.cutoff) ]
    
    def __call__(self,value):
        rval = np.zeros(self.width)
        for i in range(self.cutoff):
            rval += self.splines[i](value)*self.s[i]*self.v[i]
        rval = np.multiply(rval,self.sdevs) + self.means
        return rval