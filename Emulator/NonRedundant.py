import numpy as np
from Emulator import Emulator

class NonRedundant(Emulator):
    
    def normalize(self):
        ''' Compute normalizations of the data '''
        self.means = [ np.average(self.matrix[:,i]) for i in range(self.width) ]
        self.sdevs = [ np.std(self.matrix[:,i])     for i in range(self.width) ]
        
        ''' Mapping between emulated and total channels '''
        self.map   = {}
        self.emulc = 0
        for i,sdev in enumerate(self.sdevs):
            if sdev == 0:
                continue
            self.map[self.emulc] = i
            self.emulc+=1
        
        ''' Renormalize channels to be emulated '''
        self.renor = []
        for i in range(self.height):
            row = []
            for j in range(self.width):
                if self.sdevs[j] != 0:
                    row.append( (self.matrix[i][j]-self.means[j]) / self.sdevs[j] )
            self.renor.append(row)
        
        self.renor  = np.array(self.renor)
    
    def __call__(self,value):
        ''' Actual Magnet representation '''
        rval = np.zeros(self.width)
        
        ''' Loop over principal components '''
        emul = np.zeros(self.emulc)
        for i in range(self.cutoff):
            emul += self.splines[i](value)*self.s[i]*self.v[i]
        
        ''' Convert back to possibly constant channels '''
        for key in self.map:
            rval[self.map[key]] = emul[key]
        
        ''' Rescale back data '''
        rval = (np.multiply(rval,self.sdevs) + self.means)
        return rval
