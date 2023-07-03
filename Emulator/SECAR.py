import numpy as np
from SECARMagnet import SECARMagnet

class SECAR:
    
    ''' Digital copy of SECAR '''
    
    def __init__(self,**kwargs):
        ''' This assumes that the magnets are to be evaluated sequentially '''
        self.magnets_paths = kwargs["paths"]
        self.cutoff        = kwargs["cutoff"]
        self.magnets       = [ SECARMagnet(path=path,cutoff=self.cutoff) for path in self.magnets_paths ]
        self.nMagnets      = len(self.magnets)
        
    ''' Based on Ruchi's Main.py '''
    def transport(self,x0,map_coeff,map_power):
        x = np.zeros(x0.shape, dtype="float64", order="F")
        samples = x0.shape[0]
        n = map_coeff.shape[0]  # Beam array dimensions
        m = map_coeff.shape[1]  # Number of coefficients

        for l in range(samples):  # Loop over the samples
            for i in range(n):    # Loop over beam array dimensions
                j = 0
                while j < m and map_coeff[i][j] != 0:  # Loop over coefficients
                
                    power_temp = 1
                    for k in range(n):    # Loop over beam array dimensions
                        power_temp = power_temp*x0[l][k]**map_power[i][j][k]
                    x[l][i] = x[l][i] + map_coeff[i][j]*power_temp
                    j = j + 1
        return x 
        
    def emulate(self,**kwargs):
        
        rays     = kwargs["rays"]
        currents = kwargs["currents"]
        
        ''' Emulate Cosy '''
        self.transports = np.array([ np.array(self.magnets[i].evaluate(currents[i])) for i in range(len(currents)) ])
        
        ''' Emulate the transport of the beams '''
        self.state = rays
        for i in range(self.nMagnets):
            self.state = self.transport(self.state,self.transports[i],self.magnets[i].powers)
        
        ''' Return the ray state variable for chaining '''
        return self.state
