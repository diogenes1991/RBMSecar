import numpy as np
import matplotlib.pyplot as plt
from Hamiltonian import *

class RealHamiltonian(Hamiltonian):
    
    def __init__(self,**kwargs):
        self.matrix_dimension = lambda dim : int(dim*(dim+1)/2)
        self.type             = np.float64
        self.matrix_method    = build_pdm
        super().__init__(**kwargs)
    
    def predict_real(self,alpha):
        evl,evc = self.predict_pcs(alpha)
        rval    = np.zeros(self.width)
            
        ''' Loop over principal components '''
        for j in range(self.dim):
            rval += evc.T[0][j]*self.u[:,j]
            
        return rval
    
    def test(self):
        ''' Loop over eigenfunctions '''
        for index in range(self.dim):
            ''' Loop over known scales '''
            for scale in self.axis:
                evl,evc = self.predict_pcs(scale)
                rval    = np.zeros(self.width)
                
                ''' Loop over principal components '''
                for i in range(self.dim):
                    rval += evc.T[index][i]*self.u[:,i]
                
                ''' Plot normalized wavefunctions '''
                rval = rval / np.linalg.norm(rval)
                plt.plot(rval)
            
            plt.show()
            
    def get_error(self):
        error = []
        ''' Loop over known scales '''
        for i,scale in enumerate(self.axis):
            
            evl,evc = self.predict(scale)
            rval    = np.zeros(self.width)
            
            ''' Loop over principal components '''
            for j in range(self.dim):
                rval += evc.T[0][j]*self.u[:,j]
            
            ''' Plot normalized wavefunctions '''
            rval = rval / np.linalg.norm(rval)
            
            ''' Sign at the center '''
            s = np.sign(rval[int(len(rval)/2)])
            
            ''' Why did we transpose the data again? '''
            error.append(np.linalg.norm(s*rval - self.data[:,i]))
        
        return error
