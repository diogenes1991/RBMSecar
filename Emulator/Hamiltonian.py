''' 
    Fitting the Hamiltonian 
    Model H = H0 + alpha * H1,  with H0, H1 Symmetric
    
'''

import scipy as sp
import scipy.optimize as op
import numpy as np
import random

def build_sym_mat(coefficients,n):
    matrix = [[0] * n for _ in range(n)]  # Initialize matrix with zeros
    cc=0
    for i in range(n):
        for j in range(n):
            if j<i:
                continue
            coefficient = coefficients[cc]  # Correct indexing
            matrix[i][j] = coefficient
            matrix[j][i] = coefficient
            cc+=1
    matrix = np.array(matrix)
    return matrix

class Hamiltonian:
    def __init__(self,axis,data):
        self.axis              = axis
        self.data              = data
        self.width,self.height = self.data.shape
        self.u,self.s,self.v   = np.linalg.svd(self.data)
    
    def build(self,dimen,order):
        
        ''' Here we build an approximate to the Galerkin Equations '''
        self.dim = dimen
        self.ord = order
        
        ''' First we project every datapoint onto the first dim principal components '''
        self.low_dim_rep = np.array([ [ np.dot(self.u[:,i],self.data[:,j]) for i in range(self.dim) ] for j in range(self.height) ])
        
        ''' Normalize once '''
        for i in range(len(self.low_dim_rep)):
            self.low_dim_rep[i] = self.low_dim_rep[i]/np.linalg.norm(self.low_dim_rep[i])
        
        ''' Setup a minimization problem '''
        self.numEntries = int(self.dim*(self.dim+1)/2)
        self.vars       = np.array( [ i+1 for i in range(self.ord*self.numEntries)] )
        norm = lambda v : np.sum(np.square(v))
        
        ''' Enforce Eigenvectors '''
        def objective(variables):
            rval = 0
            
            ''' Build H0, H1, ... '''
            self.hs = [ build_sym_mat(variables[i*self.numEntries:(i+1)*self.numEntries],self.dim) for i in range(self.ord) ]
            for i,vec in enumerate(self.low_dim_rep):
                    
                ''' Loop over orders '''
                self.H = np.zeros((self.dim,self.dim))
                for j in range(self.ord):
                    self.H += (self.axis[i]**j)*self.hs[j]
                
                ''' Compute Eigenvectors and add residue as difference '''
                evl,evc = np.linalg.eigh(self.H)
                rval += norm( evc.T[0]*np.sign(evc.T[0][0]) - vec*np.sign(vec[0]) )
                
            return rval
        
        ''' Minimization call '''
        self.op_solution = op.minimize(objective,self.vars)
        self.minimum = self.op_solution['x']
    
    def predict_pcs(self,alpha):
        H = np.zeros((self.dim,self.dim))
        for j in range(self.ord):
            H += (alpha**j)*self.hs[j]
        return np.linalg.eigh(H)
    
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
                evl,evc = self.predict(scale)
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