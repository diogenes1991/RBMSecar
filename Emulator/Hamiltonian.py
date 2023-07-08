import scipy as sp
import scipy.optimize as op
import numpy as np
import random
from BlackBoxRBM import BlackBoxRBM

def build_pdm(coefficients,n,parity=1):
    ''' Builds a parity defined matrix '''
    matrix = [[0] * n for _ in range(n)]  # Initialize matrix with zeros
    cc=0
    for i in range(n):
        for j in range(n):
            if ( j<=i if parity==-1 else j<i) :
                continue
            coefficient = coefficients[cc]  # Correct indexing
            matrix[i][j] = coefficient
            matrix[j][i] = parity*coefficient
            cc+=1
    matrix = np.array(matrix)
    return matrix

def build_hermitian_matrix(coefficients,n):
    cutpoint = int(n*(n+1)/2)
    return build_pdm(coefficients[:cutpoint],n,1) + 1j*build_pdm(coefficients[cutpoint:],n,-1)
    
class Hamiltonian(BlackBoxRBM):
    def __init__(self,**kwargs):
        self.axis              = kwargs["axis"]
        self.data              = kwargs["data"]
        self.width,self.height = self.data.shape
        self.u,self.s,self.v   = np.linalg.svd(self.data)
     
    ''' Enforce Eigenvectors '''
    def objective(self,variables):
        rval = 0
            
        ''' Build H0, H1, ... '''
        self.hs = [ self.matrix_method(variables[i*self.numEntries:(i+1)*self.numEntries],self.dim) for i in range(self.ord) ]
        for i,vec in enumerate(self.low_dim_rep):
                
            ''' Loop over orders '''
            self.H = np.zeros((self.dim,self.dim),self.type)
            for j in range(self.ord):
                self.H += (self.axis[i]**j)*self.hs[j]
            
            ''' Compute Eigenvectors and add residue as difference '''
            evl,evc = np.linalg.eigh(self.H)
            rval += np.linalg.norm( evc.T[0]*np.sign(evc.T[0][0]) - vec*np.sign(vec[0]) )
         
        return rval
         
    def build(self,dimension,order):
        
        ''' Here we build an approximate to the Galerkin Equations '''
        self.dim = dimension
        self.ord = order
        
        ''' First we project every datapoint onto the first dim principal components '''
        self.low_dim_rep = np.array([ [ np.dot(self.u[:,i],self.data[:,j]) for i in range(self.dim) ] for j in range(self.height) ])
        
        ''' Normalize once '''
        for i in range(len(self.low_dim_rep)):
            self.low_dim_rep[i] = self.low_dim_rep[i]/np.linalg.norm(self.low_dim_rep[i])
        
        ''' Setup a minimization problem '''
        self.numEntries = self.matrix_dimension(self.dim)
        self.vars       = np.array( [ 2*random.random()-1 for i in range(self.ord*self.numEntries)] )
        
        ''' Minimization call '''
        self.op_solution = op.minimize(self.objective,self.vars)
        self.minimum = self.op_solution['x']
    
    def predict_pcs(self,alpha):
        H = np.zeros((self.dim,self.dim),self.type)
        for j in range(self.ord):
            H += (alpha**j)*self.hs[j]
        return np.linalg.eigh(H)
    
    def predict_real(self,alpha):
        NotImplemented
        
    def test(self):
        NotImplemented
    
    def get_error(self):
        NotImplemented   
