import os
import numpy as np
from Emulator import Emulator
from CosyIO import CosyIO
from NonRedundant import NonRedundant

class CosyMagnet(NonRedundant):
    
    ''' Wrapped emulator for the Cosy format '''
    
    def __init__(self,**kwargs):
        
        ''' Load all the datafile names and assume the scale is the name  '''
        self.path   = kwargs["path"]
        self.data   = os.listdir(self.path)
        scales      = [ float(name) for name in self.data ]
        
        ''' Format the flattened data matrix, file = row '''
        matrix = []
        
        self.cosyIO = CosyIO()
        for i in range(len(self.data)):
            ''' 
                Here we should save the magnet specs and complain if 
                some magnets do not have the same specs
            '''
            self.cosyIO.read(os.path.join(self.path,self.data[i]))
            
            ''' Flatten the matrices '''
            big_row = []
            for mat in self.cosyIO.matrices:
                for row in mat:
                    big_row.append(row)
            matrix.append(big_row)
        matrix = np.array(matrix)
        
        ''' Emulator build '''
        self.cutoff = kwargs["cutoff"]
        super().__init__(scales,matrix,self.cutoff)
        
    def write(self,value,path):
        rval     = self(value)
        mat_len  = len(self.cosyIO.lookup.keys())
        matrices = [ rval[mat_len*i:mat_len*(i+1)] for i in range(self.cosyIO.nLetters) ]
        self.cosyIO.matrices = matrices
        self.cosyIO.write(path)
        