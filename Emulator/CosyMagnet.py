import os
import numpy as np
from Emulator import Emulator
from CosyMatrix import CosyMatrix

class CosyMagnet(Emulator):
    
    ''' Wrapped emulator for the Cosy format '''
    
    def __init__(self,**kwargs):
        
        ''' Load all the datafile names and assume the scale is the name  '''
        self.path   = kwargs["path"]
        self.data   = os.listdir(self.path)
        scales      = [ float(name) for name in self.data ]
        
        ''' Format the flattened data matrix, file = row '''
        matrix = []
        for i in range(len(self.data)):
            M = CosyMatrix(os.path.join(self.path,self.data[i]))
            big_row = []
            for mat in M.matrices:
                for row in mat:
                    big_row.append(row)
            matrix.append(big_row)
        matrix = np.array(matrix)
        
        ''' Emulator build '''
        self.cutoff = kwargs["cutoff"]
        super().__init__(scales,matrix,self.cutoff)