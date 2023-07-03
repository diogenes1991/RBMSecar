import numpy as np
from CosyMagnet import CosyMagnet

class SECARMagnet(CosyMagnet):
    
    ''' Interface between cosy and secar '''
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.powers  = np.array( [ np.array([np.array(self.cosyIO.key_to_list(key,[])) for key in self.cosyIO.lookup.keys()]) for j in range(self.cosyIO.nLetters) ] ) 
    
    ''' Secar level interface call    '''
    def evaluate(self,value):
        ''' This code is repeated in cosyIO'''
        rval = super().__call__(value)
        rval     = self(value)
        mat_len  = len(self.cosyIO.lookup.keys())
        matrices = np.array([ np.array(rval[mat_len*i:mat_len*(i+1)]) for i in range(self.cosyIO.nLetters) ])
        return matrices