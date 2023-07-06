import numpy as np
from CosyIO import CosyIO

class StaticSECARMagnet:
    
    ''' Lookup table magnets '''
    
    def __init__(self,**kwargs):
        self.path     = kwargs["path"]
        
        ''' Load magnet data '''
        self.cosyIO   = CosyIO()
        self.cosyIO.read(self.path)
        
        ''' Format matrices '''
        self.matrices = np.array(self.cosyIO.matrices)
        
        ''' Powers, this is repeated and should be moved to CosyIO '''
        self.powers  = np.array( [ np.array([np.array(self.cosyIO.key_to_list(key,[])) for key in self.cosyIO.lookup.keys()]) for j in range(self.cosyIO.nLetters) ] ) 
        
    def evaluate(self,current):
        ''' Return read matrices '''
        return self.matrices