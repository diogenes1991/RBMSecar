from itertools import combinations_with_replacement
import numpy as np

def generate_lookup(**kwargs):
    
    low_pow = kwargs["lowPower"]
    hig_pow = kwargs["highPower"]
    num_let = kwargs["nLetters"]
    
    ''' Generate all nLetters words '''
    letters = [ chr(i) for i in range(num_let) ]
    words = []
    for r in range(low_pow, hig_pow+1):
        combinations_list = combinations_with_replacement(letters, r)
        words.extend([''.join(combo) for combo in combinations_list])
    
    ''' Convert from words to a lookup table '''
    word_to_list = lambda word       : [  str(word.count(letter)) for letter in letters ]
    list_to_key  = lambda base,array : list_to_key(base+array[0],array[1:]) if len(array) else base
    word_to_key  = lambda word       : list_to_key('',word_to_list(word))
    
    ''' Lookup table '''
    lookup   = dict( [ (word_to_key(words[i]),i) for i in range(len(words)) ] )
    
    return lookup


class CosyMatrix:    
    def __init__(self,path,nLetters=8,lowPower=0,highPower=4):
        
        ''' We generate a static lookup table to synchronize all matrices '''
        self.lookup = generate_lookup(lowPower=lowPower,highPower=highPower,nLetters=nLetters)
        
        file = open(path,"r")
        self.matrices = []
        local_mat = np.zeros(len(self.lookup.keys()))
        count = 1
        for line in file:
            count += 1
            vec = line.split()
            
            ''' Skip specification lines '''
            if vec[0] == "L" or vec[0] == "P" or vec[0] == "A" or vec[0] == "I":
                continue
            
            ''' Append at the end of matrix character '''
            if vec[0][0] == "-":
                self.matrices.append(local_mat)
                local_mat = np.zeros(len(self.lookup.keys()))
            
            else:
                ''' Build the key to lookup '''
                word = ""
                for i in range(3,len(vec)):
                    word+=vec[i]
                    
                try:
                    local_mat[self.lookup[word]] = float(vec[1])
                except:
                    print(word,count)
        file.close()