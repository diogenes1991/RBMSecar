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


class CosyIO:
    
    def __init__(self,nLetters=8,lowPower=0,highPower=4):
        
        self.nLetters  = nLetters
        self.lowPower  = lowPower
        self.higiPower = highPower
        
        ''' We generate a static lookup table to synchronize all matrices '''
        self.lookup = generate_lookup(lowPower=lowPower,highPower=highPower,nLetters=nLetters)
        self.inverse  = dict( [(self.lookup[key],key) for key in self.lookup] )
        
        ''' Convert keys to list of powers '''
        self.key_to_list = lambda k,l : self.key_to_list(k[1:],l+[int(k[0])]) if len(k)!=0 else l
        
    def read(self,path):
        file = open(path,"r")
        self.matrices = []
        local_mat = np.zeros(len(self.lookup.keys()))
        count = 1
        for line in file:
            count += 1
            vec = line.split()
            
            ''' Save specifications '''
            if vec[0] == "L":
                self.L = float(vec[1])
                continue
            if vec[0] == "P":
                self.P = [float(v[:-1]) for v in vec[2:]]
                continue
            if vec[0] == "A":
                self.A = [float(v[:-1]) for v in vec[2:]]
                continue
                
            ''' Skip Matrix header '''
            if vec[0] == "I":
                self.header = line
                continue
            
            ''' Append at the end of matrix character '''
            if vec[0][0] == "-":
                self.footer = line
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
        
    def write(self,path):
        ''' Comments for later '''
        f = open(path,"w")
        f.write(f" L    {self.L:.16E}\n")
        f.write(f" P ({self.P[0]:.16E}, {self.P[1]:.16E}, {self.P[2]:.16E})\n")
        f.write(f" A ({self.A[0]:.16E}, {self.A[1]:.16E}, {self.A[2]:.16E})\n")
        for i in range(len(self.matrices)):
            f.write(self.header)
            cc = 0
            for j in range(len(self.matrices[i])):
                if self.matrices[i][j] != 0:
                    powers = self.key_to_list(self.inverse[cc],[])
                    line = f"   {cc+1:3d}  {self.matrices[i][j]:20.16E}  "
                    line += f"{np.sum(powers):2d} "
                    for value in powers:
                        line += f"{value:1d} "
                    cc+=1
                    f.write(line)
                    f.write("\n")
            f.write(self.footer)
        f.close()