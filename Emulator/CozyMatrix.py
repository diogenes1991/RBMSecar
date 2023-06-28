class CozyMatrix:
    def __init__(self,path):
        file = open(path,"r")
        self.matrices = []
        local_mat = 331*[0]
        count = 1
        for line in file:
            count += 1
            vec = line.split()
            if vec[0] == "L" or vec[0] == "P" or vec[0] == "A" or vec[0][0] == "-":
                continue
            if vec[0] == "I" and len(local_mat)!=0:
                self.matrices.append(local_mat)
                #print(local_mat)
                local_mat = 331*[0]
            elif vec[0] != "I":
                word = ""
                for i in range(3,len(vec)):
                    word+=vec[i]
                try:
                    local_mat[lookup[word]-1] = float(vec[1])
                except:
                    print(    