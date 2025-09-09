import numpy as np

def norma(x,p):
    if p == "inf":
        temp = 0
        for i in x:
            if abs(i) > temp:
                temp = abs(i)
            else:
                pass
        return temp
    else:
        res: float = 0
        for i in x:
            res = res + i**p
            res = res**(1/p)
        return res

v = np.array([1,2,3])
q = "inf"
print(norma(v,q))
q = 2
print(norma(v,q))

def normaliza(X,p):    
    res: list[float] = []
    for v in X:
        res.append(norma(v,p))
    return res

q = 2
X = [np.array([1,2,3]),np.array([3,3,3]),np.array([1,1,1])]
print(normaliza(X,q))










