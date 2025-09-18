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

def normaMatMC(A, q, p, Np):
   
    n = A.shape[1]
    mejorvalor = 0
    mejorvector = " "

    for _ in range(Np):
       
        x = np.random.randn(n) #genero los vectores aleatorios (los np)

 
        normadex = norma(x, p)
        if normadex != 0:
            for i in range(n):
                x[i] = x[i] / normadex

        
            Ax = np.dot(A, x) #me calcula la norma "q" Ax
            val = norma(Ax, q)

            if val > mejorvalor:
                mejorvalor = val
                mejorvector = np.array(x) #actualizo el mejor

    return mejorvalor, mejorvector


A = np.array([[1, 2],
              [3, 4]], dtype=float)


valor, xopt = normaMatMC(A, q=2, p=2, Np=5)

print("Norma inducida estimada:", valor)
print("Vector óptimo aproximado:", xopt)










