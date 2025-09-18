# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 11:18:54 2025

@author: Julu
"""

import numpy as np


def error(x,y):
    res: int
    x64 = np.float64(x)
    y64 = np.float64(y)
    res = x64 - y64
    if res < 0:
        res = -res  
    return res


def error_relativo(x,y):
    res: int
    x64 = np.float64(x)
    err_abs: int = error(x,y)
    res = err_abs / abs(x64) if x64 != 0 else np.nan
    return res


def matricesIguales(A, B):
    res: bool = True
    for fil in range(A.shape[0]):
        for col in range(A.shape[1]):
            if not(sonIguales(A[fil][col],B[fil][col])):
                res = False
    return res

## L02 ----------------------------------

def rota(theta):
    return np.array([[np.cos(theta), -np.sin(theta)],
                     [np.sin(theta),  np.cos(theta)]])


def escala(s):
    res: list = []
    for i in range(len(s)):
        vector = []
        for j in range(len(s)):
            if i == j:
                vector.append(s[i])
            else:
                vector.append(0)
        vector = np.array(vector)
        res.append(vector)
    res = np.array(res)
    return res


def rotayescala(theta, s):
    return escala(s) @ rota(theta)


def afin(theta, s, b):
    res: list = []
    res.append(np.array([rotayescala(theta,s)[0][0], rotayescala(theta,s)[0][1], b[0]]))
    res.append(np.array([rotayescala(theta,s)[1][0], rotayescala(theta,s)[1][1], b[1]]))
    res.append(np.array([0,0,1]))
    return np.array(res)


def traspuesta(A):
    A = np.array(A)
    return A.T
    """
    filas = A.shape[0]
    columnas = A.shape[1]
    fil = 0
    col = 0
    res: list = []
    for fil in range(0, filas):
        vector = []
        for col in range(0, columnas):
            vector.append(A[col][fil])
        res.append(vector)
    res = np.array(res)
    """



def transafin(v, theta, s, b):
    xy1 = np.array([v[0],v[1],1])
    res = afin(theta, s, b) @ xy1
    return res[:2]


## L03 ----------------------------------

def norma(x, p):
   # x = np.array(x, dtype=float)  # convierte cualquier entrada a array
    if p == 'inf':
        actual = 0
        for i in x:
            if abs(i) > actual:
                actual = abs(i)
        return actual

    else:
        res: float = 0
        for i in x:
            res = res + abs(i)**p
        res = res**(1/p)
        return res


def normaliza(X, p):    
    res = []
    for v in X:
        v = np.array(v, dtype=float)
        n = norma(v, p)
        if n != 0:
            res.append(v / n)
        else:
            res.append(v)
    return res    


def normaMatMC(A, q, p, Np):
    maximoValor: float = 0
    maximoVector = None
    for it in range(Np):
        x = np.random.randn(A.shape[1])
        x = x / norma(x, p)
        Ax = A @ x
        valorActual = norma(Ax, q) / norma(x, p)
      #  normaVectEntrada = norma(x, q)
        if maximoValor < valorActual:
            maximoValor = valorActual
            maximoVector = x
    return maximoValor, maximoVector


def normaExacta(A, p):
    maximo = 0
    if p == 1:
        return normaExacta(traspuesta(A), 'inf')
    elif p == 'inf':
        for x in A:
            actual = norma(x,1)
            if maximo < actual:
                maximo = actual
        return maximo


def condMC(A, p, Np):
    return normaMatMC(A, p, p, Np)[0]*normaMatMC(np.linalg.inv(A), p, p, Np)[0]


def condExacta(A, p):
    return normaExacta(A, p)*normaExacta(np.linalg.inv(A), p)


## ----------------------------------------



    
#---------------------------------
# TESTS
#---------------------------------

def sonIguales(x, y, atol=1e-08):
    return np.allclose(error(x, y), 0, atol=atol)

assert not sonIguales(1, 1.1)
assert sonIguales(1, 1 + np.finfo('float64').eps)
assert not sonIguales(1, 1 + np.finfo('float32').eps)
assert not sonIguales(np.float16(1), np.float16(1) + np.finfo('float32').eps)
assert sonIguales(np.float16(1), np.float16(1) + np.finfo('float16').eps, atol=1e-3)

assert np.allclose(error_relativo(1, 1.1), 0.1)
assert np.allclose(error_relativo(2, 1), 0.5)
assert np.allclose(error_relativo(-1, -1), 0)
assert np.allclose(error_relativo(1, -1), 2)

assert matricesIguales(np.diag([1, 1]), np.eye(2))
assert matricesIguales(np.linalg.inv(np.array([[1, 2], [3, 4]])) @ np.array([[1, 2], [3, 4]]), np.eye(2))
assert not matricesIguales(np.array([[1, 2], [3, 4]]).T, np.array([[1, 2], [3, 4]]))

## ------------------------------------------

# Tests para rota
assert np.allclose(rota(0), np.eye(2))
assert np.allclose(rota(np.pi / 2), np.array([[0, -1], [1, 0]]))
assert np.allclose(rota(np.pi), np.array([[-1, 0], [0, -1]]))

# Tests para escala
assert np.allclose(escala([2, 3]), np.array([[2, 0], [0, 3]]))
assert np.allclose(escala([1, 1, 1]), np.eye(3))
assert np.allclose(escala([0.5, 0.25]), np.array([[0.5, 0], [0, 0.25]]))

# Tests para rota y escala
assert np.allclose(rotayescala(0, [2, 3]), np.array([[2, 0], [0, 3]]))
assert np.allclose(rotayescala(np.pi / 2, [1, 1]), np.array([[0, -1], [1, 0]]))
assert np.allclose(rotayescala(np.pi, [2, 2]), np.array([[-2, 0], [0, -2]]))

# Tests para afin
assert np.allclose(
    afin(0, [1, 1], [1, 2]),
    np.array([[1, 0, 1],
              [0, 1, 2],
              [0, 0, 1]])
)

assert np.allclose(
    afin(np.pi / 2, [1, 1], [0, 0]),
    np.array([[0, -1, 0],
              [1,  0, 0],
              [0,  0, 1]])
)

assert np.allclose(
    afin(0, [2, 3], [1, 1]),
    np.array([[2, 0, 1],
              [0, 3, 1],
              [0, 0, 1]])
)

# Tests para transafin
assert np.allclose(
    transafin(np.array([1, 0]), np.pi / 2, [1, 1], [0, 0]),
    np.array([0, 1])
)

assert np.allclose(
    transafin(np.array([1, 1]), 0, [2, 3], [0, 0]),
    np.array([2, 3])
)

assert np.allclose(
    transafin(np.array([1, 0]), np.pi / 2, [3, 2], [4, 5]),
    np.array([4, 7])
)

## ------------------------------------------

# Tests L03-Normas

# Tests norma
assert(np.allclose(norma(np.array([1,1]),2),np.sqrt(2)))
assert(np.allclose(norma(np.array([1]*10),2),np.sqrt(10)))
assert(norma(np.random.rand(10),2)<=np.sqrt(10))
assert(norma(np.random.rand(10),2)>=0)

# Tests normaliza
# Tests normaliza
for x in normaliza([np.array([1]*k) for k in range(1,11)],2):
    assert(np.allclose(norma(x,2),1))
for x in normaliza([np.array([1]*k) for k in range(2,11)],1):
    assert(not np.allclose(norma(x,2),1) )
for x in normaliza([np.random.rand(k) for k in range(1,11)],'inf'):
    assert( np.allclose(norma(x,'inf'),1) )


# Tests normaExacta

assert(np.allclose(normaExacta(np.array([[1,-1],[-1,-1]]),1),2))
assert(np.allclose(normaExacta(np.array([[1,-2],[-3,-4]]),1),6))
assert(np.allclose(normaExacta(np.array([[1,-2],[-3,-4]]),'inf'),7))
assert(normaExacta(np.array([[1,-2],[-3,-4]]),2) is None)
assert(normaExacta(np.random.random((10,10)),1)<=10)
assert(normaExacta(np.random.random((4,4)),'inf')<=4)

# Test normaMC

nMC = normaMatMC(A=np.eye(2),q=2,p=1,Np=100000)
assert(np.allclose(nMC[0],1,atol=1e-3))
assert(np.allclose(np.abs(nMC[1][0]),1,atol=1e-3) or np.allclose(np.abs(nMC[1][1]),1,atol=1e-3))
assert(np.allclose(np.abs(nMC[1][0]),0,atol=1e-3) or np.allclose(np.abs(nMC[1][1]),0,atol=1e-3))

nMC = normaMatMC(A=np.eye(2),q=2,p='inf',Np=100000)
assert(np.allclose(nMC[0],np.sqrt(2),atol=1e-3))
assert(np.allclose(np.abs(nMC[1][0]),1,atol=1e-3) and np.allclose(np.abs(nMC[1][1]),1,atol=1e-3))

A = np.array([[1,2],[3,4]])
nMC = normaMatMC(A=A,q='inf',p='inf',Np=1000000)
assert(np.allclose(nMC[0],normaExacta(A,'inf'),rtol=2e-1)) 

# Test condMC

A = np.array([[1,1],[0,1]])
A_ = np.linalg.solve(A,np.eye(A.shape[0]))
normaA = normaMatMC(A,2,2,10000)
normaA_ = normaMatMC(A_,2,2,10000)
condA = condMC(A,2,10000)
assert(np.allclose(normaA[0]*normaA_[0],condA,atol=1e-3))

A = np.array([[3,2],[4,1]])
A_ = np.linalg.solve(A,np.eye(A.shape[0]))
normaA = normaMatMC(A,2,2,10000)
normaA_ = normaMatMC(A_,2,2,10000)
condA = condMC(A,2,10000)
assert(np.allclose(normaA[0]*normaA_[0],condA,atol=1e-3))

# Test condExacta

A = np.random.rand(10,10)
A_ = np.linalg.solve(A,np.eye(A.shape[0]))
normaA = normaExacta(A,1)
normaA_ = normaExacta(A_,1)
condA = condExacta(A,1)
assert(np.allclose(normaA*normaA_,condA))

A = np.random.rand(10,10)
A_ = np.linalg.solve(A,np.eye(A.shape[0]))
normaA = normaExacta(A,'inf')
normaA_ = normaExacta(A_,'inf')
condA = condExacta(A,'inf')
assert(np.allclose(normaA*normaA_,condA))


## ------------------------------------------

# Tests L04-LU

# Tests LU

L0 = np.array([[1,0,0],[0,1,0],[1,1,1]])
U0 = np.array([[10,1,0],[0,2,1],[0,0,1]])
A =  L0 @ U0
L,U,nops = calculaLU(A)
assert(np.allclose(L,L0))
assert(np.allclose(U,U0))


L0 = np.array([[1,0,0],[1,1.001,0],[1,1,1]])
U0 = np.array([[1,1,1],[0,1,1],[0,0,1]])
A =  L0 @ U0
L,U,nops = calculaLU(A)
assert(not np.allclose(L,L0))
assert(not np.allclose(U,U0))
assert(np.allclose(L,L0,atol=1e-3))
assert(np.allclose(U,U0,atol=1e-3))
assert(nops == 13)

L0 = np.array([[1,0,0],[1,1,0],[1,1,1]])
U0 = np.array([[1,1,1],[0,0,1],[0,0,1]])
A =  L0 @ U0
L,U,nops = calculaLU(A)
assert(L is None)
assert(U is None)
assert(nops == 0)

## Tests res_tri

A = np.array([[1,0,0],[1,1,0],[1,1,1]])
b = np.array([1,1,1])
assert(np.allclose(res_tri(A,b),np.array([1,0,0])))
b = np.array([0,1,0])
assert(np.allclose(res_tri(A,b),np.array([0,1,-1])))
b = np.array([-1,1,-1])
assert(np.allclose(res_tri(A,b),np.array([-1,2,-2])))
b = np.array([-1,1,-1])
assert(np.allclose(res_tri(A,b,inferior=False),np.array([-1,1,-1])))

A = np.array([[3,2,1],[0,2,1],[0,0,1]])
b = np.array([3,2,1])
assert(np.allclose(res_tri(A,b,inferior=False),np.array([1/3,1/2,1])))

A = np.array([[1,-1,1],[0,1,-1],[0,0,1]])
b = np.array([1,0,1])
assert(np.allclose(res_tri(A,b,inferior=False),np.array([1,1,1])))

# Test inversa

ntest = 10
iter = 0
while iter < ntest:
    A = np.random.random((4,4))
    A_ = inversa(A)
    if not A_ is None:
        assert(np.allclose(np.linalg.inv(A),A_))
        iter += 1

# Matriz singular devería devolver None
A = np.array([[1,2,3],[4,5,6],[7,8,9]])
assert(inversa(A) is None)




# Test LDV:

L0 = np.array([[1,0,0],[1,1.,0],[1,1,1]])
D0 = np.diag([1,2,3])
V0 = np.array([[1,1,1],[0,1,1],[0,0,1]])
A =  L0 @ D0  @ V0
L,D,V,nops = calculaLDV(A)
assert(np.allclose(L,L0))
assert(np.allclose(D,D0))
assert(np.allclose(V,V0))

L0 = np.array([[1,0,0],[1,1.001,0],[1,1,1]])
D0 = np.diag([3,2,1])
V0 = np.array([[1,1,1],[0,1,1],[0,0,1.001]])
A =  L0 @ D0  @ V0
L,D,V,nops = calculaLDV(A)
assert(np.allclose(L,L0,1e-3))
assert(np.allclose(D,D0,1e-3))
assert(np.allclose(V,V0,1e-3))

# Tests SDP

L0 = np.array([[1,0,0],[1,1,0],[1,1,1]])
D0 = np.diag([1,1,1])
A = L0 @ D0 @ L0.T
assert(esSDP(A))

D0 = np.diag([1,-1,1])
A = L0 @ D0 @ L0.T
assert(not esSDP(A))

D0 = np.diag([1,1,1e-16])
A = L0 @ D0 @ L0.T
assert(not esSDP(A))

L0 = np.array([[1,0,0],[1,1,0],[1,1,1]])
D0 = np.diag([1,1,1])
V0 = np.array([[1,0,0],[1,1,0],[1,1+1e-10,1]]).T
A = L0 @ D0 @ V0
assert(not esSDP(A))
