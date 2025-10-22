# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 09:49:16 2025

@author: Julu
"""

import numpy as np

def calcularAx(A,x):
    columnas = A.shape[1]
    filas_A = A.shape[0]
    filas_x = x.shape[0]
    res: list = []
    for i in range(filas_A):
        suma = 0
        for j in range(min(columnas, filas_x)):  # recorrer columnas
            suma += A[i][j] * x[j]
        res.append(suma)
    res = np.array(res)    
    return res

def calcularxA(x,A):
    columnas = A.shape[1]
    res: list = []
    for col in range(columnas):
        res.append(productoInterno(x, A[col]))
    res = np.array(res)    
    return res

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

def productoInterno(u,v,tol=1e-12):
    res = 0
    n = u.shape[0]
    for fil in range(n):
        res += u[fil]*v[fil]
    return res


def signo(a):
    if a < 0:
        return -1
    else:
        return 1
    
def productoExterno(u,v):
    res = []
    for i in range(len(u)):
        vector = []
        for j in range(len(v)):
            vector.append(u[i]*v[j])
        vector = np.array(vector)
        res.append(vector)
    res = np.array(res)
    return res

def traspuesta(A):   
    filas = A.shape[0]
    columnas = A.shape[1]
    fil = 0
    col = 0
    res: list = []
    for fil in range(0, filas):
        vector = []
        for col in range(0, columnas):
            vector.append(A[col][fil])
        vector = np.array(vector)
        res.append(vector)
    res = np.array(res)
    return res

def multiplicacionMatricial(A, B, tol=1e-12):
    columnas_A = A.shape[1]
    filas_A = A.shape[0]
    filas_B = B.shape[0]
    columnas_B = B.shape[1]
    if columnas_A != filas_B:
        return None
    else:
        res = []
        for fil in range(filas_A):
            vector = []
            for col in range(columnas_B):        
                vector.append(productoInterno(A[fil,:], B[:,col]))
            res.append(vector)
        res = np.array(res)
        return res

def QR_con_HH(A,tol=1e-12):   
    n = A.shape[0]
    m = A.shape[1]
    if m < n:
        return None
    else:
        R = A.copy()
        Q = np.eye(m)
        for k in range(n):
            x = R[k:, k]
            a = -signo(x[0])*norma(x,2)
            u = x - a*(np.eye(x.shape[0])[0])
            if norma(u, 2) > tol:
                u = u/norma(u,2)
                Hk = np.eye(len(u)) - 2*(productoExterno(u, u))
                H = np.eye(m)
                H[k:, k:] = Hk
                R = multiplicacionMatricial(H, R)
                Q = multiplicacionMatricial(Q, traspuesta(H))
    return Q, R

def metpot2k(A, tol=1e-15, K=1000):
    n = A.shape[1]
    v = np.random.randn(n)
    v /= np.linalg.norm(v)
    w = calcularAx(A, v)
    e = productoInterno(w, v)
    k = 0
    while (abs(e - 1) > tol) and (k < K):
        v = w
        w = calcularAx(A, v)
        w /= np.linalg.norm(w)
        e = productoInterno(w, v)
        k += 1
    lam = productoInterno(calcularxA(w, A), w)
    #error = e - 1
    it = k
    return w, lam, it
    """
    Método de la potencia para aproximar el autovector y autovalor dominante.

    Parámetros
    ----------
    A : array-like (n x n)
        Matriz cuadrada.
    tol : float
        Tolerancia para el criterio de convergencia.
    K : int
        Número máximo de iteraciones.

    Retorna
    -------
    v : ndarray
        Autovector aproximado (normalizado con norma 1).
    lam : float
        Autovalor aproximado asociado (por cociente de Rayleigh).
    it : int
        Número de iteraciones realizadas.
    """

def diagRH(A, tol=1e-12, K=1000):
    n = A.shape[0]
    Q_total = np.eye(n)
    k: int = 0
    while k < K:
        Q, R = QR_con_HH(A)
        A = multiplicacionMatricial(R, Q)
        Q_total = multiplicacionMatricial(Q_total, Q)
        off_diag_norm = np.sqrt(np.sum(np.tril(A, -1)**2))
        if off_diag_norm < tol:
            break
        k += 1
    D = np.diag(np.diag(A))
    return Q_total, D
    
    """
    Diagonalización de una matriz simétrica real.

    Si A es simétrica, devuelve matrices S y D tales que A = S D S^T,
    donde las columnas de S son autovectores ortonormales
    y D es una matriz diagonal con los autovalores.

    Parámetros
    ----------
    A : array-like (n x n)
        Matriz cuadrada a diagonalizar.

    Retorna
    -------
    S : ndarray
        Matriz de autovectores (ortonormales, en columnas).
    D : ndarray
        Matriz diagonal de autovalores.
        Si A no es simétrica, devuelve None.
    """
    
#### TESTEOS
# Tests metpot2k

S = np.vstack([
    np.array([2,1,0])/np.sqrt(5),
    np.array([-1,2,5])/np.sqrt(30),
    np.array([1,-2,1])/np.sqrt(6)
              ]).T

# Pedimos que pase el 95% de los casos
exitos = 0
for i in range(100):
    D = np.diag(np.random.random(3)+1)*100
    A = S@D@S.T
    v,l,_ = metpot2k(A,1e-15,1e5)
    if np.abs(l - np.max(D))< 1e-8:
        exitos += 1
assert exitos > 95


#Test con HH
exitos = 0
for i in range(100):
    v = np.random.rand(9)
    #v = np.abs(v)
    #v = (-1) * v
    ixv = np.argsort(-np.abs(v))
    D = np.diag(v[ixv])
    I = np.eye(9)
    H = I - 2*np.outer(v.T, v)/(np.linalg.norm(v)**2)   #matriz de HouseHolder

    A = H@D@H.T
    v,l,_ = metpot2k(A, 1e-15, 1e5)
    #max_eigen = abs(D[0][0])
    if abs(l - D[0,0]) < 1e-8:         
        exitos +=1
assert exitos > 95



# Tests diagRH
D = np.diag([1,0.5,0.25])
S = np.vstack([
    np.array([1,-1,1])/np.sqrt(3),
    np.array([1,1,0])/np.sqrt(2),
    np.array([1,-1,-2])/np.sqrt(6)
              ]).T

A = S@D@S.T
SRH,DRH = diagRH(A,tol=1e-15,K=1e5)
assert np.allclose(D,DRH)
assert np.allclose(np.abs(S.T@SRH),np.eye(A.shape[0]),atol=1e-7)


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
# Pedimos que pase el 95% de los casos
exitos = 0
for i in range(100):
    A = np.random.random((5,5))
    A = 0.5*(A+A.T)
    S,D = diagRH(A,tol=1e-15,K=1e5)
    ARH = S@D@S.T
    e = normaExacta(ARH-A,p='inf')
    if e < 1e-5: 
        exitos += 1
assert exitos >= 95
