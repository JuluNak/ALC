# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 11:18:54 2025

@author: Julu
"""

import numpy as np

def calcularxA(x,A):
    columnas = A.shape[1]
    res: list = []
    for col in range(columnas):
        res.append(productoInterno(x, A[col]))
    res = np.array(res)    
    return res

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

## L01 ---------------------------

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
    filas = A.shape[0]
    columnas = A.shape[1]
    fil = 0
    col = 0
    res: list = []
    for fil in range(0, columnas):
        vector = []
        for col in range(0, filas):
            vector.append(A[col][fil])
        vector = np.array(vector)
        res.append(vector)
    res = np.array(res)
    return res


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


### Funciones L05-QR
def productoInterno(u,v,tol=1e-12):
    res = 0
    n = u.shape[0]
    for fil in range(n):
        res += u[fil]*v[fil]
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

def QR_con_GS(A,tol=1e-12,retorna_nops=False):
    columnas = A.shape[1]
    filas = A.shape[0]
    if columnas != filas:
        return None
    else:
        Q = []
        q1 = A[0]/norma(A[0], 2)
        qj = q1
        Q.append(qj)
        for j in range(1, columnas):
            qj = A[j]
            for k in range(0,j):
                r = productoInterno(A[j], Q[k])
                qj = qj - r*Q[k]
            r = norma(qj, 2)
            qj = qj/r
            Q.append(qj)
        Q = np.array(Q)
        R = multiplicacionMatricial(traspuesta(Q), A)
        return Q, R
    

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
"""
                Hk = np.eye(m-k+1) - 2*multiplicacionMatricial(u, traspuesta(u))
                Hk2 = np.eye(k-1)
                elementos = []
                for i in range(Hk2.shape[1]):
                    for j in range(Hk2.shape[1]):
                        elementos.append(Hk2[i][j])
                    for j in range(Hk2.shape[1],m):
                        elementos.append(0)
                for i in range(Hk.shape[1]):
                    for j in range(Hk2.shape[1]):
                        elementos.append(0)
                    for j in range(Hk2.shape[1],m):
                        elementos.append(Hk[i][j])
                Hk2 = np.array(Hk2)
                Hk2 = elementos.reshape(m,m)
"""

                    
def calculaQR(A,metodo,tol=1e-12):
    if metodo == 'RH':
        return QR_con_HH(A)
    elif metodo == 'GS':
        return QR_con_GS(A)
    else:
        return None
    """
    A una matriz de n x n 
    tol la tolerancia con la que se filtran elementos nulos en R    
    metodo = ['RH','GS'] usa reflectores de Householder (RH) o Gram Schmidt (GS) para realizar la factorizacion
    retorna matrices Q y R calculadas con Gram Schmidt (y como tercer argumento opcional, el numero de operaciones)
    Si el metodo no esta entre las opciones, retorna None
    """
    
#LABORATORIO 7)

#1
#a) QUIERO GENERAR UNA MATRIZ CUYOS VALORES ESTEN ENTRE 0 Y 1 Y LA SUMA DE LAS COLUMNAS SEA 1

def transiciones_al_azar_continuas(n):
   
   
    T = np.random.rand(n, n)  
    
    for j in range(n):
        colum=[]
        for i in range(n):
            colum.append(T[i][j])
            
        sumacolumna = sumacol(colum)
        for i in range(n):
            T[i, j] = T[i, j] / sumacolumna # dividir cada elemento por la suma

    return T


def sumacol(lista):
    res=0
    for x in lista:
        res+=x
    return res

#b) QUIERO GENERAR UNA MATRIZ CON UN CIERTO "x", tal que los aij > x  --> aij=1 , aij < x --> aij=0 , esto para todo aij, ademas luego tengo q normalizar las filas




def transicion_al_azar_uniforme(n, thres):
 
    T = np.random.rand(n, n)
    
  
    for i in range(n):
         for j in range(n):
             if T[i, j] < thres:
                 T[i, j] = 1
             else:
                 T[i, j] = 0

    
    for j in range(n):
        colum = []
        for i in range(n):
            colum.append(T[i][j])
        sumacolumna = sumacol(colum)
        
     
        
        if sumacolumna == 0: #cASO SUMA DE COL ==0, CHEQUEAR ESTO
            
             for i in range(n): 
                 T[i, j] = 1 / n
        else:
             for i in range(n):
                 T[i, j] = T[i, j] / sumacolumna

    return T
    

    
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
"""
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
"""

# Tests L05-QR:

import numpy as np

# --- Matrices de prueba ---
A2 = np.array([[1., 2.],
               [3., 4.]])

A3 = np.array([[1., 0., 1.],
               [0., 1., 1.],
               [1., 1., 0.]])

A4 = np.array([[2., 0., 1., 3.],
               [0., 1., 4., 1.],
               [1., 0., 2., 0.],
               [3., 1., 0., 2.]])

# --- Funciones auxiliares para los tests ---
def check_QR(Q,R,A,tol=1e-10):
    # Comprueba ortogonalidad y reconstrucción
    assert np.allclose(Q.T @ Q, np.eye(Q.shape[1]), atol=tol)
    assert np.allclose(Q @ R, A, atol=tol)

# --- TESTS PARA QR_by_GS2 ---
Q2,R2 = QR_con_GS(A2)
check_QR(Q2,R2,A2)

Q3,R3 = QR_con_GS(A3)
check_QR(Q3,R3,A3)

Q4,R4 = QR_con_GS(A4)
check_QR(Q4,R4,A4)

# --- TESTS PARA QR_by_HH ---
Q2h,R2h = QR_con_GS(A2)
check_QR(Q2h,R2h,A2)

Q3h,R3h = QR_con_HH(A3)
check_QR(Q3h,R3h,A3)

Q4h,R4h = QR_con_HH(A4)
check_QR(Q4h,R4h,A4)

# --- TESTS PARA calculaQR ---
Q2c,R2c = calculaQR(A2,metodo='RH')
check_QR(Q2c,R2c,A2)

Q3c,R3c = calculaQR(A3,metodo='GS')
check_QR(Q3c,R3c,A3)

Q4c,R4c = calculaQR(A4,metodo='RH')
check_QR(Q4c,R4c,A4)

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



# Tests L08
def svd_reducida(A,k="max",tol=1e-15):
    m = A.shape[0]
    n = A.shape[1]
    r = min(n,m)
    if m >= n:
        E, V = diagRH(multiplicacionMatricial(traspuesta(A), A))
        B = multiplicacionMatricial(A, V)
        if B.ndim == 1:
            B = B.reshape(-1, 1)
        U = []
        for j in range(r):
            u = multiplicacionMatricial(A, V[:, j])
            u /= norma(u, 2)
            U[:, j] = u
        """
            v = B[j]
            v /= norma(v, 2)
            U.append(v)
        U = np.array(U)
        """
    else:
        B = multiplicacionMatricial(traspuesta(A), V)
        if B.ndim == 1:
            B = B.reshape(-1, 1)
        U = []
        for j in range(r):
            u = multiplicacionMatricial(A, V[:, j])
            u /= norma(u, 2)
            U[:, j] = u
        """
            v = B[j]
            v /= norma(v, 2)
            U.append(v)
        U = np.array(U)
        """
        E, V = diagRH(multiplicacionMatricial(traspuesta(A), A))
    return U, E, V
    
#print(svd_reducida(genera_matriz_para_test(3)))
    
"""
    A la matriz de interes (de m x n)
    k el numero de valores singulares (y vectores) a retener.
    tol la tolerancia para considerar un valor singular igual a cero
    Retorna hatU (matriz de m x k), hatSig (vector de k valores singulares) y hatV (matriz de n x k)
"""
    #raise NotImplementedError("Implementar svd_reducida!")

    
# Matrices al azar
def genera_matriz_para_test(m,n=2,tam_nucleo=0):
    if tam_nucleo == 0:
        A = np.random.random((m,n))
    else:
        A = np.random.random((m,tam_nucleo))
        A = np.hstack([A,A])
    return(A)

def test_svd_reducida_mn(A,tol=1e-15):
    m,n = A.shape
    hU,hS,hV = svd_reducida(A,tol=tol)
    nU,nS,nVT = np.linalg.svd(A)
    r = len(hS)+1
    assert np.all(np.abs(np.abs(np.diag(hU.T @ nU))-1)<10**r*tol), 'Revisar calculo de hat U en ' + str((m,n))
    assert np.all(np.abs(np.abs(np.diag(nVT @ hV))-1)<10**r*tol), 'Revisar calculo de hat V en ' + str((m,n))
    assert len(hS) == len(nS[np.abs(nS)>tol]), 'Hay cantidades distintas de valores singulares en ' + str((m,n))
    assert np.all(np.abs(hS-nS[np.abs(nS)>tol])<10**r*tol), 'Hay diferencias en los valores singulares en ' + str((m,n))

for m in [2,5,10,20]:
    for n in [2,5,10,20]:
        for _ in range(10):
            A = genera_matriz_para_test(m,n)
            test_svd_reducida_mn(A)


# Matrices con nucleo

m = 12
for tam_nucleo in [2,4,6]:
    for _ in range(10):
        A = genera_matriz_para_test(m,tam_nucleo=tam_nucleo)
        test_svd_reducida_mn(A)

# Tamaños de las reducidas
A = np.random.random((8,6))
for k in [1,3,5]:
    hU,hS,hV = svd_reducida(A,k=k)
    assert hU.shape[0] == A.shape[0], 'Dimensiones de hU incorrectas (caso a)'
    assert hV.shape[0] == A.shape[1], 'Dimensiones de hV incorrectas(caso a)'
    assert hU.shape[1] == k, 'Dimensiones de hU incorrectas (caso a)'
    assert hV.shape[1] == k, 'Dimensiones de hV incorrectas(caso a)'
    assert len(hS) == k, 'Tamaño de hS incorrecto'

