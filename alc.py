# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 11:51:05 2025

@author: Julu
"""
# --- IMPORTS ---
import numpy as np
import random
import matplotlib.pyplot as plt

# --- LIBRERIAS ---
def esCuadrada(A):
    res: bool = True
    filas = A.shape[0]
    columnas = A.shape[1]
    if not(filas == columnas):
        res = False
    return res

def triangSup(A):
    filas = A.shape[0]
    columnas = A.shape[1]
    res: list = []
    for i in range(filas):
        fila = []
        for j in range(columnas):
            if i < j:
                fila.append(A[i, j])
            else:
                fila.append(0)
        res.append(fila) 
    res = np.array(res)
    return res   

def triangInf(A):
    filas = A.shape[0]
    columnas = A.shape[1]
    res: list = []
    for i in range(filas):
        fila = []
        for j in range(columnas):
            if j < i:
                fila.append(A[i, j])
            else:
                fila.append(0)
        res.append(fila) 
    res = np.array(res)
    return res

def diagonal(A):
    filas = A.shape[0]
    columnas = A.shape[1]
    res: list = []
    for i in range(filas):
        fila = []
        for j in range(columnas):
            if j == i:
                fila.append(A[i, j])
            else:
                fila.append(0)
        res.append(fila) 
    res = np.array(res)
    return res

def traza(A):
    filas = A.shape[0]
    columnas = A.shape[1]
    maximo = min([filas,columnas])
    res: int = 0
    for indice in range(0,maximo):
        res += A[indice, indice]
    return res

def traspuesta(A):
    filas = A.shape[0]
    columnas = A.shape[1]
    res: list = []
    for fil in range(filas):
        vector = []
        for col in range(columnas):
            vector.append(A[fil, col])
        res.append(vector)
    res = np.array(res)
    return res

def esSimetrica(A):
    res: bool = True
    filas = A.shape[0]
    columnas = A.shape[1]
    trasp = traspuesta(A)
    for fil in range(filas):
        for col in range(columnas):
            if not(A[fil, col] == trasp[fil, col]):
                res = False
    return res

def calcularAx(A,x):
    columnas_A = A.shape[1]
    filas_x = x.shape[0]
    res: list = []
    if filas_x == columnas_A:
        for i in range(filas_x):
            suma = 0
            for j in range(filas_x):
                prod = A[i, j] * x[j]
                suma += prod
            res.append(suma)
    res = np.array(res)    
    return res

def intercambiarFilas(A, i, j):
    columnas = A.shape[1]
    for col in range(columnas):
        valor = A[i, col]
        A[i, col] = A[j, col]
        A[j, col] = valor
    return A

def sumar_fila_multiplo(A, i, j, s):
    columnas = A.shape[1]
    for col in range(columnas):
        A[i, col] += A[j, col]*s
    return A

def esDiagonalmenteDominante(A):
    res: bool = True
    longitud = A.shape[0]
    for i in range(longitud):
        suma = 0
        for j in range(longitud):
            if i != j:
                suma += abs(A[i,j])
        print(suma)
        print(abs(A[i,i]))
        if abs(A[i,i]) < suma:
            res = False
    return res      
    
def matrizCirculante(v):
    res: list = []
    longitud = v.shape[0]
    res.append(v)
    vector = v
    for ind in range(1,longitud):
        vector = desplazarDerecha(vector)
        res.append(vector)
    res = np.array(res)
    return res
        
def desplazarDerecha(v):
    res: list = []
    ultimo: int = v[v.shape[0] - 1]
    v = v[:v.shape[0]-1]
    res.append(ultimo)
    for ind in range(0,v.shape[0]):
        res.append(v[ind])
    res = np.array(res)
    return res

def matrizVandermonde(v):
    res: list = []
    n: int = v.shape[0]
    for i in range(0,n):
        vector = v**i
        res.append(vector)
    res = np.array(res)
    return res

def numeroAureo(n):
    f0: int = 0
    f1: int = 1
    a = f0
    b = f1
    for i in range(n):
        temp = a + b   
        a = b          
        b = temp            
    return b/a

def matfibo(m):
    F = [0, 1]  
    for k in range(2, m+1):
        F.append(F[k-1] + F[k-2])
    return F

def hilbert(n):
    res=[]
    for i in range(n):
        fila=[]
        for j in range(n):
            fila.append(1/(i+j+1))
        res.append(fila)
    return res

print(hilbert(3))
        
    
    #EJERCICIO 17
    

    


# Definir polinomio
def evaluar_polinomio1(x):
    return x*5 - x*4 + x*3 - x*2 + x - 1


x = np.linspace(-1, 1, 100) #GENERA 100 PUNTOS ENTRE -1 , 1


y1 = evaluar_polinomio1(x) #EVALUA EL POL


plt.plot(x, y1, label='x^5 - x^4 + x^3 - x^2 + x - 1')
plt.xlabel('x') #muestra el nombre del ejex
plt.ylabel('Valor del polinomio') #nombre eje y
plt.title('Polinomio P1(x) entre -1 y 1') #le pone el titulo
plt.legend() #el cosito q te dice q la funcion azul por ej es esta
plt.grid(True) #cuadricula del fondo
plt.show()

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
    return np.allclose(A, B, atol=np.finfo('float64').eps)

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

def transafin(v, theta, s, b):
    xy1 = np.array([v[0],v[1],1])
    res = afin(theta, s, b) @ xy1
    return res[:2]

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

def descomplu(A):
    n = A.shape[0]
    L = np.zeros((n, n), dtype=float)
    for i in range(n):
        L[i, i] = 1.0 #lo define como uno, pues es lower
    U = A.copy()
    ops = 0  #cont
    for i in range(n):
        if U[i, i] == 0:  # pivote=0
            return None
        for j in range(i+1, n):
            L[j,i] = U[j, i] / U[i, i]
            ops += 1
            for m in range(i, n):
                U[j, m] = U[j, m] - L[j, i] * U[i, m]
                ops += 2  
    return L, U, ops

def trigangularsupyinf(L, b, lower=True):
    L = np.array(L, dtype=float)
    n= L.shape[0]
    x = []
    if lower:  # triangularizacion inferior 
        for i in range(n):
            h = 0.0
            for j in range(i):
                h += L[i, j] * x[j]
            x.append((b[i] - h) / L[i, i])   
    else:  # triangularizacion superior 
        for i in range(n-1, -1, -1):
            h = 0.0
            for j in range(n-1, i, -1):
                h += L[i, j] * x[n-1-j]
            x.append((b[i] - h) / L[i, i])
        x = x[::-1]  # invierto la lista
    return x

def cholesky(A):
    n = len(A)
    L = np.zeros_like(A, dtype=float)
    for i in range(n):
        for j in range(i + 1):
            suma = 0
            for k in range(j):
                suma += L[i][k] * L[j][k]
            if i == j:
                L[i][j] = np.sqrt(A[i][i] - suma)
            else:
                L[i][j] = (A[i][j] - suma) / L[j][j]
    return L

def calcularxA(x,A):
    columnas = A.shape[1]
    res: list = []
    for col in range(columnas):
        res.append(np.inner(x, A[col]))
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
                vector.append(np.inner(A[fil,:], B[:,col]))
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
                r = np.inner(A[j], Q[k])
                qj = qj - r*Q[k]
            r = norma(qj, 2)
            qj = qj/r
            Q.append(qj)
        Q = np.array(Q)
        R = multiplicacionMatricial(traspuesta(Q), A)
        return Q, R

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
            a = -np.sign(x[0])*norma(x,2)
            u = x - a*(np.eye(x.shape[0])[0])
            if norma(u, 2) > tol:
                u = u/norma(u,2)
                Hk = np.eye(len(u)) - 2*(np.outer(u, u))
                H = np.eye(m)
                H[k:, k:] = Hk
                R = multiplicacionMatricial(H, R)
                Q = multiplicacionMatricial(Q, traspuesta(H))
    return Q, R
                    
def calculaQR(A,metodo,tol=1e-12):
    if metodo == 'RH':
        return QR_con_HH(A)
    elif metodo == 'GS':
        return QR_con_GS(A)
    else:
        return None

def metpot2k(A, tol=1e-15, K=1000):
    n = A.shape[1]
    v = np.random.randn(n)
    v /= np.linalg.norm(v)
    w = calcularAx(A, v)
    e = np.inner(w, v)
    k = 0
    while (abs(e - 1) > tol) and (k < K):
        v = w
        w = calcularAx(A, v)
        w /= np.linalg.norm(w)
        e = np.inner(w, v)
        k += 1
    lam = np.inner(calcularxA(w, A), w)
    #error = e - 1
    it = k
    return w, lam, it
"""
def diagRH(A, tol=1e-12, K=1000):
    n = A.shape[0]
    Q_total = np.eye(n)
    A_bis = A.copy()
    k: int = 0
    while k < K:
        Q, R = QR_con_HH(A_bis)
        A_bis = multiplicacionMatricial(R, Q)
        Q_total = multiplicacionMatricial(Q_total, Q)
        off_diag_norm = np.sqrt(np.sum(np.tril(A_bis, -1)**2))
        if off_diag_norm < tol:
            break
        k += 1
    E = np.diag(A_bis)
# Ordenar por magnitud descendente
    idx = np.argsort(np.abs(E))[::-1]
    E = E[idx]
    Q_total = Q_total[:, idx]

    # Asegurar signo determinista de cada autovector
    for i in range(Q_total.shape[1]):
        if Q_total[0, i] < 0:
            Q_total[:, i] *= -1    
    return Q_total, E
"""
def diagRH(A, tol=1e-15, K=1000):
    """
    A: matriz simétrica de tamaño n×n
    tol: tolerancia para decidir convergencia
    K: número máximo de iteraciones

    Retorna (S, D) tal que A ≈ S D S^T
    Si A no es simétrica, retorna None.
    """
    # Verificar que A sea cuadrada
    if A.shape[0] != A.shape[1]:
        return None

    # 🔹 Forzar simetría numérica (evita que retorne None por redondeo)
    A = (A + A.T) / 2

    # 🔹 Comprobación (más laxa)
    if not np.allclose(A, A.T, atol=1e-10):
        return None

    n = A.shape[0]
    A_k = np.copy(A)
    S = np.eye(n)  # acumulador de autovectores

    for k in range(K):
        # Descomposición QR
        Q, R = np.linalg.qr(A_k)

        # Nueva iteración
        A_k = R @ Q
        S = S @ Q  # acumulamos los autovectores

        # Chequear convergencia
        off_diag_norm = np.sqrt(np.sum(np.tril(A_k, -1) ** 2))
        if off_diag_norm < tol:
            break

    # 🔹 Ordenar autovalores y autovectores (opcional pero recomendable)
    D_diag = np.diag(A_k)
    idx = np.argsort(-np.abs(D_diag))  # de mayor a menor
    D = np.diag(D_diag[idx])
    S = S[:, idx]

    return S, D


def transiciones_al_azar_continuas(n):
    mat = []
    # Generación de matrices
    for i in range(n):
        vector = []
        for j in range(n):
            vector.append(random.uniform(0, 1))
        vector = np.array(vector, dtype=float)
        mat.append(vector)
    mat = np.array(mat)
    # Normalizo columnas
    for j in range(n):
        suma_col = np.sum(mat[:, j])
        if suma_col > 0:
            mat[:, j] /= suma_col
        else:
        # caso borde: columna toda en ceros → reemplazo por distribución uniforme
            mat[:, j] = 1.0 / n    
    return mat

def transiciones_al_azar_uniformes(n,thres):
    # Generación de matrices 
    mat = []
    for i in range(n):
        vector = []
        for j in range(n):
            vector.append(random.uniform(0, 1))
        vector = np.array(vector, dtype=float)
        mat.append(vector)
    mat = np.array(mat)    
    # Modificación por thres
    for i in range(n):
        for j in range(n):
            if mat[i,j] < thres:
                mat[i,j] = 1.0
            else:
                mat[i,j] = 0.0
    # Normalizo columnas
    for j in range(n):
        suma_col = np.sum(mat[:, j])
        if suma_col > 0:
            mat[:, j] /= suma_col
        else:
        # caso borde: columna toda en ceros -> reemplazo por distribución uniforme
            mat[:, j] = 1.0 / n
    return mat

def nucleo(A,tol=1e-15):
    # Hago la diagonalizacion
    C, D = diagRH(multiplicacionMatricial(traspuesta(A), A))
    # Consigo los autovalores
    autovalores = []
    for i in range(D.shape[0]):
        autovalores.append(diagonal(D)[i][i])
    # Busco los indices del nucleo
    indices_nucleo = []
    n = len(autovalores)
    for i in range(n):
        if abs(autovalores[i]) <= tol:
            indices_nucleo.append(i)
    
    if not indices_nucleo:
        return np.array([])
    
    nucleo = C[:,indices_nucleo]

    if nucleo.shape[1] == 1:
        nucleo = nucleo.reshape((n, 1))

    return nucleo

def crea_rala(listado,m_filas,n_columnas,tol=1e-15):
    dims = (m_filas,n_columnas)
    if listado == []:
        return {}, dims
    
    elementos = []
    for i in range(len(listado[2])):
        if not abs(listado[2][i]) <= tol:
            elementos.append(((listado[0][i],listado[1][i]), listado[2][i]))
    diccionario = dict(elementos)

    return diccionario, dims

def multiplica_rala_vector(A,v):
    dicc, (m, n) = A  # Desempaquetamos la tupla (diccionario, dimensiones)
    w = np.zeros(m)
    for (i, j), valor in dicc.items():
        w[i] += valor * v[j]
    return w
"""
def svd_reducida(A,k="max",tol=1e-15):
    m = A.shape[0]
    n = A.shape[1]
    r = min(n,m)
    if m >= n:
        V, E = diagRH(multiplicacionMatricial(traspuesta(A), A))
        print(E, V)
        idx = np.argsort(E)[::-1]
        print(idx)
        #print(E ,V)
        B = multiplicacionMatricial(A, V)
    else:
        V, E = diagRH(multiplicacionMatricial(A, traspuesta(A)))
        B = multiplicacionMatricial(traspuesta(A), V)
    if B.ndim == 1:
        B = B.reshape(-1, 1)
    U = []
    for j in range(r):
        u = multiplicacionMatricial(A, V[:, j].reshape(-1, 1))
        u /= norma(u, 2)
        U.append(u.flatten())            # guardo como fila
    U = np.column_stack(U)
    return U, E, V
    
#print(svd_reducida(genera_matriz_para_test(3)))
"""    
def svd_reducida(A, k="max", tol=1e-15):
    m, n = A.shape

    # Paso i: diagonalizar A^T A
    ATA = multiplicacionMatricial(traspuesta(A), A)
    ATA = (ATA + traspuesta(ATA)) / 2  # asegurar simetría
    V, D = diagRH(ATA)

    if V is None or D is None:
        return None, None, None

    # Autovalores y valores singulares
    E = np.diag(D)
    S = np.sqrt(np.abs(E))

    # Filtrar valores singulares pequeños
    mask = S > tol
    S = S[mask]
    V = V[:, mask]

    # Paso ii: calcular U = A V / σ
    U_cols = []
    for j in range(len(S)):
        Avj = multiplicacionMatricial(A, V[:, j].reshape(-1, 1))
        u_j = Avj / S[j]
        # 🔹 normalización explícita
        u_j = u_j / np.linalg.norm(u_j)
        U_cols.append(u_j.flatten())

    U = np.column_stack(U_cols)

    # 🔹 re-ortogonalización por Gram-Schmidt para mayor precisión
    for i in range(U.shape[1]):
        for j in range(i):
            U[:, i] -= np.dot(U[:, j], U[:, i]) * U[:, j]
        U[:, i] /= np.linalg.norm(U[:, i])

    return U, S, V

# --------- FUNCIONES TP ------------

import os

def cargarDataset(carpeta):
    """
    Dado el path a una carpeta con subcarpetas 'train' y 'val',
    cada una conteniendo 'gatos' y 'perros' con archivos .npy de embeddings,
    retorna Xt, Yt, Xv, Yv.
    """

    # ---------- ENTRENAMIENTO ----------
    X_list_t = []
    Y_list_t = []
    path_train = os.path.join(carpeta, "train")

    for clase in ["gatos", "perros"]:
        carpeta_clase = os.path.join(path_train, clase)
        if not os.path.exists(carpeta_clase):
            continue

        for archivo in os.listdir(carpeta_clase):
            if archivo.endswith(".npy"):
                embedding = np.load(os.path.join(carpeta_clase, archivo)).reshape(-1, 1)
                X_list_t.append(embedding)
                if clase == "gatos":
                    Y_list_t.append(np.array([[1], [0]]))
                else:
                    Y_list_t.append(np.array([[0], [1]]))

    Xt = np.hstack(X_list_t)
    Yt = np.hstack(Y_list_t)

    # ---------- VALIDACIÓN ----------
    X_list_v = []
    Y_list_v = []
    path_val = os.path.join(carpeta, "val")

    for clase in ["gatos", "perros"]:
        carpeta_clase = os.path.join(path_val, clase)
        if not os.path.exists(carpeta_clase):
            continue

        for archivo in os.listdir(carpeta_clase):
            if archivo.endswith(".npy"):
                embedding = np.load(os.path.join(carpeta_clase, archivo)).reshape(-1, 1)
                X_list_v.append(embedding)
                if clase == "gatos":
                    Y_list_v.append(np.array([[1], [0]]))
                else:
                    Y_list_v.append(np.array([[0], [1]]))

    Xv = np.hstack(X_list_v)
    Yv = np.hstack(Y_list_v)
    return Xt, Yt, Xv, Yv
    

print(cargarDataset("template-alumnos"))
    


