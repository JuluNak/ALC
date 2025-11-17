# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 11:51:05 2025

@author: Julu
"""
# --- IMPORTS ---
import numpy as np
import random
import matplotlib.pyplot as plt
import zipfile
import io

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
    for col in range(columnas):
        fila = []
        for fil in range(filas):
            fila.append(A[fil, col])
        res.append(fila)
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

def calculaLU(A):
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

def res_tri(L, b, lower=True):
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

def diagRH(A, tol=1e-12, K=1000):
    n = A.shape[0]
    if n == 0:
        return np.eye(0), np.zeros((0, 0))  # caso base: matriz vacía
    
    if n == 1:
        return np.eye(1), A.copy()

    # 1. Primer autovector y autovalor dominante
    v, lam, _ = metpot2k(A, tol, K)
    v = v.reshape(-1, 1)
    v = v / norma(v, 2)

    # 2. Construcción de la reflexión de Householder
    e1 = np.zeros((n, 1))
    e1[0, 0] = 1.0
    u = e1 - np.sign(v[0, 0]) * v
    if norma(u, 2) < 1e-14:  # evitar división por cero
        H = np.eye(n)
    else:
        u = u / norma(u, 2)
        H = np.eye(n) - 2 * np.outer(u, u)

    # 3. Aplicar reflexión
    B = multiplicacionMatricial(multiplicacionMatricial(H, A), traspuesta(H))

    # 4. Caso base recursivo
    if n == 2:
        S = H
        D = B
    else:
        A2 = B[1:, 1:]
        S2, D2 = diagRH(A2, tol, K)

        D = np.zeros_like(A)
        D[0, 0] = lam
        D[1:, 1:] = D2

        S = np.eye(n)
        S[1:, 1:] = S2
        S = multiplicacionMatricial(H, S)

    # 5. Redondear errores numéricos pequeños
    D[np.abs(D) < 1e-14] = 0

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

def svd_reducida(A, k="max", tol=1e-15):
    m, n = A.shape
    r = min(m, n)

    # SVD reducida: A = U Σ Vᵗ
    if m >= n:
        # Caso rectangular alto (más filas que columnas)
        S, D = diagRH(multiplicacionMatricial(traspuesta(A), A))
        V = S
        B = multiplicacionMatricial(A, V)
    else:
        # Caso rectangular ancho
        S, D = diagRH(multiplicacionMatricial(A, traspuesta(A)))
        U = S
        B = multiplicacionMatricial(traspuesta(A), U)

    # Extraer los valores singulares (positivos)
    sigmas = np.sqrt(np.abs(np.diag(D)))

    # Filtrar según k o tolerancia
    if k == "max":
        k = r
    mask = sigmas > tol
    sigmas = sigmas[mask]
    k = min(k, len(sigmas))

    sigmas = sigmas[:k]

    # Construcción de U
    if m >= n:
        U = np.zeros((m, k))
        for j in range(k):
            if sigmas[j] > tol:
                u = B[:, j] / sigmas[j]
                U[:, j] = u
    else:
        V = np.zeros((n, k))
        for j in range(k):
            if sigmas[j] > tol:
                v = B[:, j] / sigmas[j]
                V[:, j] = v

    # Matriz diagonal Σ
    Σ = np.zeros((k, k))
    np.fill_diagonal(Σ, sigmas)
    if m >= n:
        return U, Σ, V
    else:
        return U, Σ, V
"""
# Chequeos internos
    print(f"\n--- Test interno ({m},{n}) ---")
    print("UᵀU ≈ I?\n", np.round(U.T @ U, 4))
    print("VᵀV ≈ I?\n", np.round(V.T @ V, 4))
    print("A ≈ U Σ Vᵀ ? error =", np.linalg.norm(A - U @ np.diag(sigmas) @ V.T))
"""


# ITEM 1 (versión descomprimida)
"""
def cargar_conjunto(ruta_conjunto):
    X_list, Y_list = [], []
    
    for root, _, files in os.walk(ruta_conjunto):
        for file in files:
            if file.endswith('.npy'):
                ruta_completa = os.path.join(root, file)
                
                # Cargar el embedding
                x = np.load(ruta_completa)
                X_list.append(x.reshape(-1, 1))
                
                # Etiqueta según carpeta
                if "cats" in root.lower():
                    y = np.array([[1], [0]])
                elif "dogs" in root.lower():
                    y = np.array([[0], [1]])             
                Y_list.append(y)
    
    # Convertir a matrices
    X = np.hstack(X_list)
    Y = np.hstack(Y_list)
    
    return X, Y


def cargarDataset(carpeta_base):
    # Entrenamiento
    ruta_train = os.path.join(carpeta_base, "template-alumnos", "dataset", "cats_and_dogs", "train")
    Xt, Yt = cargar_conjunto(ruta_train)
    
    # Validación
    ruta_val = os.path.join(carpeta_base, "template-alumnos", "dataset", "cats_and_dogs", "val")
    Xv, Yv = cargar_conjunto(ruta_val)
    
    return Xt, Yt, Xv, Yv

Xt, Yt, Xv, Yv = cargarDataset("template-alumnos")

print("Xt:", Xt.shape)
print("Yt:", Yt.shape)
print("Xv:", Xv.shape)
print("Yv:", Yv.shape)
"""

# ITEM 1 (versión comprimida)
def cargarDataset2(zip_base_path):

    def cargar_conjunto2(zf_inner, tipo):
        X_list = []
        Y_list = []

        nombres = zf_inner.namelist()

        for nombre in nombres:
            # Buscamos los archivos del conjunto correspondiente
            if nombre.endswith(".npy") and ("/" + tipo + "/") in nombre:

                archivo = zf_inner.open(nombre)
                datos = archivo.read()
                archivo.close()

                arr = np.load(io.BytesIO(datos))
                arr = arr.reshape(-1, 1)
                X_list.append(arr)

                # Determinar etiqueta según el nombre
                if "cats" in nombre:
                    y = np.array([[1.0], [0.0]])
                    Y_list.append(y)
                elif "dogs" in nombre:
                    y = np.array([[0.0], [1.0]])
                    Y_list.append(y)

        X = np.hstack(X_list)
        Y = np.hstack(Y_list)
        return X, Y

    # Abrimos el zip principal
    zip_externo = zipfile.ZipFile(zip_base_path, "r")
    nombres_externos = zip_externo.namelist()

    # Buscamos el dataset.zip interno
    dataset_zip_name = None
    for n in nombres_externos:
        if n.endswith("dataset.zip"):
            dataset_zip_name = n

    archivo_interno = zip_externo.open(dataset_zip_name)
    datos_zip_interno = archivo_interno.read()
    archivo_interno.close()
    zip_externo.close()

    # Abrimos el zip interno
    zip_interno = zipfile.ZipFile(io.BytesIO(datos_zip_interno), "r")

    Xt, Yt = cargar_conjunto2(zip_interno, "train")
    Xv, Yv = cargar_conjunto2(zip_interno, "val")

    zip_interno.close()

    return Xt, Yt, Xv, Yv


# Ejemplo de uso
Xt, Yt, Xv, Yv = cargarDataset2("template-alumnos.zip")

print("Xt:", Xt.shape)
print("Yt:", Yt.shape)
print("Xv:", Xv.shape)
print("Yv:", Yv.shape)



    
def pinvSVD(U, S, V, Y):
    m, n = S.shape
    tol = 1e-15
    
    # Creo la matriz S+ (n x m)
    S_plus = np.zeros((n, m))
    
    for i in range(min(m, n)):
        if S[i, i] > tol:
            S_plus[i, i] = 1 / S[i, i]
    
    # Calculo X+
    X_plus = multiplicacionMatricial(multiplicacionMatricial(V, S_plus), traspuesta(U))
    
    # Calculo W
    W = multiplicacionMatricial(X_plus, Y)
    return W  


#item2)

def pinvEcuacionesNormales(X,Y,L):
    #A = Xt . X . Vamos a factorizarlo a A = L . Lt para aplicar cholesky. 
    #IDEA PARA ESTE EJERCICIO : ASUMIMOS (ESPERO Q SEA ASI) QUE YO POR EJEMPLO RECIBO UNA X DE 3X2 ENTONCES LA L Q ME VIENE YA ES LA CHOLESKY DE (XT.X) EN ESTE CASO
    #PARA EL PRIMER CASO : LA IDEA RESOLVER EL SISTEMA L.LT.U=XT (SIENDO U LA PSEUDOINV DE X) , ENTONCES LO Q HAGO ES "DIVIDIRLO EN 2 PARTES" , PRIMER TOMO LA MULTIPLICACION
    #LT.U = B , (LO LLAMO B), PARA Q ME QUEDE L.B=XT ( ESTO ES PQ L Y LT SON MATRICES TRIANGULARES Y TENGO UNA FUNCION Q ME AYUDA CON ESTAS)
    #ENTONCES RESUELVO L.B=XT , OBTENGO B (!ACLARACION!: CUANDO CREO B PARA SABER SUS DIMENSIONES COMO L.B=XT , B TIENE Q TENER LA CANT DE COLUMNAS DE L(=CANT FILAS DE LT)(COMO FILAS) Y LA CANT DE COLUMNAS DE XT (COMO COLUMNAS) )
    #LUEGO CON B OBTENIDA , RESUELVO LT.U=B , PARA CREAR U NECESITO SABER LAS DIMENSIONES Q TENGO Q TENER 
    # POR ENUNCIADO SABEMOS Q U = (XT.X) A LA MENOS . XT , ES DECIR SI SON IGUALES <-> TIENEN LAS MISMAS DIMENSIONES  SI XT ES DE PXN Y X=ES DE NXP --> (XT.X) ES DE PXP (SU INVERSA IGUAL )
    # Y SI HAGO PXP POR (LAS DIMENSIONES DE XT(ES DECIR PXN))=ME QUEDA Q U ES DE PXN(MISMAS DIMS Q XT)
   
    XT=traspuesta(X)
    dimsX=np.shape(X)
    LT=traspuesta(L)
    dimxt=np.shape(XT)
    dimx=np.shape(X)
    filasx=dimsX[0]
    columnasx=dimx[1]
    if dimsX[0]>dimsX[1]:
    
        filasB = XT.shape[0]
        columnasB = XT.shape[1]
        B= np.zeros((filasB, columnasB))
        for j in range(columnasB):
            b = []
            for i in range(filasB):
             b.append(XT[i][j]) #donde b es cada vector columna de la matriz x traspuesta , asi puedo usar la funcion q me agarra una matriz y un vector y resuelve
            b = np.array(b)
            x_sol = resolverTriangular(L, b, "inferior") # aca es donde por cada columna uso la func
            for k in range(filasB):
                B[k][j] = x_sol[k]# aca lo agrego a la matriz
        
        
        filasU=dimxt[0]
        columnasU=dimxt[1]
        U = np.zeros((filasU , columnasU))
        for j in range(columnasU):
             b = []
             for i in range(filasU):
                  b.append(B[i][j])
             b = np.array(b)
             x_sol = resolverTriangular(LT, b, "superior")
             for k in range(filasU):
                U[k][j] = x_sol[k]
            
            
    #luego como W=U.Y
    
        W=multiplicar(Y,U) ###CHEQUEAR ESTO!!!!!!!
    
        return W

    elif dimsX[0]<dimsX[1]:

        #V x (X.XT) = XT
        #V x L.LT = XT
        #traspong todo =
        #L.LT.VT=X
        #LT.VT=B2
        #L.B2=X
        #TENIENDO B2 --> LT.VT=B2 --> HALLARIAMOS VT , LUEGO TRASPONER VT  --> CONSEGUIS V
        
        #IDEA PARA ESTE IF : EN LA CREACION DE B2 ES DECIR DE LT.VT, VAMOS A PENSARLO COMO ANTES O SINO COMO SABEMOS LA DIM DE V ,TMB SABRIAMOS LA DE VT Y POR ENDE LA DE B2
        #SI QUIERO Q L.B2 =X , ENTONCES FILASB2= COLUMNAS L(filas x) Y COLUMNAS B2 =COLUMNAS X ( EN ESTE CASO COLUMNAS X )(EN ESTE CASO L ES CHOLESKY DE X.XT , LUEGO DIM DE L ES nxn)
       
        B2= np.zeros(( filasx,columnasx))
        filasB2=filasx
        columnasB2=columnasx

        for j in range(columnasB2):
            b=[]
            for i in range(filasB2):
                b.append(X[i][j])#CHEQUEAR CON LOS TESTTTT!!!!
            b = np.array(b)
            x_sol = resolverTriangular(L, b, "inferior") # aca es donde por cada columna uso la func
            for k in range(filasB2):
                B2[k][j] = x_sol[k]

        #LT.VT=B
        filasvt=filasx
        columnasvt=columnasx
        VT=np.zeros((filasvt,columnasvt))
     
        for j in range (columnasvt):
            b2=[]
            for i in range(filasvt):
                b2.append(B2[i][j])
            b2=np.array(b2)
            xsol2=resolverTriangular(LT, b2, "superior")
            for k in range(filasvt):
                VT[k][j] = xsol2[k]###VER TEMA INDICES

        V=traspuesta(VT)
        W2=multiplicar(Y,V)
        return W2
    else:
        # pseudo(X) = inv(X)
        #Tenemos que WX = Y. 
        #Solo pasamos X al otro lado. Quedaria W = T.X^-1
        inv_X = inversa(X)
        return multiplicar(Y, inv_X)


#item4)

def pinvHouseHolder(Q,R,Y): 
    #IDEAS: PLANTEEMOS LO Q NOS DICE EL ALGORITMO 3 , YO QUIERO HALLAR PRIMERO X+,LLAMEMOSLO "V",(ALL ESTO USANDO QUE (QR= XT))
    #LUEGO V= Q.(RT)-1 
    #AHORA MULTIPLICO AMBOS LADOS POR RT Y LLEGO A QUE V.RT=Q
    #LUEGO TRASPONGO LA IGUALDAD PARA  Q ME QUEDE MAS COMODO Y PODER USAR LA FUNCION Q"RESOLVERTRIANGULAR)
    #ME QUEDA QUE R.VT=QT
    #FINALMENTE W=V.Y
    #PARENTESIS R ES TRIANG SUP
    #:)
    #ACLARACION SOBRE CREACION DE MATRIZ VT , SOBRE SUS DIMENSIONES : VT TIENE Q TENER LA SIGUIENTE DIMENSION , R ES CUADRADA LUEGO R ES nxn Q ES RECTANGULAR)? DIRIA ASI 
    #HAGAMOS Q ES DE pxn esto implicaria q QT es nxp y si yo estoy BUSCANDO LA DIMESION DE UNA MATRIZ TAL QUE VALGA R.VT=QT 
    #nxn X ?? = nxp , SABEMOS Q LA CANT DE FILAS DE VT TIENE Q SER N PUES SINO NO VALDRIA LA MULTIPLICACION Y SI COMO RESULTANTE ME QUEDA UNA MATRIZ DE PXN 
    #ENTONCES VT TIENE Q SER DE DIM NXP
    #Y POR ESTO ME CREO LA MATRIZ VT CON LAS MISMAS DIM Q QT
    
    
    QT = traspuesta(Q)
    dimsQT=np.shape(QT)
    n=dimsQT[0]
    p=dimsQT[1]
    VT = np.zeros((n, p))
    
    for j in range (p):
        b=[]
        for i in range(n):
            b.append(QT[i][j])#ESTOS INDICEES VAN BIEN??
        b=np.array(b)
        soluc=resolverTriangular(R,b,"superior")
        
        for k in range(n):
            VT[k][j]=soluc[k]#MISMO PARA ACA , ESTOS INDICES TAN BIEN??
            
        
    V=traspuesta(VT)
    
    W=multiplicar(Y,V)
    
    return W
    
def pinvGramSchmidt(Q, R, Y):
    
    QT = traspuesta(Q)
    dimsQT=np.shape(QT)
    n=dimsQT[0]
    p=dimsQT[1]
    VT = np.zeros((n, p))
    
    for j in range (p):
        b=[]
        for i in range(n):
            b.append(QT[i][j])
        b=np.array(b)
        soluc=resolverTriangular(R,b,"superior")
        
        for k in range(n):
            VT[k][j]=soluc[k]
            
        
    V=traspuesta(VT)
    
    W=multiplicar(Y,V)
    
    return W

    
    


