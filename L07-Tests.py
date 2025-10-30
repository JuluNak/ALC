import numpy as np
import random

# AUXILIARES

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

def signo(a):
    if a < 0:
        return -1
    else:
        return 1
    
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
                vector.append(np.inner(A[fil,:], B[:,col]))
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
                Hk = np.eye(len(u)) - 2*(np.outer(u, u))
                H = np.eye(m)
                H[k:, k:] = Hk
                R = multiplicacionMatricial(H, R)
                Q = multiplicacionMatricial(Q, traspuesta(H))
    return Q, R

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

def diagonal(A):
    filas = A.shape[0]
    columnas = A.shape[1]
    maximo = min([filas,columnas])
    res: list = []
    for indice in range(0,maximo):
        vector = []
        vector.append(A[indice][indice])
        res.append(vector)
    res = np.array(res)
    return res

def calcularAx(A,x):
    columnas = A.shape[1]
    filas_A = A.shape[0]
    filas_x = x.shape[0]
    col = 0
    res: list = []
    while col < columnas and col < filas_x:
        fil = 0
        suma = 0
        for fil in range(0,filas_A):
            prod = A[col][fil] * x[fil]
            suma += prod
        col += 1
        res.append(suma)
    res = np.array(res)    
    return res

# PRINCIPALES

def transiciones_al_azar_continuas(n):
    """
    n la cantidad de filas (columnas) de la matriz de transición.
    Retorna matriz T de n x n normalizada por columnas, y con entradas al azar en el intervalo [0,1]
    """
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
        
  #  raise NotImplementedError("Implementar transiciones_al_azar_continuas")

def transiciones_al_azar_uniformes(n,thres):
    """
    n la cantidad de filas (columnas) de la matriz de transición.
    thres probabilidad de que una entrada sea distinta de cero.
    Retorna matriz T de n x n normalizada por columnas. 
    El elemento i,j es distinto de cero si el número generado al azar para i,j es menor o igual a thres. 
    Todos los elementos de la columna $j$ son iguales 
    (a 1 sobre el número de elementos distintos de cero en la columna).
    """
    # Generación de matrices 
    # Opcion chat: mat = np.random.uniform(0, 1, (n, n))
    
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
   # raise NotImplementedError("Implementar transiciones_al_azar_uniformes")

def nucleo(A,tol=1e-15):
    """
    A una matriz de m x n
    tol la tolerancia para asumir que un vector esta en el nucleo.
    Calcula el nucleo de la matriz A diagonalizando la matriz traspuesta(A) * A (* la multiplicacion matricial), usando el medodo diagRH. El nucleo corresponde a los autovectores de autovalor con modulo <= tol.
    Retorna los autovectores en cuestion, como una matriz de n x k, con k el numero de autovectores en el nucleo.
    """
    # Hago la diagonalizacion
    C, D = diagRH(multiplicacionMatricial(traspuesta(A), A))
    # Consigo los autovalores
    autovalores = diagonal(D)
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
    # raise NotImplementedError("Implementar nucleo")

def crea_rala(listado,m_filas,n_columnas,tol=1e-15):
    """
    Recibe una lista listado, con tres elementos: lista con indices i, lista con indices j, y lista con valores A_ij de la matriz A. Tambien las dimensiones de la matriz a traves de m_filas y n_columnas. Los elementos menores a tol se descartan.
    Idealmente, el listado debe incluir unicamente posiciones correspondientes a valores distintos de cero. Retorna una lista con:
    - Diccionario {(i,j):A_ij} que representa los elementos no nulos de la matriz A. Los elementos con modulo menor a tol deben descartarse por default. 
    - Tupla (m_filas,n_columnas) que permita conocer las dimensiones de la matriz.
    """
    dims = (m_filas,n_columnas)
    if listado == []:
        return {}, dims
    
    elementos = []
    for i in range(len(listado[2])):
        if not listado[2][i] <= tol:
            elementos.append(((listado[0][i],listado[1][i]), listado[2][i]))
    diccionario = dict(elementos)

    return diccionario, dims
    #raise NotImplementedError("Implementar crea_rala")

def multiplica_rala_vector(A,v):
    """
    Recibe una matriz rala creada con crea_rala y un vector v. 
    Retorna un vector w resultado de multiplicar A con v
    """
    return calcularAx(A, v)
    #raise NotImplementedError("Implementar multiplica_rala_vector")



def es_markov(T,tol=1e-6):
    """
    T una matriz cuadrada.
    tol la tolerancia para asumir que una suma es igual a 1.
    Retorna True si T es una matriz de transición de Markov (entradas no negativas y columnas que suman 1 dentro de la tolerancia), False en caso contrario.
    """
    n = T.shape[0]
    for i in range(n):
        for j in range(n):
            if T[i,j]<0:
                return False
    for j in range(n):
        suma_columna = sum(T[:,j])
        if np.abs(suma_columna - 1) > tol:
            return False
    return True

def es_markov_uniforme(T,thres=1e-6):
    """
    T una matriz cuadrada.
    thres la tolerancia para asumir que una entrada es igual a cero.
    Retorna True si T es una matriz de transición de Markov uniforme (entradas iguales a cero o iguales entre si en cada columna, y columnas que suman 1 dentro de la tolerancia), False en caso contrario.
    """
    if not es_markov(T,thres):
        return False
    # cada columna debe tener entradas iguales entre si o iguales a cero
    m = T.shape[1]
    for j in range(m):
        non_zero = T[:,j][T[:,j] > thres]
        # all close
        close = all(np.abs(non_zero - non_zero[0]) < thres)
        if not close:
            return False
    return True


def esNucleo(A,S,tol=1e-5):
    """
    A una matriz m x n
    S una matriz n x k
    tol la tolerancia para asumir que un vector esta en el nucleo.
    Retorna True si las columnas de S estan en el nucleo de A (es decir, A*S = 0. Esto no chequea si es todo el nucleo
    """
    for col in S.T:
        res = A @ col
        if not np.allclose(res,np.zeros(A.shape[0]), atol=tol):
            return False
    return True

## TESTS
# transiciones_al_azar_continuas
# transiciones_al_azar_uniformes
for i in range(1,100):
    T = transiciones_al_azar_continuas(i)
    assert es_markov(T), f"transiciones_al_azar_continuas fallo para n={i}"
    
    T = transiciones_al_azar_uniformes(i,0.3)
    assert es_markov_uniforme(T), f"transiciones_al_azar_uniformes fallo para n={i}"
    # Si no atajan casos borde, pueden fallar estos tests. Recuerden que suma de columnas DEBE ser 1, no valen columnas nulas.
    T = transiciones_al_azar_uniformes(i,0.01)
    assert es_markov_uniforme(T), f"transiciones_al_azar_uniformes fallo para n={i}"
    T = transiciones_al_azar_uniformes(i,0.01)
    assert es_markov_uniforme(T), f"transiciones_al_azar_uniformes fallo para n={i}"
  
# nucleo
A = np.eye(3)
S = nucleo(A)
assert S.shape[0]==0, "nucleo fallo para matriz identidad"
A[1,1] = 0
S = nucleo(A)
msg = "nucleo fallo para matriz con un cero en diagonal"
assert esNucleo(A,S), msg
assert S.shape==(3,1), msg
assert abs(S[2,0])<1e-2, msg
assert abs(S[0,0])<1e-2, msg

v = np.random.random(5)
v = v / np.linalg.norm(v)
H = np.eye(5) - np.outer(v, v)  # proyección ortogonal
S = nucleo(H)
msg = "nucleo fallo para matriz de proyeccion ortogonal"
assert S.shape==(5,1), msg
v_gen = S[:,0]
v_gen = v_gen / np.linalg.norm(v_gen)
assert np.allclose(v, v_gen) or np.allclose(v, -v_gen), msg
  
# crea rala
listado = [[0,17],[3,4],[0.5,0.25]]
A_rala_dict, dims = crea_rala(listado,32,89)
assert dims == (32,89), "crea_rala fallo en dimensiones"
assert A_rala_dict[(0,3)] == 0.5, "crea_rala fallo"
assert A_rala_dict[(17,4)] == 0.25, "crea_rala fallo"
assert len(A_rala_dict) == 2, "crea_rala fallo en cantidad de elementos"

listado = [[32,16,5],[3,4,7],[7,0.5,0.25]]
A_rala_dict, dims = crea_rala(listado,50,50)
assert dims == (50,50), "crea_rala fallo en dimensiones con tol"
assert A_rala_dict.get((32,3)) == 7
assert A_rala_dict[(16,4)] == 0.5
assert A_rala_dict[(5,7)] == 0.25

listado = [[1,2,3],[4,5,6],[1e-20,0.5,0.25]]
A_rala_dict, dims = crea_rala(listado,10,10)
assert dims == (10,10), "crea_rala fallo en dimensiones con tol"
assert (1,4) not in A_rala_dict
assert A_rala_dict[(2,5)] == 0.5
assert A_rala_dict[(3,6)] == 0.25
assert len(A_rala_dict) == 2

# caso borde: lista vacia. Esto es una matriz de 0s
listado = []
A_rala_dict, dims = crea_rala(listado,10,10)
assert dims == (10,10), "crea_rala fallo en dimensiones con lista vacia"
assert len(A_rala_dict) == 0, "crea_rala fallo en cantidad de elementos con lista vacia"

# multiplica rala vector
listado = [[0,1,2],[0,1,2],[1,2,3]]
A_rala = crea_rala(listado,3,3)
v = np.random.random(3)
v = v / np.linalg.norm(v)
res = multiplica_rala_vector(A_rala,v)
A = np.array([[1,0,0],[0,2,0],[0,0,3]])
res_esperado = A @ v
assert np.allclose(res,res_esperado), "multiplica_rala_vector fallo"

A = np.random.random((5,5))
A = A * (A > 0.5) 
listado = [[],[],[]]
for i in range(5):
    for j in range(5):
        listado[0].append(i)
        listado[1].append(j)
        listado[2].append(A[i,j])
        
A_rala = crea_rala(listado,5,5)
v = np.random.random(5)
assert np.allclose(multiplica_rala_vector(A_rala,v), A @ v)
