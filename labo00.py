# -*- coding: utf-8 -*-
"""
Created on Sun Aug 24 14:34:10 2025

@author: Julu
"""

## > <    

import numpy as np
import matplotlib.pyplot as plt  # librería para graficar

def esCuadrada(A):
    res: bool = True
    filas = A.shape[0]
    columnas = A.shape[1]
    if not(filas == columnas):
        res = False
    return res

A = np.array([[5,3,11],[15,9,33],[20,12,44]])
print(esCuadrada(A))


def triangSup(A):
    filas = A.shape[0]
    columnas = A.shape[1]
    col = 0
    res: list = []
    for col in range(0, columnas):
        vector = []
        fil = 0
        for fil in range(0, filas):
            if fil < col:
                vector.append(A[col][fil])
            elif col < fil:
                vector.append(0)
            fil += 1
        res.append(vector) 
    res = np.array(res)
    return res

print(triangSup(A))
    

def triangInf(A):
    filas = A.shape[0]
    columnas = A.shape[1]
    col = 0
    res: list = []
    for col in range(0, columnas):
        vector = []
        fil = 0
        for fil in range(0, filas):
            if col < fil:
                vector.append(A[col][fil])
            elif fil < col:
                vector.append(0)
            fil += 1
        res.append(vector) 
    res = np.array(res)
    return res
"""
        fil = filas - 1
        while col < fil and fil > -1:
            vector.append(A[col][fil])
            fil -= 1
            """
print(triangInf(A))    


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

print(diagonal(A))  


def traza(A):
    filas = A.shape[0]
    columnas = A.shape[1]
    maximo = min([filas,columnas])
    res: int = 0
    for indice in range(0,maximo):
        res = res + A[indice][indice]
    return res

print(traza(A))


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
        res.append(vector)
    res = np.array(res)
    return res

print(traspuesta(A))


def esSimetrica(A):
    res: bool = True
    filas = A.shape[0]
    columnas = A.shape[1]
    fil = 0
    col = 0
    trasp = traspuesta(A)
    for fil in range(0, filas):
        for col in range(0, columnas):
            if not(A[col][fil] == trasp[col][fil]):
                res = False
    return res

print(esSimetrica(A))
B = np.array([[5,4,3],[4,9,4],[3,4,5]])
print(esSimetrica(B))


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
        
x = np.array([1,1,1])
print(calcularAx(A,x))


def intercambiarFilas(A, i, j):
    columnas = A.shape[1]
    for col in range(0,columnas):
        valor = A[col][i]
        A[col][i] = A[col][j]
        A[col][j] = valor
    return A

i: int = 0
j: int = 1
print(intercambiarFilas(A, i, j))


def sumar_fila_multiplo(A, i, j, s):
    columnas = A.shape[1]
    for col in range(0,columnas):
        A[col][i] += A[col][j]*s
    return A

A = np.array([[5,3,11],[15,9,33],[20,12,44]])
s = 2
print(sumar_fila_multiplo(A, i, j, s))


def esDiagonalmenteDominante(A):
    res: bool = True
    diag = diagonal(A)
    tras = traspuesta(A)
    longitud = diag.shape[0]
    for col in range(0, longitud):
        if sumaElementosAbs(tras[col]) - 2*valorAbs(diag[col]) >= 0:
            res = False
    return res      
    
def sumaElementosAbs(v):
    longitud = v.shape[0]
    res: int = 0
    for ind in range(0,longitud):
        res += valorAbs(v[ind])
    return res
    
def valorAbs(s):
    if s < 0:
        s = -s
    return s
    
C = np.array([[-36,3,11],[15,-16,-33],[20,12,45]])
print(esDiagonalmenteDominante(C))


def matrizCirculante(v):
    res: list = []
    longitud = v.shape[0]
    res.append(v)
    vector = v
    for ind in range(1,longitud):
        vector = desplazarDerecha(vector)
        res.append(vector)
    res = np.array(res)
    res = traspuesta(res)
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
    
v = np.array([1,2,3])
print(matrizCirculante(v))


def matrizVandermonde(v):
    res: list = []
    n: int = v.shape[0]
    for i in range(0,n):
        vector = v**i
        res.append(vector)
    res = np.array(res)
    res = traspuesta(res)
    return res
        
print(matrizVandermonde(v))


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

# Número de términos a calcular
N = 20

# Generamos las aproximaciones
aprox = [numeroAureo(i) for i in range(2, N+1)]  # desde n=2 para evitar división por cero
indices = list(range(2, N+1))

# Valor real del número áureo
phi = (1 + 5**0.5) / 2

# Gráfico
plt.figure(figsize=(8,5))
plt.plot(indices, aprox, marker='o', label='Aproximación $F_{n+1}/F_n$')
plt.axhline(y=phi, color='r', linestyle='--', label='Número áureo φ')
plt.title('Aproximación del número áureo con la sucesión de Fibonacci')
plt.xlabel('n')
plt.ylabel('F(n+1)/F(n)')
plt.legend()
plt.grid(True)
plt.show()

print(numeroAureo(8))


def matrizFibonacci(n):
    res: list = []
    vector: list = []
    f0: int = 0
    f1: int = 1
    a = f0
    b = f1
    for i in range(n):
        temp = a + b
        a = b    
        b = temp
        vector.append(temp)
        res.append(vector)
    res = np.array(res)
    return res
    
print(matrizFibonacci(8))
        


    