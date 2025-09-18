# -*- coding: utf-8 -*-
"""
Created on Sun Aug 31 11:54:12 2025

@author: Julu
"""

import numpy as np
import matplotlib.pyplot as plt

# 1
print(0.3 + 0.25)
print(0.3 - 0.25)


# 2
print(np.sqrt(2)**2-2)
x=4
print(np.sqrt(2*x**2 + 1) - 1)
print((2*x**2)/(np.sqrt(2*x**2 + 1) + 1))

# Valores de x
x = np.linspace(1e14, 1e16, 100)

# Expresiones equivalentes
f = np.sqrt(x**2 + 1) - x
g = 1 / (np.sqrt(x**2 + 1) + x)

# Graficar
plt.figure(figsize=(8,5))
plt.plot(x, f, label=r"$\sqrt{x^2+1} - x$", marker='o', markersize=4, linestyle="--")
plt.plot(x, g, label=r"$\frac{1}{\sqrt{x^2+1}+x}$", marker='x', markersize=4)
plt.yscale("log")  # Escala logarítmica para ver mejor
plt.xlabel("x")
plt.ylabel("Valor de la expresión")
plt.title("Comparación de estabilidad numérica")
plt.legend()
plt.grid(True)
plt.show()


# 3
def acumulacionError():
    l=[]
    x1 = np.sqrt(2)
    x = x1
    for i in range(100):
        x = x*x/np.sqrt(2)
        l.append(x)
    plt.plot(l)
    return l

print(acumulacionError())

"""
# 4
n = 7

s = np.float32(0)
for i in range(1,10**n+1):
    s = s + np.float32(1/i)
print("suma = ", s)

s = np.float32(0)
for i in range(1,5*10**n+1):
    s = s + np.float32(1/i)
print("suma = ", s)

s = np.float32(0)
for i in range(2*10**n,0,-1):
    s = s + np.float32(1/i)
print("suma = ", s)
"""

# 5
def matricesIguales(A, B):
    res: bool = True
    for fil in range(A.shape[0]):
        for col in range(A.shape[1]):
            if np.float32(A[fil][col]) != np.float32(B[fil][col]):
                res = False
    return res

A = np.array([[4,2,0],[2,7,5],[1,9,22/3]])
L = np.array([[1,0.5,0],[0,1,5/6],[0,0,1]])
U = np.array([[4,0,0],[2,6,0],[1,8.5,0.25]])
print(matricesIguales(A, L@U))


# 6

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

T = np.array([[5,4,3],[4,9,4],[3,4,5]])
A = np.array(np.random.rand(4,4))
print(esSimetrica(A.T@A))
print(esSimetrica(A.T@((A*0.25)/0.25)))
print(esSimetrica(A.T@((A*0.2)/0.2)))
