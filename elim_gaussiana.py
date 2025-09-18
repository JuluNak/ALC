#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Eliminacion Gausianna
"""
import numpy as np

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

def sumar_fila_multiplo(A, i, j, s):
    columnas = A.shape[1]
    for col in range(0,columnas):
        A[col][i] += A[col][j]*s
    return A

def matrizIdentidadCuad(n):
    res = []
    for i in range(n):
        vector = []
        for j in range(n):
            if j == i:
                vector.append(1)
            else:
                vector.append(0)
        vector = np.array(vector)
        res.append(vector)
    res = np.array(res)
    return res
    
def determinanteGauss(A):
    n = A.shape[0]
    res = 0
    for i in range(n):
        for j in range(n):
            res += ((-1)**(i+1))*A[i][j]*determinanteLaplace(A[A[:i] + A[i+1:]])

print(determinanteLaplace(np.array([[0, 4, 2, 6],[2, 6, 3, 8],[6, 1, 6, -7],[4, 5, 4, 2]])))

def permutacion(A, f1, f2):
    m = A.shape[1]
    P = matrizIdentidadCuad(m)
    for col in range(0,m):
        valor = A[col][f1]
        valor = P[col][f1]
        A[col][f1] = A[col][f2]
        P[col][f1] = P[col][f2]
        A[col][f2] = valor
        P[col][f2] = valor
    return A, P
    """
    n = A.shape[1]
    nueva_mat: list = []
    A = traspuesta(A)
    for fila in range(n):
        if fila == f1:
            nueva_mat.append(A[f2])
        elif fila == f2:
            nueva_mat.append(A[f1])
        else:
            nueva_mat.append(A[fila])
    nueva_mat = np.array(nueva_mat)
    nueva_mat = traspuesta(nueva_mat)
    
    return nueva_mat
    """
def triangular(A, piv):
    # fabricar la matriz identidad
    m=A.shape[0]
    M = matrizIdentidadCuad(m)

    
    # detectar los multiplos que hacen que la componente de la matriz vaya a cero
    # fabricar la matriz M
    for i in range(piv + 1,m):    
        if A[piv][piv] != 0:
            ri = -(A[piv][i]/A[piv][piv])
        else:
            ri = 0
        A = sumar_fila_multiplo(A, i, piv, ri)
        M[piv][i] += ri
    
    return A, M
    
    """
    et = []
    for j in range(m):
        if j == piv:
            et.append(1)
        else:
            et.append(0)
    """
#print(triangular(np.array([[0, 4, 2, 6],[2, 6, 3, 8],[6, 1, 6, -7],[4, 5, 4, 2]]), 1))
#print(permutacion(np.array([[0, 4, 2, 6],[2, 6, 3, 8],[6, 1, 6, -7],[4, 5, 4, 2]]), 0, 2))
print(triangular(permutacion(np.array([[0, 4, 2, 6],[2, 6, 3, 8],[6, 1, 6, -7],[4, 5, 4, 2]]), 0, 2)[0], 0))
    
"""
def elim_gaussiana(A):
    cant_op = 0
    m=A.shape[0]
    n=A.shape[1]
    Ac = A.copy()
    
    if m!=n:
        print('Matriz no cuadrada')
        return
    
    ## desde aqui -- CODIGO A COMPLETAR
    # claves: detectar que el triangulo inferior de U sea nulo por cada operación
    # 1- Si el elemento iterado en la diagonal es cero, permuta por fila cuya componente de misma columna es el menor. 
    # Si es otro, fabrica una matriz de modificación
    Ps = []
    Ms = []
    for i in range(n-1):
        if Ac[i][i] == 0:
            minimo = i+1
            for fila in range(i+1,n):
                if Ac[i][fila] < minimo:
                    minimo = fila
            Ac, P = permutacion(Ac, i, fila)
            Ps.append(P)
            cant_op += 1
        Ac, M = triangular(Ac, i)
        Ms.append(M)
        cant_op += 1
            
                 




                
    ## hasta aqui, calculando L, U y la cantidad de operaciones sobre 
    ## la matriz Ac
            
    
    return L, U, cant_op


def main():
    n = 7
    B = np.eye(n) - np.tril(np.ones((n,n)),-1) 
    B[:n,n-1] = 1
    print('Matriz B \n', B)
    
    L,U,cant_oper = elim_gaussiana(B)
    
    print('Matriz L \n', L)
    print('Matriz U \n', U)
    print('Cantidad de operaciones: ', cant_oper)
    print('B=LU? ' , 'Si!' if np.allclose(np.linalg.norm(B - L@U, 1), 0) else 'No!')
    print('Norma infinito de U: ', np.max(np.sum(np.abs(U), axis=1)) )

if __name__ == "__main__":
    main()
     
"""