# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 11:43:12 2025

@author: gdtsu
"""

import nunmpy as np

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

# Principales


# 2
def pinvEcuacionesNormales(X, L, Y):
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
    n = dimsX[0]
    p = dimsX[1]
    LT=traspuesta(L)
    dimxt=np.shape(XT)
    if n>p:
    
        p = XT.shape[0]
        n = XT.shape[1]
        B= np.zeros((p, n))
        for j in range(n):
            b = XT[:, j]
            x_sol = res_tri(L, b, "inferior") # aca es donde por cada columna uso la func
            for k in range(p):
                B[k][j] = x_sol[k]# aca lo agrego a la matriz
        
        
        filasU=dimxt[0]
        columnasU=dimxt[1]
        U = np.zeros((filasU , columnasU))
        for j in range(columnasU):
             b = []
             for i in range(filasU):
                  b.append(B[i][j])
             b = np.array(b)
             x_sol = res_tri(LT, b, "superior")
             for k in range(filasU):
                U[k][j] = x_sol[k]
            
            
    #luego como W=U.Y
    
        W=multiplicacionMatricial(U,Y) ###CHEQUEAR ESTO!!!!!!!
    
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
       
        B2= np.zeros((n, p))
        filasB2= n
        columnasB2=p

        for j in range(columnasB2):
            b=[]
            for i in range(filasB2):
                b.append(X[i][j])#CHEQUEAR CON LOS TESTTTT!!!!
            b = np.array(b)
            x_sol = res_tri(L, b, "inferior") # aca es donde por cada columna uso la func
            for k in range(filasB2):
                B2[k][j] = x_sol[k]

        #LT.VT=B
        filasvt=n
        columnasvt=p
        VT=np.zeros((filasvt,columnasvt))
     
        for j in range (columnasvt):
            b2=[]
            for i in range(filasvt):
                b2.append(B2[i][j])
            b2=np.array(b2)
            xsol2=res_tri(LT, b2, "superior")
            for k in range(filasvt):
                VT[k][j] = xsol2[k]###VER TEMA INDICES

        V=traspuesta(VT)
        W2=multiplicacionMatricial(V,Y)
        return W2
    else:
        # pseudo(X) = inv(X)
        #Tenemos que WX = Y. 
        #Solo pasamos X al otro lado. Quedaria W = T.X^-1
        inv_X = inversa(X)
        return multiplicacionMatricial(Y, inv_X)

# Caso n > p

X = np.array([
    [1,  2,  3],
    [0,  1,  4],
    [2, -1,  0],
    [1,  1,  1],
    [3,  0,  2]
], dtype=float)

X2 = np.array([
    [1,  2,  3],
    [0,  1,  4]
], dtype=float)

Y = np.array([
    [1, 0,2],
    [0, 1,3],
    [1, 1,4],
    [2, 1,5],
    [0, 2,6]
], dtype=float)


L = cholesky(multiplicacionMatricial(traspuesta(X), X))

print(pinvEcuacionesNormales(X, L, Y))

