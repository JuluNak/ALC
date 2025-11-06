import numpy as np
#EJERCICIO 2 (ALGORITMO 1)
#HAY 2 MANERAS PQ NO SE CUAL ES LA CORRECTA INTERPRETACION DE LA CONSIGNA
# EN LA MANERA 1 ESTAN LOS 3 CASOS
# EN LA MANERA 2 FALTA EL ULTIMO
#SE PODRA USAR LA INVERSA O ME EXIGE MUCHA COMPLEJIDAD??
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
                L[i][j] = (A[i][j] - suma) /   L[j][j]
    return L
        

 
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

def resolverTriangular(A, b, tipo):
  
    n = A.shape[0]
    x = np.zeros_like(b, dtype=float)

    if tipo == "inferior":
        for i in range(n):
            suma = 0
            for k in range(i):
                suma += A[i][k] * x[k]
            x[i] = (b[i] - suma) / A[i][i]
    else:  
        for i in range(n - 1, -1, -1):
            suma = 0
            for k in range(i + 1, n):
                suma += A[i][k] * x[k]
            x[i] = (b[i] - suma) / A[i][i]

    return x



def multiplicar(A, B):
  
    filasA = A.shape[0]
    colsA = A.shape[1]
    colsB = B.shape[1]
    C = np.zeros((filasA, colsB))
    for i in range(filasA):
        for j in range(colsB):
            suma = 0
            for k in range(colsA):
                suma += A[i][k] * B[k][j]
            C[i][j] = suma
    return C


def intento2(X,Y):
    XT=traspuesta(X)
  
    shapex=np.shape(X)
    
    A=multiplicar(XT,X)
    if shapex[0]>shapex[1]:
        
        L=cholesky(A)
        LT=traspuesta(L)
        #TENDRIA Q RESOLVER EL SISTEMA (L.Lt.U=xt)->(Lt.U=V),luego (l.v=xt)->(y como V=Lt.U)->resuelvo y obtengo U,donde U esla pseudo inversa de X
        
    
  
   
    
        filasB = XT.shape[0]
        columnasB = XT.shape[1]
        B= np.zeros((filasB, columnasB))

   #ACA RESUELVO L.V=XT
        for j in range(columnasB):
            b = []
            for i in range(filasB):
             b.append(XT[i][j]) #donde b es cada vector columna de la matriz x traspuesta , asi puedo usar la funcion q me agarra una matriz y un vector y resuelve
            b = np.array(b)
            x_sol = resolverTriangular(L, b, "inferior") # aca es donde por cada columna uso la func
            for i in range(filasB):
                B[i][j] = x_sol[i]# aca lo agrego a la matriz

  #ACA RESUELVO LT.U=V
        U = np.zeros((filasB, columnasB))
        for j in range(columnasB):
             b = []
             for i in range(filasB):
                  b.append(B[i][j])
             b = np.array(b)
             x_sol = resolverTriangular(LT, b, "superior")
             for i in range(filasB):
                U[i][j] = x_sol[i]
            
            
    #luego como W=U.Y
    
        W=multiplicar(Y,U) ###CHEQUEAR ESTO!!!!!!!
    
        return W


 #EN ESTE IF TENGO Q RESOLVER (U.L.LT=XT) ->LUEGO Traspongo todo cuestion de que  ME QUEDE POR PROPS L.LT.UT=X , esto vale pues como L es simetrica (L.LT=L.LT)
        #luego hay q hacer lo mismo q antes pero en ves de crear B con las dimensiones de XT las creo con las de X, 
        #IDEA!VOY A TOMAR LT.UT=B, LUEGO VOY A RESOLVER L.B=X , DE ACA CONSIGO B Y LUEGO RESUELVO B=LT.UT , ASI OBTENGO UT , LA TRASPONGO PARA OBTENER U Y LISTO :)
 
        
    if shapex[0]<shapex[1]:
        A2 = multiplicar(X, XT)  # XXᵀ
        L = cholesky(A2)
        LT = traspuesta(L)

        filasB2 = X.shape[0]
        columnasB2 = X.shape[1]
        B2 = np.zeros((filasB2, columnasB2))   
        for j in range(columnasB2):
            b2 = []
            for i in range(filasB2):
             b2.append(X[i][j]) #donde b es cada vector columna de la matriz x traspuesta , asi puedo usar la funcion q me agarra una matriz y un vector y resuelve
            b2 = np.array(b2)
            x_sol2 = resolverTriangular(L, b2, "inferior") # aca es donde por cada columna uso la func
            for i in range(filasB2):
                B2[i][j] = x_sol2[i]# aca lo agrego a la matriz 
                
        Ut = np.zeros((filasB2, columnasB2))
        for j in range(columnasB2):
             b2 = []
             for i in range(filasB2):
                  b2.append(B2[i][j])
             b2 = np.array(b2)
             x_sol2 = resolverTriangular(LT, b2, "superior")
             for i in range(filasB2):
                     Ut[i][j] = x_sol2[i]
        U2=traspuesta(Ut)
            
            
    #luego como W=U.Y
    
        W2=multiplicar(Y,U2) ###CHEQUEAR ESTO!!!!!!!
        return W2    
    
    
    
                
                    
    if shapex[0]==shapex[1]:
        
        W3=multiplicar(Y,traspuesta(X)) #PUES EN ESTE CASO PSEUEDOINVERSA DE X ES IGUAL A TRASPUESA DE X
        return W3
    
    
    
    #OPCION 2(NO SE SI ENTENDI EL ENUNCIADO BIEN ANTES)
    
def pinvEcuacionesNormales(L, Y):
  #IDEA! TENGO Q RESOLVER (L.LT.U=Y) , DONDE U ES LA MATRIZ ESTA PSEUDOINVERSA)?CREO, LUEGO HAGO LT.U=B, -> L.B=Y, OBTENGO B , LUEGO DESPEJO U SABIENDO QUE B=LT.U
  #LUEGO W=U.Y
#CONSIDERANDO QUE L.LT=XT.X
    LT = traspuesta(L)
    
    filas = L.shape[0]     
    columnas = L.shape[1]   #PREGUNTAR!!!
    if filas>columnas:
        B = np.zeros((filas, columnas))

        for j in range(columnas):
            b = []
            for i in range(filas):
                b.append(Y[i][j])
            b = np.array(b)
            x_sol = resolverTriangular(L, b, "inferior")
            for i in range(filas):
              B[i][j] = x_sol[i]

        LT = traspuesta(L)
        U = np.zeros((filas, columnas))

        for j in range(columnas):
            b2 = []
            for i in range(filas):
                b2.append(B[i][j])
            b2 = np.array(b)
            x_sol2 = resolverTriangular(LT, b2, "superior")
            for i in range(filas):
             U[i][j] = x_sol2[i]

    
        W=multiplicar(U,Y)
        return W

    if filas<columnas: #IDEA COMO TENGO Q RESOLVER (U.L.LT=Y) , LO Q VOY A HACER PARA Q SEA MAS COMODO ES TRASPONER TTODO EL PARENTESIS 
        #, CUESTION DE Q ME QUEDE (POR PROPS DE TRASPONER) -> (L.LT.UT=YT) , LUEGO TOMO B=LT.UT, Y RESUELVO L.B=YT , -> OBTENGO B , Y LUEGO
        # ENCUENTRO UT , DESPEJANDOLO DE LT.UT=B, FINALMENTE TRASPONGO UT ASI OBTENGO U , Y W=U.Y 
    
    #CONSIDERANDO QUE L.LT=X.XT
        
        YT = traspuesta(Y) #PQ NO ME DEJA PONERLO A LA MISMA ALTURA Q LOS DEMAS??
        filasY = YT.shape[0]
        columnasY = YT.shape[1]
        B = np.zeros((filasY, columnasY))

        for j in range(columnasY):
            b = []
            for i in range(filasY):
                 b.append(YT[i][j])
            b = np.array(b)
            x_sol = resolverTriangular(L, b, "inferior")
            for i in range(filasY):
                  B[i][j] = x_sol[i]

        LT = traspuesta(L)
        Ut = np.zeros((filasY, columnasY))
        for j in range(columnasY):
            b = []
            for i in range(filasY):
              b.append(B[i][j])
            b = np.array(b)
            x_sol = resolverTriangular(LT, b, "superior")
            for i in range(filasY):
                 Ut[i][j] = x_sol[i]

        U2 = traspuesta(Ut)
        W2 = multiplicar(Y, U2)
        return W2
    
    if filas == columnas:
        #falta esto en esta version del ejercicio!!!
        # TENDRE Q USAR LA INVERSA??
    
    
    
    
    
    
