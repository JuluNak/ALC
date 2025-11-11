import numpy as np

#IMPORTANTE!!! , ENTIENDO Q LA ESTRUCTURA DEL TP ES:
#CARGAS TUS DATOS DE X Y Y ...
# LUEGO ENTIENDO Q POR EJEMPLO EN EL DE CHOLESKY SERIA ALGO ASI:
#L=CHOLESKY (XTxX)
#W=pinv.....(L,Y)
# Y ASI CON LOS DEMAS ...
#IMPORTANTE PREGUNTAR ESTO EL MARTES!!!(PARA VER SI LO ESTOY ENTENDIENDO BIEN)
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
def inversa(A):
    return np.linalg.inv(A) #Vale usar linalg? es parte de numpy pero habria que preguntar
     

def intento3(X,Y,L):
    #A = X . Xt . Vamos a factorizarlo a A = L . Lt para aplicar cholesky. 
    XT=traspuesta(X)
    dims=np.shape(X)
    LT=traspuesta(L)
    dimxt=np.shape(XT)
    dimx=np.shape(X)
    filasx=dimx[0]
    columnasx=dimx[1]
    if dims[0]>dims[1]:
        A = multiplicar(X,XT)
        J=cholesky(A)
        filasB = LT.shape[0]
        columnasB = XT.shape[1]
        B= np.zeros((filasB, columnasB))
        for j in range(columnasB):
            b = []
            for i in range(filasB):
             b.append(XT[i][j]) #donde b es cada vector columna de la matriz x traspuesta , asi puedo usar la funcion q me agarra una matriz y un vector y resuelve
            b = np.array(b)
            x_sol = resolverTriangular(L, b, "inferior") # aca es donde por cada columna uso la func
            for i in range(filasB):
                B[i][j] = x_sol[i]# aca lo agrego a la matriz
        
        
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
                U[i][j] = x_sol[k]
            
            
    #luego como W=U.Y
    
        W=multiplicar(Y,U) ###CHEQUEAR ESTO!!!!!!!
    
        return W

    elif dims[0]<dims[1]:

        #V x (X.XT) = XT
        #V x L.LT = XT
        #traspong todo =
        #L.LT.VT=X
        #LT.VT=B
        #L.B=X
        #TENIENDO B --> LT.VT=B --> HALLARIAMOS VT , LUEGO TRASPONER VT  --> CONSEGUIS V
        A = multiplicar(X,XT)
        B2= np.zeros(( columnasx,columnasx))
        filasB2=columnasx
        columnasB2=columnasx

        for i in range(columnasB2):
            b=[]
            for j in range(filasB2):
                b.append(X[j][i])#CHEQUEAR CON LOS TESTTTT!!!!
            b = np.array(b)
            x_sol = resolverTriangular(L, b, "inferior") # aca es donde por cada columna uso la func
            for i in range(filasB2):
                B2[i][j] = x_sol[i]

        #LT.VT=B
        filasvt=columnasx
        columnasvt=columnasx
        VT=np.zeros((filasvt,columnasvt))
     
        for i in range (columnasvt):
            b2=[]
            for j in range(filasvt):
                b2.append(B2[i][j])
            b2=np.array(b2)
            xsol2=resolverTriangular(LT, b2, "superior")
            for k in range(filasvt):
                VT[i][j] = xsol2[k]###VER TEMA INDICES

        V=traspuesta(VT)
        W2=multiplicar(Y,V)
        return W2
    else:
        # pseudo(X) = inv(X)
        #Tenemos que WX = Y. 
        #Solo pasamos X al otro lado. Quedaria W = T.X^-1
        inv_X = inversa(X)
        return multiplicar(Y, inv_X)
        

"""
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
Usar desc LU para calcular la inversa ---> importante para hacer 3er if de intento3()
"""

        
    







"""def intento2(X,Y,L):
    XT=traspuesta(X)
  
  
    shapex=np.shape(X)
    


    if shapex[0]>shapex[1]:
        A=multiplicar(XT,X)
        
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
    
    
    
    #OPCION 2(NO SE SI ENTENDI EL ENUNCIADO BIEN ANTES)"""
    
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



        # ACA ABAJO HAY UNA IDEA SOBRE EL DE GRAM-SMITH Y HOUSEHOLDER:
        # PREGUNTAR Q ONDA!!!
X=......(ALGUN ARRAY)
Y=.......
XT.X=multiplicar(XT,X)
L=cholesky(XT.X)
W=pinvEcuacionesNormales(L,Y)

#ASI SERIA LA ESTRUCTURA ???
    


        def intentoejer5(Q,R,Y): 
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
    
    for i in range (p):
        b=[]
        for j in range(n):
            b.append(QT[i][j])#ESTOS INDICEES VAN BIEN??
        b=np.array(b)
        soluc=resolverTriangular(R,b,"superior")
        
        for k in range(n):
            VT[i][j]=soluc[k]#MISMO PARA ACA , ESTOS INDICES TAN BIEN??
            
        
    V=traspuesta(VT)
    
    W=multiplicar(Y,V)
    
    return W
    
def pinvGramSchmidt(Q, R, Y):
    
    QT = traspuesta(Q)
    dimsQT=np.shape(QT)
    n=dimsQT[0]
    p=dimsQT[1]
    VT = np.zeros((n, p))
    
    for i in range (p):
        b=[]
        for j in range(n):
            b.append(QT[i][j])
        b=np.array(b)
        soluc=resolverTriangular(R,b,"superior")
        
        for k in range(n):
            VT[i][j]=soluc[i]
            
        
    V=traspuesta(VT)
    
    W=multiplicar(Y,V)
    
    return W


# COMO TENGO LA DUDA DE SI ESTA BIEN LA MANERA Q ESTOY TENIENDO DE IMPLEMENTAR EL TP , TENGO Q PREGUNTAR PQ LA MISMA SOLUC ME SIRVE PARA AMBOS CASOS!!!
    
