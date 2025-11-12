def intento3(X,Y,L):
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
        A = multiplicar(XT,X)
        L=cholesky(A)
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
        A = multiplicar(X,XT)
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
    
