#COSAS Q FALTAN: HACER LA FUNCION INVERSA PARA EL EJERCICIO 2 ITEM 3
# VER Q PASE TODOS LOS TEST
# DEL EJERCICIO 4 LAS FUNCIONES ESTAN HECHAS Y DIRIA Q ESTAN BIEN COMO YA PARA PONERLAS EN EL MODULOALC 
# FALTARIA TESTEAR
#ESTO Q SERIAN LAS FUNCIONES IRIAN EN EL MODULOALC.PY
#CON TODAS LAS FUNCIONES DEL MODULO
#ENTIENDO Q EN OTRA CARPETA ES DONDE HAY Q PONER ONDA:
# CHOLESKY(A)=L, SIENDO A=X.XT(POR EJ ES UN CASO)
# Y AHI BUENO W=PINVECNORMALES(X,Y,L) ETC ETC ETC
#FALTARIA HACER ESA CARPETA

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


# COMO TENGO LA DUDA DE SI ESTA BIEN LA MANERA Q ESTOY TENIENDO DE IMPLEMENTAR EL TP , TENGO Q PREGUNTAR PQ LA MISMA SOLUC ME SIRVE PARA AMBOS CASOS!!!


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



def espseudoinv(X,Xp,tol=1e-8):

    matrizprimeracond=multiplicar(multiplicar(X,Xp),X)
    matrizsegundacond=multiplicar(multiplicar(Xp,X),Xp)
    matrizterceracond=traspuesta(multiplicar(X,Xp))
    matrizcuartacond=traspuesta(multiplicar(Xp,X))

    XXp = multiplicar(X, Xp)
    XpX = multiplicar(Xp, X)
    
    def tolerancia(y, z):
        diferencia = np.abs(y - z)
        maxerror = np.max(diferencia)
        return maxerror < tol
    
    condicion1= tolerancia(matrizprimeracond,X)
    condicion2=tolerancia(matrizsegundacond,Xp)
    condicion3=tolerancia(matrizterceracond,XXp))
    condicion4=tolerancia(matrizcuartacond,XpX)
    
    if condicion1 and condicion2 and condicion3 and condicion4 :
        return True
    else :
        return False





Q, R = QR_con_GS(traspuesta(Xt)) 
W = pinvGramSchmidt(Q, R, Yt)
    
        







    
