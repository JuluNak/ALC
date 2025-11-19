#Estos son los cambios q me acuerdo:

#EN Pinvhouseholder:
# ANTES ERA W=V.Y , AHORA ES W=Y.V
#MISMO PARA PINVGRAMSMITH 
#EN RESOLVETRIANG (O EL NOMBRE Q TENGA): ANTES RECIBIA FALSE O TRUE
#AHORA RECIBE "INFERIOR " O " SUPERIOR" (NO SE CUAL DE LAS 2 HABBRIA Q USAR , AUNQUE SEAN LITERALMENTELO MISMO
#HAY CAMBIOS EN "QR_GS" Y " QR_HH" (COMO NO ME ACUERDO LOS CAMBIO TE COPIO LAS NUEVAS FUNCIONES ABAJO )
#COMENTE EL ITEM 1 ( DESCOMENTALO EN EL PRINCIPAL PORFI)
#CREERIA Q NADA MAS


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
        soluc=res_tri(R,b,"superior")
        
        for k in range(n):
            VT[k][j]=soluc[k]#MISMO PARA ACA , ESTOS INDICES TAN BIEN??
            
        
    V=traspuesta(VT)
    
    W=multiplicacionMatricial(Y,V)
    
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
        soluc=res_tri(R,b,"superior")
        
        for k in range(n):
            VT[k][j]=soluc[k]
            
        
    V=traspuesta(VT)
 
    W=multiplicacionMatricial(Y,V)
    
    return W


def res_tri(L, b, tipo):
    L = np.array(L, dtype=float)
    n= L.shape[0]
    x = []
    if  tipo == "inferior":  # triangularizacion inferior 
        for i in range(n):
            h = 0.0
            for j in range(i):
                h += L[i, j] * x[j]
            x.append((b[i] - h) / L[i, i])   
    if tipo=="superior":  # triangularizacion superior 
        for i in range(n-1, -1, -1):
            h = 0.0
            for j in range(n-1, i, -1):
                h += L[i, j] * x[n-1-j]
            x.append((b[i] - h) / L[i, i])
        x = x[::-1]  # invierto la lista
    return x



def QR_con_GS(A,tol=1e-12,retorna_nops=False):
    
    columnas = A.shape[1]
    filas = A.shape[0]
    if columnas > filas:
       
        return None
    else:
       
        Q = []
        q1 = A[:,0]/norma(A[:,0], 2)
        qj = q1
        Q.append(qj)
       
        for j in range(1, columnas):
         
            qj = A[:,j]
            for k in range(0,j):
                r = np.inner(A[:,j], Q[k])
                qj = qj - r*Q[k]
            r = norma(qj, 2)
          
            qj = qj/r
            Q.append(qj)
           
          
        Q = np.array(Q)     
        Q = traspuesta(Q) 
         
        R = multiplicacionMatricial(traspuesta(Q), A)
       
        R[np.abs(R) < 1e-12] = 0

        return Q, R

def QR_con_HH(A,tol=1e-12):
    
    n = A.shape[0]
    m = A.shape[1]
    if n < m:
        return None
    else:
        R = A.copy()
        Q = np.eye(n)
        for k in range(m):
            x = R[k:, k]
            a = -np.sign(x[0])*norma(x,2)
            u = x - a*(np.eye(x.shape[0])[0])
            if norma(u, 2) > tol:
                u = u/norma(u,2)
                Hk = np.eye(len(u)) - 2*(np.outer(u, u))
                H = np.eye(n)
                H[k:, k:] = Hk
                R = multiplicacionMatricial(H, R)
                Q = multiplicacionMatricial(Q, H)
                
        Q = Q[:, :m]
        R = R[:m, :]  
        R[np.abs(R) < 1e-12] = 0
               
    return Q, R

