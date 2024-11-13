#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Created on Fri May 03 23:07:02 2024

@author: searCampuzano
@email: searcampuzanob@gmail.com

Implementación del algoritmo de Rosenbrock a la funcion de distribución
doble Gumbel.

"""

import pandas as pd
import numpy as np
import pylab as pl
import math

def DDG(gastos):
    """ Implementación del algoritmo de Rosenbrock a la funcion de distribución
    doble Gumbel.
    """     
    #********************************************************************************************
    # Inicia la función bisección
    #********************************************************************************************
    
    def biseccion(vector, a1, a2, b1, b2, P, m):
        
        tol = 1e-13
        gEstimados = np.zeros ((m, 1))
    
        for j in range (m):    
            
            b = 100000
            a = 0       
            fx = vector[j]
            fa = fx-((P*math.exp(-1*math.exp((b1-a)/a1))) + \
                ((1-P)*math.exp(-1*math.exp((b2-a)/a2))))
    
            for i in range (1000):  
                
                p = a + ((b - a) / 2)
                fp = fx-((P*math.exp(-1*math.exp((b1-p)/a1))) + \
                    ((1-P)*math.exp(-1*math.exp((b2-p)/a2))))      # f(p)
                
                if ((fp == 0) | (((b - a) / 2) < tol)):
                    break
                if ((fa * fp) > 0):
                    a = p
                    fa = fp
                if ((fa * fp) < 0):
                    b = p
    
            gEstimados[j] = p                 # Datos estimados
    
        return(gEstimados)
        
    #********************************************************************************************
    # Inicia la función errorC
    #********************************************************************************************
    
    def errorC(m, gastos, tRpNEx, pIn):
        E1 = 0
        E = np.zeros(m)
        for i in range (m):                 # F(qMedido)
            
            E[i] = ((pIn[4]*math.exp(-1*math.exp((pIn[2]-gastos[i])/pIn[0]))) + \
            ((1-pIn[4])*math.exp(-1*math.exp((pIn[3]-gastos[i])/pIn[1]))))
            sumaC = (tRpNEx[i] - E[i])**2  # Función de errores cuadrados pesados
            E1 = E1 + sumaC
            
        return(E1)
        
    #********************************************************************************************
    # Inicia el programa
    #********************************************************************************************

    if (type(gastos) == str):
        print ("\aPara poder realizar el análisis de datos seleccione un archivo valido")
    else:
        tamanio = gastos.shape
        if (len(tamanio) == 1):            # Si el array viene de una dimensión, este se
            gastos.resize(len(gastos), 1)  # redimensiona
        gastos01 = gastos.copy()
        m = gastos.size    
    
    matriz01 = np.zeros ((m, 9))
    
    matriz02 = np.zeros ((m, 3))
    matriz02 [:] = -999
    matriz02 [0, 0] = 2
    matriz02 [1, 0] = 5
    matriz02 [2, 0] = 10
    matriz02 [3, 0] = 20
    matriz02 [4, 0] = 50
    matriz02 [5, 0] = 100
    matriz02 [6, 0] = 200
    matriz02 [7, 0] = 500
    matriz02 [8, 0] = 1000
    matriz02 [9, 0] = 2000
    matriz02 [10, 0] = 5000
    matriz02 [11, 0] = 10000
    n = 12

    EEstandart = np.full((m, 2), -999, float)

    for i in range (n):
        matriz02 [i ,1] = 1.0 - (1.0 / matriz02 [i, 0]) #Columna de F(x)
    
    for i in range (m):
        
        matriz01 [i, 0] = i + 1           # Se crea la columna con el No de orden columna '1'
        matriz01 [i, 1] = gastos01.max()  # Se crea la columna con los gastos "columna 2"
        j = gastos01.argmax()             # registrados y ordenados
        gastos01 [j, 0] = -1              # en orden ascendente
        matriz01 [i, 2] = (m + 1) / matriz01 [i, 0]  # Se crea la columna con Tr "columna 3"
        matriz01 [i, 3] = 1 - (1 / matriz01 [i, 2])  # Se crea la columna de F(x), "a partir de Tr (columna 4)"
        matriz01 [i, 4] = 1 / matriz01[i, 2] # F(qEmpirico)=k/n+1
    
    #********************************************************************************************
    # Se decide cuantos valores son los ciclonicos
    #********************************************************************************************
    
    tR = matriz01 [:, 2]
    dReg = matriz01 [:, 1]
    
    pl.subplot(1, 1, 1)
    pl.scatter(tR, dReg, label='Datos Registrados', marker='.', color='b')
    pl.legend(loc="best")
    pl.title("Periodo de retorno VS gasto")
    pl.ylabel("Gastos(m³/s)")
    pl.xlabel("Periodo de retorno, (años)")
    pl.semilogx(True)
    pl.grid(True, which='both', color='g', linestyle='-', linewidth=0.5)
    pl.show()
    ndc = int(input("Número de datos de origen ciclonico\n (debe de ser mayor a 1): "))
    
    #********************************************************************************************
    # estimación de párametros
    #********************************************************************************************
    
    media1 = 0
    media2 = 0
    foo = 0
    ds1 = 0
    ds2 = 0
    P = (m - ndc) / m                   # Estimación probabilidad
    
    for i in range (ndc):               # Estimación media2
        media2 = media2 + matriz01[i, 1]
        
    media2 = media2 / ndc
    
    for i in range (ndc):
        
        foo = (matriz01[i, 1] - media2) ** 2
        ds2 = ds2 + foo
        
    ds2 = (ds2/(ndc-1)) ** 0.5          # Estimación desviación estandart 2
    foo = 0
    
    for i in range (m - ndc):           # Estimación media1
        media1 = media1 + matriz01[i + ndc, 1]
        
    media1 = media1 / (m-ndc)
    
    for i in range (m - ndc):
        
        foo = (matriz01[i + ndc, 1] - media1) ** 2
        ds1 = ds1 + foo
        
    ds1 = (ds1/((m-ndc) - 1)) ** 0.5    # Estimación desviación estandart 1
    a1 = ((6)**0.5/math.pi) * ds1
    a2 = ((6)**0.5/math.pi) * ds2
    b1 = media1 - (0.5772 * a1)
    b2 = media2 - (0.5772 * a2)
    
                    # Se imprimen los parametros iniciales
    
    print("\x1b[1;31;m" + "\n             Parametros Iniciales")
    print("\x1b[1;32;m" + "\nmedia 1:", "\x1b[1;35;m", media1)
    print("\x1b[1;32;m" + "desviación estandart 1: ", "\x1b[1;35;m",ds1)
    print("\x1b[1;32;m" + "\nmedia 2: ", "\x1b[1;35;m", media2)
    print("\x1b[1;32;m" + "desviación estantart 2: ", "\x1b[1;35;m",ds2)
    print("\x1b[1;32;m" + "\na1: ", "\x1b[1;35;m",a1)
    print("\x1b[1;32;m" + "b1: ", "\x1b[1;35;m",b1)
    print("\x1b[1;32;m" + "a2: ", "\x1b[1;35;m",a2)
    print("\x1b[1;32;m" + "b2: ", "\x1b[1;35;m",b2, "\n")
    print("\x1b[1;32;m" + "Probabilidad: ", "\x1b[1;35;m",P, "\n")
    
    gEstimados = biseccion(matriz01[:, 3], a1, a2, b1, b2, P, m) # Se invoca a la función "bisección", para obtener los gastos estimados.
    matriz01[:, 7] = gEstimados[:, 0]
    
    E = 0
    
    for i in range (m): 
                    # F(qMedido)
        matriz01[i, 5] = ((P*math.exp(-1*math.exp((b1-matriz01[i, 1])/a1))) + \
                ((1-P)*math.exp(-1*math.exp((b2-matriz01[i, 1])/a2))))
        matriz01[i, 6] = (matriz01[i, 3] - matriz01[i, 5])**2  # Función de errores cuadrados pesados
        E = E + matriz01[i, 6]
        
    print("\x1b[1;32;m" + "E = ", "\x1b[1;35;m", E)
    
    #********************************************************************************************

                    # Algoritmo de Rosenbrock

    #********************************************************************************************
    ''' 
                Condiciones iniciales del algoritmo de Rosenbrock
    '''             
    #********************************************************************************************
    
    pIn0 = np.array([a1, a2, b1, b2, P])
    pIn = np.array([a1, a2, b1, b2, P])    # punto inicial
    nVec = 5                               # Número de vectores
    da1 = np.array([1, 0, 0, 0, 0])        # Vectores dirección
    da2 = np.array([0, 1, 0, 0, 0])
    db1 = np.array([0, 0, 1, 0, 0])
    db2 = np.array([0, 0, 0, 1, 0])
    dP = np.array([0, 0, 0, 0, 1])
    
    delta1 = 0.01  #3 Tamaño de paso
    delta2 = 0.01
    delta3 = 0.01
    delta4 = 0.01
    delta5 = 0.01
    
    foo = np.array([0, 0, 0, 0, 0])
    pIn1 = np.array([0, 0, 0, 0, 0])
    
    #k = 1
    redireccion = 0
    i = 0
    
    #********************************************************************************************
    # Inicia el algoritmo
    #********************************************************************************************
    
    while True:
        
        i += 1
        print ("\x1b[1;36;m" + "\ni = ", i)
        
        for j in range (nVec):
            
            if (j == 0):              # Para el vector d1
                
                foo = delta1 * da1
                pIn1 = pIn + foo
                E1 = errorC(m, matriz01[:, 1], matriz01[:, 3], pIn1)
                print("\x1b[1;35;m" + "j=1")
                print ("f(yi + Djdj)=", E1)
                
                if (E1<E):
                    
                    print("\x1b[1;32;m" + "Exito")
                    delta1 = delta1 * 2
                    print ("\x1b[1;35;m" + "D1 =", delta1)
                    pIn = pIn1
                    E = E1
                    j1 = 1
                    
                else:
                    
                    print("\x1b[1;31;m" + "Fracaso")
                    delta1 = delta1 * (-0.5)
                    print ("\x1b[1;35;m" + "D1 =", delta1)
                    j1 = 0
                    
            elif (j == 1):                     # Para el vector d2
                
                foo = delta2 * da2
                pIn1 = pIn + foo
                E1 = errorC(m, matriz01[:, 1], matriz01[:, 3], pIn1)
                print("\nj=2")
                print ("f(yi + Djdj)=", E1)
                
                if (E1<E):
                    
                    print("\x1b[1;32;m" + "Exito")
                    delta2 = delta2 * 2
                    print ("\x1b[1;35;m" + "D2 =", delta2)
                    pIn = pIn1
                    E = E1
                    j2 = 1
                    
                else:
                    
                    print("\x1b[1;31;m" + "Fracaso")
                    delta2 = delta2 * (-0.5)
                    print ("\x1b[1;35;m" + "D2 =", delta2, "\n")
                    j2 = 0
                    
            elif (j == 2):                     # Para el vector d3
                
                foo = delta3 * db1
                pIn1 = pIn + foo
                E1 = errorC(m, matriz01[:, 1], matriz01[:, 3], pIn1)
                print("\nj=3")
                print ("f(yi + Djdj)=", E1)
                
                if (E1<E):
                    
                    print("\x1b[1;32;m" + "Exito")
                    delta3 = delta3 * 2
                    print ("\x1b[1;35;m" + "D3 =", delta3)
                    pIn = pIn1
                    E = E1
                    j3 = 1
                    
                else:
                    
                    print("\x1b[1;31;m" + "Fracaso")
                    delta3 = delta3 * (-0.5)
                    print ("\x1b[1;35;m" + "D3 =", delta3, "\n")
                    j3 = 0
            
            elif (j == 3):                     # Para el vector d4
                
                foo = delta4 * db2
                pIn1 = pIn + foo
                E1 = errorC(m, matriz01[:, 1], matriz01[:, 3], pIn1)
                print("\nj=4")
                print ("f(yi + Djdj)=", E1)
                
                if (E1<E):
                    
                    print("\x1b[1;32;m" + "Exito")
                    delta4 = delta4 * 2
                    print ("\x1b[1;35;m" + "D4 =", delta4)
                    pIn = pIn1
                    E = E1
                    j4 = 1
                    
                else:
                    
                    print("\x1b[1;31;m" + "Fracaso")
                    delta4 = delta4 * (-0.5)
                    print ("\x1b[1;35;m" + "D4 =", delta4, "\n")
                    j4 = 0
            
            else:
                
                foo = delta5 * dP              # Para el vector d5
                pIn1 = pIn + foo
                E1 = errorC(m, matriz01[:, 1], matriz01[:, 3], pIn1)
                print("\nj=5")
                print ("f(yi + Djdj)=", E1)
                
                if (E1<E):
                    
                    print("\x1b[1;32;m" + "Exito")
                    delta5 = delta5 * 2
                    print ("\x1b[1;35;m" + "D5 =", delta5)
                    pIn = pIn1
                    E = E1
                    j5 = 1
                    
                else:                          
                    
                    print("\x1b[1;31;m" + "Fracaso")
                    delta5 = delta5 * (-0.5)
                    print ("\x1b[1;35;m" + "D5 =", delta5, "\n")
                    j5 = 0    
                    
        if (i == 100):
            break
        if (E < 0.055):               # condición de paro
            
            print ("\x1b[1;36;m" + "<==========================================================================>")
            print("Termina la búsqueda\n")
            print("f(x1,x2) = ", E)
            print("En la coordenada: ", pIn)
            print("Número final de redirecciones: ", redireccion)
            print ("<==========================================================================>")
            break
        
    #******************************************************************************************
    # Comienza una redirección
    #******************************************************************************************
    
        if (j1 == 0) and (j2 == 0) and (j3 == 0) and (j4 == 0) and (j5 == 0):   # Comienza una nueva redirección.
            
            redireccion += 1 
            print ("<============================================================>")
            delta1 = 0.01
            delta2 = 0.01
            delta3 = 0.01
            delta4 = 0.01
            delta5 = 0.01
            print("Comienza una redirección")
            print ("redirección No: " + "\x1b[1;36;m", redireccion )
            print("\x1b[1;35;m" + "Punto actual: ", pIn)
            recorrido = pIn - pIn0
            matriz = np.array([[da1[0], da2[0], db1[0], db2[0], dP[0]], \
                                [da1[1], da2[1], db1[1], db2[1], dP[1]], \
                                [da1[2], da2[2], db1[2], db2[2], dP[2]], \
                                [da1[3], da2[3], db1[3], db2[3], dP[3]], \
                                [da1[4], da2[4], db1[4], db2[4], dP[4]]])
            lambdai = np.linalg.solve(matriz, recorrido)
            print ("Lambda: ", lambdai)
            di = np.array([da1, da2, db1, db2, dP])
            
            for j in range (5):
                
                a = np.array([0, 0, 0, 0, 0])
                
                if (lambdai[j] == 0):           # para lambda1[j] = 0
                    a = di[0, :]
                else:                           # Para lambdaj != 0
                    
                    for p in range (j, nVec):   # Se obtiene aj
                        #print("\nj: ", j)
                        #print ("p: ", p)
                        foo = lambdai[p] * di[p] 
                        a = a + foo
                        #print ('a:', a, "\n<=============================================>")
                        
                if (j == 0):                    # Vector da1
                    
                    bj1 = a
                    norma = np.linalg.norm(bj1)
                    da1 = bj1 / norma
                    
                elif (j == 1):                  # Vector da2
                    
                    bj2 = a - ((np.dot(a, da1)) * da1)
                    norma = np.linalg.norm(bj2)
                    da2 = bj2 / norma

                elif (j == 2):                  # Vector db1
                    
                    bj3 = a - (((np.dot(a, da1)) * da1) + ((np.dot(a, da2)) * da2))
                    norma = np.linalg.norm(bj3)
                    db1 = bj3 / norma
                    
                elif (j == 3):                  # Vector db2
                    
                    bj4 = a - (((np.dot(a, da1)) * da1) + ((np.dot(a, da2)) * da2) + \
                               ((np.dot(a, db1)) * db1))
                    norma = np.linalg.norm(bj4)
                    db2 = bj4 / norma
                    
                else:                           # Vector dP
                    
                    bj5 = a - (((np.dot(a, da1)) * da1) + ((np.dot(a, da2)) * da2) + \
                               ((np.dot(a, db1)) * db1) + ((np.dot(a, db2)) * db2))
                    norma = np.linalg.norm(bj5)
                    dP = bj4 / norma
                    print ("\x1b[1;35;m" + "Nuevos vectores: \nd1 = ", da1, "\nd2 = ", da2)
                    print ("d3 = ", db1, "\nd4 = ", db2, "\nd5 = ", dP)
            
        print ("\nPunto actual: ", pIn)        
        print ("\n<===================================================================>")    
            
    #******************************************************************************************
    # Ajuste de la función con los parametros optimizados y errores cuadráticos mínimos
    #******************************************************************************************
    
    a1 = pIn[0]; a2 = pIn[1]; b1 = pIn[2]; b2 = pIn[3]; P = pIn[4]
    gEstimados = biseccion(matriz01[:, 3], a1, a2, b1, b2, P, m) # Se invoca a la función "bisección" para datos estimados.
    matriz01[:, 8] = gEstimados[:, 0]
    gEstimados1 = biseccion(matriz02[:, 1], a1, a2, b1, b2, P, n) # Se invoca a la función bisección para los datos extrapolados
    matriz02[:12, 2] = gEstimados1[:, 0]
    
    foo0 = 0.0
    foo1 = 0.0
    EE0 = 0.0
    EE1 = 0.0
    
    for j in range (len(matriz01)): # Error cuadrático mínimo
        
        foo0 = (matriz01[j, 1] - matriz01[j, 7]) ** 2
        foo1 = (matriz01[j, 1] - matriz01[j, 8]) ** 2
        EE0 = EE0 + foo0
        EE1 = EE1 + foo1
    
    EE0 = (EE0 / (m - 5)) ** 0.5 # Tenía error --> EE1 = (EE1 / (58 - 5)) ** 0.5
    EE1 = (EE1 / (m - 5)) ** 0.5 # Tenía error --> EE1 = (EE1 / (58 - 5)) ** 0.5
    
    print ("Función Doble Gumbel, error cuadrático mínimo (Sin Optimizar): ", EE0)
    print ("Función Doble Gumbel, error cuadrático mínimo (Optimizado): ", EE1)
    EEstandart[0, 0] = EE0
    EEstandart[0, 1] = EE1
    #******************************************************************************************
    # Graficos
    #******************************************************************************************

    titulo0 = "Distribucion Doble Gumbel (sin optimizar)\n EE= " + str(EE0)
    titulo1 = "Distribucion Doble Gumbel (Algoritmo de Rosenbrock)\n EE= " + str(EE1)

    tR = matriz01 [:, 2]
    dReg = matriz01 [:, 1]
    dAjust = matriz01 [:, 7]
    dAjust1 = matriz01 [:, 8]
    dTrExtrap = matriz02 [:12, 0]
    dExtrap = matriz02 [:12, 2]
    
    
    pl.subplot(2, 1, 1)
    pl.subplots_adjust(hspace=0.4)
    pl.scatter(tR, dReg, label='Datos Registrados', marker='.')
    pl.plot(tR, dAjust, color="r", linewidth="1.0", linestyle="-", label ="Datos Ajustados")
    pl.legend(loc="best")
    pl.title(titulo0)
    pl.ylabel("Gastos(m³/s)")
    #pl.xlabel("Periodo de retorno") 
    pl.semilogx(True)
    pl.grid(True, which='both', color='g', linestyle='-', linewidth=0.5)
    
    pl.subplot(2, 1, 2)
    pl.scatter(tR, dReg, label='Datos Registrados', marker='.')
    pl.plot(tR, dAjust1, color="r", linewidth="1.0", linestyle="-", label ="Datos Ajustados")
    pl.legend(loc="best")
    pl.title(titulo1)
    pl.ylabel("Gastos(m³/s)")
    pl.xlabel("Periodo de retorno") 
    pl.semilogx(True)
    pl.grid(True, which='both', color='g', linestyle='-', linewidth=0.5)
    pl.savefig("salidas/DistribucionDGumbelRosen0.png", dpi=1200)
    pl.show()
    
    pl.subplot(2, 1, 1)
    pl.subplots_adjust(hspace=0.4)
    pl.scatter(tR, dReg, label='Datos Registrados', marker='.')
    pl.plot(tR, dAjust1, color="r", linewidth="1.0", linestyle="-", label ="Datos Ajustados")
    pl.legend(loc="best")
    pl.title("Ajuste Doble Gumbel (Algoritmo de Rosenbrock)")
    pl.ylabel("Gastos(m³/s)")
    #pl.xlabel("Periodo de retorno") 
    pl.semilogx(True)
    pl.grid(True, which='both', color='g', linestyle='-', linewidth=0.5)
    
    pl.subplot(2, 1, 2)
    pl.scatter(tR, dReg, label='Datos Registrados', marker='.')
    #pl.scatter(dTrExtrap, dExtrap, label='Datos Extrapolados', marker='.', color='r')
    pl.plot(dTrExtrap, dExtrap, color="r", linewidth="1.0", linestyle="-", label ="Datos Extrapolados")
    pl.legend(loc="best")
    pl.title("Gastos Extrapolados")
    pl.ylabel("Gastos (m³/s)")
    pl.xlabel("Periodo de Retorno en años (Tr)")
    pl.semilogx(True)
    pl.grid(True, which='both', color='g', linestyle='-', linewidth=0.5)
    
    pl.savefig("salidas/DistribucionDGumbelRosen1.png", dpi=1200)
    pl.show()
    
    #******************************************************************************************
    # Salidas
    # Se crea el Data Frame que servirá de salida a una sheet de excel
    #******************************************************************************************
    
    columnas = ['No Orden', 'Gastos Registrados', 'Tr (Annios)', 'F(x)', 'F(qEmpirico)=k/n+1',\
                'F(qMedido)', 'Función de errores cuadrados pesados', 'Gastos estimados',\
                'Gastos estimados optimizados (Rosenbrock)']
    cD = pd.DataFrame(matriz01, columns = columnas)
    
    cD.insert(9, 'Tr', matriz02 [:, 0])
    cD.insert(10, 'F(X)', matriz02 [:, 1])
    cD.insert(11, 'Valor Extrapolado (Rosenbrock)', matriz02 [:, 2])
    cD.insert(12, 'Error Estandart "sin optimizar"', EEstandart[:, 0])
    cD.insert(13, 'Error Estandart "algoritmo de Rosenbrock"', EEstandart[:, 1])

    #*********************** * ********************************************************************
    # se manda el data frame al principal
    #********************************************************************************************

    return(cD)

#******************************************************************************************
# El puto Fin!!  ヽ(。_°)ノ 
#******************************************************************************************
