# -*- coding: utf-8 -*-
"""
Created on Fri May 03 23:07:02 2024

@author: Sear Campuzano
em@il: searcampuzanob@gmail.com

Falta MV

"""

import math
import numpy as np
import pandas as pd
import pylab as pl

def DLPIII(gastos, xmedia, sEst, Sp):
    
    """ *********AJUSTE DISTRIBUCION Log Pearson III******************************* """
    
    if (type(gastos) == str):
        print ("\aPara poder realizar el análisis de datos seleccione un archivo valido")
    else:
        tamanio = gastos.shape
        if (len(tamanio) == 1):            # Si el array viene de una dimensión, este se
            gastos.resize(len(gastos), 1)  # redimensiona
        gastos01 = gastos.copy()
        m = gastos.size
        
        # ******************************************************************************************
        #  Se crean los arrays de salida.
        # ******************************************************************************************
        
        V_Ut = np.zeros ((m, 4))       # Contendra Ut, V. para los valores registrados y para
        V_Ut [:] = -999                # los valores extrapolados
        matriz01 = np.zeros ((m, 6))   
        matriz02 = np.zeros ((m, 4))
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
        b0, b1, b2, b3, b4, b5 = 2.2515517, 0.802853, 0.010328, 1.432788, 0.189269, 0.001308
        EEstandart = np.full((m, 1), -999, float)

        #********************************************************************************************
        #  Se llena la matriz01
        #********************************************************************************************

        for i in range (m):
            matriz01 [i, 0] = i + 1             # Se crea la columna para el No de Orden
            matriz01 [i, 1] = gastos01.max()    # Se crea la columna que contiene los gastos registrados
                                                # ordenado en forma descendente
            j = gastos01.argmax()               # Returns the indices of the maximum values along
            gastos01 [j, 0] = -1                # an axis
        
            matriz01 [i, 2] = (m + 1)/matriz01[i,0]  # Se crea la columna con los periodos de retorno
                                                    # (Tr) a partir de la columna de los gastos
                                                    # registrados ordenados
            matriz01 [i, 3] = 1-(1/matriz01[i, 2])   # Se crea la columna para los valores de F(x)
                                                    # a partir de Tr
        gsesg = 0                                    
        for i in range (m):   
            prueba = (matriz01[i, 1] - xmedia) ** 3
            gsesg = gsesg + prueba
        
        gsesg = gsesg / (m * (Sp ** 3))       # Coeficiente de asimetria sesgado
        ginsesg = (gsesg * (m ** 2)) / ((m - 1) * (m - 2))  # Coeficiente de asimetria insesgado
        
        #************************************************************************************
        #  Se estiman los parametros, método de momentos
        #************************************************************************************
        #  Método directo
        #************************************************************************************
        mu1 = 0
        mu2 = 0
        mu3 = 0
        for i in range (len(gastos)):
            mu1 = mu1 + gastos[i, 0]
            mu2 = mu2 + pow(gastos[i, 0], 2)
            mu3 = mu3 + pow(gastos[i, 0], 3)
            
        mu1 = mu1 / len(gastos)
        mu2 = mu2 / len(gastos)
        mu3 = mu3 / len(gastos)
        
        B = (math.log(mu3)-(3*math.log(mu1))) / (math.log(mu2)-(2*math.log(mu1)))
        print ('B:', B)
        C = pow(B-3, -1)
        if ((B>3.5) & (B<=6)):
            print ("Método directo")
            A = -0.23019 + (1.65262*C) + (0.20911*C**2) - (0.04557*C**3)
            alfa = pow(A+3, -1)
            beta = (math.log(mu2) - (2*math.log(mu1))) / ((math.log((1-alfa)**2)) - math.log(1-(2*alfa)))
            y0 = math.log(mu1) + (beta*math.log(1-alfa))
            
        elif ((B>3) & (B<=3.5)):
            print ("Método directo")
            A = -0.45157 + (1.99955*C)
            alfa = pow(A+3, -1)
            beta = (math.log(mu2) - (2*math.log(mu1))) / ((math.log((1-alfa)**2)) - math.log(1-(2*alfa)))
            y0 = math.log(mu1) + (beta*math.log(1-alfa))
            
        elif ((B<3) | (B>6)):  # Método indirecto
            print ("Método indirecto")
            y = 0
            for i in range (len(gastos)):       
                y = y + math.log(gastos[i, 0])
                
            y = y / len(gastos)                # y = 1/n * suma (lnxi)
            print('y:', y)
            
            S = 0
            for i in range (len(gastos)):
                S = S + pow(math.log(gastos[i, 0])-y, 2) 
            
            S = pow(S/(len(gastos)-1), 0.5)    # S = raiz((suma(xi-y)^2)/(n-1))
            print('S:', S)
            
            g = 0
            for i in range (len(gastos)):
                g = g + pow((math.log(gastos[i, 0]))-y, 3)
            m = len(gastos)
            g = m*g / ((m-1)*(m-2)*pow(S, 3))  # g = (n*suma(ln(xi-y)^3)) / ((n-1)*(n*2)*S^3)
            print('g:', g)
            
            beta = (2 / g)**2
            alfa = S / beta**0.5
            y0 = y - (S*beta**0.5)
            
            print('\nalfa:', alfa)
            print('beta:', beta)
            print('y0:', y0)
        
        EE_Mom = 0
        
        for i in range (m):
            if (matriz01[i, 3] > 0.5):
                V_Ut[i, 0] = (math.log(1 / pow(1-matriz01[i, 3], 2))) ** 0.5
                V_Ut[i, 1] = V_Ut[i, 0] - ((b0 + (b1 * V_Ut[i, 0]) + (b2 * pow(V_Ut[i, 0], 2))) /\
                (1 + (b3 * V_Ut[i, 0]) + (b4 * pow(V_Ut[i, 0], 2)) + (b5 * pow(V_Ut[i, 0], 3))))
        
            if (matriz01[i, 3] <= 0.5):
                V_Ut[i, 0] = (math.log(1 / pow(matriz01[i, 3], 2))) ** 0.5
                V_Ut[i, 1] = V_Ut[i, 0] - ((b0 + (b1 * V_Ut[i, 0]) + (b2 * pow(V_Ut[i, 0], 2))) /\
                (1 + (b3 * V_Ut[i, 0]) + (b4 * pow(V_Ut[i, 0], 2)) + (b5 * pow(V_Ut[i, 0], 3))))
                V_Ut[i, 1] = -V_Ut[i, 1]       # Cambio de signo

            matriz01[i, 4] = math.exp(y0 + ((alfa * beta) * ((1 - pow(9 * beta, -1) + (V_Ut[i,1] * \
            (pow(9 * beta, -1)) ** 0.5)) ** 3)))                # Se calcula x
            prueba = (matriz01 [i, 1] - matriz01 [i, 4]) ** 2  # Error Estandart momentos
            EE_Mom = EE_Mom + prueba
        
        EE_Mom = (EE_Mom/(m-3)) ** 0.5  # Se estima el Error Estandart para 3 parametros momentos
        print ("\nDistribución LogPearson III, error Estandart (Momentos): ", EE_Mom)
        EEstandart[0, 0] = EE_Mom

        #************************************************************************************
        #  Se estiman los parametros, método de máxima verosimilitud
        #************************************************************************************
        
        
        
        
        
        
        
        
        
        
        
        
        #********************************************************************************************
        # Se llena la matriz02. Datos extrapolados
        #********************************************************************************************
        
        n = 12
        
        for i in range (n):
            matriz02 [i, 1] = 1 - pow(matriz02[i,0], -1)  # F(x) valores extrapolados
            if (matriz02[i, 1] > 0.5):
                V_Ut[i, 2] = (math.log(1 / pow(1 - matriz02[i, 1], 2))) ** 0.5
                V_Ut[i, 3] = V_Ut[i, 2] - ((b0 + (b1 * V_Ut[i, 2]) + (b2 * pow(V_Ut[i, 2], 2))) /\
                (1 + (b3 * V_Ut[i, 2]) + (b4 * pow(V_Ut[i, 2], 2)) + (b5 * pow(V_Ut[i, 2], 3))))
        
            if (matriz02[i, 1] <= 0.5):
                V_Ut[i, 2] = (math.log(1 / pow(matriz02[i, 1], 2))) ** 0.5
                V_Ut[i, 3] = V_Ut[i, 2] - ((b0 + (b1 * V_Ut[i, 2]) + (b2 * pow(V_Ut[i, 2], 2))) /\
                (1 + (b3 * V_Ut[i, 2]) + (b4 * pow(V_Ut[i, 2], 2)) + (b5 * pow(V_Ut[i, 2], 3))))
                V_Ut[i, 3] = -V_Ut[i, 3]
                
            matriz02[i, 2] = math.exp(y0 + (alfa * beta) * ((1 - pow(9 * beta, -1) + (V_Ut[i,3] * \
            (pow(9 * beta, -1)) ** 0.5)) ** 3))         # Se calcula x Extrapolados
        
        #********************************************************************************************
        # Se crea el DataFrame final
        #********************************************************************************************
        
        columnas = ['No Orden', 'Valor Registrado','Tr (Anios)', 'F(x)', 'Valor Ajustado (momentos)',\
                    'Valor Ajustado (max ver)']
        cD = pd.DataFrame(matriz01, columns = columnas)
        
        cD.insert(6, 'Tr', matriz02 [:, 0])
        cD.insert(7, 'F(X)', matriz02 [:, 1])
        cD.insert(8, 'Valor Extrapolado (momentos)', matriz02 [:, 2])
        cD.insert(9, 'Valor Extrapolado (max ver)', matriz02 [:, 3])
        cD.insert(10, 'Error Estandart "Momentos"', EEstandart[:, 0])

        #*********************** * ********************************************************************
        # Graficos
        #********************************************************************************************

        titulo0 = "Distribucion LogPearson III (Momentos)\n EE= " + str(EE_Mom)

        tR = matriz01 [:, 2]
        dReg = matriz01 [:, 1]
        dAjustMom = matriz01 [:, 4]
        # AjustMax = matriz01 [:, 5]
        dExtrapMom = matriz02 [:12, 2]
        # dExtrapMax = matriz02 [:12, 3]
        dTrExtrap = matriz02 [:12, 0]
        
        #*********************** * ********************************************************************
        
        pl.subplot(2,1,1)
        pl.subplots_adjust(hspace=0.3)
        pl.scatter(tR, dReg, label='Datos Registrados', marker='.')
        pl.plot(tR, dAjustMom, color="r", linewidth="1.0", linestyle="-", label ="Datos Ajustados")
        pl.legend(loc="best")
        pl.title(titulo0)
        pl.ylabel("Gastos (m^3/s)")
        pl.semilogx(True)
        pl.grid(True, which='both', color='g', linestyle='-', linewidth=0.5)
        
        pl.subplot(2,1,2)
        pl.scatter(tR, dReg, label='Datos Registrados', marker='.')
        pl.plot(dTrExtrap, dExtrapMom, color="r", linewidth="1.0", linestyle="-", label ="Datos Extrapolados")
        pl.legend(loc="best")
        pl.ylabel("Gastos (m^3/s)")
        pl.xlabel("(Tr) Periodos de Retorno")
        pl.semilogx(True)
        pl.grid(True, which='both', color='g', linestyle='-', linewidth=0.5)
        
        pl.savefig("salidas/LogPearsonIII_Momentos.png", dpi=1200)
        pl.show()

        #*********************** * ********************************************************************
        # se manda el data frame al principal
        #********************************************************************************************

        return(cD)

        #*********************** * ********************************************************************
        # Fucking end!!
        #********************************************************************************************
