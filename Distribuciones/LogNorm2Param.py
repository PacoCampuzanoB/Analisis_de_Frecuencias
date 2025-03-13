# -*- coding: utf-8 -*-
"""

Created on Fri May 03 23:07:02 2024

@author: searCampuzano
searcampuzanob@gmail.com

******************************************************************************************
Estimación de los parámetros por momentos y máxima
verosimilitud (es el mismo)
******************************************************************************************

*********AJUSTE DISTRIBUCION LOGNORMAL DE 2 PARAMETROS*******************************

"""
import pandas as pd
import numpy as np
import scipy.stats as st
import pylab as pl

def DLN2P(gastos):
    if (type(gastos) == str):
        print ("\aPara poder realizar el análisis de datos seleccione un archivo valido")
    else:
        tamanio = gastos.shape
        if (len(tamanio) == 1): # Si el array viene de una dimensión, este se
            gastos.resize(len(gastos), 1) # redimensiona
        gastos01 = gastos.copy()
        m = gastos01.size
        rv1 = st.norm
        print(f"\n{m} datos")
        #******************************************************************************************
        #Se crean los arrays de salida.
        #******************************************************************************************
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
        EEstandart = np.full((m, 1), -999, float)

        #********************************************************************************************
        # Se llena la matriz01
        #********************************************************************************************
        
        for i in range (m):
            matriz01 [i, 0] = i + 1             # Se crea la columna para el No de Orden
            matriz01 [i, 1] = gastos01.max()    # Se crea la columna que contiene los gastos registrados
                                                # ordenado en forma descendente
            j = gastos01.argmax()               # Returns the indices of the maximum values along
            gastos01 [j, 0] = -1                # an axis

            matriz01 [i, 2] = (m + 1)/matriz01[i,0]     # Se crea la columna con los periodos de retorno
                                                        # (Tr) a partir de la columna de los gastos
                                                        # registrados ordenados
            matriz01 [i, 3] = 1-(1/matriz01[i, 2])      # Se crea la columna para los valores de F(x)
                                                        # a partir de Tr
            matriz01 [i, 4] = -rv1.isf(matriz01[i, 3])  #Se hace el ajuste (z)

        #************************************************************************************
        #Se Calculan los parametros "Seria bueno poner este codigo en una funcion"
        #************************************************************************************
        EE = 0
        mu = 0
        sigma = 0
        for i in range (m):
            prueba = np.log(matriz01[i, 1])
            mu = mu + prueba
        mu = mu/m                                       # Calculo para obtener mu
        for i in range (m):
            prueba = (np.log(matriz01[i, 1])-mu)**2
            sigma = sigma + prueba
        sigma = np.sqrt(sigma/m)                        # Calculo para obtener sigma

        #************************************************************************************
        #Se Calcula la columna de los gastos ajustados
        #************************************************************************************

        for i in range (m):
            matriz01 [i, 5] = np.exp((matriz01[i, 4] * sigma) + mu)
            prueba = (matriz01 [i, 1] - matriz01 [i, 5]) ** 2  # Error Estandart
            EE = EE + prueba
        EE = (EE/(m-2)) ** 0.5                                 # Se estima el Error Estandart para 2 parametros
        print ("Error Estandart (Momentos y Máxima verosimilitud), Distribución LogNormal 2 parámetros: ", EE)
        EEstandart[0, 0] = EE

        #********************************************************************************************
        # Se llena la matriz02
        #********************************************************************************************

        n=12
        for j in range (n):
            matriz02 [j, 1] = 1.0 - (1.0/matriz02 [j, 0])            # se crea la columna de F(x)
            matriz02 [j, 2] = -rv1.isf(matriz02[j, 1])               # Se realiza el ajuste (z)
            matriz02 [j, 3] = np.exp((matriz02[j, 2] * sigma) + mu)  # Valores Extrapolados

        #********************************************************************************************
        # Se crea el DataFrame final
        #********************************************************************************************

        columnas = ['No Orden', 'Valor Registrado','Tr (Anios)', 'F(x)', 'z', 'Valor Ajustado']
        cD = pd.DataFrame(matriz01, columns = columnas)
        cD.insert(6, 'Tr', matriz02 [:, 0])
        cD.insert(7, 'F(X)', matriz02 [:, 1])
        cD.insert(8, 'Z', matriz02 [:, 2])
        cD.insert(9, 'Valor Extrapolado', matriz02 [:, 3])
        cD.insert(10, 'Error Estandart "Momentos y Máxima verosimilitud"', EEstandart[:, 0])

        #*********************** * ********************************************************************
        # Graficos
        #********************************************************************************************

        titulo = "Distribucion LogNormal 2 Parámetros\n EE= " + str(EE)
        tR = matriz01 [:, 2]
        dReg = matriz01 [:, 1]
        dAjust = matriz01 [:, 5]
        dExtrap = matriz02 [:12, 3]
        dTrExtrap = matriz02 [:12, 0]

        pl.subplot(2,1,1)
        pl.subplots_adjust(hspace=0.3)
        pl.scatter(tR, dReg, label='Datos Registrados', marker='.')
        pl.plot(tR, dAjust, color="r", linewidth="1.0", linestyle="-", label ="Datos Ajustados")
        pl.legend(loc="best")
        pl.title(titulo)
        pl.ylabel("Gastos (m^3/s)")
        pl.semilogx(True)
        pl.grid(True, which='both', color='g', linestyle='-', linewidth=0.5)

        pl.subplot(2,1,2)
        pl.scatter(tR, dReg, label='Datos Registrados', marker='.')
        pl.plot(dTrExtrap, dExtrap, color="r", linewidth="1.0", linestyle="-", label ="Datos Extrapolados")
        pl.legend(loc="best")
        pl.ylabel("Gastos (m^3/s)")
        pl.xlabel("(Tr) Periodos de Retorno")
        pl.semilogx(True)
        pl.grid(True, which='both', color='g', linestyle='-', linewidth=0.5)

        pl.savefig("salidas/DistribucionLogNormal2P.png", dpi=1200)
        pl.show()

        #*********************** * ********************************************************************
        # se manda el data frame al principal
        #********************************************************************************************

        return(cD)

        #*********************** * ********************************************************************
        # Fin del programa
        #********************************************************************************************
