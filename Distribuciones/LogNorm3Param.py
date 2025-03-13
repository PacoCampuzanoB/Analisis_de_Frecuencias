# -*- coding: utf-8 -*-
"""
Created on Fri May 03 23:07:02 2024

@author: searCampuzano
searcampuzanob@gmail.com

******************************************************************************************
Estimación de los parámetros por momentos y máxima
verosimilitud.
Para encontrar los parámetros por máxima verosimilitud se resolverá la función a 
través del metodo de bisección, modulo paramMaxVerLN3P
******************************************************************************************

*********AJUSTE DISTRIBUCION LOGNORMAL DE 3 PARAMETROS*******************************
"""

import pandas as pd
import numpy as np
import scipy.stats as st
import pylab as pl
import Distribuciones.paramMaxVerLN3P

def DLN3P(gastos, mediaDr, desvEst, S):
    if (type(gastos) == str):
        print ("\aPara poder realizar el análisis de datos seleccione un archivo valido")
    else:
        tamanio = gastos.shape
        if (len(tamanio) == 1):            # Si el array viene de una dimensión, este se
            gastos.resize(len(gastos), 1)  # redimensiona
        gastos01 = gastos.copy()
        m = gastos.size
        rv1 = st.norm

        #******************************************************************************************
        #Se crean los arrays de salida.
        #******************************************************************************************

        matriz01 = np.zeros ((m, 7))
        matriz02 = np.zeros ((m, 5))
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
        EEstandart = np.full((m, 2), -999, float)

        #********************************************************************************************
        # Se llena la matriz01
        #********************************************************************************************
        
        for i in range (m):
            matriz01 [i, 0] = i + 1                      # Se crea la columna para el No de Orden
            matriz01 [i, 1] = gastos01.max()             # Se crea la columna que contiene los gastos registrados
            j = gastos01.argmax()                        # ordenado en forma descendente
            gastos01 [j, 0] = -1
            matriz01 [i, 2] = (m + 1)/matriz01[i,0]      # Se crea la columna con los periodos de retorno
                                                        # (Tr) a partir de la columna de los gastos registrados ordenados
            matriz01 [i, 3] = 1-(1/matriz01[i, 2])       # Se crea la columna para los valores de F(x) a partir de Tr
            matriz01 [i, 4] = -rv1.isf(matriz01[i, 3])   # Se hace el ajuste (z)

        gsesg = 0
        for i in range (m):                              # Se calcula la desviacion estandart
            
            prueba = (matriz01[i, 1] - mediaDr) ** 3
            gsesg = gsesg + prueba

        gsesg = gsesg / (m * (S**3))                     # Coeficiente de asimetria sesgado
        ginsesg = (gsesg * (m**2)) / ((m-1) * (m-2))     # Coeficiente de asimetria insesgado
        w = ((((ginsesg**2)+4)**0.5) -ginsesg) / 2       # Omega
        etaX = desvEst / mediaDr                         # Eta X
        etaZ = (1 - (w**(2.0/3))) / w**(1.0/3)           # Eta Z
        print(f'gsesg = {gsesg}')
        print(f'ginsesg = {ginsesg}')
        print(f'omega = {w}')
        print(f'etaX = {etaX}')
        print(f'etaZ = {etaZ}')

        #************************************************************************************
        #Se Calculan los parámetros por el método de momentos
        #************************************************************************************

        x0_Mom = mediaDr * (1 - (etaX/etaZ))
        muY_Mom = np.log(desvEst / etaZ) - (0.5 * np.log(etaZ**2 + 1))
        sigmaY_Mom = (np.log(etaZ**2 + 1)) ** 0.5
        print(f'x0_Mom = {x0_Mom}')
        print(f'muY_Mom = {muY_Mom}')
        print(f'sigmaY_mom = {sigmaY_Mom}')

        #************************************************************************************
        #Se calculan los parámetros por el método Máxima Verosimilitud
        #************************************************************************************
        
        parametrosMV = Distribuciones.paramMaxVerLN3P.parametros(gastos)  # Se estiman los parámetros 
        x0_MV = parametrosMV[0]     # Parámetro (X0) por el metodo de máxima verosimilitud
        paramMV = parametrosMV[1]
        muY_MV =paramMV[0]          # Parametro (muy) por el metodo de máxima verosimilitud
        sigmaY_MV = paramMV[1]      # Parámetro (sigmay) por el metodo de máxima verosimilitud
        Fx = paramMV[2]             # Resultado de la ecuación F(x) descrita en la función paramMaxVerLN3P
        noIteraciones = paramMV[3]  # No de iteraciones que se llevó la función para estimar los parametros
        print(f'parametrosMv = {parametrosMV}')
        print(f'x0_Mv = {x0_MV}')
        print(f'paramMV = {paramMV}')
        print(f'muY_MV = {muY_MV}')
        print(f'sigmaY_MV = {sigmaY_MV}')
        print(f'Fx = {Fx}')
        print(f'NoIteraciones = {noIteraciones}')
        #************************************************************************************
        #Se Calcula la columna de los gastos ajustados y Error Estandart
        #************************************************************************************

        EE_Mom = 0
        EE_MV = 0
        for i in range (m):
            matriz01 [i, 5] = x0_Mom + np.exp(muY_Mom + (matriz01[i, 4] * sigmaY_Mom))  # Gastos ajustados momentos
            matriz01 [i, 6] = x0_MV + np.exp(muY_MV + (matriz01[i, 4] * sigmaY_MV))     # Gastos ajustados máxVer
            prueba = (matriz01 [i, 1] - matriz01 [i, 5]) ** 2   # Error Estandart momentos
            EE_Mom = EE_Mom + prueba
            prueba1 = (matriz01 [i, 1] - matriz01 [i, 6]) ** 2  # Error Estandart máxVer
            EE_MV = EE_MV + prueba1
            
        EE_Mom = (EE_Mom/(m-3)) ** 0.5  # Se estima el Error Estandart para 3 parametros
        EE_MV = (EE_MV/(m-3)) ** 0.5    # Se estima el Error Estandart para 3 parametros máxVer
        EEstandart[0, 0] = EE_Mom
        EEstandart[0, 1] = EE_MV

        #********************************************************************************************
        # Se llena la matriz02
        #********************************************************************************************

        n=12
        for j in range (n):
            matriz02 [j, 1] = 1.0 - (1.0/matriz02 [j, 0])  # Se crea la columna de F(x)
            matriz02 [j, 2] = -rv1.isf(matriz02[j, 1])     # Se realiza el ajuste (z)
            matriz02 [j, 3] = x0_Mom + np.exp(muY_Mom + (matriz02[j, 2] * sigmaY_Mom))  # Valores Extrapolados momentos
            matriz02 [j, 4] = x0_MV + np.exp(muY_MV + (matriz02[j, 2] * sigmaY_MV))     # Valores Extrapolados máxVer

        #********************************************************************************************
        # Se crea el DataFrame final
        #********************************************************************************************

        columnas = ['No Orden', 'Valor Registrado','Tr (Anios)', 'F(x)', 'z', 'Valor Ajustado Momentos', 'Valor Ajustado MaxVer']
        cD = pd.DataFrame(matriz01, columns = columnas)

        cD.insert(7, 'Tr', matriz02 [:, 0])
        cD.insert(8, 'F(X)', matriz02 [:, 1])
        cD.insert(9, 'Z', matriz02 [:, 2])
        cD.insert(10, 'Valor Extrapolado', matriz02 [:, 3])
        cD.insert(11, 'Valor Extrapolado MaxVer', matriz02 [:, 4])
        cD.insert(12, 'Error Estandart "Momentos"', EEstandart[:, 0])
        cD.insert(13, 'Error Estandart "Máxima verosimilitud"', EEstandart[:, 1])

        #*********************** * ********************************************************************
        # Graficos
        #********************************************************************************************

        titulo0 = "Distribucion LogNormal 3 Parámetros (Momentos)\n EE= " + str(EE_Mom)
        titulo1 = "Distribucion LogNormal 3 Parámetros (MaxVer)\n EE= " + str(EE_MV)

        tR = matriz01 [:, 2]
        dReg = matriz01 [:, 1]
        dAjustMom = matriz01 [:, 5]
        dAjustMax = matriz01 [:, 6]
        dExtrapMom = matriz02 [:12, 3]
        dExtrapMax = matriz02 [:12, 4]
        dTrExtrap = matriz02 [:12, 0]
                
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
        
        pl.savefig("salidas/LogNormal3PMomentos.png", dpi=1200)
        pl.show()
        
        #*********************** * ********************************************************************
        
        pl.subplot(2,1,1)
        pl.subplots_adjust(hspace=0.3)
        pl.scatter(tR, dReg, label='Datos Registrados', marker='.')
        pl.plot(tR, dAjustMax, color="r", linewidth="1.0", linestyle="-", label ="Datos Ajustados")
        pl.legend(loc="best")
        pl.title(titulo1)
        pl.ylabel("Gastos (m^3/s)")
        pl.semilogx(True)
        pl.grid(True, which='both', color='g', linestyle='-', linewidth=0.5)
        
        pl.subplot(2,1,2)
        pl.scatter(tR, dReg, label='Datos Registrados', marker='.')
        pl.plot(dTrExtrap, dExtrapMax, color="r", linewidth="1.0", linestyle="-", label ="Datos Extrapolados")
        pl.legend(loc="best")
        pl.ylabel("Gastos (m^3/s)")
        pl.xlabel("(Tr) Periodos de Retorno")
        pl.semilogx(True)
        pl.grid(True, which='both', color='g', linestyle='-', linewidth=0.5)
        
        pl.savefig("salidas/LogNormal3PMaxV.png", dpi=1200)
        pl.show()
        
        #*********************** * ********************************************************************
        
        pl.subplot(2,1,1)
        pl.subplots_adjust(hspace=0.3)
        pl.scatter(tR, dReg, label='Datos Registrados', marker='.')
        pl.plot(tR, dAjustMom, color="r", linewidth="1.0", linestyle="-", label ="Datos Ajustados Mom")
        pl.plot(tR, dAjustMax, color="g", linewidth="1.0", linestyle="-", label ="D Ajustados Max Ver")
        pl.legend(loc="best")
        pl.title("Distribucion LogNormal 3 Parametros Metodos Momentos y \nMaxima Verosimilitud")
        pl.ylabel("Gastos (m^3/s)")
        pl.semilogx(True)
        pl.grid(True, which='both', color='g', linestyle='-', linewidth=0.5)
        
        pl.subplot(2,1,2)
        pl.scatter(tR, dReg, label='Datos Registrados', marker='.')
        pl.plot(dTrExtrap, dExtrapMom, color="r", linewidth="1.0", linestyle="-", label ="D Extrap Mom")
        pl.plot(dTrExtrap, dExtrapMax, color="g", linewidth="1.0", linestyle="-", label ="D Extrap Max Ver")
        pl.legend(loc="best")
        pl.ylabel("Gastos (m^3/s)")
        pl.xlabel("(Tr) Periodos de Retorno")
        pl.semilogx(True)
        pl.grid(True, which='both', color='g', linestyle='-', linewidth=0.5)
        
        pl.savefig("salidas/LogNormal3PMomMaxV.png", dpi=1200)
        pl.show()
        
        #*********************** * ********************************************************************
        
        print ("Error estandart método de momentos, ajuste LogNormal 3 parámetros: ", EE_Mom)
        print ("Error estandart método de máxima verosimilitud, ajuste LogNormal 3 parámetros: ", EE_MV)
        print ("\a")

        #*********************** * ********************************************************************
        # se manda el data frame al principal
        #********************************************************************************************

        return(cD)

        #*********************** * ********************************************************************
        # Fin del programa
        #********************************************************************************************
