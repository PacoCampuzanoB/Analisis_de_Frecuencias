# -*- coding: utf-8 -*-
"""
Created on Fri May 03 23:07:02 2024

@author: searCampuzano
searcampuzanob@gmail.com

******************************************************************************************
Estimación de los parámetros por momentos y máxima
verosimilitud.
******************************************************************************************

*********AJUSTE DISTRIBUCION EXPONENCIAL DE 2 PARAMETROS*******************************

"""

import pandas as pd
import numpy as np
import pylab as pl

def DE2P(gastos, mediaDr, desvEst):
    
    if (type(gastos) == str):
        print ("\aPara poder realizar el análisis de datos seleccione un archivo valido")
    else:
        tamanio = gastos.shape
        if (len(tamanio) == 1):            # Si el array viene de una dimensión, este se
            gastos.resize(len(gastos), 1)  # redimensiona
        gastos01 = gastos.copy()
        m = gastos01.size
    
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
        EEstandart = np.full((m, 2), -999, float)

        #********************************************************************************************
        # Se llena la matriz01
        #********************************************************************************************
        
        for i in range (m):
            matriz01 [i, 0] = i + 1                  # Se crea la columna para el No de Orden
            matriz01 [i, 1] = gastos01.max()         # Se crea la columna que contiene los gastos 
            j = gastos01.argmax()                    # registrados ordenado en forma descendente
            gastos01 [j, 0] = -1
            matriz01 [i, 2] = (m + 1)/matriz01[i,0]  # Se crea la columna con los periodos de retorno  
                                                    # (Tr) a partir de la columna de los gastos registrados ordenados
            matriz01 [i, 3] = 1-(1/matriz01[i, 2])   # Se crea la columna para los valores de F(x) a partir de Tr

        #********************************************************************************************
        # Se calcula los parametros alfa y beta para momentos
        #********************************************************************************************

        betaMom = desvEst
        alfaMom = mediaDr - desvEst
        
        #********************************************************************************************
        # Se calcula los parametros alfa y beta para maxima verosimilitud
        #********************************************************************************************
        
        betaMaxV = 0
        for i in range (m):
            
            prueba = matriz01 [i, 1] - matriz01 [m - 1, 1]
            betaMaxV = betaMaxV + prueba
        
        betaMaxV = betaMaxV / (m - 1)
        alfaMaxV = matriz01 [m-1, 1] - (betaMaxV / (m - 1))
        
        #************************************************************************************
        #Se Calculan la columnas de los gastos ajustados por el metodo de momentos y
        #metodo de maxima verosimilitud.
        #************************************************************************************
        
        EEMom = 0
        EEMaxV = 0 
        for i in range (m):
        
            # Momentos
            matriz01 [i, 4] = alfaMom - (betaMom * np.log(1 - matriz01 [i, 3]))
            prueba = (matriz01 [i, 1] - matriz01 [i, 4]) ** 2
            EEMom = EEMom + prueba

            # Maxima Verosimilitud
            matriz01 [i, 5] = alfaMaxV - (betaMaxV * np.log(1 - matriz01 [i, 3]))
            prueba = (matriz01 [i, 1] - matriz01 [i, 5]) ** 2
            EEMaxV = EEMaxV + prueba
            
        # Se estiman los errores estandart para 2 parametros    
        EEMom = (EEMom / (m - 2)) ** 0.5
        EEMaxV = (EEMaxV /(m - 2)) ** 0.5
        EEstandart[0, 0] = EEMom
        EEstandart[0, 1] = EEMaxV

        #********************************************************************************************
        # Se llena la matriz02
        #********************************************************************************************
        
        n = 12
        for i in range (n):
                
            matriz02 [i, 1] = 1.0 - (1.0 / matriz02 [i, 0])  # se crea la columna de F(x)
                                                            # Valores Extrapolados
                                # Momentos
            matriz02 [i, 2] = alfaMom - (betaMom * np.log(1 - matriz02 [i, 1]))
                                # Maxima Verosimilitud
            matriz02 [i, 3] = alfaMaxV - (betaMaxV * np.log(1 - matriz02 [i, 1]))  
            
        #********************************************************************************************
        # Se crea el DataFrame final
        #********************************************************************************************
        
        columnas = ['No Orden', 'Valor Registrado','Tr (Anios)', 'F(x)', 'Valor Ajustado (Momentos)', 'Valor Ajustado (Maxima Verosimilitud)']
        cD = pd.DataFrame(matriz01, columns = columnas)
        cD.insert(6, 'Tr', matriz02 [:, 0])
        cD.insert(7, 'F(X)', matriz02 [:, 1])
        cD.insert(8, 'Valor Extrapolado (Momentos)', matriz02 [:, 2])
        cD.insert(9, 'Valor Extrapolado (Maxima Verosimilitud', matriz02 [:, 3])
        cD.insert(10, 'Error Estandart "Momentos"', EEstandart[:, 0])
        cD.insert(11, 'Error Estandart "Máxima verosimilitud"', EEstandart[:, 1])

        #*********************** * ********************************************************************
        # Graficos
        #********************************************************************************************

        titulo0 = "Distribucion Exponencial 2 Parámetros (Momentos)\n EE= " + str(EEMom)
        titulo1 = "Distribucion Exponencial 2 Parámetros (MaxVer)\n EE= " + str(EEMaxV)

        tR = matriz01 [:, 2]
        dReg = matriz01 [:, 1]
        dAjustMom = matriz01 [:, 4]
        dAjustMax = matriz01 [:, 5]
        dExtrapMom = matriz02 [:12, 2]
        dExtrapMax = matriz02 [:12, 3]
        dTrExtrap = matriz02 [:12, 0]
        
        #********************************************************************************************
        
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
        
        pl.savefig("salidas/Exponencial2PMom.png", dpi=1200)
        pl.show()
        
        #********************************************************************************************
        
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
        
        pl.savefig("salidas/Exponencial2PMaxV.png", dpi=1200)
        pl.show()
        
        #********************************************************************************************
        
        pl.subplot(2,1,1)
        pl.subplots_adjust(hspace=0.3)
        pl.scatter(tR, dReg, label='Datos Registrados', marker='.')
        pl.plot(tR, dAjustMom, color="r", linewidth="1.0", linestyle="-", label ="D Ajustados Mom")
        pl.plot(tR, dAjustMax, color="g", linewidth="1.0", linestyle="-", label ="D Ajustados Max Ver")
        pl.legend(loc="best")
        pl.title("Distribucion Exponencial 2 Parametros\n (Metodos de Momentos y Máxima Verosimilitud)")
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
        
        pl.savefig("salidas/Exponencial2PMomMaxV.png", dpi=1200)
        pl.show()
        
        #********************************************************************************************
        
        print ("Distribucion Exponencial 2 parámetros, error estandart método de momentos: ", EEMom)
        print ("Distribucion Exponencial 2 parámetros, error estandart método de maxima verosimilitud es: ", EEMaxV)

        #*********************** * ********************************************************************
        # se manda el data frame al principal
        #********************************************************************************************

        return(cD)

#*********************** * ********************************************************************
# Fin del programa
#********************************************************************************************
