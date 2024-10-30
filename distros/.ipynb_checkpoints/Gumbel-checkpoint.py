# -*- coding: utf-8 -*-
"""
Created on Fri May 03 23:07:02 2024

@author: searCampuzano
searcampuzanob@gmail.com

*********AJUSTE DISTRIBUCION GUMBEL (MOMENTOS), (MOMENTOS L)*******************************

"""

import pandas as pd
import numpy as np
import pylab as pl
import math

def DG(gastos, mediaDr, desvEst, estacion):
    if (type(gastos) == str):
        print ("\aPara poder realizar el análisis de datos seleccione un archivo valido")
    else:
        tamanio = gastos.shape
        if (len(tamanio) == 1):            # Si el array viene de una dimensión, este se
            gastos.resize(len(gastos), 1)  # redimensiona
        gastos01 = gastos.copy()
        m = gastos01.size                     # tamaño de gastos "renglones"

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
            matriz01 [i, 0] = i + 1           # Se crea la columna con el No de orden columna '1'
            matriz01 [i, 1] = gastos01.max()  # Se crea la columna con los gastos "columna 2"
            j = gastos01.argmax()             # registrados y ordenados
            gastos01 [j, 0] = -1              # en orden ascendente
            matriz01 [i, 2] = (m + 1) / matriz01 [i, 0]  # Se crea la columna con Tr
            matriz01 [i, 3] = 1 - (1 / matriz01 [i, 2])  # Se crea la columna de F(x), "a partir de Tr"
        
        #********************************************************************************************
        # Se calcula los parametros alfa y beta para  el método de momentos
        #********************************************************************************************
        
        beta = mediaDr - (0.45 * desvEst)
        alfa = desvEst * ((6 ** 0.5) / math.pi)
        
        #********************************************************************************************
        # Se calcula los parametros alfa y beta para momentos-L
        #********************************************************************************************
        
        lambda1 = mediaDr      # alfa = media = M0
        suma = 0.0
        for i in range (m):                   # M1 = (1/n(n-1))*suma(xi*n-i)
            M1 = matriz01[i, 1] * (m-(i+1))
            suma = suma + M1  
        M1 = suma/(m*(m-1))
        lambda2 = (2*M1)-mediaDr
        
        alfaMl = lambda2/math.log(2)          # alfa momentos-L
        betaMl = lambda1 - (0.577216*alfaMl)  # beta momentos-L

        #*****************************************************************************************
        #Se Calculan la columnas de los gastos ajustados por el metodo de momentos y
        #metodo de momentos-L.
        #*****************************************************************************************
                        
                        # Para el método de momentos
                        
        EEMom = 0.0
        for i in range (m):
            matriz01 [i, 4] = beta - (alfa * np.log(-np.log(matriz01 [i, 3])))
            prueba = (matriz01 [i, 1] - matriz01 [i, 4]) ** 2
            EEMom = EEMom + prueba

        EEMom = (EEMom / (m - 2)) ** 0.5  # Se estima el error estandart para dos parametros
        print ('Distribución Gumbel, error estandart (momentos): ', EEMom)
        
                        # Para el método de momentos-L
                        
        EEMomL = 0.0
        for i in range (m):
            matriz01 [i, 5] = betaMl - (alfaMl * np.log(-np.log(matriz01 [i, 3])))
            prueba = (matriz01 [i, 1] - matriz01 [i, 5]) ** 2
            EEMomL = EEMomL + prueba
        
        EEMomL = (EEMomL / (m - 2)) ** 0.5  # Se estima el error estandart para dos parametros
        print ('Distribución Gumbel, error estandart (momentos_L): ', EEMomL)
        EEstandart[0, 0] = EEMom
        EEstandart[0, 1] = EEMomL

        #********************************************************************************************
        # Se llena la matriz02
        #********************************************************************************************
        
        n = 12
        for i in range (n):
            matriz02 [i ,1] = 1.0 - (1.0 / matriz02 [i, 0]) #Columna de F(x)
            
                        # Valores extrapolados
                        # Momentos
                        
            matriz02 [i, 2] = beta - (alfa * np.log(-np.log(matriz02 [i, 1])))
            
                        # Momentos-L
                        
            matriz02 [i, 3] = betaMl - (alfaMl * np.log(-np.log(matriz02 [i, 1])))
    
        #********************************************************************************************
        # Se crea el DataFrame final
        #********************************************************************************************
        
        columnas = ['No Orden', 'Gastos Registrados', 'Tr (Anios)', 'F(x)', 'Valor Ajustado (Momentos)', 'Valor Ajustado (Momentos-L)']
        cD = pd.DataFrame(matriz01, columns = columnas)
        
        cD.insert(6, 'Tr', matriz02 [:, 0])
        cD.insert(7, 'F(X)', matriz02 [:, 1])
        cD.insert(8, 'Valor Extrapolado (Momentos)', matriz02 [:, 2])
        cD.insert(9, 'Valor Extrapolado (Momentos-L)', matriz02[:, 3])
        cD.insert(10, 'Error Estandart "Momentos"', EEstandart[:, 0])
        cD.insert(11, 'Error Estandart "Momentos-L"', EEstandart[:, 1])


        #*********************** * ********************************************************************
        # Salidas
        #********************************************************************************************
        
        # writer = pd.ExcelWriter('DistGumbel.xlsx', engine = 'xlsxwriter')
        # cD.to_excel(writer, 'DistGumbel')
        # writer.save()

        #*********************** * ********************************************************************
        # Graficos
        #********************************************************************************************
        titulo0 = "Distribucion Gumbel (Momentos)\n EE= " + str(EEMom)
        titulo1 = "Distribucion Gumbel (Momentos-L)\n EE= " + str(EEMomL)

        tR = matriz01 [:, 2]
        dReg = matriz01 [:, 1]
        dAjustMom = matriz01 [:, 4]
        dAjustMax = matriz01 [:, 5]
        dExtrapMom = matriz02 [:12, 2]
        dExtrapMax = matriz02 [:12, 3]
        dTrExtrap = matriz02 [:12, 0]
        titulo_grafico_1 = "salidas/" + estacion + "_Gumbel_Mom.png"
        titulo_grafico_2 = "salidas/" + estacion + "_Gumbel_Mom-L.png"
        titulo_grafico_3 = "salidas/" + estacion + "_Gumbel_Mom_Mom-L.png"
        
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
        
        pl.savefig(titulo_grafico_1, dpi=1200)
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
        
        pl.savefig(titulo_grafico_2, dpi=1200)
        pl.show()
        
        #*********************** * ********************************************************************
        
        pl.subplot(2,1,1)
        pl.subplots_adjust(hspace=0.3)
        pl.scatter(tR, dReg, label='Datos Registrados', marker='.')
        pl.plot(tR, dAjustMom, color="r", linewidth="1.0", linestyle="-", label ="Datos Ajustados Mom")
        pl.plot(tR, dAjustMax, color="g", linewidth="1.0", linestyle="-", label ="Datos Ajustados Mom-L")
        pl.legend(loc="best")
        pl.title("Distribucion Gumbel\n Momentos y Momentos-L")
        pl.ylabel("Gastos (m^3/s)")
        pl.semilogx(True)
        pl.grid(True, which='both', color='g', linestyle='-', linewidth=0.5)
        
        pl.subplot(2,1,2)
        pl.scatter(tR, dReg, label='Datos Registrados', marker='.')
        pl.plot(dTrExtrap, dExtrapMom, color="r", linewidth="1.0", linestyle="-", label ="D Extrap Mom")
        pl.plot(dTrExtrap, dExtrapMax, color="g", linewidth="1.0", linestyle="-", label ="D Extrap Mom-L")
        pl.legend(loc="best")
        pl.ylabel("Gastos (m^3/s)")
        pl.xlabel("(Tr) Periodos de Retorno")
        pl.semilogx(True)
        pl.grid(True, which='both', color='g', linestyle='-', linewidth=0.5)
        
        pl.savefig(titulo_grafico_3, dpi=1200)
        pl.show()

        #*********************** * ********************************************************************
        # se manda el data frame al principal
        #********************************************************************************************

        return(cD)

#*********************** * ********************************************************************
# Fin del script! (>‿◠)✌
#********************************************************************************************
