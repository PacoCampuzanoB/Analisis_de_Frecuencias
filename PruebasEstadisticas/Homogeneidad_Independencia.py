# -*- coding: utf-8 -*-
"""
Created on Fri May 03 22:07:02 2024

@author: searCampuzano
searcampuzanob@gmail.com
****************************PRUEBAS DE HOMOGENEIDAD*************************************
************************PRUEBA ESTADISTICA DE HELMERT "Ln 44-88"************************
************************PRUEBA ESTADISTICA T DE STUDENT "Ln 92-143"**********************
************************PRUEBA ESTADISTICA DE CRAMER "Ln "147-198"**********************


****************************PRUEBA INDEPENDENCIA DE EVENTOS ****************************
****************************Ln "202-254"*************************************
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as st

def hI(datos1):
    
    gastos1 = datos1

    if (type(gastos1) == str):

        print ("\aPara poder realizar el análisis de datos seleccione un archivo valido")

    else:

        gastos = gastos1.copy()
        tamanio = gastos.shape
        if (len(tamanio) == 1):            # Si el array viene de una dimensión, este se
            
            gastos.resize(len(gastos), 1)  # redimensiona
            
        m = gastos.size          
        media = gastos.mean()

    #*********************** * ********************************************************************
    # Prueba estadística de Helmert
    #*********************** * ********************************************************************
        
        diferencias = np.zeros ((m, 1))    # Se crea la matriz de diferencias (xi - media)
        s = 0 
        c = 0

        for i in range (m):

            diferencias [i, 0] = gastos [i, 0] - media
            if i > 0:
                if ((diferencias[i-1, 0] > 0) and (diferencias[i, 0] > 0)) or ((diferencias[i-1, 0] < 0) and (diferencias[i, 0] < 0)):
                    s = s + 1
                else:
                    c = c + 1
                    
        h = (m - 1) ** 0.5
        
        #*********************** * ********************************************************************
                    # Salidas de Helmert
        #*********************** * ********************************************************************

        print ('\n \n')
        print ('*********************Prueba estadistica de Helmert *************************')
        print ('Total de secuencias S = ', s)
        print ('Total de cambios C = ', c)
        
        if ((-h <= (s-c)) and ((s-c) <= h)):

            condicionH = 'La condición "-(nj-1)<=(S-C)<=(nj-1)" se cumple.'
            resMh = 'La Muestra es Homogénea'                      #resultado muestra Helmert
            print ('S - C =', s-c ,', y (n - 1)^0.5 = ', h)
            print ('\n', condicionH)
            print (-h,"<=", s - c, "<=", h)
            print (resMh)
            
        else:
            
            condicionH = 'La condición "-(nj-1)<=(S-C)<=(nj-1)" no se cumple.'
            resMh = 'La Muestra no es Homogénea'
            print ('La muestra no es homogénea:S - C =', s-c ,', y (n - 1)^0.5 = ', h)
            print (condicionH)
            print (-h,"<=", s - c, "<=", h)

        print ('***************************************************************************')
        print ('\n')

    #*********************** * ***************************************************************    
    # PRUEBA ESTADISTICA T DE STUDENT
    #*********************** * ******************************************************************** 

        n = m // 2
        muestra1 = gastos[0:n, 0]
        muestra2 = gastos[n:m, 0]          # Checar en el caso n impar
        media1 = muestra1.mean() 
        media2 = muestra2.mean()                            
        desvEst1 = 0                       # Desviación estandart
        desvEst2 = 0
        
        for i in range (n):

            prueba = (muestra1 [i,] - media1) ** 2
            desvEst1 = desvEst1 + prueba
            prueba = (muestra2 [i,] - media2) ** 2
            desvEst2 = desvEst2 + prueba
            
        desvEst1 = (desvEst1 / (n-1)) ** 0.5    
        desvEst2 = (desvEst2 / (n-1)) ** 0.5    
        
        d = (n ** -1) + (n ** -1)
        f = ((n * desvEst1**2) + (n * desvEst2 ** 2)) / (n + n -2)
        
        td = (media1 - media2) / ((f * d) ** 0.5)
        
        v = (2 * n) - 2                    # grados de libertad
        a = st.t.interval(0.95, v)         # t de student para una probabilidad de 95% y grados de libertad v
        Tv = a[-1]                         # O nivel de significancia del 5%
        td = abs (td)
        
        #*********************** * ********************************************************************
                    # Salidas T de Student
        #*********************** * ********************************************************************

        print ('*****************Prueba estadistica de T de Student*************************')

        if (td < Tv):
            
            resMTs = 'La Muestra es Homogénea'
            condicionT = 'td<Tv, se cumple'
            print (resMTs)
            print ('La condición: td = ', td,'<','Tv = ', Tv, ", se cumple")
    
        else:
            
            resMTs = 'La Muestra no es Homogénea'
            condicionT = 'td>Tv, no se cumple'
            print (resMTs)
            print ('td = ', td,'>','Tv = ', Tv, ', no se cumple')
    
        print ('***************************************************************************')   
        print ('\n')    
        
    #*********************** * ****************************************************************    
    #PRUEBA ESTADISTICA DE CRAMER
    #*********************** * ********************************************************************

        n60 = (m * 60) // 100
        n30 = (m * 30) // 100
        
        bloque1 = gastos [-n60:m, 0]       # se crean los arrays para los dos
        bloque2 = gastos [-n30:m, 0]       # diferentes bloques
        
        desvEst = 0
        for i in range (m):

            prueba = (gastos [i, 0] - media) ** 2
            desvEst = desvEst + prueba

        desvEst = np.sqrt(desvEst / (m - 1))     
        
        x30 = bloque2.mean()
        x60 = bloque1.mean()
        
        tau30 = (x30 - media) / desvEst
        tau60 = (x60 - media) / desvEst
        
        t60 = ((n60 * (m - 2)) / ((m - n60) * (1 + tau60 ** 2))) * abs(tau60)
        t30 = ((n30 * (m - 2)) / ((m - n30) * (1 + tau30 ** 2))) * abs(tau30)
        
        #*********************** * ********************************************************************
                    # Salidas de pruebas estadistica de Cramer
        #*********************** * ********************************************************************

        print ('************************PRUEBA ESTADISTICA DE CRAMER**********************')

        if ((t30 < Tv) and (t60 < Tv)):
            
            resMc = 'La Muestra es Homogénea'
            condicionC = 'Se cumple la condicion (t30 < Tv) y (t60 < Tv)'
            print (resMc)
            print ('t60=', t60)
            print ('t30=', t30)
            print ('Tv=', Tv)
            print (condicionC)
    
        else:
            
            resMc = 'La Muestra no es Homogénea'
            condicionC = 'No se cumple la condicion (t30 < Tv) y (t60 < Tv)'
            print (resMc)
            print ('t60=', t60)
            print ('t30=', t30)
            print ('Tv=', Tv)
            print (condicionC)
            
        print ('****************************************************************************')
        
    #*********************** * ***************************************************************    
    # PRUEBA INDEPENDENCIA DE EVENTOS 
    #*********************** * ********************************************************************
    
        diferencias2 = np.zeros ((m, 1))   # Se crea la matriz de diferencias (xi - media)^2
        
        for i in range (m):

            diferencias2 [i, 0] = (gastos [i, 0] - media) ** 2

        suma = diferencias2.sum()
        k = m // 3
        matrizK = np.zeros ((m, k))
        matrizRk = np.zeros ((k, 4))
        n = 0
        
        for i in range (k):

            n = n + 1
            prueba = 0
            
            for j in range (m - n):
                
                matrizK [j, i] = (gastos [j, 0] - media) * (gastos [j + n, 0] - media)
                prueba = prueba + matrizK [j, i]
            
            matrizRk [i, 2] = prueba / suma
            matrizRk [i, 0] = n
            matrizRk [i, 1] = (-1 - (1.96 * (m - i - 1) ** 0.5)) / (m - i)
            matrizRk [i, 3] = (-1 + (1.96 * (m - i - 1) ** 0.5)) / (m - i)
        
    #*********************** * ********************************************************************
    # Correlograma
    #********************************************************************************************
  
        numeracion = matrizRk [:, 0]
        limiteInf = matrizRk [:, 1]
        limiteMax = matrizRk [:, 3]
        rK = matrizRk [:, 2]
        fig, axes =  plt.subplots()
        axes.plot(numeracion, limiteInf, 'r-')
        plt.legend(['Limite superior e inferior'], loc='lower left')
        axes.plot(numeracion, limiteMax, 'r-')
        axes.plot(numeracion, rK, 'b.')
        plt.xlabel("k")
        plt.ylabel("")
        plt.title("Correlograma")
        plt.grid(True, which='both', color='g', linestyle='-', linewidth=0.5)
        plt.show()
        fig.savefig("salidas/Correlograma.png", dpi=1500) 
        print ('\n')
        print ('****************************************************************************')
        print ('     Tabla. Cálculo de los valores de rk y de los limites de confianza.\
                \n \n ', matrizRk)
        print ('****************************************************************************')
        
        #*********************** * ********************************************************************
        # Se genera el DataFrame
        #********************************************************************************************

        diccPest = {'----------': ['**********', resMh, 'Total de secuencias. S: ', 'Total de \
        cambios C:', 'S-C:', '(n-1)^0.5:', condicionH, -h, '**********', resMTs, 'td:', \
        'Tv:', condicionT, '**********', resMc, 't60', 't30', 'Tv', condicionC], '~~~~~~\
        ~~~~': ['Prueba estadística de Helmert', '----------', s, c, s-c, h, '----------',\
        '<= 5 <=', 'Prueba estadística de T de student', '-----------', td, Tv, '--------\
        --', 'Prueba estadística de Cramer', '----------', t60, t30, Tv, '----------'], \
        '__________': ['**********', '----------', '----------', '----------', '----------'\
        , '----------', '----------', h, '**********', '----------', '----------', '----\
        ------', '----------', '**********', '----------', '----------', '----------', \
        '----------', '----------']}

        dF = pd.DataFrame(diccPest)
        return(dF)
#*********************** * ********************************************************************
# Fin _Homogeneidad_Independencia (◕‿◕)
#********************************************************************************************