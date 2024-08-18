# -*- coding: utf-8 -*-
"""
Created on Fri May 03 23:07:02 2024

@author: searCampuzano
searcampuzanob@gmail.com

********************************************************************************************
Estimación de parametros por el método de momentos y máxima verosimilitud (es el mismo en este caso).
********************************************************************************************

**************************AJUSTE DISTRIBUCION NORMAL**************************************

* Salida pseudoaleatoria (random variates, "rvs")
* Funcion de densidad de prob. (probability density function, "pdf")
* Funcion de distribucion (cumulative distribution functio, "cdf")


Methods

rvs(loc=0, scale=1, size=1, random_state=None)	Random variates.
pdf(x, loc=0, scale=1)	    Probability density function.
logpdf(x, loc=0, scale=1)	 Log of the probability density function.
cdf(x, loc=0, scale=1)	    Cumulative density function.
logcdf(x, loc=0, scale=1)	Log of the cumulative density function.
sf(x, loc=0, scale=1)	   Survival function (1 - cdf — sometimes more accurate).
logsf(x, loc=0, scale=1)	Log of the survival function.
ppf(q, loc=0, scale=1)	   Percent point function (inverse of cdf — percentiles).

isf(q, loc=0, scale=1)	   Inverse survival function (inverse of sf).

moment(n, loc=0, scale=1)	Non-central moment of order n
stats(loc=0, scale=1,      moments='mv')	Mean(‘m’), variance(‘v’), skew(‘s’), and/or kurtosis(‘k’).
entropy(loc=0, scale=1)	   (Differential) entropy of the RV.
fit(data, loc=0, scale=1)	Parameter estimates for generic data.
expect(func, loc=0, scale=1, lb=None, ub=None, conditional=False, **kwds)	Expected value of a function (of one argument) with respect to the distribution.
median(loc=0, scale=1)	   Median of the distribution.
mean(loc=0, scale=1)	      Mean of the distribution.
var(loc=0, scale=1)	      Variance of the distribution.
std(loc=0, scale=1)	      Standard deviation of the distribution.
interval(alpha, loc=0, scale=1)	Endpoints of the range that contains alpha percent of the distribution

https://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.stats.norm.html

"""

import pandas as pd
import numpy as np
import scipy.stats as st
import pylab as pl


def DN(gastos, mediaDr, desvEst):

    if (type(gastos) == str):
        print ("\aPara poder realizar el análisis de datos seleccione un archivo valido")
    else:
        tamanio = gastos.shape
        if (len(tamanio) == 1):            # Si el array viene de una dimensión, este se
            gastos.resize(len(gastos), 1)  # redimensiona
        gastos01 = gastos.copy()
        m = gastos01.size
        rv1 = st.norm
                                        # Se crean los arrays de salida.
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
            matriz01 [i, 0] = i + 1                     # Se crea la columna pra el No de Orden
            matriz01 [i, 1] = gastos01.max()            # Se crea la columna que contiene los gastos registrados ordenado en forma descendente
            j = gastos01.argmax()                       # Returns the indices of the maximum values along an axis
            gastos01 [j, 0] = -1
            matriz01 [i, 2] = (m + 1)/matriz01[i,0]     # Se crea la columna con los periodos de retorno (Tr)
                                                        # a partir de la columna de los gastos registrados ordenados
            matriz01 [i, 3] = 1-(1/matriz01[i, 2])      # Se crea la columna para los valores de F(x)
            matriz01 [i, 4] = -rv1.isf(matriz01[i, 3])  # Se hace el ajuste (z)
    
        #********************************************************************************************
        # Error Standart
        #********************************************************************************************

        EE = 0
        for i in range (m):
            matriz01 [i, 5] = mediaDr + (desvEst * matriz01 [i, 4])  # Datos ajustados (x)
            gastos01 [i, 0] = (matriz01 [i, 1] - matriz01 [i, 5]) ** 2
            EE = EE + gastos01 [i,0]
                                                        # Se estima el Error Estandart para 2 parametros

        EE = (EE/(m-2)) ** 0.5
        print ("Distribución Normal, Error Estandart (Momentos y Máxima verosimilitud): ", EE)
        EEstandart[0, 0] = EE

        #********************************************************************************************
        # Se llena la matriz02
        #********************************************************************************************

        n=12
        for j in range (n):
            matriz02 [j, 1] = 1.0 - (1.0/matriz02 [j, 0])            # se crea la columna de F(x)
            matriz02 [j, 2] = -rv1.isf(matriz02[j, 1])               # Se realiza el ajuste (z)
            matriz02 [j, 3] = mediaDr + (desvEst * matriz02 [j, 2])  # Valores Extrapolados

        #********************************************************************************************
        # Se crea el DataFrame final
        #********************************************************************************************

        columnas = ['No Orden', 'Valor Registrado','Tr (Anios)', 'F(x)', 'z', 'Valor Ajustado']
        cD = pd.DataFrame(matriz01, columns = columnas)
        cD.insert(6, 'Tr', matriz02 [:, 0])
        cD.insert(7, 'F(X)', matriz02 [:, 1])
        cD.insert(8, 'Z', matriz02 [:, 2])
        cD.insert(9, 'Valor Extrapolado', matriz02 [:, 3])
        cD.insert(10, 'Error Estandart "Momentos y Máxima verosimilitud"', EEstandart [:, 0])

        #*********************** * ********************************************************************
        # Graficos
        #********************************************************************************************

        error0 = str(EE)
        error1 = "Distribucion Normal\n EE="
        titulo = error1 + error0

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
        pl.ylabel("Gastos (m³/s)")
        pl.semilogx(True)
        pl.grid(True, which='both', color='g', linestyle='-', linewidth=0.5)

        pl.subplot(2,1,2)
        pl.scatter(tR, dReg, label='Datos Registrados', marker='.')
        pl.plot(dTrExtrap, dExtrap, color="r", linewidth="1.0", linestyle="-", label ="Datos Extrapolados")
        pl.legend(loc="best")
        pl.ylabel("Gastos (m³/s)")
        pl.xlabel("(Tr) Periodos de Retorno")
        pl.semilogx(True)
        pl.grid(True, which='both', color='g', linestyle='-', linewidth=0.5)

        pl.savefig("salidas/DistribucionNormalb.png", dpi=1000)
        pl.show()
        #*********************** * ********************************************************************
        # se manda el data frame al principal
        #********************************************************************************************

        return(cD)

        #*********************** * ********************************************************************
        # Fin!ヽ(。_°)ノ !
        #********************************************************************************************

# https://matplotlib.org/api/markers_api.html#module-matplotlib.markers
"""
  marker description
    "."	point
    ","	pixel
    "o"	circle
    "v"	triangle_down
    "^"	triangle_up
    "<"	triangle_left
    ">"	triangle_right
    "1"	tri_down
    "2"	tri_up
    "3"	tri_left
    "4"	tri_right
    "8"	octagon
    "s"	square
    "p"	pentagon
    "P"	plus (filled)
    "*"	star
    "h"	hexagon1
    "H"	hexagon2
    "+"	plus
    "x"	x
    "X"	x (filled)
    "D"	diamond
    "d"	thin_diamond
    "|"	vline
    "_"	hline
    TICKLEFT	tickleft
    TICKRIGHT	tickright
    TICKUP	tickup
    TICKDOWN	tickdown
    CARETLEFT	caretleft (centered at tip)
    CARETRIGHT	caretright (centered at tip)
    CARETUP	caretup (centered at tip)
    CARETDOWN	caretdown (centered at tip)
    CARETLEFTBASE	caretleft (centered at base)
    CARETRIGHTBASE	caretright (centered at base)
    CARETUPBASE	caretup (centered at base)
    "None", " " or ""	nothing
    '$...$'	render the string using mathtext.
    verts	a list of (x, y) pairs used for Path vertices. The center of the marker is located at (0,0) and the size is normalized.
    path	a Path instance.
"""
