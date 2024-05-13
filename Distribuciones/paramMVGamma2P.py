# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 12:28:19 2018

@author: searCampuzano
searcampuzanob@gmail.com

Estimadores de parámetros para la distribución Gamma de 2 parámetros, por el método de 
máxima verosimilitud.
Los parametros son:
    alfa = muX/beta
    beta
    
* beta se obtiene al resolver:
    
    F(beta) = muY - ln(muX) + ln(beta) - psi(beta) = 0
    
* empleando como valor inicial:
    
    beta = (1+raiz(1+(4/3*C)))/(4*C)    
    
* donde:
    
    C = ln(muX) - muY
    
    muY = 1/n Sumatoria|i=1, n|(ln(muX))
    muX = 1/n Sumatoria|i=1, n| xi
    
Donde psi(beta) es la función digamma de beta.

* Para resolver la ecuación F(beta) se utiliza el método de Newton-Raphson.

    p = beta - F(beta)/F'(beta)
    
Donde F'(beta) = 1/beta - psi'(beta)

psi'(beta) es la aproximación de la función trigamma de beta.

"""

import math

def parametros(datos):
    muY = 0
    muX = 0
    for i in range (len(datos)):
        muY = muY + math.log(datos[i, 0])
        muX = muX + datos[i, 0]
    
    muY = muY/len(datos)
    muX = muX/len(datos)
    C = math.log(muX) - muY
    beta = (1 + pow(1 + ((4*C)/3), 0.5)) / (4*C)
    tol = 1e-13
    
    for i in range (1000):
        psi = math.log(beta+2) - pow(2*(beta+2),-1) - pow(12*(beta+2),-2) + pow(120*(beta+2),-4)\
        - pow(252*(beta+2),-6) - pow(beta+1,-1) - pow(beta, -1)  # Psi(beta), aproximación de la función bigamma
        psiP = pow(beta+2,-1) + pow(2*(beta+2),-2) + pow(6*(beta+2),-3) - pow(30*(beta+2),-5)\
        + pow(42*(beta+2),-7) - pow(30*(beta+2),-9) + pow(beta+1,-2) + pow(beta, -2)  # Psi'(beta), aproximación de la función trigamma
        
        Fb = muY - math.log(muX) + math.log(beta) - psi
        p = beta - (Fb/(pow(beta,-1)-psiP))
        diferencia = abs(p-beta)
        if diferencia <= tol:
            beta = p
            break
        else:
            beta = p

    alfa = muX/beta
    return(alfa, beta)
#********************************************************************************************
#  Fin
#********************************************************************************************