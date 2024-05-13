# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 21:25:08 2018

@author: Sear Campuzano 
    
email: searcampuzanob@gmail.com

Estimadores de parámetros para la distribución Gamma de 3 parámetros, por el método de 
máxima verosimilitud.
Los parametros son:


"""


import math

def parametros(datos):
    
#******************************************************************************************** 
    
    def estimadores(x0, datos):
        sum1 = 0.0      # suma(xi-xo)
        sum2 = 0.0      # suma(xi-xo)^-1
        sum3 = 0.0      # suma (Ln(xi-xo))
        n = len(datos)
        
        #print ("X0: ", x0)
        for i in range (n):
            suma1 = datos[i, 0] - x0
            sum1 = sum1 + suma1
            suma2 = 1 / (datos[i, 0] - x0)
            sum2 = sum2 + suma2
            suma3 = math.log(datos[i, 0] - x0)
            sum3 = sum3 + suma3
        #print ("suma1: ", sum1, "suma2: ", sum2, "suma3: ", sum3)
        #********************************************************************************************
        beta = 1 / (1-(n**2/(sum1*sum2))) 
        #print ("Beta: ", beta)
        #********************************************************************************************
        alfa = sum1/n-n/sum2
        #print ("Alfa: ", alfa)
        #********************************************************************************************
        psi = math.log(beta+2) - pow(2*(beta+2),-1) - pow(12*(beta+2),-2) + pow(120*(beta+2),-4)\
        - pow(252*(beta+2),-6) - pow(beta+1,-1) - pow(beta, -1)  # Psi, aproximación de la función bigamma
        #print ("Psi: ", psi)
        #******************************************************************************************** 
        
        Fb = sum3 - (n*math.log(alfa)) - (n*psi)
        #print ("F(b)", Fb)
        return (beta, alfa, Fb, psi)
#********************************************************************************************  
        
    b = 99 * (datos[datos.argmin(), 0] / 100)  # 90% del dato minimo
    a = -10000
    tol = 1e-13
    ra = estimadores (a, datos)  # resultados de f(a)
    fa = ra[2]                   # f(a)
    
    for i in range (1000):  
        p = a + ((b - a) / 2)
        #print ("i: ", i)
        rp = estimadores(p, datos)   # resultados de f(p)
        fp = rp[2]                   # f(p)
        
        
        if ((fp == 0) | (((b - a) / 2) < tol)):
            #print ("El valor de X0 es:", p)
            #print ("i =", i)
            break
        if ((fa * fp) > 0):
            a = p
            fa = fp
        if ((fa * fp) < 0):
            b = p
        
    #print ("\n\nXo: ", p)
    #print ("Beta: ", rp[0])
    #print ("Alfa: ", rp[1])
    #print ("F(Xo): ", rp[2])
    #print ("Psi(beta): ", rp[3])
    return (p, rp[0], rp[1])
#********************************************************************************************
#  Fin
#********************************************************************************************