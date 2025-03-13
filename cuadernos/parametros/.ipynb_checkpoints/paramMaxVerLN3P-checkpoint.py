# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 07:11:19 2018

@author: searCampuzano
searcampuzanob@gmail.com

estimación de parametros por máxima verosimilitud para la distribución
LogNormal de tres parametros.

Los parámetros por Máxima Verosimilitud son:

            muY=sumatoria(ln(xi-xo)/n)
            sigmaY=(1/n)Sumatoria(ln(xi-xo)-MuY)^2

    El estimador xo se obtiene al resolver la siguiente ecuación

    F(xo)=Sumatoria((1/(xi-xo))(My-Sy2))-Sumatoria(ln(xi-xo)/(xi-xo))
    
Para encontrar los parámetros se estimará la función a por medio del metodo
de bisección, se usará una tolerancia de 1*10-13. el valor inicial será 
el el 95% del valor más pequeño de la muestra.

"""


import math

def parametros(datos):
    
#********************************************************************************************
    # Función esimadores
#********************************************************************************************
    
    def estimadores(x0, datos):
        Mu_y = 0.0
        sigma_y = 0.0    
        
        for i in range (len(datos)):
            mu = math.log(datos[i, 0] - x0)     # ln (xi - x0)
            Mu_y = Mu_y + mu                    # Mu_y = suma ln(xi - x0)
                                                
        mu_y = Mu_y / len(datos)                # mu_y = suma [ln (xi - x0) / n]
        #print ("mu_y: ", mu_y)
        
        for i in range (len(datos)):
            sigma = (math.log(datos[i, 0] - x0) - mu_y) ** 2
            sigma_y = sigma_y + sigma
        
        sigma_y = (sigma_y / len(datos)) ** 0.5
        #print ("sigma_y: ", sigma_y)
        
        # F(xo)
        
        b = 0.0
        a = 0.0
        for i in range (len(datos)):
            d = (mu_y - pow(sigma_y, 2)) / (datos[i, 0] - x0)  
            c = (math.log(datos[i, 0] - x0)) / (datos[i, 0] - x0)
            a = a + d
            b = b + c
            
        F1 = a - b
        #print ("F1 = ", F1)
        return (mu_y, sigma_y, F1, i)
    
#********************************************************************************************
#********************************************************************************************
#********************************************************************************************

    b = 95 * (datos[datos.argmin(), 0] / 100)  # 95% del dato minimo
    a = -10000.0
    tol = 1e-13
    
    for i in range (1000):  
        p = a + ((b - a) / 2)
        rp = estimadores(p, datos)   # resultados de f(p)
        fp = rp[2]                   # f(p)
        ra = estimadores (a, datos)  # resultados de f(a)
        fa = ra[2]                   # f(a)
    
        if ((fp == 0) | (((b - a) / 2) < tol)):
            #print "El valor de X0 es:", p
            #print "i =", i
            break
        if ((fa * fp) > 0):
            a = p
            fa = fp
        if ((fa * fp) < 0):
            b = p
    print ("\n\nx0:", p)
    print ("mu_y:", rp[0]) 
    print ("sigma_y:", rp[1])
    print ("F(xo):", rp[2])
    return (p, rp)
#********************************************************************************************
#  Fin! ヽ(。_°)ノ 
#********************************************************************************************