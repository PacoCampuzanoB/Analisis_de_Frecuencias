# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 09:01:26 2018

Comienza el proyecto de análisis de frecuencias

@author: SearCampuzano
searcampuzanob@gmail.com

* sript que carga un archivo con formato txt, csv, xlsx, y lo importa como
    un array numpy

* el cometido de este script es que este se integre en un modulo que 
    será parte de el proyecto de analisis de frecuencias.

Notas:
+ Falta aplicar una excepción para validar arrays que no tengan datos tipo 
    float.
"""
import numpy as np
import pandas as pd
from tkinter import filedialog

def dImport():

##########################################################################################
# Formato
# Se incluyen dos excepciones, una en caso de que no esté cargado el archivo con el que
# se va a trabajar y la otra para el caso de que se seleccione cancelar al importar el
# archivo.
# la referencia de las excepciones la dejaré en un commit de git.
##########################################################################################
    
    while True:
        
        try:   
            rutaU = filedialog.askopenfilename(initialdir = "/home/sear/Documentos/",
                                        title = "Abrir solo archivos .txt, .csv, .xlsx",
                                        filetypes = (("Text files",
                                                    "*.txt*"),
                                                    ("csv files",
                                                    "*.csv*"),
                                                    ("Excel files",
                                                    "*.xlsx*"),
                                                    ("all files",
                                                    "*.*")))

        except TypeError:
            print ("\nSeleccionó cancelar, no se importó el archivo.")
            respuesta = -1
            return(respuesta)
            break
        
        ruta = str(rutaU)
        formato = ruta[ruta.find("."): ] 
        #print ("Ruta:", ruta)
        
##########################################################################################
# Nombre del archivo
##########################################################################################
        
        if (ruta != ""):
            listaBarras = []
            for i in range (len(ruta)):
                if (ruta[i] == "/"):
                    listaBarras.append(i + 1)
            archivo = ruta[listaBarras[-1]: ]    
            print ("\nLa ruta es:", ruta)      
            print("\nEl archivo es de formato:", formato) 
            print("El nombre del archivo es:", archivo)
            
##########################################################################################
# Se abre el archivo según sea su formato
##########################################################################################

            if (formato == ".txt"):
                datos = np.loadtxt (ruta, dtype=float)
                return(datos)
                break
            elif (formato == ".csv"):    
                datos = pd.read_csv(ruta, header=None).values
                datos.astype(float)
                return(datos)
                break
            elif (formato == ".xlsx"):
                xlsx = pd.ExcelFile(ruta)
                df = pd.read_excel(xlsx)
                datos = df.values
                datos.astype(float)
                return(datos)
                break
            else:
                print("\nSeleccionó un archivo con formato diferente a los formatos\
                '.txt', 'csv', o 'xlsx', por favor seleccione otra vez: ")

##########################################################################################
# Fin datosImport
##########################################################################################