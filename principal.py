#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: searCampuzano
searcampuzanob@gmail.com

* Script principal del proyecto análisis de frecuencias

"""
#########################################################################################

from tkinter import *
from tkinter import messagebox
import Distribuciones.Exp1Param
import Distribuciones.Exp2Param
import Distribuciones.Gamma2Par
import Distribuciones.Gamma3Par
import Distribuciones.GumbelDouble 
import Distribuciones.Gumbel
import Distribuciones.LogNorm2Param
import Distribuciones.LogNorm3Param
import Distribuciones.LogPearsonIII
import Distribuciones.Normal
import Importador.datos_txt_csv_xlsx
import PruebasEstadisticas.Homogeneidad_Independencia
import numpy as np
import pandas as pd

#########################################################################################

def iniciar():
    global datos1
    datos1 = Importador.datos_txt_csv_xlsx.dImport()
    return(datos1)

#########################################################################################

def todo(datos1):

    media = np.mean(datos1)
    sP = np.std(datos1)
    sM =  np.std(datos1, ddof=1)
    with pd.ExcelWriter("salidas/todosLosAjustes.xlsx") as writer:
        dfPe = PruebasEstadisticas.Homogeneidad_Independencia.hI(datos1)
        dfPe.to_excel(writer, sheet_name='PruebasEstadisticas', index=False)
        
        dfDn = Distribuciones.Normal.DN(datos1, media, sM)
        dfDn.to_excel(writer, "Normal")
        
        dfDln2p = Distribuciones.LogNorm2Param.DLN2P(datos1)
        dfDln2p.to_excel(writer, sheet_name="LogNormal2Param", index=False)
        
        dfDln3p = Distribuciones.LogNorm3Param.DLN3P(datos1, media, sM, sP)
        dfDln3p.to_excel(writer, "LogNormal3Param")
        
        dfDe1p = Distribuciones.Exp1Param.DE1P(datos1, media)
        dfDe1p.to_excel(writer, sheet_name="Exponencial_1_Param", index=False)
        
        dfDe2p = Distribuciones.Exp2Param.DE2P(datos1, media, sM)
        dfDe2p.to_excel(writer, sheet_name="Exponencial_2_Param", index=False)
        
        dfDG = Distribuciones.Gumbel.DG(datos1, media, sM)
        dfDG.to_excel(writer, sheet_name="Gumbel", index=False)
        
        dfDg2p = Distribuciones.Gamma2Par.DG2P(datos1, media, sM)
        dfDg2p.to_excel(writer, sheet_name="Gamma_2_Param", index=False)
        
        dfDg3p = Distribuciones.Gamma3Par.DG3P(datos1, media, sM, sP)
        dfDg3p.to_excel(writer, sheet_name="Gamma_3_Param", index=False)
        
        dfDlpIII = Distribuciones.LogPearsonIII.DLPIII(datos1, media, sM, sP)
        dfDlpIII.to_excel(writer, sheet_name="LogPearsonIII", index=False)
        
        dfDdg = Distribuciones.GumbelDouble.DDG(datos1)
        dfDdg.to_excel(writer, sheet_name="DobleGumbel", index=False)
    messagebox.showinfo(title="Info!", message="proceso terminado! ☜(ˆ▿ˆ)")
#########################################################################################

def pedidos(num):
    if ((num in lista) == True):
        pass    
    else:    
        lista.append(num)
    lista.sort()
    print("Petición: ", lista)

#########################################################################################

def ejecutar(datos1):
    
    media = np.mean(datos1)
    sP = np.std(datos1)
    sM =  np.std(datos1, ddof=1)
    with pd.ExcelWriter("salidas/ajustes_pedidos.xlsx") as writer:
        if ("PruebasEstadisticas" in lista):
            dfPe = PruebasEstadisticas.Homogeneidad_Independencia.hI(datos1)
            dfPe.to_excel(writer, sheet_name='PruebasEstadisticas', index=False)
        if ("Normal" in lista):
            dfDn = Distribuciones.Normal.DN(datos1, media, sM)
            dfDn.to_excel(writer, sheet_name="Normal", index=False)
        if ("LogNormal2P" in lista):
            dfDln2p = Distribuciones.LogNorm2Param.DLN2P(datos1)
            dfDln2p.to_excel(writer, sheet_name="LogNormal2Param", index=False)
        if ("LogNormal3P" in lista):
            dfDln3p = Distribuciones.LogNorm3Param.DLN3P(datos1, media, sM, sP)
            dfDln3p.to_excel(writer, sheet_name="LogNormal3Param", index=False)
        if ("Exponencial1P" in lista):
            dfDe1p = Distribuciones.Exp1Param.DE1P(datos1, media)
            dfDe1p.to_excel(writer, sheet_name="Exponencial_1_Param", index=False)
        if ("Exponencial2P" in lista):
            dfDe2p = Distribuciones.Exp2Param.DE2P(datos1, media, sM)
            dfDe2p.to_excel(writer, sheet_name="Exponencial_2_Param", index=False)
        if ("Gumbel" in lista):
            dfDG = Distribuciones.Gumbel.DG(datos1, media, sM)
            dfDG.to_excel(writer, sheet_name="Gumbel", index=False)
        if ("Gamma2P" in lista):
            dfDg2p = Distribuciones.Gamma2Par.DG2P(datos1, media, sM)
            dfDg2p.to_excel(writer, sheet_name="Gamma_2_Param", index=False)
        if ("Gamma3P" in lista):
            dfDg3p = Distribuciones.Gamma3Par.DG3P(datos1, media, sM, sP)
            dfDg3p.to_excel(writer, sheet_name="Gamma_3_Param", index=False)
        if ("LogPearsonIII" in lista):
            dfDlpIII = Distribuciones.LogPearsonIII.DLPIII(datos1, media, sM, sP)
            dfDlpIII.to_excel(writer, sheet_name="LogPearsonIII", index=False)
        if ("DobleGumbel" in lista):
            dfDdg = Distribuciones.GumbelDouble.DDG(datos1)
            dfDdg.to_excel(writer, sheet_name="DobleGumbel", index=False)
    messagebox.showinfo(title="Info!", message="proceso terminado! ☜(ˆ▿ˆ)")
#########################################################################################

def limpiar():
    global lista
    lista = []
    print("Sin Peticiones: ", lista)

#########################################################################################

ventana = Tk()
ventana.geometry("900x506")
ventana.title("FDP's")
imagen = PhotoImage(file="image/mar.png")
fondo = Label(ventana, image=imagen).place(x=0, y=0)
lista = []

#########################################################################################

botonTodo = Button(ventana, bg="#059", fg="white", text="Ajuste de todas las FDP's \
'ALV :V'", command=lambda: todo(datos1)).place(x=250, y=50)

botonImport = Button(ventana, text="Importar archivo",\
command=lambda: iniciar()).place(x=50, y=50)

botonSalida = Button(ventana, text = "Salir", command=ventana.destroy).place(x=450, y=400)

botonEjecutar = Button(ventana,bg="#006", fg="white", text = "Ejecutar",\
command=lambda: ejecutar(datos1)).place(x=50, y=400)

botonLimpiar = Button(ventana,bg="#006", fg="white", text = "Limpiar",\
command=lambda: limpiar()).place(x=250, y=400)

#########################################################################################
# Botones distribuciones
#########################################################################################

botonPest = Button(ventana, bg="#059", fg="white", text="Pruebas estadísticas",\
command=lambda: pedidos("PruebasEstadisticas")).place(x=50, y=80)

botonNormal = Button(ventana, bg="#059", fg="white", text="Normal",\
command=lambda: pedidos("Normal")).place(x=50, y=110)

botonImport = Button(ventana, bg="#059", fg="white", text="LogNormal 2P",\
command=lambda: pedidos("LogNormal2P")).place(x=50, y=140)

botonImport = Button(ventana, bg="#059", fg="white", text="LogNormal 3P",\
command=lambda: pedidos("LogNormal3P")).place(x=50, y=170)

botonImport = Button(ventana, bg="#059", fg="white", text="Exponencial 1P",\
command=lambda: pedidos("Exponencial1P")).place(x=50, y=200)

botonImport = Button(ventana, bg="#059", fg="white", text="Exponencial 2P",\
command=lambda: pedidos("Exponencial2P")).place(x=250, y=80)

botonImport = Button(ventana, bg="#059", fg="white", text="Gumbel",\
command=lambda: pedidos("Gumbel")).place(x=250, y=110)

botonImport = Button(ventana, bg="#059", fg="white", text="Gamma 2P",\
command=lambda: pedidos("Gamma2P")).place(x=250, y=140)

botonImport = Button(ventana, bg="#059", fg="white", text="Gamma 3P",\
command=lambda: pedidos("Gamma3P")).place(x=250, y=170)

botonImport = Button(ventana, bg="#059", fg="white", text="LogPearson III",\
command=lambda: pedidos("LogPearsonIII")).place(x=250, y=200)

botonImport = Button(ventana, bg="#059", fg="white", text="Doble Gumbel",\
command=lambda: pedidos("DobleGumbel")).place(x=50, y=230)

#########################################################################################

ventana.mainloop()

#########################################################################################
# （っ＾▿＾）っ FIN!!
#########################################################################################
