#-*- coding: utf-8 -*-

import os
import pandas as pd

direccion = r'C:\\Users\\sear2\\Music\\estaciones_clim\\bcs\\'
archivos = os.listdir(direccion)
archivos.sort()
direcciones = []
estaciones = []
multiplos = []

for i in range(len(archivos)):
    direcciones.append(direccion + archivos[i])
    estacion = archivos[i][3:8]
    estaciones.append('annio_' + estacion)
    estaciones.append('maximo_' + estacion)
    estaciones.append('num_datos_' + estacion)
    print('Estación --> ', estacion)
    if i == 0:
        multiplos.append(i)
    else:
        multiplos.append(i * 3)
df = pd.DataFrame(columns=estaciones,
                index=range(200))

for i in range(len(direcciones)):
    
    with open(direcciones[i], encoding='utf-8') as mensual:
        lista_mensual = mensual.readlines()
    
    variable = 23
    renglon = lista_mensual[23].split('\t')
    inicio = renglon[0]
    annio = []
    record = []
    no_reg = []
    while inicio != "MÍNIMA":
        annio.append(inicio)                           # Año
        lista = []                                     # Lista de los datos maximos mensuales (12)
        for j in range(1, 13):
            if (renglon[j] != ''):
                lista.append(float(renglon[j]))
            else:
                lista.append(-999.0)
        record.append(max(lista))
        no_reg.append(int(renglon[-1].split('\n')[0])) # Número de registros mensuales
        variable += 1
        renglon = lista_mensual[variable].split('\t')
        inicio = renglon[0]
        columna = multiplos[i]
    for k in range(len(annio)):                        # Se llena el dataframe
        df.iat[k, columna] = annio[k]
        df.iat[k, columna + 1] = record[k]
        df.iat[k, columna + 2] = no_reg[k]
df.to_excel( r'C:\\Users\\sear2\\Music\\estaciones_clim\\salidas\\salida_pruebas0_bcs.xlsx', 
            header=True, 
            sheet_name='salida', 
            index=False)
