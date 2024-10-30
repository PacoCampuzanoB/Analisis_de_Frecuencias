import pandas as pd
import numpy as np
import distros.Gumbel
import distros.GumbelDouble
import os

df = pd.read_excel('salida1_pruebas.xlsx', sheet_name=0, header=0, usecols=None)

columnas = list(df.columns)

maximos = []
estaciones = []
for i in range(1, len(columnas), 3):
    maximos.append(columnas[i])
    estaciones.append(columnas[i].split('_')[1])

for i in range (len(estaciones)):
    print('Estación --> ', estaciones[i])
    estacion1 = df[maximos[i]]
    estacion1 = estacion1.dropna()
    print('La estación ', estaciones[i], ' tiene', len(estacion1), ' datos\n\n')
    if len(estacion1) > 12:
        datos = estacion1.to_numpy()
        media = np.mean(datos)
        sm = np.std(datos, ddof=1)
        excel_salida = 'salidas/' + estaciones[i] + '_gumbel_dgumbel.xlsx'
        hoja_nombre_g = estaciones[i] + '_Gumbel'
        hoja_nombre_dg = estaciones[i] + '_Doble_Gumbel'
        estacion = estaciones[i]
        with pd.ExcelWriter(excel_salida) as writer:
            dfDG = distros.Gumbel.DG(datos, media, sm, estacion)
            dfDG.to_excel(writer, sheet_name=hoja_nombre_g, index=False)
        
            dfDdg = distros.GumbelDouble.DDG(datos, estacion)
            dfDdg.to_excel(writer, sheet_name=hoja_nombre_dg, index=False)
    else:
        print("La estación ", estaciones[i], 'no se procesó por falta de datos\n\n')
