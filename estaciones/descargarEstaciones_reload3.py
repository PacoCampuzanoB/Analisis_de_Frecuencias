#-*- coding: utf-8 -*-
"""
Created on Tu jul 05 10:53:41 2022

@author: searcampuzanob@gmail.com

--- Para descargar la información de estaciones climatológicas hay que disponer de un archivo xlsx
    con una columna con las claves de cada estación por renglón, como en el ejemplo "estaciones.xlsx"---
    
Nota: si se quiere decargar los archivos de datos diarios, entonces hay que comentar las líneas 33 y 34 que
corresponden a los mensuales y extremos, si se quieren descargar los datos mensuales hay que comentar las líneas 
32 y 34 y así :V, para descargar datos diarios, mensuales y extremos no hay que comentar las líneas.
"""
###################################################################################
import requests
import pandas as pd
import os

#output_folder = "/home/sear/Documentos/estaciones_clim/chiapas2/"
output_folder = r"C:\Users\sear2\Music\estaciones_clim\bcs"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

df_estaciones = pd.read_excel("bcs.xlsx",
                dtype={"Estaciones":object})

df_estaciones['Estaciones'] = df_estaciones['Estaciones'].astype(str)
df_estaciones["Estaciones"] = df_estaciones["Estaciones"].str.zfill(5)
dicc = {'01': 'ags', '02': 'bc', '03': 'bcs', '04': 'camp', '05': 'coah', '06': 'col', '07': 'chis', '08': 'chih',
        '09': 'df', '10': 'dgo', '11': 'gto', '12': 'gro', '13': 'hgo', '14': 'jal', '15': 'mex', '16': 'mich', '17': 'mor',
        '18': 'nay', '19': 'nl', '20': 'oax', '21': 'pue', '22': 'qro', '23': 'qroo', '24': 'slp', '25': 'sin', '26': 'son',
        '27': 'tab', '28': 'tamps', '29': 'tlax', '30': 'ver', '31': 'yuc', '32': 'zac'}

urls = []
for i in range(len(df_estaciones)):
    estado = df_estaciones.iloc[i, 0][0:2]
    estado = dicc[estado]
    # urls.append('https://smn.conagua.gob.mx/tools/RESOURCES/Normales_Climatologicas/Diarios/' + estado + '/dia' + df_estaciones.iloc[i, 0] + '.txt')
    urls.append('https://smn.conagua.gob.mx/tools/RESOURCES/Normales_Climatologicas/Mensuales/' + estado + '/mes' + df_estaciones.iloc[i, 0] + '.txt')
    # urls.append('https://smn.conagua.gob.mx/tools/RESOURCES/Normales_Climatologicas/Med-Extr/' + estado + '/medex' + df_estaciones.iloc[i, 0] + '.txt')
print(urls)

files=[]
for i in range(len(urls)):
    derivados = urls[i].split("/")
    file = derivados[-1]
    myfile = requests.get(urls[i], verify=False)
    print(myfile)
    file_path = os.path.join(output_folder, file) 
    #open(file, 'wb').write(myfile.content)
    with open(file_path, 'wb') as f:
        f.write(myfile.content)
###################################################################################
# (☞ﾟヮﾟ)☞ Fin
###################################################################################
