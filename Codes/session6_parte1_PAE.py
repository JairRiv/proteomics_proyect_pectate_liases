import json, glob, os, re
import numpy as np
import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict 

#Buscamos el path de bionnfo en donde esta la informacion de alpha fold
BASE = os.path.expanduser("~/LCG/3rdyear/bioinfo/alpha_fold")
PROTEIN_CHAIN = "A" #Especificamos que sea la de la proteina, no importa mucho porque no hay ligando

#Creamos una lista de carpetas que contienen los resultados de AlphaFold, con niveles "low" y "high" y 5 repeticiones cada uno
carpetas = [f"fold_{nivel}_{i}" for nivel in ("low", "high") for i in range(1, 6)]

#Esta funcion busca el primer archivo que coincida con el patron dentro de la carpeta
def buscar(carpeta, patron):
    hits = sorted(glob.glob(os.path.join(BASE, carpeta, "**", patron), recursive=True))
    return hits[0] if hits else None

#Este vector filas va a contener diccionarios con los resultados de cada carpeta
filas = []
#El for recorre cada carpeta, extrae el nombre de la secuencia, y busca los archivos relevantes para calcular PAE promedio
for c in carpetas:
    nombre = c.replace("fold_", "").replace("_", "")   #Solo para mostrar en la tabla, no es necesario para el calculo
    #Lit, solo borra el _


    #Busca el  vector del PAE y encuentra el promedio para cada prediccion
    f_pae = buscar(c, "*full_data_0.json")
    if f_pae is None:
        cands = [f for f in glob.glob(os.path.join(BASE, c, "**", "*confidences*.json"), recursive=True)
                 if "summary" not in os.path.basename(f)]
        f_pae = sorted(cands)[0] if cands else None

    if f_pae:
        d = json.load(open(f_pae))
        pae = np.array(d["pae"])
        mask = np.array(d["token_chain_ids"]) == PROTEIN_CHAIN
        pae_prom = pae[np.ix_(mask, mask)].mean()
    else:
        pae_prom = np.nan
    
    filas.append({"secuencia": nombre, "PAE_promedio": pae_prom})

#Creamos un dataframe con la prediccion y el PAE primedio
df = pd.DataFrame(filas)
#Crea un csv con el pae
df.to_csv("resultados_alphafold.csv", index=False)