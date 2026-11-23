#!/usr/bin/env python3
#Codigo para buscar secuencias naturales y generadas por LLM, analizar su longitud y perplejidad, identificar residuos de unión en un PDB, ejecutar ProteinMPNN y analizar los resultados.
#Creado por: Diana Barrientos, Leonardo Gonzales y Jair Rivera
#Date: 22/Nov/2026

#Primero importamos las librerías necesarias para el análisis de secuencias y estructuras de proteínas

#Os permite verificas si es que existen archivos, crear carpaetas y unir rutas
#Re permite trabajar con expresiones regulares, importantes para extraer información de los headers de las secuencias
#sys permite trabajar con argumentos de línea de comandos
#random permite generar números aleatorios, útil para seleccionar secuencias naturales de manera aleatoria
#glob permite buscar archivos en directorios usando patrones
#shutil permite copiar y mover archivos y directorios
#subprocess permite ejecutar comandos del sistema desde python
import os, re, sys, glob, random, shutil, subprocess

# Requests permite hacer solicitudes HTTP, es importante para descargar secuencias naturales de UniProt
import requests

#Permite manipular dataframes, armar csv, agrupar, etc.
import pandas as pd

#Hace figuras lindas
import matplotlib.pyplot as plt

#Tiene modulos de biopython para leer secuencias y estructuras de proteínas
from Bio import SeqIO



#Archivos a usar
LLM_FASTA = "team_pl_llms.fasta" #Secuencias generadas por LLM con perplejidad en el header
NAT_FASTA = "secuencias_naturales.fasta" # Secuencias naturales, si no existe se descarga de UniProt
SEED = 37 #Semilla para extraer secuencias naturales (extraemos 1000 de forma aleatoria)

filas = []

#Este for lee el FASTA del LLM y extraer la perplejidad del header
for r in SeqIO.parse(LLM_FASTA, "fasta"):
    m = re.search(r"PPL:\s*([0-9.]+)", r.description) #Buscamos en el header de cada secuencia la perplejidad, que está en el formato "PPL: 0.1234"
    if m is None:
        raise ValueError(f"Sin PPL en: {r.description}") #Si no encuentra la perplejidad, lanza un error
    filas.append({"Identificador": r.id,
                  "Secuencia": str(r.seq), #En dado caso, guardamos el identificador, la secuencia y la perplejidad en un diccionario
                  "Perplejidad": float(m.group(1))})

df = pd.DataFrame(filas)
print(f"{len(df)} secuencias leídas (esperado: 1000)")
assert df["Identificador"].is_unique, "Hay identificadores repetidos"

#Creamos un dataframe con las filas y verificamos que no haya identificadores repetidos
df = pd.DataFrame(filas)

#Posterior a extraerlos, seleccionamos los valores que tengan mayor y menor perplejidad
sel = pd.concat([df.nsmallest(5, "Perplejidad"), df.nlargest(5, "Perplejidad")])
#Posterior a esto, guardamos los resultados en un CSV
sel.to_csv("parte1_perplejidad_top_bottom.csv", index=False)

#Verificación de que los identificadores y secuencias coinciden con el FASTA original
orig = {r.id: str(r.seq) for r in SeqIO.parse(LLM_FASTA, "fasta")}
assert all(orig[i] == s for i, s in zip(sel["Identificador"], sel["Secuencia"]))
print("Si paso el check de verificación de secuencias:)")

# Para las secuencias naturales, si no existe el archivo, lo descargamos de UniProt usando la API REST
if not os.path.exists(NAT_FASTA):
    print("No existen secuencias naturales, las vamos a descargar de UniProt")
    #Mandamos la solicitud a UniProt para obtener las secuencias naturales de la familia Pfam PF00544 en formato FASTA
    r = requests.get("https://rest.uniprot.org/uniprotkb/stream", #Mandamos la solicitud a uniprot, solicitamos las secuencias de la familia Pfam PF00544 en formato FASTA
                     params={"query": "xref:pfam-PF00544", "format": "fasta"},
                     timeout=300) #Si tarda más de 5 minutos, lanza un error
    r.raise_for_status()
    open(NAT_FASTA, "w").write(r.text) #Guardamos las secuencias naturales en un archivo FASTA

nat = [str(r.seq) for r in SeqIO.parse(NAT_FASTA, "fasta")] #Extraemos las secuencias naturales del archivo FASTA
print(f"{len(nat)} naturales disponibles") #Contamos cuántas secuencias naturales tenemos disponibles
random.seed(SEED) #Fijamos la semilla para que la selección aleatoria sea reproducible
nat = random.sample(nat, min(1000, len(nat))) #Escogemos 1000 secuencias de forma aleatoria

len_llm = df["Secuencia"].str.len() #Creamos una serie con la longitud de las secuencias generadas por LLM
len_nat = pd.Series([len(s) for s in nat]) #Creamos una serie con la longitud de las secuencias naturales
fig, ax = plt.subplots(figsize=(8, 5)) #Creamos una figura y un eje para graficar la distribución de longitudes de las secuencias generadas por LLM y las naturales
maxlen = int(max(len_llm.max(), len_nat.max())) #Calculamos la longitud máxima de las secuencias para definir los bins del histograma
bins = range(0, maxlen + 25, 25) #Definimos los bins del histograma, de 0 a la longitud máxima, con un paso de 25 en 25
ax.hist(len_llm, bins=bins, alpha=0.6, color="tab:blue", #Colores lindos para las secuencias generadas por LLM y las naturales
    label=f"Generadas por LLM (n={len(len_llm)})")
ax.hist(len_nat, bins=bins, alpha=0.6, color="tab:orange", #Para las secuencias naturales, usamos un color diferente y una etiqueta diferente
label=f"Naturales (n={len(len_nat)})")
ax.set_xlabel("Longitud (número de residuos)")
ax.set_ylabel("Número de secuencias")
ax.set_title("Distribución de longitudes: secuencias LLM vs. naturales")
ax.legend()
plt.tight_layout()
plt.savefig("parte1_longitudes.png", dpi=300) #Guardamos la figura en un archivo PNG
plt.close()
print(f"LLM: media={len_llm.mean():.0f}, mediana={len_llm.median():.0f} | "f"Naturales: media={len_nat.mean():.0f}, mediana={len_nat.median():.0f}") #Imprimimos las estadísticas de longitud de las secuencias

