#!/usr/bin/env python3
"""
Calcula el pLDDT promedio de cada modelo de AlphaFold3 obtenido de las
5 secuencias con mejor perplejidad y las 5 secuencias con peor perplejidad 
obtenidas de la Tarea 1. Después calcula el promedio global de todos los modelos.

Para cada archivo *_model_N.cif (ya que así los devuelve AlphaFold) 
toma el valor B_iso_or_equiv de los átomos de la cadena proteica 
(registros ATOM) y calcula la media.

Excluye HETATM (ligandos, iones y agua) y los templates de templates/.

Uso:
    python session6_parte1_pLDDT.py carpeta_predicciones/ -o plddt_promedio.csv
"""

import argparse
import csv
import math
import re
import sys
from pathlib import Path

# En los archivos mmCIF se puede encontrar texto entre comillas y seguido de espacio o fin de línea,
# también entre comillas dobles o como normalmente sería, sin espacios.
TOKEN_RE = re.compile(r"'[^']*'(?=\s|$)|\"[^\"]*\"(?=\s|$)|\S+")

def leer_atom_site(texto: str):
    """
    Input:
    - Contenido de un archivo mmCIF como string.

    Output:
    - Lista de diccionarios uno por átomo del bucle _atom_site.
    """
    lineas = texto.splitlines() #Separa el texto en líneas
    i = 0
    n = len(lineas)

    while i < n: # Recorre las líneas hasta encontrar el bucle _atom_site,
                 # ósea la línea que diga _loop
        if lineas[i].strip() != "loop_":
            i += 1
            continue

        j = i + 1
        encabezados = [] # Va a contener los encabezados de las columnas
        while j < n and lineas[j].strip().startswith("_"): 
            # Se asegura que corresponda a una columna de datos
            encabezados.append(lineas[j].strip()) # Agrega el encabezado a la lista
            j += 1 # Va a leer todos los enabezados

        if not encabezados or not encabezados[0].startswith("_atom_site."):
            #Si no hay encabezados o no se encuentran columnas que comiencen con
            # _atom_site. entonces no es la matriz que buscamos
            i = j
            continue

        columnas = [h.split(".", 1)[1] for h in encabezados]
        # Va a contener el nombre de las columnas sin el prefijo _atom_site
        # porque ya sabemos que estamos en la matriz de átomos.
        tokens = []
        k = j
        while k < n:
            s = lineas[k].strip() # Separa la línea en tokens sin espacios al inicio ni al final
            if s == "" or s.startswith("#"): #Ignora comentarios y líneas vacías
                k += 1
                continue
            if s.startswith(("loop_", "data_", "_")): #Se termina si encuentra otra matriz de datos
                break
            tokens.extend(TOKEN_RE.findall(s)) # Agrega los tokens encontrados en la línea a la lista de tokens
            k += 1

        ncol = len(columnas) # Número de columnas que se esperan en la matriz de átomos
        if len(tokens) % ncol != 0: # Si el número de tokens no es múltiplo del número de columnas, 
                                    # entonces hay un error en el archivo
            raise ValueError(
                f"Número de valores ({len(tokens)}) no es múltiplo de columnas ({ncol})"
            )

        return [
            # Se crea un diccionario para cada átomo, donde las claves del diccionario
            # son los nombres de las columnas y los valores son los tokens correspondientes
            dict(zip(columnas, tokens[r:r + ncol]))
            for r in range(0, len(tokens), ncol)
        ]

    raise ValueError("No se encontró el bucle _atom_site en el archivo")


def plddt_promedio(texto: str):
    """
    Input:
    - Contenido de un archivo mmCIF como string.

    Output:
    - Tupla (media de pLDDT, número de átomos) para la cadena proteica.
    """
    atomos = leer_atom_site(texto) #Contiene la lista de diccionarios con los átomos del archivo mmCIF
    valores = [
        # Revisa que el átomo sea de la cadena proteica (group_PDB == "ATOM") y toma su valor B_iso_or_equiv
        float(a["B_iso_or_equiv"])
        for a in atomos
        if a["group_PDB"] == "ATOM"
    ]
    if not valores:
        raise ValueError("No hay átomos ATOM con B_iso_or_equiv")
    # Devuelve la media de los valores B_iso_or_equiv y el número de átomos ATOM encontrados
    return math.fsum(valores) / len(valores), len(valores) 


def recolectar_cifs(carpeta: Path):
    """
    Input:
    - Carpeta donde se encuentran los archivos *_model_*.cif de la predicción de AlphaFold3.

    Output:
    - Generador de tuplas (nombre_modelo, ruta) para cada archivo *_model_*.cif encontrado en la carpeta.
    """
    for cif in sorted(carpeta.glob("*_model_*.cif")):
        # Solo busca dentro de la carpeta especificada, no en subcarpetas por el .glob
        # Como los archivos de AlphaFold3 de interés terminan en _model_#_.cif, así se filtra
        yield cif.stem, cif
        # Un return que pausa para continuar en la siguiente iteración y no devolver una lista con
        # todos los modelos que se descargaron de AlphaFold3 y los devuelve sin la terminal .cif
        

def main():
    # Pedir los argumentos de la línea de comandos
    parser = argparse.ArgumentParser(description="pLDDT promedio de los modelos de AF3")
    parser.add_argument("carpeta", type=Path, help="Carpeta con las predicciones")
    parser.add_argument("-o", "--salida", type=Path, default=Path("plddt_promedio.csv"))
    args = parser.parse_args()

    if not args.carpeta.is_dir():
        sys.exit(f"Error: {args.carpeta} no es una carpeta")

    resultados = [] # Almacena: nombre del modelo, número de átomos y pLDDT promedio
    for nombre, ruta in recolectar_cifs(args.carpeta): # Se itera en todos los archivos *_model_*.cif
        media, natomos = plddt_promedio(ruta.read_text(encoding="utf-8"))
        resultados.append((nombre, natomos, media))
        print(f"{nombre}: pLDDT promedio = {media:.2f} ({natomos} átomos)")

    if not resultados:
        sys.exit("No se encontraron archivos *_model_*.cif en la carpeta")

    medias = [m for _, _, m in resultados] # Lista de los pLDDT promedio de cada modelo
    promedio_global = math.fsum(medias) / len(medias) # Calcula el promedio global de todos los modelos

    # Crea el archivo CSV con los resultados
    with open(args.salida, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["modelo", "n_atomos", "plddt_promedio"])
        for nombre, natomos, media in resultados:
            writer.writerow([nombre, natomos, round(media, 2)])
        writer.writerow(["promedio_global", "", round(promedio_global, 2)])

    print(f"\nPromedio de {len(medias)} modelos: {promedio_global:.2f}")
    print(f"Resultados guardados en {args.salida}")


if __name__ == "__main__":
    main()
