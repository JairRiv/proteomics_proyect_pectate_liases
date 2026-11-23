#!/usr/bin/env python3
"""
Calcula el pLDDT promedio de cada modelo de AlphaFold3 y el promedio global.

Para cada archivo *_model_N.cif toma el valor B_iso_or_equiv de los átomos
de la cadena proteica (registros ATOM) y calcula la media aritmética.
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

# Tokens de mmCIF: valores simples o entre comillas simples/dobles
TOKEN_RE = re.compile(r"'[^']*'(?=\s|$)|\"[^\"]*\"(?=\s|$)|\S+")


def leer_atom_site(texto: str):
    """Devuelve una lista de diccionarios, uno por átomo, del bucle _atom_site."""
    lineas = texto.splitlines()
    i, n = 0, len(lineas)

    while i < n:
        if lineas[i].strip() != "loop_":
            i += 1
            continue

        j = i + 1
        encabezados = []
        while j < n and lineas[j].strip().startswith("_"):
            encabezados.append(lineas[j].strip())
            j += 1

        if not encabezados or not encabezados[0].startswith("_atom_site."):
            i = j
            continue

        columnas = [h.split(".", 1)[1] for h in encabezados]
        tokens = []
        k = j
        while k < n:
            s = lineas[k].strip()
            if s == "" or s.startswith("#"):
                k += 1
                continue
            if s.startswith(("loop_", "data_", "_")):
                break
            tokens.extend(TOKEN_RE.findall(s))
            k += 1

        ncol = len(columnas)
        if len(tokens) % ncol != 0:
            raise ValueError(
                f"Número de valores ({len(tokens)}) no es múltiplo de columnas ({ncol})"
            )

        return [
            dict(zip(columnas, tokens[r:r + ncol]))
            for r in range(0, len(tokens), ncol)
        ]

    raise ValueError("No se encontró el bucle _atom_site en el archivo")


def plddt_promedio(texto: str):
    """Calcula (media de pLDDT, número de átomos) para la cadena proteica."""
    atomos = leer_atom_site(texto)
    valores = [
        float(a["B_iso_or_equiv"])
        for a in atomos
        if a["group_PDB"] == "ATOM"
    ]
    if not valores:
        raise ValueError("No hay átomos ATOM con B_iso_or_equiv")
    return math.fsum(valores) / len(valores), len(valores)


def recolectar_cifs(carpeta: Path):
    """Yield (nombre_modelo, ruta) solo de los modelos de la predicción.

    glob (no rglob) evita entrar a templates/, y el patrón *_model_*.cif
    excluye los templates, cuyos nombres contienen 'template_hit'.
    """
    for cif in sorted(carpeta.glob("*_model_*.cif")):
        yield cif.stem, cif


def main():
    parser = argparse.ArgumentParser(description="pLDDT promedio de los modelos de AF3")
    parser.add_argument("carpeta", type=Path, help="Carpeta con las predicciones")
    parser.add_argument("-o", "--salida", type=Path, default=Path("plddt_promedio.csv"))
    args = parser.parse_args()

    if not args.carpeta.is_dir():
        sys.exit(f"Error: {args.carpeta} no es una carpeta")

    resultados = []
    for nombre, ruta in recolectar_cifs(args.carpeta):
        media, natomos = plddt_promedio(ruta.read_text(encoding="utf-8"))
        resultados.append((nombre, natomos, media))
        print(f"{nombre}: pLDDT promedio = {media:.2f} ({natomos} átomos)")

    if not resultados:
        sys.exit("No se encontraron archivos *_model_*.cif en la carpeta")

    medias = [m for _, _, m in resultados]
    promedio_global = math.fsum(medias) / len(medias)

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
