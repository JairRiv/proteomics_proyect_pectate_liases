import json, glob, os, re
import numpy as np
import pandas as pd

BASE = "/home/dbarrientos/Desktop/genomica/bioinfo_5to/Parte2/Parte2"
PROTEIN_CHAIN = "A"


def pae_promedio(ruta_json, chain=PROTEIN_CHAIN):
    """PAE promedio entre todos los pares de residuos de la cadena indicada."""
    d = json.load(open(ruta_json))
    pae = np.array(d["pae"])
    ids = np.array(d["token_chain_ids"])
    mask = ids == chain
    if not mask.any():
        # Si la cadena no existe, usa la primera que aparezca
        mask = ids == ids[0]
    return pae[np.ix_(mask, mask)].mean()


def nombre_del_job(carpeta):
    """Lee el nombre del job desde job_request.json (lista o diccionario)."""
    archivos = glob.glob(os.path.join(BASE, carpeta, "*job_request.json"))
    if not archivos:
        return carpeta
    try:
        req = json.load(open(archivos[0]))
        if isinstance(req, list):
            req = req[0]
        return req.get("name", carpeta)
    except Exception:
        return carpeta


def grupo_del_nombre(nombre):
    n = nombre.lower()
    if "low" in n:
        return "low"
    if "high" in n:
        return "high"
    return "sin_grupo"


# Detecta todas las carpetas fold_num*, ordenadas por número
carpetas = sorted(
    [c for c in os.listdir(BASE)
     if c.startswith("fold_num") and os.path.isdir(os.path.join(BASE, c))],
    key=lambda c: int(re.search(r"\d+", c).group())
)

filas = []
for c in carpetas:
    archivos = sorted(glob.glob(os.path.join(BASE, c, "*full_data_*.json")))
    if not archivos:
        print(f"[!] Sin archivos full_data en {c}")
        continue

    paes = [pae_promedio(f) for f in archivos]
    nombre = nombre_del_job(c)

    filas.append({
        "carpeta": c,
        "secuencia": nombre,
        "grupo": grupo_del_nombre(nombre),
        "PAE_modelo0": paes[0],
        "PAE_promedio_5_modelos": np.mean(paes),
        "PAE_min": np.min(paes),
        "PAE_max": np.max(paes),
        "n_modelos": len(paes),
    })

df = pd.DataFrame(filas)
pd.set_option("display.width", 200)
print(df.to_string(index=False))
df.to_csv("resultados_alphafold.csv", index=False)

# Resumen por grupo (low vs high)
if df["grupo"].nunique() > 1 or df["grupo"].iloc[0] != "sin_grupo":
    resumen = df.groupby("grupo")["PAE_promedio_5_modelos"].agg(["count", "mean", "std"])
    print("\nResumen por grupo:")
    print(resumen)
    resumen.to_csv("resumen_por_grupo.csv")
