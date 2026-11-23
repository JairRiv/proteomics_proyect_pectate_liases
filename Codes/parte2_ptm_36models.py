import json, glob, os, re
import pandas as pd

# Carpeta donde están fold_num1 ... fold_num36
BASE = os.path.expanduser("~/LCG/3rdyear/bioinfo/proteomics_proyect_pectate_liases/alpha_fold/Parte2")

# Busca el primer archivo que coincida con el patrón dentro de una carpeta
def buscar(carpeta, patron):
    hits = sorted(glob.glob(os.path.join(BASE, carpeta, "**", patron), recursive=True))
    return hits[0] if hits else None

# Orden NUMÉRICO: con sorted() normal quedaría fold_num1, fold_num10, fold_num11, ..., fold_num2
carpetas = sorted(
    [d for d in os.listdir(BASE) if re.fullmatch(r"fold_num\d+", d)],
    key=lambda d: int(re.search(r"\d+", d).group()),
)
print(f"{len(carpetas)} carpetas encontradas")   # esperado: 36

filas = []
for c in carpetas:
    # Modelo 0 = el de mayor ranking que entrega AlphaFold
    f_sum = buscar(c, "*summary_confidences_0.json")

    if f_sum is None:
        print(f"[AVISO] {c}: no encontré summary_confidences_0.json")
        ptm = float("nan")
    else:
        with open(f_sum) as f:
            ptm = json.load(f)["ptm"]

    filas.append({"fold": c, "pTM": ptm})

df = pd.DataFrame(filas)

# Chequeos rápidos
assert len(df) == 36, f"Esperaba 36 folds, encontré {len(df)}"
print("pTM faltantes:", df["pTM"].isna().sum())
print(df["pTM"].describe().round(3))

df.to_csv("resultados_alphafold_ptm.csv", index=False)
df