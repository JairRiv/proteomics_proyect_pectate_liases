# Analisis proteomico de las proteínas de las pectato liasas

**Por: Diana Barrientos, Leonardo González y Jair Rivera**

**22/Nov/2026**

Proyecto del módulo de Biología Estructural (LCG, ENES Juriquilla, UNAM) sobre la
familia enzimática de las **pectato liasas**. Combina tres líneas de trabajo:

1. **Exploración de estructuras en el PDB** con la API de RCSB (R).
2. **Análisis de secuencias generadas por un modelo de lenguaje (LLM)** frente a
   secuencias naturales (Python).
3. **Diseño de un *binder* de patulina con ProteinMPNN** y evaluación de las
   secuencias mediante `global_score` (Python).



## Parte A. Exploración de pectato liasas en el PDB (R)

**Archivo:** `pectate_lyases_initial_exploration.Rmd`

La familia se definió con una búsqueda *full text* de la frase **"pectate lyase"**
en la RCSB Search API, lo que recupera entradas cuyas anotaciones mencionan la
frase y no solo aquellas con ese nombre. Los datos adicionales (método, fecha de
depósito, longitudes, ligandos) se obtienen con la API GraphQL de RCSB.

Preguntas que responde:

1. Número de estructuras de la familia.
2. Distribución por método experimental.
3. Organismos más representados.
4. Distribución de resoluciones (rayos X).
5. Línea de tiempo de depósitos.
6. Distribución de longitudes de las proteínas.
7. Cinco ligandos más frecuentes.
8. Exploración estructural de **7BBV** (pectato liasa B de *Verticillium dahliae*)
   en Mol\*: plegamiento, sitio activo (Arg208), sitio de unión a Ca²⁺ e
   interpretación.

### Resultados principales

- **413** estructuras en el PDB.
- Método: **356** rayos X, **51** criomicroscopía electrónica, **6** RMN en solución.
- Organismos más frecuentes: *Enterobacteria phage HK620* (31), *Escherichia coli*
  (26), *Dickeya chrysanthemi* (25), *Lederbergvirus P22* (23) y *Homo sapiens* (20).
  Los 10 principales concentran cerca del 42 % de las entradas.
- La mayoría de las estructuras de rayos X tienen resolución entre 1.5 y 2.5 Å.
- Depósitos por año: primer depósito en 1993, con un máximo de 39 en 2015.
- Longitud de las proteínas: mediana de 373 residuos y media de 446 (rango 18 a 3,250).
- Ligandos más frecuentes: `CA` (82), `GOL` (44), `NA` (44), `SO4` (37) y `TRS` (30).

### Cómo ejecutarlo

Requiere R y los paquetes `httr2`, `jsonlite` y `rmarkdown`:

```r
install.packages(c("httr2", "jsonlite", "rmarkdown"))
rmarkdown::render("pectate_lyases_initial_exploration.Rmd")
```

El documento consulta RCSB en vivo, así que necesita conexión a internet y los
números pueden variar si el PDB se actualiza.

---

## Parte B. Secuencias generadas por un LLM vs. naturales (Python)

**Archivo:** `codigo.py`

El FASTA del LLM tiene el formato:

```
>seq_88_high_3 PPL: 7.39
MLSRIPVTASALAGVLGLCQVVNPGTSSSTFSEKAVYVTSAWAKKAVTNNQAMNITQKAS...
```

El script:

1. Extrae identificador, secuencia y **perplejidad** de cada registro.
2. Selecciona las **5 secuencias de menor** y las **5 de mayor** perplejidad y las
   guarda en `parte1_perplejidad_top_bottom.csv` (columnas: `Identificador`,
   `Secuencia`, `Perplejidad`).
3. Verifica que cada secuencia del CSV coincida con la del FASTA original.
4. Descarga secuencias naturales de UniProt (dominio Pfam PF00544) con la API REST
   y muestrea 1,000 al azar (`random.seed(42)`).
5. Compara las longitudes de ambos conjuntos en un histograma
   (`parte1_longitudes.png`).

> Las etiquetas `low`/`high` del identificador no son confiables: el ordenamiento
> se hace siempre con el valor numérico de `PPL`.

### Resultados

- 1,000 secuencias del LLM; 20,663 naturales disponibles en UniProt (consulta del
  22 de noviembre de 2026), de las que se muestrean 1,000.
- Longitud de las secuencias del LLM: media 354, mediana 382.
- Longitud de las naturales: media 429, mediana 404.
- Las secuencias generadas son, en promedio, más cortas que las naturales.

### Dependencias

```bash
pip install biopython pandas matplotlib requests
```

| Librería                  | Uso                                       |
| ------------------------- | ----------------------------------------- |
| `biopython` (`Bio.SeqIO`) | Leer FASTA.                               |
| `pandas`                  | Tablas, ordenamiento y exportación a CSV. |
| `matplotlib`              | Histograma de longitudes.                 |
| `requests`                | Descarga desde la API REST de UniProt.    |

### Uso

```bash
python codigo.py
```

---

## Parte C. Binder de patulina con ProteinMPNN

A partir de `binder_patulin_315.pdb` (diseñado con RFdiffusion; ligando `I8L` =
patulina en la cadena B) se generan secuencias con [ProteinMPNN](https://github.com/dauparas/ProteinMPNN):

- **100 secuencias** por condición.
- **3 temperaturas** de muestreo: 0.01, 0.1 y 1.
- **2 condiciones:** residuos de unión fijados (`fixed_residues_binder_315.txt`)
  o sin fijar.
- Total: **6 archivos `.fasta`** (3 temperaturas × 2 condiciones).

Análisis posterior:

- Para cada archivo, las 3 secuencias de menor y las 3 de mayor `global_score`
  se reúnen en un CSV con las columnas `Identificador`, `Secuencia`,
  `Global_score`, `T_de_generación` y `residuos_fijados`.
- Dos *boxplots* de `global_score` por temperatura (uno para residuos fijados y
  otro para residuos sin fijar).

En ProteinMPNN, un `global_score` **menor** indica que la secuencia es más
compatible con el esqueleto. El primer registro de cada salida (`.fa`) es la
secuencia nativa y se excluye del análisis.

**Estado:** en progreso. Fecha límite de la entrega: 24 de noviembre de 2026.

```bash
git clone https://github.com/dauparas/ProteinMPNN
```

---

## Reproducibilidad

| Elemento                          | Valor                               |
| --------------------------------- | ----------------------------------- |
| Semilla del muestreo de naturales | 42                                  |
| Semilla de ProteinMPNN            | 37                                  |
| Consulta de UniProt               | `xref:pfam-PF00544` (formato FASTA) |
| Consulta de RCSB                  | *full text*: `"pectate lyase"`      |
| Fecha de las consultas            | Noviembre de 2026                   |

Las bases de datos (PDB, UniProt) se actualizan constantemente, por lo que los
conteos pueden cambiar con el tiempo.

## Limitaciones

- La búsqueda *full text* en el PDB también recupera proteínas con plegamiento
  similar que no son pectato liasas canónicas (por ejemplo, proteínas de fagos).
- El número de estructuras depositadas no refleja la diversidad real de la familia.
- La perplejidad baja no garantiza que una secuencia sea funcional; las secuencias
  del LLM requieren validación estructural y experimental.
- 7BBV no incluye el sustrato ni la conformación activa; para estudiarlos haría
  falta simulación (por ejemplo, BioEmu).
