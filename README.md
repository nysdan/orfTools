# orfTools

Paquete R para encontrar, traducir y visualizar ORFs (open reading frames) en
genomas virales pequeños, desarrollado para el proyecto
[SGLS](../../../README.md) (identificación de proteínas Sgl putativas en
fagos ssRNA(+) Fiersviridae/Leviviridae).

## Instalación

Dentro del entorno pixi del proyecto (o cualquier R con `devtools`):

```r
devtools::load_all("workflow/R/orfTools")   # para desarrollo, sin instalar
# o
devtools::install("workflow/R/orfTools")
```

Dependencias: `Biostrings` (Bioconductor), `dplyr`, `ggplot2`, `gggenes`,
`tibble`.

## Funciones

### `find_all_orfs(fasta_file, min_length = 30, keep_stop_codon = FALSE)`

Escanea un archivo FASTA (una o más secuencias) en los 6 marcos de lectura
(3 en la hebra +, 3 en la reversa complementaria), extrae los ORFs con
[`find_orfs_in_seq()`] y traduce cada uno a proteína. Devuelve un tibble con
una fila por ORF: coordenadas genómicas (siempre expresadas sobre la hebra
+, incluso para ORFs de la hebra -), marco de lectura, longitud y secuencia
de aminoácidos.

```r
orfs <- find_all_orfs("genoma.fasta", min_length = 30)
```

Por defecto `keep_stop_codon = FALSE` elimina el `*` final que deja
`Biostrings::translate()`, dejando una secuencia proteica lista para
herramientas río abajo (p. ej. predicción de dominios TM con DeepTMHMM).

### `find_orfs_in_seq(seq, strand, frame_offset, min_length = 30)`

Función de bajo nivel usada por `find_all_orfs()`: escanea **una** hebra en
**un** marco de lectura fijo y devuelve la lista de ORFs encontrados (start
codon ATG/GTG/TTG hasta el primer stop en fase). Rara vez se llama
directamente — se expone porque es útil para tests o para escanear
secuencias que no vienen de un archivo FASTA.

### `detect_and_classify_overlaps(orf_table)`

Toma una tabla de ORFs (la salida de `find_all_orfs()`) y devuelve todos los
pares que se solapan en coordenadas genómicas, clasificando el solapamiento
como `"Embedded overlap"` (uno anidado dentro del otro) o `"Terminal
overlap"` (parcial), e indicando si comparten hebra y/o marco de lectura.
Relevante para candidatos Sgl, que típicamente solapan en marco +1 con el
gen mayor que los contiene (ver `docs/methods.md` del proyecto).

```r
overlaps <- detect_and_classify_overlaps(orfs)
```

### `plot_orfs(orf_table)`

Grafica el mapa génico de una tabla de ORFs con `gggenes`, separando por
hebra (arriba/abajo del eje) y por marco de lectura.

```r
plot_orfs(orfs)
```

## Notas de diseño

- **Definición de ORF:** dado un start codon, se reporta el primer stop en
  fase que le sigue; tras reportar un ORF la búsqueda continúa después de
  ese stop. Esto significa que ORFs anidados dentro de uno ya reportado
  (mismo marco, start más interno) **no** se listan por separado — se asume
  traducción desde el primer start disponible en cada marco.
- **Coordenadas:** en la tabla final, los ORFs de la hebra `-` se reportan
  con `start`/`end` ya convertidos a coordenadas de la hebra `+`, igual que
  los de la hebra `+`, para poder graficar/comparar ambas hebras en el mismo
  eje.
- **Codones de inicio:** se aceptan ATG, GTG y TTG (no solo ATG), ya que
  varios ORFs pequeños (candidatos a Sgl) usan starts no canónicos.
