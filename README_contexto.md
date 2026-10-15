# workflow/R/orfTools — nota de contexto

Este es el paquete R propio para predicción exhaustiva de ORFs (todos los marcos), pieza del **paso 2** del pipeline (rama Levi). Tiene su propio `README.md`, `DESCRIPTION` y control de versiones git independiente.

## Rol en el pipeline
`find_all_orfs()` predice ORFs en los 6 marcos sobre un genoma ya confirmado como fago ssRNA, sin depender de señales de inicio canónicas (RBS) — necesario porque las Sgls carecen de RBS propio y Prodigal las saltea.

## Ejemplo de uso
```r
devtools::load_all(".")
orfs <- find_all_orfs(genoma_fasta, min_aa = 30, max_aa = 150)
```

## Nota de mantenimiento
Antes de consolidar scripts sueltos de `../scripts/*.R` acá, verificar solapamiento de funciones (mismo nombre = decidir versión vigente). Pendiente registrado en la bitácora de migración.
