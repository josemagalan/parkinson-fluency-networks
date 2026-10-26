# Datos de entrada

El repositorio **no incluye datos**. La lista exacta de ficheros que leen los scripts, con su tamaño y huella md5, está en `datos_entrada_md5.csv`. Ese fichero lo genera `analisis/python/99_exportar_repositorio.py` a partir de los `entradas_md5.csv` de cada output y de las entradas del script 17. Las rutas son relativas a `code/`.

## Disposición esperada

```
<raíz del repositorio>/
├── datos sexo_edad_parkinson_control.csv        orden ON/OFF por paciente (columna «primera evaluacion»)
├── parkinson/                                   repositorio histórico de la R2 (solo lectura)
│   ├── data/                                    datos brutos y diccionarios psicolingüísticos
│   │   ├── FLUIDEZ PARKINSON.xlsx               respuestas de fluidez (6 hojas: 3 tareas × PD/HC)
│   │   ├── Lexesp data.xlsx                     imaginabilidad, concreción, familiaridad
│   │   ├── AoA_verbos.xlsx, AoA_sustantivos.xlsx
│   │   ├── Buscador frecuencia oral Subtlex-Esp.xlsx
│   │   ├── buscador motor content.xlsx
│   │   ├── datos sexo_edad_parkinson_control.csv, datos_sexo_edad_controles.xlsx
│   ├── networks/  results/                      GraphML, ARFF y tablas de métricas de la R2
│   ├── hub and spoke model/{directed,undirected}/csvs/   particiones core-periphery R2
│   ├── resultsCorePeriphery/                    Excel de las Tablas 2/3 R2
│   └── corePeriphery.R, Sequences.R, processing.R, naiveBayes.R,
│       mainParkinson.R, mainParkinson_Dirigidas.R   (se auditan; solo se evalúan funciones puras)
└── code/
    └── historial/elena_20261022_v3/datos sexo_edad_parkinson_control.csv   copia del adjunto (opcional*)
```

\* `03_orden_sesiones.R` comprueba, si existe esta copia, que coincide byte a byte con el CSV de la raíz (procedencia: correo de E. Herrera del 22-10-2026). El correo original (`.msg`) es privado. Solo se registra su huella si está presente y no es necesario para ejecutar.

## Comprobación

Desde `code/`:

```bash
python -I analisis/python/99_exportar_repositorio.py --comprobar   # código sincronizado
```

y comparar las huellas de los datos con `datos_entrada_md5.csv` (por ejemplo, con `md5sum` o `Get-FileHash -Algorithm MD5`). Cada output guarda además su propio `entradas_md5.csv`.

## Codificación

Las 18 particiones CSV históricas tienen codificación mixta (16 UTF-8 y 2 Latin-1). Una fila está entrecomillada entera: `n202`, la celda `partir<TAB>llegar` de controles en acciones. `02_reproducir_R2.R` y `17_replica_SBM_historico.py` lo tratan explícitamente. No hay que convertir los ficheros.

## Restricciones

Los datos de fluidez y la demografía proceden de 21 pacientes y 24 controles evaluados por el equipo de la Universidad de Oviedo. Los diccionarios psicolingüísticos son de terceros. Su publicación depende de las autorizaciones éticas y de las licencias correspondientes; está pendiente de decidir con los coautores.
