# Entorno de ejecución

## Referencia (outputs vigentes del proyecto)

- **R 4.5.2 (2025-10-31 ucrt), Windows 11 x64.** Paquetes y versiones en `paquetes_R.csv`, tomados de los `sessionInfo.txt` de los outputs 01–16. `randomForest` 4.7-1.2 se instaló en una biblioteca local (`analisis/R/library/`, no versionada); `15_ML_agrupado.R` la añade a `.libPaths()` si existe.
- **Python 3.13, Linux, para el script 17:** numpy 2.5.3, scipy 1.18.1 y networkx 3.7 (`requirements-python.txt`). La librería `core_periphery_sbm` va incluida en `analisis/python/vendor/` (commit 90fb289 de https://github.com/ryanjgallagher/core_periphery_sbm, sin modificar). No hay que instalarla: no está en PyPI.
- Las huellas de los outputs de referencia están en `referencia/outputs_md5_windows.csv` (01–16, Windows) y `referencia/outputs_md5_linux.csv` (17–20, Linux). Los outputs 18–20 se generaron en Linux con R 4.3.3 y las mismas versiones de paquetes, a partir de los outputs vigentes de Windows; 18 y 20 se ejecutaron dos veces (sobre los outputs de Windows y sobre los de la prueba en limpio) con resultados idénticos byte a byte.

## Instalación

R, opción A (recomendada): `renv`.

```r
install.packages("renv")
renv::restore()          # lee renv.lock (versiones de paquetes_R.csv)
```

R, opción B: instalar a mano las versiones de `paquetes_R.csv`, por ejemplo con `remotes::install_version("igraph", "2.2.1")`.

**Mínimos imprescindibles:** igraph ≥ 2.0 (los scripts usan `as_undirected()` y `get_edge_ids()`) y dplyr ≥ 1.1 (`relationship =` en las uniones).

Python:

```bash
python -m pip install -r requirements-python.txt
```

Linux: hace falta una configuración regional UTF-8 con orden alfabético C (`LANG=C.UTF-8`, `LC_COLLATE=C`, sin `LC_ALL`). Sin UTF-8, `source(..., encoding = "UTF-8")` falla con «invalid input found on input connection». Sin colación C cambian los órdenes (ver más abajo). `ejecutar_todo.sh` fija ambas. En Windows, los scripts ajustan `LC_CTYPE` ellos mismos.

## Prueba de reproducción en limpio (22-10-2026)

La cadena 01–16 se ejecutó desde cero en un contenedor Linux (Ubuntu 24.04, 2 núcleos), con la disposición de `DATOS.md` y sin el correo privado. Se usó **R 4.3.3** con las **mismas versiones de paquetes** que la referencia, instaladas desde las fuentes de CRAN. Los resultados se compararon con `98_comparar_outputs.py`:

- Los **16 scripts terminan sin error**, y `16_verificacion_integrada.R` supera sus **60 comprobaciones**.
- Ficheros de referencia: **131 idénticos** (CSV y GraphML, estos salvo el fin de línea CRLF/LF de Windows) y **11 CSV numéricamente iguales**, con una diferencia relativa máxima de 3,2·10⁻¹⁴ (redondeo en coma flotante: 02, 04, 11 y 15).
- Hay 17 binarios (`.rds`, `.png`, `.svg`) con huella distinta, como es esperable. Los `.rds` clave tienen contenido idéntico con `all.equal()`: tokens de 03, modelos de 06, 594.000 métricas de 07 y permutaciones de 12.
- Detalle por fichero en `referencia/prueba_limpia_linux_20261022.csv`.

Problemas detectados y resueltos en la prueba:

1. **Locale.** En Linux hace falta UTF-8 (`LANG=C.UTF-8`) **y orden alfabético C** (`LC_COLLATE=C`), como en la sesión de referencia de Windows. Con la colación ICU de `C.UTF-8`, `sort()`/`merge()` ordenan «agregar» antes que «Pensar». Eso cambia el orden de filas de 08 y 14 y, sobre todo, la auditoría 09 de `Sequences.R`, cuyo error depende precisamente de ese orden. `ejecutar_todo.sh` fija ambas variables.
2. **igraph < 2.0** (igraph 1.6.0 de Ubuntu 24.04): fallan 05, 07, 12, 14 y 16. Con dplyr 1.1.4, purrr 1.0.2 y ggplot2 3.4.4 (Ubuntu) el resto de scripts funciona.
3. `16_verificacion_integrada.R` leía la carpeta `07_redes_HC_bootstrap_v2`, que solo existe en la máquina de referencia. Ahora usa `_v2` si existe y, si no, la carpeta sin sufijo. Además, un fichero opcional ausente en ambos momentos (el correo privado) ya no invalida la comprobación de huellas (cambio del 22-10-2026; la versión anterior está en `_obsoleto/`).

La réplica SBM (17) es determinista: con las mismas semillas, una segunda ejecución reproduce nodo a nodo las particiones y las MDL.

## Tiempos aproximados (2 núcleos, Linux)

| Script | Tiempo |
|---|---|
| 01–04, 06, 08–11, 13, 14, 16 | < 10 s cada uno |
| 05 (12.000 grafos) | ≈ 2 min |
| 07 (11.000 extracciones, 99.000 grafos) | ≈ 6,5 min |
| 12 (999 permutaciones) | ≈ 3 min |
| 15 (LOPO, 2.000 bootstrap, 199 nulos) | ≈ 4,5 min |
| **Total 01–16** | **≈ 17 min** |
| 17 (270 + 162 cadenas, 2 procesos) | ≈ 61 min |
| 17_explorar_configuracion.sh | ≈ 36 min |
| 18 (210 redes, 211 comprobaciones) | < 30 s |
| 19 (840 cadenas hub-and-spoke + 270 de capas, 2 procesos) | ≈ 2 h 15 min |
| 20 (brechas y sensibilidades) | < 10 s |
| 19 complementario (`--solo igualada`, particiones de igualadas) | ≈ 46 min |
| 21 (Tabla 2 fiel a la R2, 352 modelos mixtos con sensibilidades) | ≈ 15 s |
| 22 (niveles en 10 extracciones igualadas) | ≈ 40 s |

## Entorno separado TabPFN v2 (opcional; D-026, Windows22-10-2026)

Python3.12.14, PyTorch2.8.0+cu128, tabpfn2.2.1, numpy2.2.6, pandas2.2.3, scipy1.15.3 y scikit-learn1.6.1. Lista completa fijada en `requirements-tabpfn-v2-lock.txt`. RTX5070Ti16GB;4 hilosTorch,1 proceso de preprocesamiento, float32 y algoritmos deterministas. La preparación40 y figura43 usan el entorno R existente. El checkpoint v2 está fijado por revisión y SHA256 en42; su licencia se conserva al descargarlo. No se usa el modelo por defecto de una biblioteca más nueva.

Crear un venv separado (`python -m venv analisis/python/.venv_tabpfn_v2`, desde `code/`), activarlo y, con ese Python, instalar primero `python -m pip install torch==2.8.0 --index-url https://download.pytorch.org/whl/cu128`, después `python -m pip install -r ../requirements-tabpfn-v2-lock.txt`. En Linux puede requerirse una distribuciónCUDA compatible distinta; registrar cualquier diferencia y no prometer identidad numérica entre GPU/sistemas diferentes.41 exige CUDA para esta configuración. No instalar el lock sobre el entorno SBM ni sobre los paquetes históricos.

41 trabaja offline con pesos locales y desactiva telemetría antes de importar TabPFN. TEMP/cachés se sitúan dentro de `analisis/python/.runtime_tabpfn_v2`.42 descarga solo pesos/licencia públicos. No se envían registros del estudio a una API. Los entornos, pesos y cachés no se incluyen en el borrador de repositorio.

## Diferencias esperables

- Los `.rds` no son idénticos byte a byte entre versiones de R o sistemas: incluyen la versión de serialización y atributos del entorno. El contenido se compara cargándolos.
- Las figuras PNG/SVG dependen de las fuentes y del dispositivo gráfico del sistema.
- `entradas_md5.csv` y `sessionInfo.txt` registran rutas y versiones locales, y se excluyen de la comparación.
- **clubSandwich** 0.5.10 (solo `45_tabla2_robustez_residuos.R`, D-029). En Windows se instala en `analisis/R/library/`, como randomForest; no se actualizan paquetes históricos.
