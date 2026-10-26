# Verbal fluency networks in Parkinson's disease on and off medication

Analysis code for the article *Network-based analysis of semantic, phonological, and action fluency in Parkinson's disease on/off medication

What the code does:

1. **Data.** It builds a token-level table from the raw verbal fluency responses (semantic, phonological and action fluency), with psycholinguistic annotation and each patient's on/off assessment order.
2. **Participant-level analyses.** It compares word production and psycholinguistic means with the participant as the unit of analysis, adjusting the paired on/off comparisons for assessment session.
3. **Networks.** It builds co-occurrence and sequence networks for each group and task, both raw and matched for word output (participant-level bootstrap and permutations).
4. **Core–periphery structure.** It fits Bayesian core–periphery stochastic block models (SBM; Gallagher et al., 2021) and estimates the lexical comparisons of Table 2 with linear mixed models, plus robustness and sensitivity checks.
5. **Classification.** It evaluates penalised multinomial regression, random forest and TabPFN v2 with leave-one-participant-out cross-validation.
6. **Serial position.** It produces the psycholinguistic profiles by response position (Figure 4).

The raw data are read from the authors' original processing pipeline. Some scripts also reproduce and audit that pipeline. The data are **not** included .

## Repository layout

```
.
├── README.md  LICENSE  DATOS.md  ENTORNO.md     (DATOS/ENTORNO: detailed notes, in Spanish)
├── ejecutar_todo.sh  ejecutar_todo.ps1         run the main pipeline (01–22) in order
├── paquetes_R.csv  renv.lock                   R package versions
├── requirements-python.txt                     Python packages for the SBM (17, 19, 46)
├── requirements-tabpfn-v2-lock.txt             TabPFN v2 environment (41)
├── MANIFEST.csv                                md5 of every code file (generated)
├── datos_entrada_md5.csv                       input files the scripts expect, with md5
├── referencia/                                 md5 fingerprints of the reference outputs
├── parkinson/                                  NOT included: original pipeline (raw data, networks, partitions)
└── code/
    ├── analisis/R/                             R scripts (00_* helpers, 01–48)
    ├── analisis/python/                        Python scripts (00a–00c, 17, 19, 41, 42, 46, 49, 98, 99)
    │   └── vendor/core_periphery_sbm/          Gallagher et al. library, commit 90fb289 (Apache-2.0)
    ├── notas/                                  pre-registered TabPFN protocol read by script 41
    └── analisis/output/                        created when the scripts run (not versioned)
```

All scripts use paths relative to `code/` and read the input data from `../parkinson/`. **Run them from `code/`.** Some script messages and docstrings refer to `code/` by its working name, `revision_R3/`. Every script starts with a header in Spanish and English that gives its purpose, inputs and outputs.

## Requirements

- **R 4.5.2** (reference run on Windows 11; Linux checks used R 4.3.3) with the versions in `paquetes_R.csv` (`renv::restore()`). The code needs igraph ≥ 2.0 and dplyr ≥ 1.1, plus `nlme` (21, 22, 45, 47, 48), `clubSandwich` ≥ 0.5.6 (45), `randomForest` and `nnet` (15), and `ggplot2` and `svglite` (figures).
- **Python 3.13** with numpy, scipy and networkx (`requirements-python.txt`) for the SBM. The SBM library is vendored because it is not on PyPI.
- Linux needs `LANG=C.UTF-8` and `LC_COLLATE=C`; `ejecutar_todo.sh` sets both. On Windows the scripts set `LC_CTYPE` themselves.
- **TabPFN v2** (41): a separate Python 3.12 environment with a CUDA GPU. See `ENTORNO.md`.

## Data

The repository contains no data. The input data are verbal fluency responses from 21 patients (each assessed on and off medication) and 24 controls, the assessment order, demographics, and third-party psycholinguistic norms (LEXESP, SUBTLEX-ESP, age-of-acquisition norms, motor content). The clinical data are not publicly available; access conditions are given in the article's Data availability statement. The lexical norms are subject to their own licences. To reproduce the analyses:

- Place the files as described in `DATOS.md`.
- Check them against `datos_entrada_md5.csv`.

## How to run

```bash
bash ejecutar_todo.sh                  # Linux/macOS; --sin-replica skips 17; --desde NN resumes at step NN
./ejecutar_todo.ps1                    # Windows PowerShell 7.2+ recommended (-SinReplica, -Desde NN)
```

The launchers run 01–22. Scripts 01–16 take about 17 min, 17 about 1.5 h, and 19 about 2 h 15 min (2 cores). At the end the launchers compare the outputs with `referencia/` using `98_comparar_outputs.py`.

Each script writes to `code/analisis/output/<script>/`, together with `sessionInfo.txt` (R) or `entorno.txt` (Python) and `entradas_md5.csv`. The R scripts never overwrite an output: if the folder already exists, they create `<script>_v2`, `_v3`, and so on. Python scripts 17, 19 and 46 overwrite their folder.

Then run the remaining scripts by hand:

1. **Second run of 04, 20 and 21.** Scripts 45, 47, 48 and 49 read the `_v2` output folders of 04, 20 and 21 (the reference runs). Running 04, 20 and 21 a second time creates them.
2. **46 → 47 → 48 → 45 → 49.**
3. **30** (Figure 4).
4. **40 → 42 → 41 → 43** (TabPFN v2; needs the separate environment above).
5. Optional: 31, 32 and 44 compare runs across platforms and with earlier tables. They need folders that are not distributed.

## Scripts

| Script | Purpose | In the article |
|---|---|---|
| `00_*.R` | Shared helpers: I/O, dictionaries, network metrics, data scenarios, Table 2 model | — |
| `01_datos_tidy` | Token table (2,965 tokens), minimal cleaning, dictionary coverage | Methods |
| `03_orden_sesiones` | Individual on/off order (11 off→on, 10 on→off) and assessment session | Methods |
| `04_produccion_sesion` | Paired on−off word production adjusted for session | Results |
| `05_redes_igualadas_preparacion` | 18 group networks (3 groups × 3 tasks × co-occurrence/sequence) and raw metrics | Results |
| `06_psicolinguistica_sesion` | Participant psycholinguistic means; paired on−off models with session | Supplementary material |
| `07_redes_HC_bootstrap` | Output-matched networks, participant bootstrap | Results; supplementary material |
| `10_perfiles_posicion_corregidos` | Serial-position profiles with participant-bootstrap bands | Sensitivity |
| `11_HC_primera_sesion_influencia` | Patients' first assessment vs controls; leave-one-out influence | Supplementary material |
| `12_permutacion_secuencia_sensibilidad` | Permutation of on/off sequence labels for network metrics | Sensitivity |
| `14_preparar_entrada_SBM` | SBM input package | — |
| `15_ML_agrupado` | Leave-one-participant-out classification, bootstrap, permutation nulls | Results; supplementary material |
| `16_verificacion_integrada` | Cross-checks of 03–15 and of the recorded input md5s | — |
| `18_exportar_redes_SBM` | 210 networks for the SBM (own-group, joint, output-matched) | — |
| `python/19_SBM_redes_limpias` | SBM fits: chain consensus, P(core), description-length comparisons | Figures 1–3; Tables S1–S2 |
| `20_tabla2_core_periferia` | Alternative specification: within-participant core−periphery gaps on a joint partition | Sensitivity |
| `21_tabla2_niveles_R2` | Group × partition mixed models with the participant as unit | Tables 2 and 3 |
| `22_tabla2_niveles_igualadas` | Table 2 model on 10 output-matched subsamples | Results |
| `45_tabla2_robustez_residuos` | Table 2 robustness: CR2 cluster-robust SEs and participant bootstrap | Results; supplementary material |
| `python/46_SBM_F_secuencia_cadenas` | 50 SBM chains for the patients' phonological sequence networks | Supplementary material |
| `47_tabla2_F_secuencia_particiones` | Table 2 model across the partitions from 46 (exploratory) | Supplementary material |
| `48_gradiente_por_grupo` | Periphery−core gradient within each group | Results |
| `python/49_multiplicidad_holm_bonferroni` | Holm and Bonferroni adjustment compared with Benjamini–Hochberg | Sensitivity |
| `30_figura4_estilo_R2` | Psycholinguistic ratings by response position | Figure 4 |
| `40`, `python/42`, `python/41`, `43` | TabPFN v2 classification (fixed configuration, no tuning) | Results; supplementary material |
| `python/98`, `python/99` | Output comparison with the reference; repository export | — |

**Audit of the original pipeline.** These scripts reproduce the analyses in `parkinson/` and document how the current analysis differs from them. They are not needed for the results in the article.

| Script | Purpose |
|---|---|
| `python/00a`–`00c` | Token coverage, node merging and replication of the original Table 2 |
| `02_reproducir_R2` | Reproduces the original Tables 2/3 and quantifies merge and encoding errors |
| `08_auditoria_ML` | Checks the original Weka exports; builds leakage-free participant folds |
| `09_auditoria_secuencia` | Checks whether the original sequence code kept the production order (it did not) |
| `13_auditoria_metricas_historicas` | Recomputes the original network metrics |
| `python/17_replica_SBM_historico` | Replicates the original core–periphery SBM fits (`17_explorar_configuracion.sh`: sampler sensitivity) |
| `31`, `32`, `44` | Cross-checks against earlier tables and between Linux and Windows runs |

## Reproducibility notes

- **Seeds.** The R scripts use `set.seed(20261022)`. The SBM scripts seed each chain (1…N). With the package versions in `paquetes_R.csv` and `requirements-python.txt`, results are identical. A clean Linux run of 01–16 reproduced all reference CSV/GraphML files, exactly or within 3·10⁻¹⁴ (`referencia/prueba_limpia_linux_20261022.csv`).
- **Code version.** Comments in the scripts were edited after the reference run. The executable code is unchanged; we checked this by comparing R parse trees and Python ASTs. The outputs record the md5 of the scripts as they were run, so a fresh run records different md5s. Script 98 ignores them; only re-running 16 or 44 against the reference outputs fails on that check.
- **Known limitations** (they do not affect the reported results):
  - Downstream scripts read fixed output folders (the un-suffixed folder, or `_v2` where stated above), so they ignore later `_vN` re-runs.
  - `00_utilidades.R` adds a fixed Windows user library to `.libPaths()` if it exists.
  - Some README/LEEME texts written by the scripts contain fixed figures.

## Licence

Code: [MIT](LICENSE). `vendor/core_periphery_sbm` keeps its Apache-2.0 licence. The data are not covered by this licence.

## Citation

Herrera, E., Suárez Coalla, M. P., Ahedo, V. & Galán, J. M. ([year]). Network-based analysis of semantic, phonological, and action fluency in Parkinson's disease on/off medication.

Gallagher, R. J., Young, J.-G., & Welles, B. F. (2021). A clarified typology of core–periphery structure in networks. *Science Advances*, 7(12), eabc9800. https://doi.org/10.1126/sciadv.abc9800

## Authors

Elena Herrera (Universidad de Oviedo, corresponding author), José Manuel Galán (Universidad de Burgos), Virginia Ahedo (Universidad de Burgos) and María Paz Suárez Coalla (Universidad de Oviedo).
