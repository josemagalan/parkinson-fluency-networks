# Propósito: contrastar tablas anteriores aportadas por Josema con la réplica 17.
# Entradas: historial/20261022_tablas_core_periphery_previas/transcripcion.csv,
# output 17 (redes, resumen_hubspoke, mdl_modelos) y texto R2; fuentes preservadas.
# Salidas: output/31_contrastar_tablas_historicas_SBM[_vN]/; no ejecuta el SBM.
# E5/E6, P5/P15: procedencia, modelos ganadores, representación y recuentos core.
# Paquetes: 00_utilidades.R (readxl/readr/dplyr/tidyr/purrr/broom).
# --- English ---
# Purpose (optional): compare transcribed historical SBM description-length tables with the script 17
#   replication (E5, E6; P5/P15). No SBM is run.
# Inputs: historial/ transcription (not public), 17 outputs. Outputs: analisis/output/31_contrastar_tablas_historicas_SBM[_vN]/.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
set.seed(20261022)
base <- "historial/20261022_tablas_core_periphery_previas"
f <- file.path(base, "transcripcion.csv")
f17 <- "analisis/output/17_replica_SBM_historico"
h <- read_csv(f, show_col_types = FALSE)
redes <- read_csv(file.path(f17, "redes.csv"), show_col_types = FALSE)
part <- read_csv(file.path(f17, "resumen_hubspoke.csv"), show_col_types = FALSE)
mdl <- read_csv(file.path(f17, "mdl_modelos.csv"), show_col_types = FALSE)
stopifnot(nrow(h) == 18L, !anyDuplicated(h$red),
          all(h$HS_0 + h$HS_1 == h$N),
          all(h$Ly_0 + h$Ly_1 + h$Ly_2 == h$N))
no_dir <- redes |> filter(variante == "no_dirigida") |>
  select(red, N_graph = nodos, L_no_dirigida = aristas, bucles, core_CSV_R2 = core_R2)
dir <- redes |> filter(variante == "dirigida_tal_cual") |>
  select(red, L_dirigida = aristas)
best <- mdl |> filter(variante == "no_dirigida") |>
  mutate(familia = ifelse(modelo == "hubspoke", "HS", "Ly")) |>
  group_by(red, familia) |> summarise(mdl_min = min(mdl), .groups = "drop") |>
  pivot_wider(names_from = familia, values_from = mdl_min, names_prefix = "replica_")
res <- h |> left_join(no_dir, by = "red", relationship = "one-to-one") |>
  left_join(dir, by = "red", relationship = "one-to-one") |>
  left_join(best, by = "red", relationship = "one-to-one") |>
  left_join(part |> filter(variante == "no_dirigida") |>
    select(red, mdl_CSV_R2 = mdl_R2, conc_consenso_R2), by = "red", relationship = "one-to-one") |>
  mutate(delta_historica_Ly_menos_HS = Ly_fit - HS_fit,
    ganador_historico = ifelse(delta_historica_Ly_menos_HS > 0, "HS", "Ly"),
    delta_replica_Ly_menos_HS = replica_Ly - replica_HS,
    ganador_replica_numerico = ifelse(delta_replica_Ly_menos_HS > 0, "HS", "Ly"),
    empate_replica_aproximado = abs(delta_replica_Ly_menos_HS) < 1,
    ganador_R2_texto = ifelse(red %in% c("graphAccionesoff", "graphFpalabraoff", "graphFpalabraon"), "Ly", "HS"),
    coincide_ganador_replica = ganador_historico == ganador_replica_numerico,
    coincide_texto_R2 = ganador_historico == ganador_R2_texto,
    coincide_N = N == N_graph,
    coincide_L_no_dirigida = L == L_no_dirigida,
    coincide_L_dirigida = ifelse(is.na(L_dirigida), NA, L == L_dirigida),
    diferencia_core_tabla_menos_CSV = HS_0 - core_CSV_R2,
    coincide_recuento_core_CSV = diferencia_core_tabla_menos_CSV == 0)
stopifnot(all(res$coincide_N), all(res$coincide_L_no_dirigida),
          !anyNA(res$replica_HS), !anyNA(res$replica_Ly))
out <- nuevo_directorio("31_contrastar_tablas_historicas_SBM")
guardar_csv(res, file.path(out, "comparacion_18_redes.csv"))
guardar_csv(res |> filter(!coincide_texto_R2 | !coincide_ganador_replica | !coincide_recuento_core_CSV),
  file.path(out, "discrepancias.csv"))
resumen <- res |> group_by(tipo) |> summarise(n = n(), HS_gana = sum(ganador_historico == "HS"),
  Ly_gana = sum(ganador_historico == "Ly"), N_coinciden = sum(coincide_N),
  L_no_dirigida_coinciden = sum(coincide_L_no_dirigida),
  recuentos_core_CSV_distintos = sum(!coincide_recuento_core_CSV),
  ganadores_replica_coinciden = sum(coincide_ganador_replica),
  ganadores_texto_R2_coinciden = sum(coincide_texto_R2), .groups = "drop")
guardar_csv(resumen, file.path(out, "resumen.csv"))
registrar_entorno(out, c(f, file.path(base, c("tabla1_secuencia_usuario.png", "tabla2_coocurrencia_usuario.png", "texto_aportado.md")),
  file.path(f17, c("redes.csv", "resumen_hubspoke.csv", "mdl_modelos.csv")),
  "manuscrito/R2_texto.md", "analisis/R/31_contrastar_tablas_historicas_SBM.R"))
writeLines(c(
  "# Contraste de tablas históricas aportadas por Josema (versiones anteriores)",
  "Las capturas aportadas no se identifican automáticamente con S1/S2 de la R2 enviada.",
  "18 filas transcritas; sumas de bloques, N y L validados frente a output17.",
  "Delta = longitud de descripción Ly menos HS: positivo favorece HS.",
  "Históricas: HS gana9/9 secuencia y7/9 co-ocurrencia; Ly gana solo F_ON/F_OFF co-ocurrencia.",
  "La frase R2 que añade Acciones_OFF co-ocurrencia contradice estos valores.",
  "L de las9 redes de secuencia coincide con grafos convertidos a no dirigidos;",
  "no coincide con el número de aristas de los grafos dirigidos guardados.",
  "Los recuentos de core difieren de CSV_R2 en13/18 redes: documentan salidas distintas,",
  "no identidad de particiones ni una causa demostrada de la discrepancia.",
  "Ganador numérico réplica17 coincide16/18: F_secuencia_OFF cambia; Animales_secuencia_ON",
  "es casi empate en17. Mínimos de réplicas vs valor histórico único no son réplica exacta.",
  "mdl_CSV_R2 es una reevaluación en17 de la partición CSV, no la MDL histórica de la captura.",
  "No se ha vuelto a ejecutar el SBM; no se cambian outputs17 ni particiones.",
  "Los bloques no vacíos Ly de la captura no revelan capas permitidas/configuración original.",
  "Comparar dos modelos core-periphery no constituye por sí solo contraste frente a ausencia de estructura."
), file.path(out, "LEEME.md"), useBytes = TRUE)
print(resumen)
print(res |> filter(!coincide_ganador_replica | !coincide_texto_R2) |>
  select(red, HS_fit, Ly_fit, delta_historica_Ly_menos_HS, delta_replica_Ly_menos_HS))
cat("Recuentos core distintos:", sum(!res$coincide_recuento_core_CSV), "de18\nSalida:", out, "\n")
