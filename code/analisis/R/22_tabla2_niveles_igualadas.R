# Propósito: comprobación E1 de la Tabla 2 fiel a la R2 (D-024, sensibilidad 6): repetir el modelo de niveles en
#   las 10 extracciones igualadas HC/ON/OFF de 07 (mismo número de listas y tokens por grupo), con la partición
#   propia de cada red igualada (SBM, script 19 con particiones de igualadas) y medias de los tokens retenidos.
# Entradas: output/03_orden_sesiones/tokens.rds; output/07_redes_HC_bootstrap[_v2]/plan_personas.rds y
#   output/18_exportar_redes_SBM/nodes.csv (verificación);
#   output/19_SBM_redes_limpias_igualadas/{particiones_consenso,resumen_redes}.csv; output/21_tabla2_niveles_R2/resultados.csv.
# Salidas: output/22_tabla2_niveles_igualadas[_vN]/ (resultados por extracción, resumen, tokens retenidos, LEEME).
# E1/E2/E5. Descriptivo: no es un contraste entre extracciones. Paquetes: 00_utilidades.R, 00_redes_metricas.R, nlme.
# Ejecutar desde revision_R3/: Rscript --vanilla analisis/R/22_tabla2_niveles_igualadas.R
# --- English ---
# Purpose: E1 check of the script 21 model on the 10 size-matched draws of script 07, each with its own
#   SBM partition. Descriptive, not a test across draws (E1, E2, E5).
# Inputs: 03, 07, 18, 19 (size-matched partitions), 21. Outputs: analisis/output/22_tabla2_niveles_igualadas[_vN]/.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
source("analisis/R/00_sensibilidad_datos.R", encoding = "UTF-8")
source("analisis/R/00_redes_metricas.R", encoding = "UTF-8")
source("analisis/R/00_niveles_tabla2.R", encoding = "UTF-8")
D <- 10L
entrada <- "analisis/output/03_orden_sesiones/tokens.rds"
d07 <- if (dir.exists("analisis/output/07_redes_HC_bootstrap_v2")) "analisis/output/07_redes_HC_bootstrap_v2" else "analisis/output/07_redes_HC_bootstrap"
# Particiones de las igualadas: salida complementaria de 19 (ejecución con --solo igualada, 22-10-2026) o, en una
# ejecución completa con la versión actual de 19, su propio particiones_consenso.csv (que ya incluye las igualadas).
d19i <- if (dir.exists("analisis/output/19_SBM_redes_limpias_igualadas")) "analisis/output/19_SBM_redes_limpias_igualadas" else "analisis/output/19_SBM_redes_limpias"
f21 <- "analisis/output/21_tabla2_niveles_R2/resultados.csv"
base <- readRDS(entrada)
redes <- c(undirected = "cooccurrence", directed = "sequence")
# ---- Reproducir las selecciones de tokens de 07 (misma secuencia aleatoria que 07 y 18) --------------------
set.seed(20261022)
B <- 1000L
pd_global <- base |> filter(cohort == "PD") |> distinct(participant, order) |> arrange(participant)
hc_global <- sort(unique(base$participant[base$cohort == "HC"]))
plan_pd_boot <- t(replicate(B, unlist(lapply(c("OFF_ON", "ON_OFF"), function(s) {
  q <- which(pd_global$order == s); sample(q, length(q), replace = TRUE)
}))))
plan_hc_boot <- t(replicate(B, sample.int(length(hc_global), 21L, replace = TRUE)))
plan_hc_cond <- t(replicate(B, sample.int(length(hc_global), 21L, replace = FALSE)))
plan07 <- readRDS(file.path(d07, "plan_personas.rds"))
stopifnot(identical(plan_hc_cond, plan07$HC_conditional), identical(plan_pd_boot, plan07$PD_boot))
datos <- escenario_datos(base, "recorded")
seleccion <- list(); si <- 0L
for (tarea in c("F", "actions", "animals")) {
  dd <- datos |> filter(task == tarea) |> arrange(participant, condition, position, subposition)
  ll <- dd |> group_by(participant, condition) |> summarise(rows = list(token_id), n = n(), .groups = "drop")
  off <- ll |> filter(condition == "OFF") |> arrange(participant); on <- ll |> filter(condition == "ON") |> arrange(participant)
  hc <- ll |> filter(condition == "HC") |> arrange(participant)
  stopifnot(identical(off$participant, pd_global$participant), identical(hc$participant, hc_global))
  for (fase in c("conditional", "participant_bootstrap")) for (b in seq_len(B)) {
    id_pd <- if (fase == "participant_bootstrap") plan_pd_boot[b, ] else seq_len(nrow(off))
    id_hc <- if (fase == "participant_bootstrap") plan_hc_boot[b, ] else plan_hc_cond[b, ]
    o <- off$rows[id_pd]; n_ <- on$rows[id_pd]; h <- hc$rows[id_hc]
    rp <- order(pmin(lengths(o), lengths(n_)), id_pd); rh <- order(lengths(h), id_hc)
    o <- o[rp]; n_ <- n_[rp]; h <- h[rh]
    k <- pmin(lengths(o), lengths(n_), lengths(h))
    io <- seleccionar_posiciones(o, k); ion <- seleccionar_posiciones(n_, k); ih <- seleccionar_posiciones(h, k)
    if (fase != "conditional" || b > D) next
    si <- si + 1L
    seleccion[[si]] <- tibble(task = tarea, draw = b,
      token_id = c(unlist(Map(`[`, o, io)), unlist(Map(`[`, n_, ion)), unlist(Map(`[`, h, ih))))
  }
}
seleccion <- bind_rows(seleccion)
# Verificación: por extracción y grupo, mismo número de tokens en HC, OFF y ON (igualación de 07).
chk <- seleccion |> left_join(base |> select(token_id, condition), by = "token_id") |> count(task, draw, condition) |>
  pivot_wider(names_from = condition, values_from = n)
stopifnot(nrow(chk) == 3L * D, all(chk$HC == chk$OFF), all(chk$ON == chk$OFF))
# Verificación: las palabras retenidas son exactamente los nodos de las redes igualadas exportadas por 18.
n18 <- read_csv("analisis/output/18_exportar_redes_SBM/nodes.csv", show_col_types = FALSE, trim_ws = FALSE) |>
  filter(grepl("_undirected_igualada", graph)) |>
  mutate(task = sub("_.*", "", graph), condition = sub("^[^_]+_([^_]+)_.*", "\\1", graph),
    draw = as.integer(sub(".*igualada", "", graph)))
mias <- seleccion |> left_join(base |> select(token_id, condition, word_norm), by = "token_id") |>
  distinct(task, condition, draw, word_norm)
stopifnot(nrow(anti_join(mias, n18, by = c("task", "condition", "draw", "word_norm"))) == 0L,
  nrow(anti_join(n18, mias, by = c("task", "condition", "draw", "word_norm"))) == 0L)
# ---- Niveles con la partición propia de cada red igualada ---------------------------------------------------
part <- read_csv(file.path(d19i, "particiones_consenso.csv"), show_col_types = FALSE, trim_ws = FALSE) |> filter(family == "igualada")
res <- read_csv(file.path(d19i, "resumen_redes.csv"), show_col_types = FALSE) |> filter(family == "igualada")
soporte <- res |> group_by(task, network, draw) |> summarise(todos_con_soporte = all(soporte_core), .groups = "drop") |>
  mutate(network = redes[network])
filas <- list(); fi <- 0L
for (b in seq_len(D)) for (red in names(redes)) {
  p <- part |> filter(draw == b, network == red) |> select(task, group, word_norm, hs_consenso)
  d <- datos |> semi_join(seleccion |> filter(draw == b), by = c("task", "token_id")) |>
    left_join(p, by = c("task", "condition" = "group", "word_norm"), relationship = "many-to-one")
  stopifnot(!anyNA(d$hs_consenso))     # toda palabra retenida está en su red igualada
  nv <- d |> pivot_longer(all_of(variables_psico), names_to = "variable", values_to = "value") |>
    filter(variable != "motor" | task == "actions", !is.na(value)) |>
    mutate(partition = if_else(hs_consenso == 0, "core", "periphery")) |>
    group_by(participant, cohort, condition, session, order, task, variable, partition) |>
    summarise(n_tokens = n(), mean = mean(value), .groups = "drop")
  fi <- fi + 1L
  filas[[fi]] <- nv |> group_by(task, variable) |> group_modify(~ ajustar_uno(.x)) |> ungroup() |>
    mutate(network = redes[[red]], draw = b)
}
por_extraccion <- bind_rows(filas) |> left_join(soporte, by = c("task", "network", "draw")) |>
  mutate(no_identificable = !todos_con_soporte | random_effects == "no estimable")
principal <- read_csv(f21, show_col_types = FALSE) |>
  select(task, network, variable, term, estimate_main = estimate, p_BH_main = p_BH, no_identificable_main = no_identificable)
resumen <- por_extraccion |> filter(term != "Intercept_HC_core") |> group_by(task, network, variable, term) |>
  summarise(draws_estimables = sum(!no_identificable),
    mean_estimate = if (any(!no_identificable)) mean(estimate[!no_identificable]) else NA_real_,
    min_estimate = if (any(!no_identificable)) min(estimate[!no_identificable]) else NA_real_,
    max_estimate = if (any(!no_identificable)) max(estimate[!no_identificable]) else NA_real_,
    draws_p05 = sum(!no_identificable & p < .05), .groups = "drop") |>
  left_join(principal, by = c("task", "network", "variable", "term")) |>
  mutate(same_sign_mean = sign(mean_estimate) == sign(estimate_main),
    main_sig_BH = !no_identificable_main & !is.na(p_BH_main) & p_BH_main < .05)
out <- nuevo_directorio("22_tabla2_niveles_igualadas")
guardar_csv(por_extraccion, file.path(out, "resultados_por_extraccion.csv"))
guardar_csv(resumen, file.path(out, "resumen.csv"))
guardar_csv(seleccion, file.path(out, "tokens_retenidos.csv"))
registrar_entorno(out, c(entrada, file.path(d07, "plan_personas.rds"), file.path(d19i, "particiones_consenso.csv"),
  file.path(d19i, "resumen_redes.csv"), f21, "analisis/output/18_exportar_redes_SBM/nodes.csv", "analisis/R/22_tabla2_niveles_igualadas.R", "analisis/R/00_niveles_tabla2.R",
  "analisis/R/00_redes_metricas.R", "analisis/R/00_sensibilidad_datos.R", "analisis/R/00_utilidades.R"))
writeLines(c("# Tabla 2 (niveles, D-024) en las redes igualadas — comprobación E1",
  sprintf("Extracciones condicionales 1–%d de 07: 21 listas por grupo y el mismo número de tokens en HC, OFF y ON.", D),
  "Selecciones reproducidas con la misma secuencia aleatoria que 07/18 (planes verificados) y tokens iguales por grupo.",
  "Partición propia de cada red igualada (19 con particiones de igualadas, 3 cadenas). Medias con los tokens retenidos.",
  "Mismo modelo que 21 (00_niveles_tabla2.R). No identificable si alguna partición de la extracción no tiene soporte.",
  "resumen.csv: estimaciones medias/rango y nº de extracciones con p < 0,05 sin ajustar, frente a la estimación principal (21).",
  "Descriptivo: 10 extracciones condicionadas a estas personas; no es un intervalo ni un contraste."), file.path(out, "LEEME.md"), useBytes = TRUE)
print(resumen |> filter(variable != "AoA") |> group_by(term) |>
  summarise(modelos = n(), estimables = sum(draws_estimables > 0), principal_sig = sum(main_sig_BH, na.rm = TRUE),
    mismo_signo = sum(same_sign_mean, na.rm = TRUE), mayoria_p05 = sum(draws_p05 >= 6), .groups = "drop"))
cat("Salida:", out, "\n")
