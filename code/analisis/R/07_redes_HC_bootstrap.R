# Propósito: igualación HC/ON/OFF, bootstrap de personas pareado y sensibilidad de adyacencias.
# Entradas: output/03_orden_sesiones/tokens.rds; escenarios hipotéticos P3/P11.
# Salidas: output/07_redes_HC_bootstrap[/_vN]/, métricas, presupuestos, sesgos e intervalos exploratorios.
# E1/E2/E4: no SBM ni p causal; bootstrap estratificado por secuencia ON/OFF.
# Paquetes: 00_utilidades.R, igraph; funciones de 00_redes_metricas.R.
# --- English ---
# Purpose: size-matched HC/ON/OFF networks, person bootstrap stratified by ON/OFF sequence (B = 1000) and
#   edge-definition sensitivity; intervals are exploratory (E1, E2, E4).
# Inputs: 03 tokens.rds. Outputs: analisis/output/07_redes_HC_bootstrap[_vN]/ (metrics, plans, simulations).
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
source("analisis/R/00_sensibilidad_datos.R", encoding = "UTF-8")
source("analisis/R/00_redes_metricas.R", encoding = "UTF-8")
suppressPackageStartupMessages(library(igraph))
set.seed(20261022)
B <- 1000L
entrada <- "analisis/output/03_orden_sesiones/tokens.rds"
base <- readRDS(entrada)
stopifnot(nrow(base) == 2965L)
pd_global <- base |> filter(cohort == "PD") |> distinct(participant, order) |> arrange(participant)
hc_global <- sort(unique(base$participant[base$cohort == "HC"]))
# Plan global: misma persona remuestreada en todas sus condiciones Y tareas.
plan_pd_boot <- t(replicate(B, unlist(lapply(c("OFF_ON", "ON_OFF"), function(s) {
  q <- which(pd_global$order == s); sample(q, length(q), replace = TRUE)
}))))
plan_hc_boot <- t(replicate(B, sample.int(length(hc_global), 21L, replace = TRUE)))
plan_hc_cond <- t(replicate(B, sample.int(length(hc_global), 21L, replace = FALSE)))
modos <- c("cooccurrence", "sequence_rebuilt", "sequence_original")
metricas <- c("nodes", "edges", "density", "transitivity", "mean_degree_total", "largest_weak_component_pct")
# Misma selección, pero distintas adyacencias; en el segundo caso quedan los nodos 1 y 3 (sin aristas).
test <- list(c(1L, 2L, 3L)); ii <- list(c(1L, 3L))
stopifnot(red_metricas(test, ii, 3L, "sequence_rebuilt")["edges"] == 1,
  red_metricas(test, ii, 3L, "sequence_original")["edges"] == 0,
  red_metricas(test, ii, 3L, "sequence_original")["nodes"] == 2)
preparar <- function(d, tarea) {
  d <- d |> filter(task == tarea) |> arrange(participant, condition, position, subposition)
  vocab <- unique(d$word_norm); d$wid <- match(d$word_norm, vocab)
  ll <- d |> group_by(participant, cohort, condition, order) |>
    summarise(words = list(wid), .groups = "drop")
  off <- ll |> filter(condition == "OFF") |> arrange(participant)
  on <- ll |> filter(condition == "ON") |> arrange(participant)
  hc <- ll |> filter(condition == "HC") |> arrange(participant)
  stopifnot(identical(off$participant, on$participant), nrow(off) == 21L, nrow(hc) == 24L)
  list(OFF = off$words, ON = on$words, HC = hc$words, PD_id = off$participant, HC_id = hc$participant,
    order = off$order, V = length(vocab))
}
extraer <- function(p, boot, b) {
  stopifnot(identical(p$PD_id, pd_global$participant), identical(p$HC_id, hc_global))
  id_pd <- if (boot) plan_pd_boot[b, ] else seq_along(p$OFF)
  id_hc <- if (boot) plan_hc_boot[b, ] else plan_hc_cond[b, ]
  off <- p$OFF[id_pd]; on <- p$ON[id_pd]; hc <- p$HC[id_hc]
  # Emparejamiento por rango de longitud: misma distribución de k en los tres grupos.
  rp <- order(pmin(lengths(off), lengths(on)), id_pd)
  rh <- order(lengths(hc), id_hc)
  off <- off[rp]; on <- on[rp]; hc <- hc[rh]
  k <- pmin(lengths(off), lengths(on), lengths(hc))
  ioff <- seleccionar_posiciones(off, k); ion <- seleccionar_posiciones(on, k); ihc <- seleccionar_posiciones(hc, k)
  stopifnot(sum(lengths(ioff)) == sum(lengths(ion)), sum(lengths(ioff)) == sum(lengths(ihc)))
  sets <- list(OFF = list(w = off, ix = ioff), ON = list(w = on, ix = ion), HC = list(w = hc, ix = ihc))
  values <- lapply(modos, function(m) do.call(rbind, lapply(sets, function(s) red_metricas(s$w, s$ix, p$V, m))))
  names(values) <- modos
  list(values = values, budget = c(tokens_each = sum(k), unique_PD = length(unique(id_pd)), unique_HC = length(unique(id_hc)),
    expected_original_edges_OFF = sum(k * (k - 1) / lengths(off)),
    expected_original_edges_ON = sum(k * (k - 1) / lengths(on)),
    expected_original_edges_HC = sum(k * (k - 1) / lengths(hc))))
}
resultados <- list(); presupuestos <- list(); cursor <- 0L; pc <- 0L
for (scenario in c("recorded", "deduplicated", "ambiguous_drop", "ambiguous_split")) {
  datos <- escenario_datos(base, scenario)
  for (tarea in c("F", "actions", "animals")) {
    # Los escenarios de la celda ambigua solo cambian acciones; no repetir las otras tareas.
    if (scenario %in% c("ambiguous_drop", "ambiguous_split") && tarea != "actions") next
    p <- preparar(datos, tarea)
    fases <- if (scenario == "recorded") c("conditional", "participant_bootstrap") else "conditional"
    for (fase in fases) {
      arr <- array(NA_real_, dim = c(B, 3L, 3L, 6L),
        dimnames = list(draw = as.character(seq_len(B)), mode = modos, condition = c("OFF", "ON", "HC"), metric = metricas))
      budget <- matrix(NA_real_, nrow = B, ncol = 6L)
      for (b in seq_len(B)) {
        x <- extraer(p, fase == "participant_bootstrap", b)
        for (j in seq_along(modos)) arr[b, j, , ] <- x$values[[j]]
        budget[b, ] <- x$budget
      }
      colnames(budget) <- names(x$budget)
      cursor <- cursor + 1L
      resultados[[cursor]] <- as.data.frame.table(arr, responseName = "value", stringsAsFactors = FALSE) |>
        as_tibble() |> mutate(draw = as.integer(draw), scenario = scenario, task = tarea, phase = fase)
      pc <- pc + 1L
      presupuestos[[pc]] <- as_tibble(budget) |> mutate(draw = seq_len(B), scenario = scenario, task = tarea, phase = fase)
      cat("Completado:", scenario, tarea, fase, B, "réplicas\n")
      flush.console()
    }
  }
}
resultados <- bind_rows(resultados); presupuestos <- bind_rows(presupuestos)
stopifnot(nrow(presupuestos) == 11000L, all(presupuestos$tokens_each > 0),
  all(presupuestos$unique_PD[presupuestos$phase == "conditional"] == 21L),
  all(presupuestos$unique_HC[presupuestos$phase == "conditional"] == 21L))
ancho <- resultados |> pivot_wider(names_from = condition, values_from = value)
diferencias <- bind_rows(
  ancho |> mutate(contrast = "ON_minus_OFF", delta = ON - OFF),
  ancho |> mutate(contrast = "HC_minus_OFF", delta = HC - OFF),
  ancho |> mutate(contrast = "HC_minus_ON", delta = HC - ON))
resumen <- diferencias |> group_by(scenario, task, phase, mode, metric, contrast) |>
  summarise(B = n(), mean = mean(delta), sd = sd(delta), q025 = unname(quantile(delta, .025)),
    q975 = unname(quantile(delta, .975)), .groups = "drop")
sesgo <- resumen |> filter(scenario == "recorded") |>
  select(-B) |> pivot_wider(names_from = phase, values_from = c(mean, sd, q025, q975)) |>
  mutate(bootstrap_bias = mean_participant_bootstrap - mean_conditional,
    point_outside_boot_percentiles = mean_conditional < q025_participant_bootstrap | mean_conditional > q975_participant_bootstrap,
    normal_SE_low_approx = mean_conditional - 1.96 * sd_participant_bootstrap,
    normal_SE_high_approx = mean_conditional + 1.96 * sd_participant_bootstrap)
sesgo_grupos <- resultados |> filter(scenario == "recorded") |>
  group_by(task, mode, condition, metric, phase) |> summarise(mean = mean(value), .groups = "drop") |>
  pivot_wider(names_from = phase, values_from = mean) |>
  mutate(bootstrap_bias = participant_bootstrap - conditional)
out <- nuevo_directorio("07_redes_HC_bootstrap")
for (nm in c("presupuestos", "resumen", "sesgo", "sesgo_grupos")) guardar_csv(get(nm), file.path(out, paste0(nm, ".csv")))
saveRDS(list(metrics = resultados, differences = diferencias), file.path(out, "simulaciones.rds"), compress = TRUE)
saveRDS(list(PD = pd_global, HC = hc_global, PD_boot = plan_pd_boot, HC_boot = plan_hc_boot,
  HC_conditional = plan_hc_cond), file.path(out, "plan_personas.rds"))
registrar_entorno(out, c(entrada, "analisis/R/07_redes_HC_bootstrap.R", "analisis/R/00_redes_metricas.R",
  "analisis/R/00_sensibilidad_datos.R", "analisis/R/00_utilidades.R",
  list.files("../parkinson/data", pattern = "\\.xlsx$", full.names = TRUE)))
writeLines(c("# Igualación de tres grupos y bootstrap — 22-10-2026",
  "21 personas por grupo y misma distribución de k/tokens en HC/ON/OFF; parejas PD intactas.",
  "1.000 extracciones condicionales y 1.000 bootstrap de personas por tarea, con estratos de secuencia PD 11/10.",
  "Bootstrap conserva ambas condiciones y todas las palabras de cada paciente elegido; HC independiente.",
  "Cada réplica incorpora una extracción de tokens: incertidumbre combinada, no bootstrap aislado de pacientes.",
  "Tres representaciones: co-ocurrencia, secuencia reconstruida, solo adyacencias originales. Esta última no iguala el número esperado de transiciones: ver presupuestos.",
  "Escenarios clínicos hipotéticos: retirar repeticiones y excluir/separar celda ambigua; no cambian el dato principal.",
  "Percentiles bootstrap y normal_SE son exploratorios, no IC calibrados ni p confirmatorios; revisar sesgo_grupos/sesgo.",
  "El remuestreo con reemplazo reduce personas distintas y soporte de palabras/aristas: no asumir que el bootstrap de grafos es insesgado.",
  "No se ajustan métricas por sesión ni se infiere reorganización causal; el estrato conserva composición, no elimina periodo.",
  "Fuentes metodológicas de cautela: https://pmc.ncbi.nlm.nih.gov/articles/PMC3958753/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC5809547/ .",
  "Sin SBM ni cambios en manuscrito."), file.path(out, "LEEME.md"), useBytes = TRUE)
print(presupuestos |> group_by(phase, task) |> summarise(mean_unique_PD = mean(unique_PD), mean_unique_HC = mean(unique_HC), .groups = "drop"))
cat("Salida:", out, "\n")
