# Propósito: redes normalizadas y rarefacción pareada ON/OFF al mínimo individual.
# Entradas: output/03_orden_sesiones/tokens.rds; sin diccionarios ni particiones.
# Salidas: output/05_redes_igualadas_preparacion[/_vN]/, GraphML, métricas, RDS y tablas.
# E1/E2: igualar n de participantes/tokens; no SBM, bootstrap poblacional ni p de redes.
# Paquetes: 00_utilidades.R + igraph. Semilla y B declarados antes de calcular métricas.
# Ejecutar desde revision_R3/: Rscript --vanilla analisis/R/05_redes_igualadas_preparacion.R
# --- English ---
# Purpose: 18 normalised networks (3 groups x 3 tasks x directed/undirected) and paired ON/OFF
#   rarefaction to each patient's minimum (B = 1000) (E1, E2; no SBM, no inference).
# Inputs: 03 tokens.rds. Outputs: analisis/output/05_redes_igualadas_preparacion[_vN]/ (metrics, GraphML, RDS).
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
suppressPackageStartupMessages(library(igraph))
set.seed(20261022)
B <- 1000L
entrada <- "analisis/output/03_orden_sesiones/tokens.rds"
tokens <- readRDS(entrada) |> arrange(participant, task, condition, position)
stopifnot(nrow(tokens) == 2965L, !anyDuplicated(tokens$token_id))
listas <- tokens |> group_by(participant, cohort, condition, task, session, order) |>
  summarise(words = list(word_norm), positions = list(position), n_words = n(), .groups = "drop")
stopifnot(nrow(listas) == 198L)
repeticiones <- tokens |> count(participant, task, condition, word_norm) |> filter(n > 1)
presupuesto <- listas |> filter(cohort == "PD") |>
  select(participant, task, order, condition, n_words) |>
  pivot_wider(names_from = condition, values_from = n_words) |>
  mutate(k = pmin(ON, OFF), retained_ON = k / ON, retained_OFF = k / OFF)
stopifnot(nrow(presupuesto) == 63L, !anyNA(presupuesto), all(presupuesto$k >= 2))
resumen_presupuesto <- presupuesto |> group_by(task) |>
  summarise(n_patients = n(), tokens_ON_full = sum(ON), tokens_OFF_full = sum(OFF),
    tokens_each_equalized = sum(k), min_k = min(k), max_k = max(k),
    retained_ON_pct = 100 * sum(k) / sum(ON), retained_OFF_pct = 100 * sum(k) / sum(OFF), .groups = "drop")
construir <- function(words_list, directed) {
  vertices <- unique(unlist(words_list, use.names = FALSE))
  edges <- lapply(words_list, function(w) {
    if (length(w) < 2L) return(matrix(character(), ncol = 2))
    if (directed) cbind(head(w, -1), tail(w, -1)) else t(combn(w, 2))
  })
  ed <- do.call(rbind, edges)
  g <- graph_from_data_frame(as.data.frame(ed, stringsAsFactors = FALSE), directed = directed,
    vertices = data.frame(name = vertices))
  E(g)$weight <- 1
  simplify(g, remove.multiple = TRUE, remove.loops = FALSE, edge.attr.comb = list(weight = "sum"))
}
medir <- function(g) {
  loops <- which_loop(g)
  n_loops <- sum(loops)
  loop_weight <- sum(E(g)$weight[loops])
  s <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list(weight = "sum"))
  n <- vcount(s); e <- ecount(s)
  # Directed: clustering en proyección no dirigida; no confundir con tríadas dirigidas.
  cl <- transitivity(as_undirected(s, mode = "collapse"), type = "global", isolates = "zero")
  tibble(nodes = n, edges = e, density = edge_density(s, loops = FALSE),
    transitivity_undirected = cl, mean_degree_total = mean(degree(s, mode = "all")),
    largest_weak_component_pct = 100 * max(components(s, mode = "weak")$csize) / n,
    loops_excluded = n_loops, loop_weight_excluded = loop_weight)
}
# Validaciones con grafos pequeños conocidos: pesos y transición tras seleccionar posiciones.
test_d <- construir(list(c("a", "b", "a")), TRUE)
test_u <- construir(list(c("a", "b", "a")), FALSE)
stopifnot(ecount(test_d) == 2L, all(E(test_d)$weight == 1),
  sum(which_loop(test_u)) == 1L, E(test_u)$weight[get_edge_ids(test_u, c("a", "b"))] == 2,
  medir(test_u)$edges == 1L, medir(test_u)$nodes == 2L)
full_graphs <- list(); full_metrics <- list()
for (tarea in c("F", "actions", "animals")) {
  for (cond in c("HC", "OFF", "ON")) {
    ws <- listas$words[listas$task == tarea & listas$condition == cond]
    for (dirigida in c(FALSE, TRUE)) {
      net <- if (dirigida) "directed" else "undirected"
      clave <- paste(tarea, cond, net, sep = "_")
      full_graphs[[clave]] <- construir(ws, dirigida)
      full_metrics[[clave]] <- medir(full_graphs[[clave]]) |> mutate(task = tarea, condition = cond, network = net,
        n_participants = length(ws), n_tokens = sum(lengths(ws)))
    }
  }
}
full_metrics <- bind_rows(full_metrics)
simulaciones <- vector("list", 3L * B * 4L); cursor <- 0L
ejemplos <- list(); configuracion <- list()
for (tarea in c("F", "actions", "animals")) {
  p <- presupuesto |> filter(task == tarea) |> arrange(participant)
  off <- listas |> filter(task == tarea, cohort == "PD", condition == "OFF") |> arrange(participant)
  on <- listas |> filter(task == tarea, cohort == "PD", condition == "ON") |> arrange(participant)
  stopifnot(identical(p$participant, off$participant), identical(p$participant, on$participant))
  configuracion[[tarea]] <- list(participants = p$participant, k = p$k, OFF = off$words, ON = on$words)
  for (b in seq_len(B)) {
    # Mantener el orden de producción tras seleccionar posiciones, incluso si la palabra se repite.
    idx_off <- lapply(seq_len(nrow(p)), function(i) sort(sample.int(length(off$words[[i]]), p$k[i])))
    idx_on <- lapply(seq_len(nrow(p)), function(i) sort(sample.int(length(on$words[[i]]), p$k[i])))
    w_off <- lapply(seq_len(nrow(p)), function(i) off$words[[i]][idx_off[[i]]])
    w_on <- lapply(seq_len(nrow(p)), function(i) on$words[[i]][idx_on[[i]]])
    stopifnot(identical(lengths(w_off), as.integer(p$k)), identical(lengths(w_on), as.integer(p$k)))
    for (dirigida in c(FALSE, TRUE)) {
      net <- if (dirigida) "directed" else "undirected"
      for (cond in c("OFF", "ON")) {
        g <- construir(if (cond == "OFF") w_off else w_on, dirigida)
        cursor <- cursor + 1L
        simulaciones[[cursor]] <- medir(g) |> mutate(task = tarea, network = net, condition = cond, draw = b,
          n_participants = nrow(p), n_tokens = sum(p$k))
        if (b == 1L) ejemplos[[paste(tarea, cond, net, sep = "_")]] <- g
      }
    }
  }
  cat("Rarefacción terminada:", tarea, B, "extracciones\n")
}
simulaciones <- bind_rows(simulaciones)
stopifnot(nrow(simulaciones) == 12000L)
metricas <- c("nodes", "edges", "density", "transitivity_undirected", "mean_degree_total", "largest_weak_component_pct")
deltas_mc <- simulaciones |> select(task, network, condition, draw, all_of(metricas)) |>
  pivot_longer(all_of(metricas), names_to = "metric", values_to = "value") |>
  pivot_wider(names_from = condition, values_from = value) |> mutate(delta_ON_OFF = ON - OFF)
resumen_mc <- deltas_mc |> group_by(task, network, metric) |>
  summarise(B = n(), mean_delta = mean(delta_ON_OFF), median_delta = median(delta_ON_OFF),
    q025_MC = unname(quantile(delta_ON_OFF, .025)), q975_MC = unname(quantile(delta_ON_OFF, .975)),
    .groups = "drop")
deltas_full <- full_metrics |> filter(condition %in% c("ON", "OFF")) |>
  select(task, network, condition, all_of(metricas)) |> pivot_longer(all_of(metricas), names_to = "metric", values_to = "value") |>
  pivot_wider(names_from = condition, values_from = value) |> mutate(full_delta_ON_OFF = ON - OFF)
comparacion <- left_join(deltas_full, resumen_mc, by = c("task", "network", "metric"))
out <- nuevo_directorio("05_redes_igualadas_preparacion")
for (nm in c("presupuesto", "resumen_presupuesto", "repeticiones", "full_metrics", "simulaciones", "deltas_mc", "resumen_mc", "comparacion")) {
  guardar_csv(get(nm), file.path(out, paste0(nm, ".csv")))
}
dir.create(file.path(out, "graphml_completos")); dir.create(file.path(out, "graphml_igualados_ejemplo_1"))
for (nm in names(full_graphs)) write_graph(full_graphs[[nm]], file.path(out, "graphml_completos", paste0(nm, ".graphml")), format = "graphml")
for (nm in names(ejemplos)) write_graph(ejemplos[[nm]], file.path(out, "graphml_igualados_ejemplo_1", paste0(nm, ".graphml")), format = "graphml")
saveRDS(list(seed = 20261022, B = B, input = entrada, matched_lists = configuracion), file.path(out, "configuracion_rarefaccion.rds"))
registrar_entorno(out, c(entrada, "analisis/R/05_redes_igualadas_preparacion.R", "analisis/R/00_utilidades.R"))
writeLines(c(
  "# Preparación de redes y rarefacción pareada — 22-10-2026",
  "18 redes completas normalizadas: tres grupos × tres tareas × dos tipos; sin ajustar particiones.",
  "Rarefacción solo en pacientes: 21 pares, k=min(n_ON,n_OFF) por paciente/tarea, 1.000 extracciones, semilla 20261022.",
  "Submuestreo de posiciones sin reemplazo, manteniendo orden original. Los mismos k en ambas condiciones.",
  "La secuencia conecta respuestas retenidas consecutivas: esto puede crear enlaces que no eran adyacentes en la lista completa.",
  "Co-ocurrencia suma todos los pares de posiciones; si una palabra se repite, se cuentan sus multiplicidades. No deduplicar listas sin validación clínica P3.",
  "GraphML conserva bucles/pesos; métricas eliminan bucles y los cuantifican. Clustering de redes dirigidas usa proyección no dirigida.",
  "mean_degree_total cuenta entrada+salida en redes dirigidas; densidad dirigida usa n*(n-1). Pesos no se tratan como distancias.",
  "q025_MC/q975_MC son percentiles condicionales del submuestreo, NO IC de pacientes ni intervalos de efecto poblacional.",
  "No hay bootstrap por paciente, p de permutación ni ajuste de sesión en las métricas; no se infiere reorganización dopaminérgica.",
  "El ejemplo GraphML 1 sirve para verificar intercambio con SBM, no representa la media de 1.000 redes.",
  "Las redes HC son provisionales por partir<TAB>llegar (P11); HC no se incluye en la rarefacción pareada ON/OFF.",
  "Pendientes: fijar sensibilidad clínica, igualación con HC, bootstrap poblacional, inferencia que respete sesión/pareado y ajuste SBM.",
  "Fuentes y manuscrito intactos; este bloque prepara E1 pero no lo cierra."
), file.path(out, "LEEME.md"), useBytes = TRUE)
print(resumen_presupuesto); cat("Salida:", out, "\n")
