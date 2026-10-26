# Propósito: exportar, como listas de aristas binarias, las redes sobre las que se ajusta el SBM R3 (D-021).
#   (a) 18 redes propias completas y normalizadas (tarea × grupo × tipo), verificadas contra el output 14;
#   (b) redes conjuntas por tarea y tipo: HC+ON+OFF («conjunta_3») y ON+OFF («conjunta_PD»), que dan una
#       partición de referencia común a los grupos que se comparan (evita la circularidad de la R2);
#   (c) redes igualadas HC/ON/OFF: las D primeras extracciones condicionales del script 07 (escenario
#       registrado), reproducidas con la misma semilla y verificadas contra simulaciones.rds de 07.
# Entradas: output/03_orden_sesiones/tokens.rds; output/14_preparar_entrada_SBM/edges.csv;
#   output/07_redes_HC_bootstrap[_v2]/{simulaciones,plan_personas}.rds (verificación); funciones de 00_redes_metricas.R.
# Salidas: output/18_exportar_redes_SBM[_vN]/ (index.csv, edges.csv, nodes.csv, verificaciones.csv, LEEME.md).
# E1/E2/E5. No ajusta el SBM (lo hace analisis/python/19_SBM_redes_limpias.py). Paquetes: 00_utilidades.R + igraph.
# Ejecutar desde revision_R3/: Rscript --vanilla analisis/R/18_exportar_redes_SBM.R
# --- English ---
# Purpose: export the R3 SBM networks as binary edge lists (D-021): 18 own-group, 12 joint (HC+ON+OFF,
#   ON+OFF) and 180 size-matched networks re-drawn with script 07's random stream (E1, E2, E5).
# Inputs: 03 tokens.rds, 14 edges.csv, 07 plans/simulations. Outputs: analisis/output/18_exportar_redes_SBM[_vN]/.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
source("analisis/R/00_sensibilidad_datos.R", encoding = "UTF-8")
source("analisis/R/00_redes_metricas.R", encoding = "UTF-8")
suppressPackageStartupMessages(library(igraph))
D <- 10L                      # extracciones igualadas exportadas (de las 1.000 de 07)
entrada <- "analisis/output/03_orden_sesiones/tokens.rds"
d14 <- "analisis/output/14_preparar_entrada_SBM"
d07 <- if (dir.exists("analisis/output/07_redes_HC_bootstrap_v2")) "analisis/output/07_redes_HC_bootstrap_v2" else "analisis/output/07_redes_HC_bootstrap"
base <- readRDS(entrada)
stopifnot(nrow(base) == 2965L)
tareas <- c("F", "actions", "animals")
modos <- c(cooccurrence = "undirected", sequence_rebuilt = "directed")
verif <- list(); vcur <- 0L
anota <- function(nombre, ok, detalle = "") {
  stopifnot(isTRUE(ok)); vcur <<- vcur + 1L
  verif[[vcur]] <<- tibble(check = nombre, passed = ok, detail = detalle)
}
# Mismo grafo que red_metricas() (sin bucles ni multiaristas, inducido en las palabras retenidas), con nombres.
red_grafo <- function(ws, ix, vocab, modo) {
  selected <- unique(unlist(Map(function(w, i) w[i], ws, ix), use.names = FALSE))
  ed <- Map(function(w, i) {
    if (length(i) < 2L) return(matrix(integer(), ncol = 2))
    if (modo == "cooccurrence") return(t(combn(w[i], 2)))
    e <- cbind(w[head(i, -1L)], w[tail(i, -1L)])
    if (modo == "sequence_original") e <- e[diff(i) == 1L, , drop = FALSE]
    e
  }, ws, ix)
  ed <- do.call(rbind, ed)
  g <- make_empty_graph(length(vocab), directed = modo != "cooccurrence")
  V(g)$name <- vocab
  if (nrow(ed)) g <- add_edges(g, as.vector(t(ed)))
  g <- induced_subgraph(g, selected)
  simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
}
filas_index <- list(); filas_edges <- list(); filas_nodes <- list(); gi <- 0L
guardar_red <- function(id, g, familia, tarea, grupo, tipo, draw = NA_integer_, n_listas = NA_integer_, n_tokens = NA_integer_) {
  gi <<- gi + 1L
  el <- as_edgelist(g, names = TRUE)
  filas_index[[gi]] <<- tibble(graph = id, family = familia, task = tarea, group = grupo, network = tipo,
    draw = draw, nodes = vcount(g), edges = ecount(g), n_lists = n_listas, n_tokens = n_tokens)
  filas_edges[[gi]] <<- tibble(graph = id, from = el[, 1], to = el[, 2])
  filas_nodes[[gi]] <<- tibble(graph = id, word_norm = V(g)$name)
}
# ---- (a) y (b): redes completas propias y conjuntas -------------------------------------------------
for (tarea in tareas) {
  d <- base |> filter(task == tarea) |> arrange(participant, condition, position)
  vocab <- unique(d$word_norm); d$wid <- match(d$word_norm, vocab)
  ll <- d |> group_by(participant, condition) |> summarise(words = list(wid), .groups = "drop")
  conjuntos <- list(HC = "HC", OFF = "OFF", ON = "ON", conjunta_3 = c("HC", "OFF", "ON"), conjunta_PD = c("OFF", "ON"))
  for (nm in names(conjuntos)) {
    ws <- ll$words[ll$condition %in% conjuntos[[nm]]]
    for (m in names(modos)) {
      g <- red_grafo(ws, lapply(ws, seq_along), vocab, m)
      fam <- if (nm %in% c("HC", "OFF", "ON")) "propia" else nm
      guardar_red(paste(tarea, nm, modos[[m]], sep = "_"), g, fam, tarea, nm, modos[[m]],
        n_listas = length(ws), n_tokens = sum(lengths(ws)))
    }
  }
}
idx_full <- bind_rows(filas_index); edges_full <- bind_rows(filas_edges)
# Verificación (a): mismas aristas binarias (sin bucles) que el output 14, que es la entrada preparada para el SBM.
e14 <- read_csv(file.path(d14, "edges.csv"), show_col_types = FALSE) |> filter(!is_loop)
for (id in idx_full$graph[idx_full$family == "propia"]) {
  dirigida <- grepl("_directed$", id)
  norm_par <- function(a, b) if (dirigida) paste(a, b, sep = "\r") else paste(pmin(a, b), pmax(a, b), sep = "\r")
  mia <- edges_full |> filter(graph == id); suya <- e14 |> filter(graph == id)
  anota(paste("Aristas = output 14:", id), setequal(norm_par(mia$from, mia$to), unique(norm_par(suya$from, suya$to))),
    paste(nrow(mia), "aristas"))
}
# Verificación (b): la conjunta contiene exactamente la unión de las aristas de sus grupos.
for (tarea in tareas) for (tipo in modos) {
  key <- function(g) { e <- edges_full |> filter(graph == g)
    if (tipo == "directed") paste(e$from, e$to, sep = "\r") else paste(pmin(e$from, e$to), pmax(e$from, e$to), sep = "\r") }
  u3 <- unique(c(key(paste(tarea, "HC", tipo, sep = "_")), key(paste(tarea, "OFF", tipo, sep = "_")), key(paste(tarea, "ON", tipo, sep = "_"))))
  uP <- unique(c(key(paste(tarea, "OFF", tipo, sep = "_")), key(paste(tarea, "ON", tipo, sep = "_"))))
  anota(paste("Conjunta_3 = unión", tarea, tipo), setequal(key(paste(tarea, "conjunta_3", tipo, sep = "_")), u3))
  anota(paste("Conjunta_PD = unión", tarea, tipo), setequal(key(paste(tarea, "conjunta_PD", tipo, sep = "_")), uP))
}
# ---- (c): reproducir las extracciones condicionales de 07 ----------------------------------------
# Se repite exactamente la secuencia de números aleatorios de 07 (misma semilla, mismos planes y mismo
# orden de bucles); las extracciones que no se exportan se consumen sin calcular métricas.
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
anota("Planes de personas = 07", identical(plan_pd_boot, plan07$PD_boot) && identical(plan_hc_boot, plan07$HC_boot) &&
  identical(plan_hc_cond, plan07$HC_conditional))
preparar <- function(dd, tarea) {
  dd <- dd |> filter(task == tarea) |> arrange(participant, condition, position, subposition)
  vocab <- unique(dd$word_norm); dd$wid <- match(dd$word_norm, vocab)
  ll <- dd |> group_by(participant, cohort, condition, order) |> summarise(words = list(wid), .groups = "drop")
  off <- ll |> filter(condition == "OFF") |> arrange(participant)
  on <- ll |> filter(condition == "ON") |> arrange(participant)
  hc <- ll |> filter(condition == "HC") |> arrange(participant)
  list(OFF = off$words, ON = on$words, HC = hc$words, PD_id = off$participant, HC_id = hc$participant, vocab = vocab)
}
extraer_sel <- function(p, boot, b) {
  id_pd <- if (boot) plan_pd_boot[b, ] else seq_along(p$OFF)
  id_hc <- if (boot) plan_hc_boot[b, ] else plan_hc_cond[b, ]
  off <- p$OFF[id_pd]; on <- p$ON[id_pd]; hc <- p$HC[id_hc]
  rp <- order(pmin(lengths(off), lengths(on)), id_pd); rh <- order(lengths(hc), id_hc)
  off <- off[rp]; on <- on[rp]; hc <- hc[rh]
  k <- pmin(lengths(off), lengths(on), lengths(hc))
  list(OFF = list(w = off, ix = seleccionar_posiciones(off, k)), ON = list(w = on, ix = seleccionar_posiciones(on, k)),
       HC = list(w = hc, ix = seleccionar_posiciones(hc, k)), tokens = sum(k))
}
sim07 <- readRDS(file.path(d07, "simulaciones.rds"))$metrics |>
  filter(scenario == "recorded", phase == "conditional", draw <= D, metric %in% c("nodes", "edges"))
datos <- escenario_datos(base, "recorded")
for (tarea in tareas) {
  p <- preparar(datos, tarea)
  for (fase in c("conditional", "participant_bootstrap")) for (b in seq_len(B)) {
    x <- extraer_sel(p, fase == "participant_bootstrap", b)
    if (fase != "conditional" || b > D) next
    for (cond in c("OFF", "ON", "HC")) for (m in names(modos)) {
      g <- red_grafo(x[[cond]]$w, x[[cond]]$ix, p$vocab, m)
      ref <- sim07 |> filter(task == tarea, draw == b, condition == cond, mode == m)
      anota(sprintf("Igualada = 07: %s %s %s extracción %d", tarea, cond, m, b),
        vcount(g) == ref$value[ref$metric == "nodes"] && ecount(g) == ref$value[ref$metric == "edges"])
      guardar_red(sprintf("%s_%s_%s_igualada%02d", tarea, cond, modos[[m]], b), g, "igualada", tarea, cond, modos[[m]],
        draw = b, n_listas = 21L, n_tokens = x$tokens)
    }
  }
}
index <- bind_rows(filas_index); edges <- bind_rows(filas_edges); nodes <- bind_rows(filas_nodes)
verificaciones <- bind_rows(verif)
stopifnot(nrow(index) == 30L + 3L * 3L * 2L * D, !anyDuplicated(index$graph), all(verificaciones$passed))
out <- nuevo_directorio("18_exportar_redes_SBM")
for (nm in c("index", "edges", "nodes", "verificaciones")) guardar_csv(get(nm), file.path(out, paste0(nm, ".csv")))
registrar_entorno(out, c(entrada, file.path(d14, "edges.csv"), file.path(d07, "simulaciones.rds"), file.path(d07, "plan_personas.rds"),
  "analisis/R/18_exportar_redes_SBM.R", "analisis/R/00_redes_metricas.R", "analisis/R/00_sensibilidad_datos.R", "analisis/R/00_utilidades.R"))
writeLines(c("# Redes exportadas para el SBM R3 — D-021",
  sprintf("%d redes: 18 propias completas, 12 conjuntas (HC+ON+OFF y ON+OFF por tarea y tipo) y %d igualadas (07, extracciones condicionales 1–%d).", nrow(index), 18L * D, D),
  "Aristas binarias, sin bucles ni multiaristas. Secuencia exportada dirigida; el SBM la usa como no dirigida (D-021).",
  "Propias: idénticas al output 14 (sin bucles). Conjuntas: unión exacta de las aristas de los grupos. Igualadas: mismos nodos y aristas que 07.",
  "Las igualadas usan la secuencia reconstruida entre posiciones retenidas (sequence_rebuilt de 07), no solo adyacencias originales.",
  "HC sigue siendo provisional por la celda partir<TAB>llegar (P11); las repeticiones se conservan (P3).",
  "No contiene particiones: ver analisis/python/19_SBM_redes_limpias.py."), file.path(out, "LEEME.md"), useBytes = TRUE)
print(index |> count(family, network)); cat("Comprobaciones:", nrow(verificaciones), "superadas\nSalida:", out, "\n")
