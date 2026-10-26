# Propósito: sensibilidad de red a permutación de secuencia manteniendo periodo y pareado.
# Entradas: tokens con session/order. Supuesto de intercambiabilidad NO verificado (P14).
# Salidas: output/12_permutacion_secuencia_sensibilidad[/_vN]/, p condicionales BH54 y estadísticos.
# E1/E2/E4/E5: exploratorio; sin afirmación de aleatorización, causalidad ni SBM.
# Paquetes: 00_utilidades.R + igraph, 00_redes_metricas.R. B=999, K=3 selecciones fijas.
# --- English ---
# Purpose: exploratory permutation of ON/OFF sequence labels (B = 999, K = 3 fixed selections, BH54)
#   for network metric differences; exchangeability assumed, not verified (P14) (E1, E2, E4, E5).
# Inputs: 03 tokens.rds. Outputs: analisis/output/12_permutacion_secuencia_sensibilidad[_vN]/.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
source("analisis/R/00_redes_metricas.R", encoding = "UTF-8")
suppressPackageStartupMessages(library(igraph))
set.seed(20261022); B <- 999L; K <- 3L
entrada <- "analisis/output/03_orden_sesiones/tokens.rds"
t <- readRDS(entrada) |> filter(cohort == "PD") |> arrange(participant, task, session, position)
personas <- t |> distinct(participant, order) |> arrange(participant)
observed <- personas$order == "OFF_ON"
stopifnot(sum(observed) == 11L, length(observed) == 21L)
assignments <- rbind(observed, t(replicate(B, {x <- rep(FALSE, 21); x[sample.int(21, 11)] <- TRUE; x})))
stopifnot(all(rowSums(assignments) == 11L))
modos <- c("cooccurrence", "sequence_rebuilt", "sequence_original")
stats <- list(); cursor <- 0L
for (tarea in c("F", "actions", "animals")) {
  d <- t |> filter(task == tarea); vocab <- unique(d$word_norm); d$wid <- match(d$word_norm, vocab)
  ll <- d |> group_by(participant, session) |> summarise(w = list(wid), .groups = "drop")
  s1 <- ll$w[match(personas$participant, ll$participant[ll$session == 1])]
  # Seleccionar primero cada sesión para evitar índices del data frame combinado.
  a <- ll |> filter(session == 1L) |> arrange(participant)
  b <- ll |> filter(session == 2L) |> arrange(participant)
  stopifnot(identical(a$participant, personas$participant), identical(b$participant, personas$participant))
  ws1 <- a$w; ws2 <- b$w; k <- pmin(lengths(ws1), lengths(ws2)); V <- length(vocab)
  selections <- lapply(seq_len(K), function(j) list(i1 = seleccionar_posiciones(ws1, k), i2 = seleccionar_posiciones(ws2, k)))
  arr <- array(NA_real_, c(B + 1L, length(modos), 6L), dimnames = list(NULL, modos,
    c("nodes", "edges", "density", "transitivity", "mean_degree_total", "largest_weak_component_pct")))
  for (r in seq_len(B + 1L)) {
    z <- assignments[r, ]; on <- Map(function(x, y, second) if (second) y else x, ws1, ws2, z)
    off <- Map(function(x, y, second) if (second) x else y, ws1, ws2, z)
    value <- matrix(0, nrow = length(modos), ncol = 6L)
    for (sel in selections) {
      ion <- Map(function(x, y, second) if (second) y else x, sel$i1, sel$i2, z)
      ioff <- Map(function(x, y, second) if (second) x else y, sel$i1, sel$i2, z)
      for (m in seq_along(modos)) value[m, ] <- value[m, ] +
        red_metricas(on, ion, V, modos[m]) - red_metricas(off, ioff, V, modos[m])
    }
    arr[r, , ] <- value / K
  }
  cursor <- cursor + 1L; stats[[tarea]] <- arr
  cat("Permutación terminada:", tarea, B, "reasignaciones ×", K, "selecciones\n"); flush.console()
}
resultados <- bind_rows(lapply(names(stats), function(task) {
  a <- stats[[task]]
  bind_rows(lapply(seq_along(modos), function(m) {
    obs <- a[1, m, ]; perm <- a[-1, m, ]
    tibble(task, mode = modos[m], metric = names(obs), observed_delta = as.numeric(obs),
      null_mean = colMeans(perm), p = vapply(seq_along(obs), function(j) (1 + sum(abs(perm[, j]) >= abs(obs[j]) - 1e-12)) / (B + 1), numeric(1)))
  }))
})) |> mutate(p_BH54 = p.adjust(p, "BH"), assumption = "sequence_exchangeability_unverified")
stopifnot(nrow(resultados) == 54L, !anyNA(resultados$p), all(resultados$p > 0 & resultados$p <= 1))
out <- nuevo_directorio("12_permutacion_secuencia_sensibilidad")
guardar_csv(resultados, file.path(out, "resultados.csv")); saveRDS(list(statistics = stats, assignments = assignments,
  persons = personas, K = K, B = B, seed = 20261022), file.path(out, "permutaciones.rds"))
registrar_entorno(out, c(entrada, "analisis/R/12_permutacion_secuencia_sensibilidad.R", "analisis/R/00_redes_metricas.R"))
writeLines(c("# Permutación condicional de secuencia — sensibilidad",
  "Se mantienen listas de primera/segunda evaluación, personas y k mínimo; solo se reasigna qué secuencia de condición recibe cada persona.",
  "11 OFF→ON y 10 ON→OFF en cada reasignación; las mismas asignaciones se usan en todas las tareas.",
  "Tres selecciones de posiciones fijadas antes de permutar, idéntico estadístico observado/nulo (media de tres redes por representación).",
  "999 reasignaciones; p bilateral (1+conteo)/(B+1), BH54 para todos los contrastes, resolución mínima 0,001.",
  "Supone intercambiabilidad de secuencia y ausencia de arrastre diferencial. No se ha confirmado aleatorización: P14.",
  "Resultados exploratorios condicionados a ese supuesto; no p confirmatorios ni prueba causal de reorganización.",
  "No incluye HC ni particiones; no sustituye la futura inferencia core-periphery."), file.path(out, "LEEME.md"))
print(resultados |> filter(p_BH54 < .05), n = 54); cat("Salida:", out, "\n")
