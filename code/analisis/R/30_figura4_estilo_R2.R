# Propósito: Figura 4 con estilo R2, posición original y bandas loess por persona.
# Entradas: output/03_orden_sesiones/tokens.rds; estilo de ../parkinson/Sequences.R.
# Salidas: output/30_figura4_estilo_R2[_vN], figuras, curvas, N, diagnóstico y sesión.
# E2/E5/E6. Descriptivo, sin tests por posición ni ajuste de condición por sesión.
# Paquetes: 00_utilidades.R (readxl/readr/dplyr/tidyr/purrr/broom), ggplot2, stats.
# Ejecutar desde revision_R3/: Rscript --vanilla analisis/R/30_figura4_estilo_R2.R
# --- English ---
# Purpose (optional): corrected Figure 4 in the R2 style: ratings by original position, LOESS and
#   person-bootstrap bands (B = 1000) (E2, E5, E6).
# Inputs: 03 tokens.rds. Outputs: analisis/output/30_figura4_estilo_R2[_vN]/ (PNG/PDF/SVG, curves, caption draft).
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
suppressPackageStartupMessages(library(ggplot2))
set.seed(20261022)
B <- 1000L; span <- .75
entrada <- "analisis/output/03_orden_sesiones/tokens.rds"
t <- readRDS(entrada)
stopifnot(nrow(t) == 2965L, !anyDuplicated(t$token_id),
          !anyDuplicated(t[c("participant", "task", "condition", "position")]))
huella_inicial <- tools::md5sum(entrada)
pd <- t |> filter(cohort == "PD") |> distinct(participant, order) |> arrange(participant)
hc <- sort(unique(t$participant[t$cohort == "HC"]))
Wpd <- matrix(0L, nrow = nrow(pd), ncol = B)
for (s in c("OFF_ON", "ON_OFF")) {
  ii <- which(pd$order == s)
  Wpd[ii, ] <- rmultinom(B, length(ii), rep(1 / length(ii), length(ii)))
}
Whc <- rmultinom(B, length(hc), rep(1 / length(hc), length(hc)))
stopifnot(all(colSums(Wpd) == 21L), all(colSums(Whc) == 24L),
          all(colSums(Wpd[pd$order == "OFF_ON", , drop = FALSE]) == 11L))

# Misma estadística de R2: loess de medias disponibles por posición (peso igual
# por posición). Se vuelve a ajustar en CADA bootstrap, sin usar el SE de loess.
# Solo se predice dentro del soporte original N>=5 y >=95% réplicas válidas.
avisos_loess <- 0L
suavizar <- function(x, y, grid) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4L) return(rep(NA_real_, length(grid)))
  tryCatch(withCallingHandlers({
    fit <- stats::loess(y ~ x, data = data.frame(x = x[ok], y = y[ok]),
      span = span, degree = 2, family = "gaussian",
      control = stats::loess.control(surface = "direct"))
    pred <- as.numeric(predict(fit, newdata = data.frame(x = grid)))
    pred[grid < min(x[ok]) | grid > max(x[ok])] <- NA_real_
    pred
  }, warning = function(w) {
    avisos_loess <<- avisos_loess + 1L
    invokeRestart("muffleWarning")
  }), error = function(e) rep(NA_real_, length(grid)))
}
curvas <- list(); perfiles <- list(); diagnostico <- list(); idx <- 0L
for (task_i in c("actions", "F", "animals")) for (v in
    c("AoA", "concreteness", "familiarity", "imageability", "logFreq", "motor")) {
  if (v == "motor" && task_i != "actions") next
  for (cond in c("HC", "ON", "OFF")) {
    dd <- t[t$task == task_i & t$condition == cond, ]
    ids <- if (cond == "HC") hc else pd$participant
    W <- if (cond == "HC") Whc else Wpd
    posiciones <- seq_len(max(t$position[t$task == task_i]))
    X <- matrix(NA_real_, nrow = length(ids), ncol = length(posiciones))
    X[cbind(match(dd$participant, ids), dd$position)] <- dd[[v]]
    available <- !is.na(X); Xzero <- X; Xzero[!available] <- 0
    n <- colSums(available)
    mu <- colSums(Xzero) / n; mu[n == 0] <- NA_real_
    den <- crossprod(W, available * 1)
    bootmeans <- crossprod(W, Xzero) / den
    bootmeans[den == 0] <- NA_real_
    valid <- colMeans(is.finite(bootmeans))
    show <- n >= 5L & valid >= .95
    stopifnot(sum(show) >= 4L)
    x <- posiciones[show]
    grid <- sort(unique(c(seq(min(x), max(x), length.out = 160L), x)))
    mean_smooth <- suavizar(x, mu[show], grid)
    boots <- vapply(seq_len(B), function(b) suavizar(x, bootmeans[b, show], grid),
                    numeric(length(grid)))
    valid_smooth <- rowMeans(is.finite(boots))
    percentil <- function(q) apply(boots, 1, function(z) {
      if (all(!is.finite(z))) NA_real_ else unname(quantile(z[is.finite(z)], q))
    })
    idx <- idx + 1L
    curvas[[idx]] <- tibble(task = task_i, variable = v, condition = cond,
      position = grid, estimate = mean_smooth, low = percentil(.025), high = percentil(.975),
      bootstrap_valid_fraction = valid_smooth)
    nproducen <- tabulate(dd$position, nbins = length(posiciones))
    perfiles[[idx]] <- tibble(task = task_i, variable = v, condition = cond,
      position = posiciones, mean = mu, n_with_score = n, n_produced = nproducen,
      coverage = ifelse(nproducen > 0, n / nproducen, NA_real_),
      bootstrap_valid_fraction = valid, show = show)
    diagnostico[[idx]] <- tibble(task = task_i, variable = v, condition = cond,
      first_position = min(x), last_position = max(x), n_positions = length(x),
      min_bootstrap_valid = min(valid_smooth),
      finite_curve = all(is.finite(mean_smooth)))
    cat(task_i, v, cond, "\n")
    flush.console()
  }
}
curvas <- bind_rows(curvas); perfiles <- bind_rows(perfiles)
diagnostico <- bind_rows(diagnostico)
stopifnot(nrow(diagnostico) == 48L, all(diagnostico$finite_curve),
          all(curvas$low <= curvas$high, na.rm = TRUE),
          identical(huella_inicial, tools::md5sum(entrada)))

raw <- t |> pivot_longer(all_of(variables_psico), names_to = "variable", values_to = "value") |>
  filter(is.finite(value), variable != "motor" | task == "actions")
# Factores y paleta del gráfico histórico, con etiquetas completas de la R2.
formato <- function(d) d |>
  mutate(task = factor(task, levels = c("actions", "F", "animals"),
    labels = c("Action fluency", "Phonological fluency", "Semantic fluency")),
    variable = factor(variable, levels = c("AoA", "concreteness", "familiarity", "imageability", "logFreq", "motor")),
    condition = factor(condition, levels = c("HC", "ON", "OFF")))
raw <- formato(raw)
plot_curves <- formato(curvas |> filter(bootstrap_valid_fraction >= .95))
palette <- c(HC = "#35b779", ON = "#781c6d", OFF = "#ed6925")
labels <- c(HC = "Control", ON = "PD on medication", OFF = "PD off medication")
fig <- ggplot(raw, aes(position, value, colour = condition)) +
  geom_point(alpha = .60, size = 1.0, stroke = 0) +
  geom_ribbon(data = plot_curves, aes(x = position, ymin = low, ymax = high,
    group = condition), inherit.aes = FALSE, fill = "grey55", alpha = .22) +
  geom_line(data = plot_curves, aes(position, estimate, colour = condition),
    inherit.aes = FALSE, linewidth = .65) +
  facet_grid(rows = vars(variable), cols = vars(task), scales = "free", drop = FALSE) +
  scale_colour_manual(values = palette, labels = labels, name = NULL) +
  labs(x = "Position in production sequence", y = "Value") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4))) +
  theme_gray(base_size = 12) +
  theme(legend.position = "bottom", legend.text = element_text(size = 14, face = "bold"),
    strip.text.x = element_text(size = 13, face = "bold"),
    strip.text.y = element_text(size = 11),
    axis.title = element_text(size = 12),
    panel.spacing = grid::unit(.12, "cm"),
    plot.margin = margin(6, 6, 6, 6))
out <- nuevo_directorio("30_figura4_estilo_R2")
guardar_csv(curvas, file.path(out, "curvas_loess_bootstrap.csv"))
guardar_csv(perfiles, file.path(out, "medias_N_cobertura.csv"))
guardar_csv(diagnostico, file.path(out, "diagnostico_curvas.csv"))
saveRDS(fig, file.path(out, "figura4_estilo_R2.rds"))
ggsave(file.path(out, "figura4_estilo_R2.png"), fig, width = 12, height = 16, dpi = 300)
ggsave(file.path(out, "figura4_estilo_R2_preview.png"), fig, width = 12, height = 16, dpi = 150)
ggsave(file.path(out, "figura4_estilo_R2.pdf"), fig, width = 12, height = 16)
ggsave(file.path(out, "figura4_estilo_R2.svg"), fig, width = 12, height = 16)
writeLines(c(
  "# Figura 4 corregida con estructura y estilo R2",
  "Preferencia de Josema, 22-10-2026: conservar el estilo R de la R2.",
  "Posición original desde dataset03. Nubes = puntuaciones individuales disponibles.",
  "Columnas: acciones/F/animales; seis variables; contenido motor solo acciones.",
  "Paleta de Sequences.R, paneles grises, puntos alfa0.6, leyenda inferior.",
  "Curvas loess sobre medias por posición, span0.75, degree2; posiciones con N>=5",
  "y al menos95% de réplicas de media válidas. Sin extrapolar fuera de ese soporte.",
  "Bandas grises: percentiles puntuales2.5/97.5 de1000 bootstrap de PERSONAS,",
  "volviendo a ajustar el mismo loess en cada réplica; no SE de geom_smooth.",
  "PD: ON/OFF y tareas juntos, estratos de secuencia11/10; HC independiente.",
  "Estos intervalos son descriptivos para la curva suavizada, no simultáneos ni",
  "tests de condición o pendientes, y no incorporan sesgo del suavizado.",
  "Se mantienen todos los puntos observados; las curvas acaban antes que las nubes",
  "cuando falta soporte. Posiciones tardías condicionan a quienes llegan a ellas.",
  "Sin ajuste por sesión ni interpretación causal. AoA/P3/P11 pendientes.",
  paste("Avisos numéricos loess capturados:", avisos_loess),
  "Archivos PNG300dpi/PDF/SVG; preview150dpi y RDS del gráfico.",
  "Ejecución opcional desde revision_R3/: Rscript --vanilla analisis/R/30_figura4_estilo_R2.R",
  "No cambia los modelos ni los outputs01–17; script30 fuera del lanzador01–19."
), file.path(out, "LEEME.md"), useBytes = TRUE)
writeLines(c(
  "Figure 4. Psycholinguistic characteristics by position in the verbal fluency response sequence.",
  "Columns show action, phonological and semantic fluency; rows show age of acquisition (AoA), concreteness, familiarity, imageability, lexical frequency (logFreq), and motor content (action fluency only). Points represent available word ratings in the original production order. Coloured curves are LOESS fits to mean ratings by position (span = 0.75). Grey bands show pointwise 95% participant-bootstrap percentile intervals for the fitted curves (1,000 resamples), retaining paired ON/OFF assessments and stratifying patients by assessment sequence. Curves are shown only over positions supported by at least five available ratings and adequate bootstrap coverage. All available individual points are retained. Later positions represent fewer participants who produced longer lists. Profiles are descriptive and are not adjusted for assessment period; the bands do not constitute tests of condition differences. AoA assignment remains provisional in this working version."
), file.path(out, "figure4_caption_draft.txt"), useBytes = TRUE)
registrar_entorno(out, c(entrada, "analisis/R/30_figura4_estilo_R2.R",
  "analisis/R/00_utilidades.R", "../parkinson/Sequences.R"))
cat("Salida:", out, "\nAvisos loess:", avisos_loess, "\n")
