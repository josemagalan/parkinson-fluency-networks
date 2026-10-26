# Propósito: regenerar Figura 4 conservando posición y usando bootstrap por persona.
# Entradas: output/03_orden_sesiones/tokens.rds; puntuaciones/cobertura actuales.
# Salidas: output/10_perfiles_posicion_corregidos[/_vN]/, CSV y figura PNG/SVG.
# E2/E5/E6: perfiles descriptivos; no pruebas por posición ni ajuste causal de condición.
# Paquetes: 00_utilidades.R + ggplot2. Bootstrap B=1000, secuencia PD estratificada.
# --- English ---
# Purpose: corrected Figure 4 profiles by true production position with person-bootstrap bands
#   (B = 1000, PD stratified by sequence). Descriptive (E2, E5, E6).
# Inputs: 03 tokens.rds. Outputs: analisis/output/10_perfiles_posicion_corregidos[_vN]/ (CSV, PNG/SVG).
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
suppressPackageStartupMessages(library(ggplot2))
set.seed(20261022); B <- 1000L
entrada <- "analisis/output/03_orden_sesiones/tokens.rds"
t <- readRDS(entrada)
pd <- t |> filter(cohort == "PD") |> distinct(participant, order) |> arrange(participant)
hc <- sort(unique(t$participant[t$cohort == "HC"]))
Wpd <- matrix(0L, nrow = nrow(pd), ncol = B)
for (s in c("OFF_ON", "ON_OFF")) {
  ii <- which(pd$order == s); Wpd[ii, ] <- rmultinom(B, length(ii), rep(1 / length(ii), length(ii)))
}
Whc <- rmultinom(B, length(hc), rep(1 / length(hc), length(hc)))
stopifnot(all(colSums(Wpd) == 21L), all(colSums(Whc) == 24L))
perfiles <- list(); cursor <- 0L
for (task in c("F", "actions", "animals")) for (v in variables_psico) {
  if (v == "motor" && task != "actions") next
  for (cond in c("OFF", "ON", "HC")) {
    dd <- t[t$task == task & t$condition == cond, ]
    ids <- if (cond == "HC") hc else pd$participant
    W <- if (cond == "HC") Whc else Wpd
    maxpos <- max(t$position[t$task == task])
    X <- matrix(NA_real_, nrow = length(ids), ncol = maxpos)
    stopifnot(!anyDuplicated(dd[c("participant", "position")]))
    X[cbind(match(dd$participant, ids), dd$position)] <- dd[[v]]
    available <- !is.na(X); Xzero <- X; Xzero[!available] <- 0
    # Ceros solo para suma matricial; los denominadores excluyen todos los valores NA.
    den <- crossprod(W, available * 1); num <- crossprod(W, Xzero)
    bootmeans <- num / den; bootmeans[den == 0] <- NA_real_
    n <- colSums(available); mu <- colSums(Xzero) / n; mu[n == 0] <- NA_real_
    lo <- apply(bootmeans, 2, function(x) if (all(is.na(x))) NA_real_ else unname(quantile(x, .025, na.rm = TRUE)))
    hi <- apply(bootmeans, 2, function(x) if (all(is.na(x))) NA_real_ else unname(quantile(x, .975, na.rm = TRUE)))
    cursor <- cursor + 1L
    perfiles[[cursor]] <- tibble(task = task, variable = v, condition = cond, position = seq_len(maxpos),
      n_with_score = n, mean = mu, bootstrap_low = lo, bootstrap_high = hi,
      valid_bootstrap_fraction = colMeans(!is.na(bootmeans)),
      show = n >= 5L & valid_bootstrap_fraction >= .95)
  }
}
perfiles <- bind_rows(perfiles)
n_producen <- t |> group_by(task, condition, position) |> summarise(n_produced = n_distinct(participant), .groups = "drop")
perfiles <- left_join(perfiles, n_producen, by = c("task", "condition", "position")) |>
  mutate(n_produced = coalesce(n_produced, 0L), coverage_at_position = if_else(n_produced > 0, n_with_score / n_produced, NA_real_))
stopifnot(all(perfiles$n_with_score <= perfiles$n_produced))
mostrar <- perfiles |> filter(show) |>
  mutate(task = factor(task, levels = c("F", "actions", "animals")),
    variable = factor(variable, levels = c("AoA", "concreteness", "familiarity", "imageability", "logFreq", "motor")),
    condition = factor(condition, levels = c("HC", "OFF", "ON")))
fig <- ggplot(mostrar, aes(position, mean, colour = condition, fill = condition)) +
  geom_ribbon(aes(ymin = bootstrap_low, ymax = bootstrap_high), alpha = .12, colour = NA) +
  geom_line(linewidth = .45) + geom_point(size = .65) +
  facet_grid(variable ~ task, scales = "free_y", drop = FALSE) +
  scale_colour_manual(values = c(HC = "#2A8C6D", OFF = "#D36A37", ON = "#80558C")) +
  scale_fill_manual(values = c(HC = "#2A8C6D", OFF = "#D36A37", ON = "#80558C")) +
  labs(x = "Position in the produced list", y = "Mean available word rating", colour = "Condition", fill = "Condition",
    caption = "Original production order; ≥5 contributors per point. Pointwise participant-bootstrap bands (1,000 draws).\nDescriptive profiles; late positions include fewer participants. AoA remains provisional.") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", panel.grid.minor = element_blank(),
    plot.caption = element_text(size = 8, hjust = 0), strip.text = element_text(face = "bold"))
out <- nuevo_directorio("10_perfiles_posicion_corregidos")
guardar_csv(perfiles, file.path(out, "perfiles.csv")); guardar_csv(n_producen, file.path(out, "n_por_posicion.csv"))
ggsave(file.path(out, "figura4_posicion_corregida.png"), fig, width = 12, height = 12, dpi = 180)
ggsave(file.path(out, "figura4_posicion_corregida.svg"), fig, width = 12, height = 12)
registrar_entorno(out, c(entrada, "analisis/R/10_perfiles_posicion_corregidos.R", "analisis/R/00_utilidades.R"))
writeLines(c("# Perfiles por posición corregidos — 22-10-2026",
  "Se usa la posición original, no la ordenación alfabética de Sequences.R.",
  "Bootstrap por persona B=1000: pacientes con condiciones/tareas juntas, estratos de secuencia 11/10; HC independiente.",
  "Medias de puntuaciones disponibles y bandas percentiles puntuales; no IC simultáneos ni pruebas de condición.",
  "Mostrar solo N con puntuación >=5 y al menos 95% réplicas con valor; CSV conserva todas las posiciones/N/coberturas.",
  "No se ajusta aquí por sesión: usar 06 para contrastes generales. Bandas superpuestas/no superpuestas no prueban diferencias.",
  "Posiciones tardías condicionan a quienes produjeron más palabras; no hay tiempo real de emisión.",
  "AoA/P3/P11 siguen pendientes; esta es una figura de trabajo, no una edición del manuscrito."), file.path(out, "LEEME.md"), useBytes = TRUE)
cat("Salida:", out, "\n")
