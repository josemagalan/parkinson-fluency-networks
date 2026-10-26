# Propósito: contrastes de la Tabla 2 (modelo de 21) en F secuencia, cuyas particiones propias ON/OFF son inestables
#   (script 19: cadenas discordantes; MDL hub-and-spoke ≈ «sin core»). En lugar de usar una partición de consenso
#   arbitraria, se repite el modelo con muchas particiones plausibles del SBM (script 46) y se combina la incertidumbre
#   de la partición con las reglas de Rubin (cada partición posterior se trata como una imputación del bloque latente).
#   Responde a qué ocurre con las dos diferencias de grupo de la R2 en F secuencia (familiaridad del core OFF−HC y
#   OFF × Periphery en logFreq) y si incluir F secuencia cambiaría las decisiones BH de la Tabla 2. Sensibilidad D-030,
#   añadida el 24-10-2026 tras conocer 19/21; no sustituye la regla «n.i.» de la Tabla 2 principal.
# Entradas: output/03_orden_sesiones/tokens.rds; output/19_SBM_redes_limpias/particiones_consenso.csv (HC y referencia);
#   output/46_SBM_F_secuencia_cadenas/particiones.csv; output/21_tabla2_niveles_R2_v2/resultados.csv (familias BH y verificación).
# Salidas: output/47_tabla2_F_secuencia_particiones[_vN]/ (estimaciones por partición, combinación de Rubin, distribución,
#   BH con 26 modelos, verificación frente a 21, LEEME).
# E2/E5. Paquetes: 00_utilidades.R + nlme (vía 00_niveles_tabla2.R).
# Ejecutar desde revision_R3/: Rscript --vanilla analisis/R/47_tabla2_F_secuencia_particiones.R
# --- English ---
# Purpose: script 21 model for F-sequence over many plausible SBM partitions (script 46), pooled with
#   Rubin's rules; exploratory sensitivity (D-030, D-032) (E2, E5).
# Inputs: 03, 19, 46 and 21_v2 outputs. Outputs: analisis/output/47_tabla2_F_secuencia_particiones[_vN]/.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
source("analisis/R/00_sensibilidad_datos.R", encoding = "UTF-8")
source("analisis/R/00_niveles_tabla2.R", encoding = "UTF-8")
set.seed(20261024)
entrada <- "analisis/output/03_orden_sesiones/tokens.rds"
f19 <- "analisis/output/19_SBM_redes_limpias/particiones_consenso.csv"
f46 <- "analisis/output/46_SBM_F_secuencia_cadenas/particiones.csv"
f21 <- "analisis/output/21_tabla2_niveles_R2_v2/resultados.csv"
base <- readRDS(entrada)
stopifnot(nrow(base) == 2965L)
datos <- escenario_datos(base) |> filter(task == "F")
p19 <- read_csv(f19, show_col_types = FALSE, trim_ws = FALSE) |>
  filter(family == "propia", task == "F", network == "directed") |> select(group, word_norm, hs = hs_consenso)
p46 <- read_csv(f46, show_col_types = FALSE, trim_ws = FALSE) |>
  mutate(group = sub("^F_(ON|OFF)_directed$", "\\1", graph))
r21 <- read_csv(f21, show_col_types = FALSE)
vars_F <- c("logFreq", "imageability", "concreteness", "familiarity", "AoA")
# ---- Medias por persona/lista/partición con una partición dada (misma definición que 21) ------------------
niveles_F <- function(part) {
  d <- datos |> left_join(part, by = c("condition" = "group", "word_norm"), relationship = "many-to-one")
  stopifnot(!anyNA(d$hs))
  d |> pivot_longer(all_of(vars_F), names_to = "variable", values_to = "value") |> filter(!is.na(value)) |>
    mutate(partition = if_else(hs == 0, "core", "periphery")) |>
    group_by(participant, cohort, condition, session, order, variable, partition) |>
    summarise(n_tokens = n(), mean = mean(value), .groups = "drop")
}
ajustar_F <- function(part) {
  nv <- niveles_F(part)
  nv |> group_by(variable) |> group_modify(~ suppressWarnings(ajustar_uno(.x))) |> ungroup()
}
hc <- p19 |> filter(group == "HC")
# ---- Verificación: con las particiones de consenso de 19 se reproduce 21 ----------------------------------
ver <- ajustar_F(p19) |> inner_join(r21 |> filter(task == "F", network == "sequence") |>
  select(variable, term, estimate_21 = estimate, p_21 = p), by = c("variable", "term")) |>
  mutate(dif = abs(estimate - estimate_21))
stopifnot(max(ver$dif, na.rm = TRUE) < 1e-6)
# ---- Repetición con las particiones de 46 ------------------------------------------------------------
claves <- p46 |> distinct(tipo, semilla, muestra) |> arrange(tipo, semilla, muestra)
por_particion <- bind_rows(lapply(seq_len(nrow(claves)), function(i) {
  k <- claves[i, ]
  pk <- p46 |> filter(tipo == k$tipo, semilla == k$semilla, muestra == k$muestra) |> select(group, word_norm, hs)
  core <- pk |> group_by(group) |> summarise(core = sum(hs == 0), .groups = "drop")
  ajustar_F(bind_rows(hc, pk)) |> mutate(tipo = k$tipo, semilla = k$semilla, muestra = k$muestra,
    core_OFF = core$core[core$group == "OFF"], core_ON = core$core[core$group == "ON"])
}))
# ---- Combinación de Rubin (gl de Barnard y Rubin, 1999) ----------------------------------------------
rubin <- function(q, se, nu_com) {
  ok <- is.finite(q) & is.finite(se); q <- q[ok]; se <- se[ok]; nu_com <- min(nu_com[ok]); m <- length(q)
  if (m < 2) return(tibble(m = m, estimate = NA_real_, SE = NA_real_, df = NA_real_, p = NA_real_, fmi = NA_real_))
  qb <- mean(q); ub <- mean(se^2); b <- var(q); tt <- ub + (1 + 1 / m) * b
  lam <- (1 + 1 / m) * b / tt
  nu_old <- if (lam > 0) (m - 1) / lam^2 else Inf
  nu_obs <- (nu_com + 1) / (nu_com + 3) * nu_com * (1 - lam)
  nu <- 1 / (1 / nu_old + 1 / nu_obs)
  tibble(m = m, estimate = qb, SE = sqrt(tt), df = nu, p = 2 * pt(-abs(qb / sqrt(tt)), nu), fmi = lam,
    conf_low = qb - qt(.975, nu) * sqrt(tt), conf_high = qb + qt(.975, nu) * sqrt(tt))
}
terminos_inf <- setdiff(unique(por_particion$term), "Intercept_HC_core")
combinado <- por_particion |> filter(term %in% terminos_inf) |> group_by(tipo, variable, term) |>
  group_modify(~ bind_cols(rubin(.x$estimate, .x$SE, .x$df),
    tibble(n_particiones = nrow(.x), no_estimables = sum(!is.finite(.x$estimate)),
      est_p025 = quantile(.x$estimate, .025, na.rm = TRUE), est_mediana = median(.x$estimate, na.rm = TRUE),
      est_p975 = quantile(.x$estimate, .975, na.rm = TRUE), prop_p05_nominal = mean(.x$p < .05, na.rm = TRUE)))) |> ungroup()
# ---- ¿Cambiaría la Tabla 2 si F secuencia entrara en las familias BH (26 modelos)? -----------------------
principal <- r21 |> filter(!no_identificable, term %in% terminos_inf) |>
  mutate(family = paste0(if_else(variable == "AoA", "AoA_", ""), term)) |> select(family, task, network, variable, term, p, p_BH_21 = p_BH)
bh26 <- bind_rows(principal |> mutate(origen = "21"),
  combinado |> filter(tipo == "muestra") |> transmute(family = paste0(if_else(variable == "AoA", "AoA_", ""), term),
    task = "F", network = "sequence", variable, term, p, origen = "47_rubin")) |>
  group_by(family) |> mutate(m_familia = n(), p_BH_26 = p.adjust(p, "BH")) |> ungroup()
cambios <- bh26 |> filter(origen == "21") |> filter((p_BH_21 < .05) != (p_BH_26 < .05))
resumen_core <- por_particion |> distinct(tipo, semilla, muestra, core_OFF, core_ON) |> group_by(tipo) |>
  summarise(particiones = n(), core_OFF_min = min(core_OFF), core_OFF_mediana = median(core_OFF), core_OFF_max = max(core_OFF),
    core_ON_min = min(core_ON), core_ON_mediana = median(core_ON), core_ON_max = max(core_ON), .groups = "drop")
# ---- Salidas ---------------------------------------------------------------------------------------
out <- nuevo_directorio("47_tabla2_F_secuencia_particiones")
guardar_csv(por_particion, file.path(out, "estimaciones_por_particion.csv"))
guardar_csv(combinado, file.path(out, "combinacion_rubin.csv"))
guardar_csv(bh26, file.path(out, "BH_26_modelos.csv"))
guardar_csv(cambios, file.path(out, "cambios_decision_BH26.csv"))
guardar_csv(resumen_core, file.path(out, "tamano_core.csv"))
guardar_csv(ver, file.path(out, "verificacion_21.csv"))
registrar_entorno(out, c(entrada, f19, f46, f21, "analisis/R/47_tabla2_F_secuencia_particiones.R", "analisis/R/00_niveles_tabla2.R",
  "analisis/R/00_sensibilidad_datos.R", "analisis/R/00_utilidades.R"))
r2 <- combinado |> filter(tipo == "muestra", (variable == "familiarity" & term == "OFF") | (variable == "logFreq" & term == "OFF_x_Periphery"))
writeLines(c("# 47: Tabla 2 en F secuencia con incertidumbre de partición (D-030)",
  "Modelo idéntico a 21 (media ~ grupo × partición + sesión; 1 | participante/lista). HC con su partición de consenso (19, estable);",
  "OFF y ON con particiones del SBM de 46: 'muestra' = muestras posteriores individuales (principal, combinación de Rubin con gl de",
  "Barnard–Rubin); 'consenso_cadena' = consenso de cada cadena (descriptivo). La cadena s de OFF se empareja con la cadena s de ON.",
  sprintf("Verificación: con las particiones de 19 se reproduce 21 en F secuencia (diferencia máxima %.1e).", max(ver$dif, na.rm = TRUE)),
  "Diferencias de grupo de la R2 en F secuencia (combinación de Rubin, muestras posteriores):",
  sprintf("  %s %s: %.2f [%.2f; %.2f], p = %.3f; p < 0,05 nominal en %.0f %% de las particiones.", r2$variable, r2$term, r2$estimate,
    r2$conf_low, r2$conf_high, r2$p, 100 * r2$prop_p05_nominal),
  sprintf("BH con 26 modelos (22 de 21 + 4 de F secuencia combinados): %d cambios de decisión en los 22 modelos de 21.", nrow(cambios)),
  "Es una sensibilidad; la Tabla 2 principal mantiene F secuencia como n.i. Las reglas de Rubin suponen que las particiones son",
  "muestras de la posterior del bloque latente; con cadenas que no se mezclan entre sí es una aproximación."),
  file.path(out, "LEEME.md"), useBytes = TRUE)
print(combinado |> filter(tipo == "muestra", variable != "AoA") |> select(variable, term, estimate, conf_low, conf_high, p, fmi, prop_p05_nominal), n = 60)
print(resumen_core); print(cambios); cat("Salida:", out, "\n")
