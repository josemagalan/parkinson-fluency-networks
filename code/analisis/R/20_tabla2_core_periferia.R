# Propósito: nueva Tabla 2 (core-periphery) con la persona como unidad, según el plan preespecificado
#   notas/plan_tabla2_R3_20261022.md (D-023). Brecha core−periferia por persona/condición/tarea/variable
#   con partición de referencia común (red conjunta del SBM, script 19). Familias:
#   G gradiente (45 personas, t de una muestra), C ON−OFF y S periodo (D~z, 21 pares), H HC−OFF/HC−ON (Welch).
# Entradas: output/03_orden_sesiones/tokens.rds; output/19_SBM_redes_limpias/{particiones_consenso,resumen_redes}.csv.
# Salidas: output/20_tabla2_core_periferia[_vN]/ (resultados, tabla2_borrador, sensibilidades, brechas, descriptivos, cobertura, figura PNG/SVG, LEEME).
# E2/E4/E5/E6. Paquetes: 00_utilidades.R (readr, dplyr, tidyr, purrr, broom); gráficos base.
# Ejecutar desde revision_R3/: Rscript --vanilla analisis/R/20_tabla2_core_periferia.R
# --- English ---
# Purpose: alternative Table 2 (D-023): within-person core-periphery gap on a common joint SBM partition;
#   families G, C, S, H with BH per family (E2, E4, E5, E6).
# Inputs: 03 tokens.rds; 19 partitions and summary. Outputs: analisis/output/20_tabla2_core_periferia[_vN]/.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
source("analisis/R/00_sensibilidad_datos.R", encoding = "UTF-8")
set.seed(20261022)
entrada <- "analisis/output/03_orden_sesiones/tokens.rds"
d19 <- "analisis/output/19_SBM_redes_limpias"
base <- readRDS(entrada)
stopifnot(nrow(base) == 2965L)
part <- read_csv(file.path(d19, "particiones_consenso.csv"), show_col_types = FALSE, trim_ws = FALSE)
res19 <- read_csv(file.path(d19, "resumen_redes.csv"), show_col_types = FALSE)
redes <- c(undirected = "cooccurrence", directed = "sequence")
principales <- c("logFreq", "imageability", "concreteness", "familiarity", "motor")
# ---- Asignar a cada token su partición en la red de referencia -----------------------------------
# referencia: "conjunta_3", "conjunta_PD" o "propia" (la red del propio grupo, como en la R2).
anotar_particion <- function(d, referencia, red) {
  p <- part |> filter(network == red, if (referencia == "propia") family == "propia" else family == referencia) |>
    select(task, group, word_norm, p_core, hs_consenso, estable)
  if (referencia == "propia") {
    d |> left_join(p, by = c("task", "condition" = "group", "word_norm"), relationship = "many-to-one")
  } else {
    d |> left_join(p |> select(-group), by = c("task", "word_norm"), relationship = "many-to-one")
  }
}
# modo: "binaria" (consenso), "estables" (excluye nodos inestables), "continua" (pesos P(core)).
brechas <- function(datos, referencia, modo = "binaria") {
  bind_rows(lapply(names(redes), function(red) {
    d <- anotar_particion(datos, referencia, red)
    # Toda palabra producida está en la red de referencia que la contiene (la partición no puede faltar).
    if (referencia != "conjunta_PD") stopifnot(!anyNA(d$p_core))
    else stopifnot(!anyNA(d$p_core[d$cohort == "PD"]))
    if (modo == "estables") d <- d |> filter(estable)
    d |> pivot_longer(all_of(variables_psico), names_to = "variable", values_to = "value") |>
      filter(variable != "motor" | task == "actions", !is.na(value), !is.na(p_core)) |>
      group_by(participant, cohort, condition, session, order, task, variable) |>
      summarise(n_core = sum(hs_consenso == 0), n_periphery = sum(hs_consenso == 1),
        mean_core = if (modo == "continua") sum(p_core * value) / sum(p_core) else mean(value[hs_consenso == 0]),
        mean_periphery = if (modo == "continua") sum((1 - p_core) * value) / sum(1 - p_core) else mean(value[hs_consenso == 1]),
        .groups = "drop") |>
      mutate(gap = mean_core - mean_periphery, network = redes[[red]],
        gap = if_else(is.finite(gap) & (modo == "continua" | (n_core > 0 & n_periphery > 0)), gap, NA_real_))
  })) |> mutate(reference = referencia, mode = modo)
}
familia <- function(var, base_fam) if_else(var == "AoA", paste0("AoA_", base_fam), base_fam)
hc3 <- function(m) {
  X <- model.matrix(m); bread <- solve(crossprod(X))
  sqrt(diag(bread %*% crossprod(X, X * (resid(m)^2 / (1 - hatvalues(m))^2)) %*% bread))
}
# ---- G: gradiente general (una brecha por persona) ------------------------------------------------
contraste_G <- function(b) {
  b |> filter(!is.na(gap)) |> group_by(participant, cohort, task, variable, network) |>
    summarise(gap = mean(gap), n_conditions = n(), .groups = "drop") |>
    group_by(task, variable, network) |> group_modify(~ {
      tt <- t.test(.x$gap)
      tibble(effect = "core_minus_periphery", n = nrow(.x), n_HC = sum(.x$cohort == "HC"), n_PD = sum(.x$cohort == "PD"),
        estimate = unname(tt$estimate), conf_low = tt$conf.int[1], conf_high = tt$conf.int[2], df = unname(tt$parameter), p = tt$p.value)
    }) |> ungroup() |> mutate(family = familia(variable, "G_gradiente"))
}
# ---- C y S: ON−OFF y periodo en la brecha (21 pares, D ~ z) -----------------------------------------
contraste_CS <- function(b) {
  pares <- b |> filter(cohort == "PD") |> select(participant, order, task, variable, network, condition, gap) |>
    pivot_wider(names_from = condition, values_from = gap) |>
    mutate(D = ON - OFF, z = if_else(order == "OFF_ON", 1, -1)) |> filter(!is.na(D))
  pares |> group_by(task, variable, network) |> group_modify(~ {
    m <- lm(D ~ z, data = .x); b0 <- coef(m); se <- sqrt(diag(vcov(m))); df <- df.residual(m); s3 <- hc3(m)
    tibble(effect = c("ON_minus_OFF", "session2_minus_session1"), n = nrow(.x), estimate = unname(b0), SE = unname(se),
      conf_low = unname(b0 - qt(.975, df) * se), conf_high = unname(b0 + qt(.975, df) * se), df = df,
      p = unname(2 * pt(-abs(b0 / se), df)), p_HC3 = unname(2 * pt(-abs(b0 / s3), df)))
  }) |> ungroup() |>
    mutate(family = familia(variable, if_else(effect == "ON_minus_OFF", "C_ON_OFF", "S_periodo")))
}
# ---- H: HC frente a PD-OFF y PD-ON (Welch, una observación por persona) ------------------------------
contraste_H <- function(b, solo_primera = FALSE) {
  b <- b |> filter(!is.na(gap))
  if (solo_primera) b <- b |> filter(cohort == "HC" | session == 1L)
  bind_rows(lapply(c("OFF", "ON"), function(cond) {
    b |> filter(condition %in% c("HC", cond)) |> group_by(task, variable, network) |> group_modify(~ {
      x <- .x$gap[.x$condition == "HC"]; y <- .x$gap[.x$condition == cond]
      if (length(y) < 3L) return(tibble())
      tt <- t.test(x, y)
      tibble(effect = paste0("HC_minus_", cond), n = length(x) + length(y), n_HC = length(x), n_PD = length(y),
        estimate = unname(diff(rev(tt$estimate))), conf_low = tt$conf.int[1], conf_high = tt$conf.int[2],
        df = unname(tt$parameter), p = tt$p.value)
    }) |> ungroup()
  })) |> mutate(family = familia(variable, "H_HC_PD"))
}
bh <- function(r) r |> group_by(family) |> mutate(p_BH = p.adjust(p, "BH"), m_family = n()) |> ungroup()
analizar <- function(datos, ref_G = "conjunta_3", ref_C = "conjunta_PD", modo = "binaria", solo_primera = FALSE) {
  bG <- brechas(datos, ref_G, modo); bC <- if (ref_C == ref_G) bG else brechas(datos, ref_C, modo)
  r <- bind_rows(contraste_G(bG) |> mutate(reference = ref_G), contraste_CS(bC) |> mutate(reference = ref_C),
    contraste_H(bG, solo_primera) |> mutate(reference = ref_G))
  list(resultados = bh(r), brechas = bind_rows(bG, if (ref_C != ref_G) bC))
}
datos <- escenario_datos(base)
main <- analizar(datos)
resultados <- main$resultados
n_fam <- resultados |> count(family)
esperado <- c(G_gradiente = 26L, C_ON_OFF = 26L, S_periodo = 26L, H_HC_PD = 52L,
  AoA_G_gradiente = 6L, AoA_C_ON_OFF = 6L, AoA_S_periodo = 6L, AoA_H_HC_PD = 12L)
stopifnot(setequal(n_fam$family, names(esperado)), all(n_fam$n[match(names(esperado), n_fam$family)] == esperado))
# Calidad de la partición de referencia (19): soporte frente a «sin core» y estabilidad entre cadenas.
calidad <- res19 |> filter(family %in% c("conjunta_3", "conjunta_PD")) |>
  transmute(reference = family, task, network = redes[network], ref_nodes = nodos, ref_core = core_consenso,
    ref_conc_entre_cadenas = conc_entre_cadenas, ref_nodos_inestables = nodos_inestables, ref_soporte_core = soporte_core,
    ref_modelo_ganador = modelo_ganador)
resultados <- resultados |> left_join(calidad, by = c("reference", "task", "network")) |>
  mutate(flag_particion_debil = !ref_soporte_core | ref_conc_entre_cadenas < 0.9)
# ---- Sensibilidades preespecificadas ------------------------------------------------------------------
sens <- list(
  propia = analizar(datos, "propia", "propia")$resultados,
  estables = analizar(datos, modo = "estables")$resultados,
  continua = analizar(datos, modo = "continua")$resultados,
  sin_repeticiones = analizar(escenario_datos(base, "deduplicated"))$resultados,
  H_primera_sesion = analizar(datos, solo_primera = TRUE)$resultados |> filter(grepl("H_HC_PD", family)))
sensibilidades <- bind_rows(lapply(names(sens), function(s) sens[[s]] |> mutate(sensitivity = s)))
comparacion <- sensibilidades |>
  select(sensitivity, family, task, variable, network, effect, n_sens = n, estimate_sens = estimate, p_BH_sens = p_BH) |>
  left_join(resultados |> select(family, task, variable, network, effect, n_main = n, estimate_main = estimate, p_BH_main = p_BH),
    by = c("family", "task", "variable", "network", "effect")) |>
  mutate(same_sign = sign(estimate_sens) == sign(estimate_main),
    same_decision_BH05 = (p_BH_sens < .05) == (p_BH_main < .05))
# ---- Descriptivos por partición (sustituyen a la Tabla 3): palabra, token y persona -------------------
desc <- bind_rows(lapply(names(redes), function(red) {
  d <- anotar_particion(datos, "conjunta_3", red) |>
    pivot_longer(all_of(variables_psico), names_to = "variable", values_to = "value") |>
    filter(variable != "motor" | task == "actions", !is.na(value)) |>
    mutate(partition = if_else(hs_consenso == 0, "core", "periphery"), network = redes[[red]])
  bind_rows(
    d |> distinct(task, network, variable, partition, word_norm, value) |> group_by(task, network, variable, partition) |>
      summarise(level = "word_type", n = n(), mean = mean(value), sd = sd(value), .groups = "drop"),
    d |> group_by(task, network, variable, partition) |>
      summarise(level = "token", n = n(), mean = mean(value), sd = sd(value), .groups = "drop"),
    d |> group_by(task, network, variable, partition, condition, participant) |> summarise(m = mean(value), .groups = "drop") |>
      group_by(task, network, variable, partition, condition) |>
      summarise(level = "person_mean", n = n(), mean = mean(m), sd = sd(m), .groups = "drop"))
}))
cobertura <- main$brechas |> group_by(reference, network, task, variable, condition) |>
  summarise(n_people = n(), n_with_gap = sum(!is.na(gap)), median_core_tokens = median(n_core),
    median_periphery_tokens = median(n_periphery), .groups = "drop")
# ---- Figura: gradiente G (co-ocurrencia y secuencia) --------------------------------------------------
out <- nuevo_directorio("20_tabla2_core_periferia")
dibujar <- function() {
  par(mfrow = c(1, 3), mar = c(4, 8, 3, 1), oma = c(2.5, 0, 0, 0))
  g <- resultados |> filter(family == "G_gradiente")
  for (tarea in c("F", "actions", "animals")) {
    d <- g |> filter(task == tarea) |> arrange(variable, network)
    # Escala común: brecha en unidades de DE entre personas de la propia brecha no tiene sentido; se estandariza por
    # la DE de la variable entre tipos de palabra de la red conjunta para poder comparar variables.
    sdw <- desc |> filter(level == "word_type", task == tarea) |> group_by(variable, network) |>
      summarise(s = sqrt(mean(sd^2)), .groups = "drop")
    d <- d |> left_join(sdw, by = c("variable", "network")) |>
      mutate(est = estimate / s, lo = conf_low / s, hi = conf_high / s)
    y <- seq_len(nrow(d))
    titulo <- c(F = "Phonological (F)", actions = "Actions", animals = "Animals")[[tarea]]
    plot(range(c(0, d$lo, d$hi)), c(.5, nrow(d) + .5), type = "n", yaxt = "n", xlab = "Core − periphery (SD of word values)",
      ylab = "", main = titulo, bty = "l")
    etiqueta <- c(logFreq = "Log frequency", imageability = "Imageability", concreteness = "Concreteness",
      familiarity = "Familiarity", motor = "Motor content")[d$variable]
    axis(2, y, paste(etiqueta, ifelse(d$network == "cooccurrence", "(co-occ.)", "(sequence)")), las = 1, cex.axis = .75)
    abline(v = 0, lty = 3, col = "grey60")
    col <- ifelse(d$flag_particion_debil, "#B0B0B0", "#244B5A")
    segments(d$lo, y, d$hi, y, col = col, lwd = 2); points(d$est, y, pch = ifelse(d$p_BH < .05, 19, 1), col = col)
  }
  mtext("Within-person gap (45 participants), joint HC+ON+OFF partition · 95% CI · filled: BH-adjusted p < .05 (26 tests) · grey: weak partition",
    side = 1, outer = TRUE, cex = .75)
}
png(file.path(out, "gradiente_core_periferia.png"), width = 1800, height = 750, res = 140); dibujar(); dev.off()
svg(file.path(out, "gradiente_core_periferia.svg"), width = 12.85, height = 5.35); dibujar(); dev.off()
# Borrador de Tabla 2 (formato ancho): estimación [IC95 %] y p BH de cada familia; AoA marcada como exploratoria.
fmt <- function(e, l, h, p) sprintf("%.2f [%.2f, %.2f]; p_BH=%s", e, l, h, formatC(p, format = "g", digits = 2))
tabla2 <- resultados |> filter(effect != "session2_minus_session1") |>
  mutate(celda = fmt(estimate, conf_low, conf_high, p_BH),
    columna = recode(effect, core_minus_periphery = "Core_minus_periphery_all45", ON_minus_OFF = "ON_minus_OFF_gap_21pairs",
      HC_minus_OFF = "HC_minus_OFF_gap", HC_minus_ON = "HC_minus_ON_gap")) |>
  select(task, variable, network, columna, celda) |>
  pivot_wider(names_from = columna, values_from = celda) |>
  mutate(exploratory = variable == "AoA") |>
  arrange(exploratory, factor(task, c("F", "actions", "animals")), variable, network)
guardar_csv(tabla2, file.path(out, "tabla2_borrador.csv"))
guardar_csv(resultados, file.path(out, "resultados.csv"))
guardar_csv(sensibilidades, file.path(out, "sensibilidades.csv"))
guardar_csv(comparacion, file.path(out, "comparacion_sensibilidades.csv"))
guardar_csv(main$brechas, file.path(out, "brechas_persona.csv"))
guardar_csv(desc, file.path(out, "descriptivos_particion.csv"))
guardar_csv(cobertura, file.path(out, "cobertura.csv"))
registrar_entorno(out, c(entrada, file.path(d19, "particiones_consenso.csv"), file.path(d19, "resumen_redes.csv"),
  "analisis/R/20_tabla2_core_periferia.R", "analisis/R/00_sensibilidad_datos.R", "analisis/R/00_utilidades.R",
  "notas/plan_tabla2_R3_20261022.md"))
writeLines(c("# Nueva Tabla 2 core-periphery — plan preespecificado D-023 (notas/plan_tabla2_R3_20261022.md)",
  "Unidad: persona. Brecha = media de la variable en tokens de palabras core − en tokens periféricos (valores disponibles).",
  "Partición de referencia común: SBM conjunta HC+ON+OFF (G, H) y ON+OFF (C, S); consenso de 10 cadenas (script 19).",
  "G: t de una muestra sobre 45 personas (PD: media de sus brechas ON/OFF). C/S: D=ON−OFF ~ z, 21 pares. H: Welch HC frente a OFF/ON.",
  "BH dentro de cada familia (G26, C26, S26, H52; AoA aparte 6/6/6/12). flag_particion_debil: sin soporte frente a «sin core» o concordancia entre cadenas < 0,9.",
  "Sensibilidades: partición propia (circular, estilo R2), solo nodos estables, P(core) continua, sin repeticiones y H con primera sesión.",
  "Descriptivos por partición en tres niveles (tipo de palabra, token y media por persona); solo la persona se usa para inferencia.",
  "C no es efecto causal de la medicación (orden no aleatorizado documentado, P14). AoA provisional (P10). HC provisional (P11)."),
  file.path(out, "LEEME.md"), useBytes = TRUE)
print(resultados |> filter(!grepl("^AoA", family)) |> group_by(family) |>
  summarise(m = n(), sig_BH05 = sum(p_BH < .05), positivos = sum(estimate > 0), debiles = sum(flag_particion_debil)))
cat("Salida:", out, "\n")
