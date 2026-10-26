# Propósito: Tabla 2 fiel a la R2 con la persona como unidad (D-024; plan notas/plan_tabla2_niveles_D024_20261022.md).
#   Misma pregunta y términos que la R2 (partición propia de cada grupo; OFF/ON frente a HC en el core, periferia−core
#   en HC e interacciones), estimados con medias por persona/lista/partición y un modelo mixto
#   media ~ grupo × partición + sesión, aleatorio 1 | participante/lista (nlme). Contrastes derivados (periferia
#   frente a HC, ON−OFF pareado), BH por columna, regla de particiones no identificables y sensibilidades.
# Entradas: output/03_orden_sesiones/tokens.rds; output/19_SBM_redes_limpias/{particiones_consenso,resumen_redes}.csv;
#   output/02_reproducir_R2/tabla2_antes_despues.csv (comparación con la R2).
# Salidas: output/21_tabla2_niveles_R2[_vN]/ (resultados, tabla formato R2, sensibilidades, comparación R2, niveles, LEEME).
# E2/E4/E5/E6. Paquetes: 00_utilidades.R + nlme (paquete recomendado de R), vía 00_niveles_tabla2.R. Análisis añadido tras D-023.
# Ejecutar desde revision_R3/: Rscript --vanilla analisis/R/21_tabla2_niveles_R2.R
# --- English ---
# Purpose: Table 2 with the R2 structure and the person as unit (D-024): lme mean ~ group x partition +
#   session, random participant/list, own-group SBM partitions, derived contrasts, BH per column (E2, E4, E5, E6).
# Inputs: 03 tokens.rds; 19 partitions; 02 R2 table. Outputs: analisis/output/21_tabla2_niveles_R2[_vN]/.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
source("analisis/R/00_sensibilidad_datos.R", encoding = "UTF-8")
set.seed(20261022)
entrada <- "analisis/output/03_orden_sesiones/tokens.rds"
d19 <- "analisis/output/19_SBM_redes_limpias"
f02 <- "analisis/output/02_reproducir_R2/tabla2_antes_despues.csv"
base <- readRDS(entrada)
stopifnot(nrow(base) == 2965L)
part <- read_csv(file.path(d19, "particiones_consenso.csv"), show_col_types = FALSE, trim_ws = FALSE) |> filter(family == "propia")
res19 <- read_csv(file.path(d19, "resumen_redes.csv"), show_col_types = FALSE) |> filter(family == "propia")
redes <- c(undirected = "cooccurrence", directed = "sequence")
# Regla D-024: modelo no identificable si la partición propia de algún grupo no tiene soporte frente a «sin core».
soporte <- res19 |> group_by(task, network) |> summarise(todos_con_soporte = all(soporte_core),
  grupos_sin_soporte = paste(group[!soporte_core], collapse = "+"), .groups = "drop") |> mutate(network = redes[network])
# ---- Medias por persona, lista y partición -----------------------------------------------------------
niveles <- function(datos, solo_estables = FALSE) {
  bind_rows(lapply(names(redes), function(red) {
    p <- part |> filter(network == red) |> select(task, group, word_norm, hs_consenso, estable)
    d <- datos |> left_join(p, by = c("task", "condition" = "group", "word_norm"), relationship = "many-to-one")
    stopifnot(!anyNA(d$hs_consenso))          # toda palabra producida está en la red propia de su grupo
    if (solo_estables) d <- d |> filter(estable)
    d |> pivot_longer(all_of(variables_psico), names_to = "variable", values_to = "value") |>
      filter(variable != "motor" | task == "actions", !is.na(value)) |>
      mutate(partition = if_else(hs_consenso == 0, "core", "periphery")) |>
      group_by(participant, cohort, condition, session, order, task, variable, partition) |>
      summarise(n_tokens = n(), mean = mean(value), .groups = "drop") |>
      mutate(network = redes[[red]])
  }))
}
source("analisis/R/00_niveles_tabla2.R", encoding = "UTF-8")
ajustar_todo <- function(nv, aleatorio = "anidado") {
  nv |> group_by(task, network, variable) |> group_modify(~ ajustar_uno(.x, aleatorio)) |> ungroup() |>
    left_join(soporte, by = c("task", "network")) |>
    mutate(no_identificable = !todos_con_soporte,
      family = paste0(if_else(variable == "AoA", "AoA_", ""), term)) |>
    group_by(family) |>
    mutate(p_BH = if_else(no_identificable | term == "Intercept_HC_core", NA_real_,
      p.adjust(if_else(no_identificable, NA_real_, p), "BH"))) |> ungroup()
}
datos <- escenario_datos(base)
nv <- niveles(datos)
resultados <- ajustar_todo(nv)
fam <- resultados |> filter(term != "Intercept_HC_core") |> group_by(family) |> summarise(m = sum(!no_identificable), .groups = "drop")
stopifnot(all(fam$m[!grepl("^AoA_", fam$family)] == 22L), all(fam$m[grepl("^AoA_", fam$family)] == 5L),
  nrow(fam) == 20L)
# ---- Sensibilidades -----------------------------------------------------------------------------------
hc3 <- function(m) { X <- model.matrix(m); br <- solve(crossprod(X))
  sqrt(diag(br %*% crossprod(X, X * (resid(m)^2 / (1 - hatvalues(m))^2)) %*% br)) }
simples <- function(nv, solo_primera = FALSE) {
  welch <- bind_rows(lapply(c("OFF", "ON"), function(cond) {
    nv |> filter(condition %in% c("HC", cond), cohort == "HC" | !solo_primera | session == 1L) |>
      group_by(task, network, variable, partition) |> group_modify(~ {
        x <- .x$mean[.x$condition == "HC"]; y <- .x$mean[.x$condition == cond]
        if (length(y) < 3L || length(x) < 3L) return(tibble())
        tt <- t.test(y, x)
        tibble(term = paste0(cond, "_minus_HC_", .y$partition), estimate = unname(diff(rev(tt$estimate))),
          conf_low = tt$conf.int[1], conf_high = tt$conf.int[2], p = tt$p.value, n = length(x) + length(y))
      }) |> ungroup()
  }))
  if (solo_primera) return(welch)
  pares <- nv |> filter(cohort == "PD") |> select(participant, order, task, network, variable, partition, condition, mean) |>
    pivot_wider(names_from = condition, values_from = mean) |> mutate(D = ON - OFF, z = if_else(order == "OFF_ON", 1, -1)) |>
    filter(!is.na(D))
  pareado <- pares |> group_by(task, network, variable, partition) |> group_modify(~ {
    m <- lm(D ~ z, data = .x); b <- coef(m)[1]; se <- sqrt(vcov(m)[1, 1]); df <- df.residual(m)
    tibble(term = paste0("ON_minus_OFF_", .y$partition), estimate = unname(b), conf_low = unname(b - qt(.975, df) * se),
      conf_high = unname(b + qt(.975, df) * se), p = unname(2 * pt(-abs(b / se), df)),
      p_HC3 = unname(2 * pt(-abs(b / hc3(m)[1]), df)), n = nrow(.x))
  }) |> ungroup()
  bind_rows(welch, pareado)
}
marcar <- function(r) r |> left_join(soporte, by = c("task", "network")) |>
  mutate(no_identificable = !todos_con_soporte, family = paste0(if_else(variable == "AoA", "AoA_", ""), term)) |>
  group_by(family) |> mutate(p_BH = if_else(no_identificable, NA_real_, p.adjust(if_else(no_identificable, NA_real_, p), "BH"))) |> ungroup()
sens <- list(
  contrastes_simples = marcar(simples(nv)),
  HC_vs_PD_primera_sesion = marcar(simples(nv, solo_primera = TRUE)),
  solo_nodos_estables = ajustar_todo(niveles(datos, solo_estables = TRUE)),
  sin_repeticiones = ajustar_todo(niveles(escenario_datos(base, "deduplicated"))),
  aleatorio_solo_participante = ajustar_todo(nv, "participante"))
sensibilidades <- bind_rows(lapply(names(sens), function(s) sens[[s]] |> mutate(sensitivity = s)))
# Equivalencias para comparar: contrastes simples frente a los términos/contrastes del modelo.
equiv <- c(OFF_minus_HC_core = "OFF", ON_minus_HC_core = "ON", OFF_minus_HC_periphery = "OFF_minus_HC_periphery",
  ON_minus_HC_periphery = "ON_minus_HC_periphery", ON_minus_OFF_core = "ON_minus_OFF_core", ON_minus_OFF_periphery = "ON_minus_OFF_periphery")
principal <- resultados |> filter(term != "Intercept_HC_core") |>
  select(task, network, variable, term, estimate_main = estimate, p_main = p, p_BH_main = p_BH, no_identificable)
comparacion <- sensibilidades |> filter(term != "Intercept_HC_core") |>
  mutate(term_main = if_else(term %in% names(equiv), unname(equiv[term]), term)) |>
  select(sensitivity, task, network, variable, term, term_main, estimate_sens = estimate, p_BH_sens = p_BH) |>
  inner_join(principal, by = c("task", "network", "variable", "term_main" = "term")) |>
  filter(!no_identificable) |>
  mutate(same_sign = sign(estimate_sens) == sign(estimate_main), same_decision_BH05 = (p_BH_sens < .05) == (p_BH_main < .05))
# ---- Comparación con la R2 (coeficientes históricos y con uniones corregidas, unidad palabra-tipo) -------
r2 <- read_csv(f02, show_col_types = FALSE) |>
  mutate(task = recode(test, Acciones = "actions", Animales = "animals", Fpalabra = "F"), network = redes[network],
    term = recode(term, groupoff = "OFF", groupon = "ON", hs_partition1 = "Periphery",
      `groupoff:hs_partition1` = "OFF_x_Periphery", `groupon:hs_partition1` = "ON_x_Periphery")) |>
  filter(term != "(Intercept)") |>
  select(task, network, variable = DependentVar, term, publicado_manuscrito, coef_R2, p_R2, coef_corregido, p_corregido)
comparacion_R2 <- r2 |> inner_join(resultados |> select(task, network, variable, term, estimate, conf_low, conf_high, p, p_BH, no_identificable),
    by = c("task", "network", "variable", "term")) |>
  mutate(sig_R2_nominal = p_R2 < .05, sig_R3_BH = !is.na(p_BH) & p_BH < .05, mismo_signo = sign(coef_R2) == sign(estimate),
    estado = case_when(no_identificable ~ "no identificable en R3 (F secuencia)",
      sig_R2_nominal & sig_R3_BH & mismo_signo ~ "se mantiene",
      sig_R2_nominal & !sig_R3_BH ~ "no se mantiene",
      !sig_R2_nominal & sig_R3_BH ~ "nuevo en R3",
      TRUE ~ "no significativo en ambas"))
resumen_R2 <- comparacion_R2 |> filter(publicado_manuscrito) |> count(term, estado) |>
  pivot_wider(names_from = estado, values_from = n, values_fill = 0)
# ---- Tabla en formato R2 ------------------------------------------------------------------------------
estrellas <- function(p) case_when(is.na(p) ~ "", p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", TRUE ~ "")
tabla2 <- resultados |>
  mutate(celda = if_else(term == "Intercept_HC_core", sprintf("%.2f", estimate),
      sprintf("%.2f [%.2f, %.2f]%s", estimate, conf_low, conf_high, estrellas(p_BH))),
    celda = if_else(no_identificable & term != "Intercept_HC_core", paste0(celda, " (n.i.)"), celda)) |>
  select(network, task, variable, term, celda, n_HC, n_OFF, n_ON) |>
  pivot_wider(names_from = term, values_from = celda) |>
  mutate(exploratory_AoA = variable == "AoA") |>
  arrange(exploratory_AoA, network, factor(task, c("actions", "F", "animals")), variable)
# Diagnósticos del modelo principal (revisión Codex de D-024): componentes de varianza, residuos y convergencia.
diagnosticos <- nv |> group_by(task, network, variable) |> group_modify(~ {
  d <- .x |> mutate(group = factor(condition, levels = c("HC", "OFF", "ON")), partition = factor(partition, c("core", "periphery")),
    ses = if_else(cohort == "HC", 0, session - 1), lista = paste(participant, condition, sep = ":"))
  m <- tryCatch(lme(mean ~ group * partition + ses, random = ~ 1 | participant / lista, data = d, method = "REML",
    control = lmeControl(maxIter = 200, msMaxIter = 200, niterEM = 50, returnObject = FALSE)), error = function(e) NULL)
  if (is.null(m)) return(tibble(converged_nested = FALSE))
  vc <- suppressWarnings(as.numeric(VarCorr(m)[, "StdDev"])); r <- resid(m, type = "normalized")
  tibble(converged_nested = TRUE, n_obs = nrow(d), n_participants = n_distinct(d$participant), n_lists = n_distinct(d$lista),
    sd_participant = vc[2], sd_list = vc[4], sd_residual = m$sigma, max_abs_norm_resid = max(abs(r)),
    n_abs_resid_gt3 = sum(abs(r) > 3), shapiro_p_resid = if (length(r) >= 3 && length(r) <= 5000) shapiro.test(r)$p.value else NA_real_)
}) |> ungroup()
cobertura <- nv |> group_by(network, task, variable, condition, partition) |>
  summarise(n_people = n(), median_tokens = median(n_tokens), min_tokens = min(n_tokens), .groups = "drop")
out <- nuevo_directorio("21_tabla2_niveles_R2")
for (nm in c("resultados", "tabla2", "sensibilidades", "comparacion", "comparacion_R2", "resumen_R2", "cobertura", "soporte", "diagnosticos")) {
  guardar_csv(get(nm), file.path(out, paste0(c(resultados = "resultados", tabla2 = "tabla2_formato_R2",
    sensibilidades = "sensibilidades", comparacion = "comparacion_sensibilidades", comparacion_R2 = "comparacion_R2",
    resumen_R2 = "resumen_comparacion_R2", cobertura = "cobertura", soporte = "soporte_particiones",
    diagnosticos = "diagnosticos_modelos")[[nm]], ".csv")))
}
guardar_csv(nv, file.path(out, "niveles_persona.csv"))
registrar_entorno(out, c(entrada, file.path(d19, "particiones_consenso.csv"), file.path(d19, "resumen_redes.csv"), f02,
  "analisis/R/21_tabla2_niveles_R2.R", "analisis/R/00_niveles_tabla2.R", "analisis/R/00_sensibilidad_datos.R", "analisis/R/00_utilidades.R",
  "notas/plan_tabla2_niveles_D024_20261022.md"))
writeLines(c("# Tabla 2 fiel a la R2 con la persona como unidad — D-024 (notas/plan_tabla2_niveles_D024_20261022.md)",
  "Análisis añadido tras conocer las brechas de D-023 (script 20), adoptado por Josema; detalles fijados antes de calcular.",
  "Partición propia de cada grupo (script 19, consenso de 10 cadenas, hub-and-spoke). Medias por persona/lista/partición.",
  "Modelo por tarea × red × variable: media ~ grupo × partición + sesión, aleatorio 1 | participante/lista (nlme, REML).",
  "Términos como la R2 (OFF, ON = diferencia con HC en el core; Periphery = periferia−core en HC; interacciones) más",
  "contrastes derivados: OFF/ON − HC en la periferia y ON − OFF en core y periferia (pareado dentro de paciente).",
  "BH por columna (22 modelos estimables; AoA aparte, 5). F secuencia: particiones ON/OFF sin soporte → no identificable (n.i.), fuera de BH.",
  "Estrellas de tabla2_formato_R2.csv con p BH; intercepto sin p. Comparación con la R2 en comparacion_R2.csv (unidad palabra-tipo en la R2).",
  "Sensibilidades: contrastes simples (Welch, D~z), HC frente a primera sesión, nodos estables, sin repeticiones, aleatorio solo participante.",
  "Contrastes derivados con la matriz de covarianzas completa; sus gl son el mínimo de los gl de nlme de los términos implicados.",
  "Gl de nlme por contención: grupo y sesión con los del nivel lista (conservadores); partición e interacciones, los del nivel observación.",
  "diagnosticos_modelos.csv: convergencia anidada, DE de participante/lista/residual y residuos normalizados (máximo, >3 y Shapiro).",
  "BH controla la FDR dentro de cada columna (familia), no globalmente entre columnas.",
  "E1 (niveles en redes igualadas): script 22. Brechas con referencia común: script 20 (complemento)."), file.path(out, "LEEME.md"), useBytes = TRUE)
print(resultados |> filter(term != "Intercept_HC_core", !no_identificable, !grepl("^AoA", family)) |>
  group_by(term) |> summarise(m = n(), sig_BH = sum(p_BH < .05), positivos = sum(estimate > 0)))
print(resumen_R2); cat("Salida:", out, "\n")
