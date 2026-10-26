# Propósito: sensibilidad HC a igual primera sesión y exclusión individual de modelos psicolingüísticos.
# Entradas: tokens con orden, modelos 06 y escenarios P3/P11 (solo hipótesis).
# Salidas: output/11_HC_primera_sesion_influencia[/_vN]/, CSV (hc, loo, loo_resumen) y LEEME.
# E2/E4/E5: Welch HC24 vs OFF11/ON10, BH32; AoA exploratoria BH6.
# Paquetes: 00_utilidades.R; sin cambios en datos ni manuscrito.
# --- English ---
# Purpose: HC vs first-session PD comparisons (Welch, BH32; AoA BH6) and leave-one-out influence
#   on the 06 paired models (E2, E4, E5).
# Inputs: 03 tokens.rds, 06 modelos.rds. Outputs: analisis/output/11_HC_primera_sesion_influencia[_vN]/.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
source("analisis/R/00_sensibilidad_datos.R", encoding = "UTF-8")
set.seed(20261022)
entrada <- "analisis/output/03_orden_sesiones/tokens.rds"; base <- readRDS(entrada)
crear <- function(d) {
  a <- d |> pivot_longer(all_of(variables_psico), names_to = "variable", values_to = "value") |>
    filter(variable != "motor" | task == "actions") |> group_by(participant, condition, session, task, variable) |>
    summarise(value = if (all(is.na(value))) NA_real_ else mean(value, na.rm = TRUE), .groups = "drop")
  b <- d |> count(participant, condition, session, task, name = "value") |> mutate(variable = "n_tokens")
  bind_rows(a, b)
}
hc <- list()
for (scenario in c("recorded", "deduplicated", "ambiguous_drop", "ambiguous_split")) {
  d <- crear(escenario_datos(base, scenario)) |> filter(session == 1L)
  hh <- list()
  for (task in c("F", "actions", "animals")) for (v in c(variables_psico, "n_tokens")) {
    if (v == "motor" && task != "actions") next
    for (cond in c("OFF", "ON")) {
      x <- d$value[d$task == task & d$variable == v & d$condition == "HC"]
      y <- d$value[d$task == task & d$variable == v & d$condition == cond]
      m <- t.test(x, y)
      hh[[paste(task, v, cond)]] <- tibble(task, variable = v, contrast = paste0("HC_minus_", cond),
        n_HC = sum(!is.na(x)), n_PD = sum(!is.na(y)), estimate = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE),
        conf_low = m$conf.int[1], conf_high = m$conf.int[2], p = m$p.value,
        family = if_else(v == "AoA", "AoA_exploratory_6", "first_session_HC_32"))
    }
  }
  hc[[scenario]] <- bind_rows(hh) |> group_by(family) |> mutate(p_BH = p.adjust(p, "BH")) |> ungroup() |> mutate(scenario = scenario)
}
hc <- bind_rows(hc)
modelos <- readRDS("analisis/output/06_psicolinguistica_sesion/modelos.rds")
loo <- bind_rows(lapply(names(modelos), function(nm) {
  m <- modelos[[nm]]; d <- model.frame(m)
  bind_rows(lapply(seq_len(nrow(d)), function(i) {
    b <- coef(lm(D ~ z, data = d[-i, ]))
    tibble(model = nm, excluded_row_in_model_frame = i, effect = c("ON_minus_OFF", "session2_minus_session1"), estimate = b)
  }))
}))
loo_resumen <- loo |> group_by(model, effect) |> summarise(min_estimate = min(estimate), max_estimate = max(estimate), .groups = "drop")
out <- nuevo_directorio("11_HC_primera_sesion_influencia")
for (nm in c("hc", "loo", "loo_resumen")) guardar_csv(get(nm), file.path(out, paste0(nm, ".csv")))
registrar_entorno(out, c(entrada, "analisis/output/06_psicolinguistica_sesion/modelos.rds",
  "analisis/R/11_HC_primera_sesion_influencia.R", "analisis/R/00_sensibilidad_datos.R"))
writeLines(c("# HC en primera sesión y sensibilidad de influencia",
  "HC24 vs PD_OFF11/PD_ON10: una observación por persona, solo primera evaluación. No es comparación pareada ON/OFF.",
  "Welch, BH32 para conteos+medias no AoA; AoA exploratoria BH6. IC95 marginales.",
  "Comparaciones entre cohortes, sin asumir asignación aleatoria ni ajustar covariables demográficas faltantes.",
  "No reemplaza modelos core-periphery. P3/P10/P11 pendientes; escenarios etiquetados como hipótesis.",
  "Exclusión de una persona por vez de modelos 06; las filas corresponden al frame completo guardado de cada modelo.",
  "No se selecciona un modelo principal según resultado ni se descartan personas automáticamente."), file.path(out, "LEEME.md"))
print(hc |> filter(scenario == "recorded", p_BH < .05), n = 38, width = 150); cat("Salida:", out, "\n")
