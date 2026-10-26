# Propósito: medias por persona ON/OFF con periodo, cobertura y sensibilidades P3/P10.
# Entradas: output/03_orden_sesiones/tokens.rds y diccionarios históricos.
# Salidas: output/06_psicolinguistica_sesion[/_vN]/, CSV/RDS y figura PNG/SVG.
# E2/E4/E5; 13 resultados principales × 2 efectos: BH26. AoA exploratoria BH6.
# Paquetes: 00_utilidades.R y 00_sensibilidad_datos.R; sin modelos de partición ni modificación de manuscrito.
# --- English ---
# Purpose: per-person psycholinguistic means, paired ON-OFF models with period, BH26 (AoA exploratory,
#   BH6) and sensitivities P3/P10 (E2, E4, E5).
# Inputs: 03 tokens.rds; dictionaries. Outputs: analisis/output/06_psicolinguistica_sesion[_vN]/ (CSV, RDS, figure).
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
source("analisis/R/00_sensibilidad_datos.R", encoding = "UTF-8")
set.seed(20261022)
entrada <- "analisis/output/03_orden_sesiones/tokens.rds"
base <- readRDS(entrada)
definir_medias <- function(d) d |>
  pivot_longer(all_of(variables_psico), names_to = "variable", values_to = "value") |>
  filter(variable != "motor" | task == "actions") |>
  group_by(participant, cohort, condition, session, order, task, variable) |>
  summarise(n_tokens = n(), n_available = sum(!is.na(value)), coverage = mean(!is.na(value)),
    mean = if (all(is.na(value))) NA_real_ else mean(value, na.rm = TRUE), .groups = "drop")
ajustar <- function(medias) {
  pares <- medias |> filter(cohort == "PD") |> select(participant, task, order, variable, condition, mean) |>
    pivot_wider(names_from = condition, values_from = mean) |>
    mutate(complete = !is.na(ON) & !is.na(OFF), D = ON - OFF, z = if_else(order == "OFF_ON", 1, -1))
  resultados <- robustos <- list(); modelos <- list()
  for (tarea in c("F", "actions", "animals")) for (var in variables_psico) {
    if (var == "motor" && tarea != "actions") next
    d <- pares |> filter(task == tarea, variable == var, complete)
    stopifnot(nrow(d) >= 10L, n_distinct(d$z) == 2L)
    m <- lm(D ~ z, data = d); b <- coef(m); se <- sqrt(diag(vcov(m))); df <- df.residual(m)
    X <- model.matrix(m); bread <- solve(crossprod(X))
    V3 <- bread %*% crossprod(X, X * (resid(m)^2 / (1 - hatvalues(m))^2)) %*% bread
    s3 <- sqrt(diag(V3)); key <- paste(tarea, var, sep = "_")
    for (j in 1:2) {
      r <- tibble(task = tarea, variable = var, effect = c("ON_minus_OFF", "session2_minus_session1")[j],
        n_pairs = nrow(d), df = df, estimate = b[j], SE = se[j],
        conf_low = b[j] - qt(.975, df) * se[j], conf_high = b[j] + qt(.975, df) * se[j],
        p = 2 * pt(-abs(b[j] / se[j]), df), SE_HC3 = s3[j], p_HC3 = 2 * pt(-abs(b[j] / s3[j]), df),
        family = if_else(var == "AoA", "AoA_exploratory_6", "psycholinguistic_26"))
      resultados[[paste(key, j)]] <- r
    }
    modelos[[key]] <- m
  }
  result <- bind_rows(resultados) |> group_by(family) |> mutate(p_BH = p.adjust(p, "BH"), p_BH_HC3 = p.adjust(p_HC3, "BH")) |> ungroup()
  stopifnot(sum(result$family == "psycholinguistic_26") == 26L, sum(result$family == "AoA_exploratory_6") == 6L)
  list(pares = pares, resultados = result, modelos = modelos)
}
medias <- definir_medias(escenario_datos(base))
main <- ajustar(medias)
descriptivos <- medias |> group_by(cohort, condition, task, variable) |>
  summarise(n_people = n(), n_with_mean = sum(!is.na(mean)), mean_of_means = mean(mean, na.rm = TRUE),
    sd_of_means = sd(mean, na.rm = TRUE), mean_coverage = mean(coverage), min_coverage = min(coverage), .groups = "drop")
incompletos <- main$pares |> filter(!complete)
sens_unique <- ajustar(definir_medias(escenario_datos(base, "deduplicated")))$resultados |>
  mutate(scenario = "within_list_deduplicated")
# AoA alternativa: una puntuación por palabra, promedio no ponderado entre diccionarios disponibles.
dic <- leer_diccionarios()
alt <- bind_rows(dic$aoa_verbs, dic$aoa_nouns) |> mutate(palabra = normalizar(palabra)) |>
  group_by(palabra) |> summarise(AoA = mean(AoA, na.rm = TRUE), .groups = "drop")
alternativa <- base; alternativa$AoA <- alt$AoA[match(alternativa$word_norm, alt$palabra)]
sens_aoa <- ajustar(definir_medias(escenario_datos(alternativa)))$resultados |> filter(variable == "AoA") |>
  mutate(scenario = "AoA_mean_across_sources")
sensibilidades <- bind_rows(sens_unique, sens_aoa)
resultados <- main$resultados
comparacion <- sensibilidades |> select(task, variable, effect, scenario, estimate_sensitivity = estimate, p_BH_sensitivity = p_BH) |>
  left_join(resultados |> select(task, variable, effect, estimate_main = estimate, p_BH_main = p_BH), by = c("task", "variable", "effect"))
# Producción tras retirar repeticiones: comparar con 04, sin cambiar el resultado vigente.
cuentas <- bind_rows(lapply(c("recorded", "deduplicated"), function(s) escenario_datos(base, s) |>
  filter(cohort == "PD") |> count(participant, task, order, condition, name = "n") |>
  pivot_wider(names_from = condition, values_from = n) |>
  group_by(task) |> group_modify(~ {
    m <- lm(I(ON - OFF) ~ I(ifelse(order == "OFF_ON", 1, -1)), data = .x)
    b <- coef(m); se <- sqrt(diag(vcov(m))); df <- df.residual(m)
    tibble(effect = c("ON_minus_OFF", "session2_minus_session1"), estimate = b,
      conf_low = b - qt(.975, df) * se, conf_high = b + qt(.975, df) * se, p = 2 * pt(-abs(b / se), df))
  }) |> ungroup() |> mutate(scenario = s, p_BH6 = p.adjust(p, "BH"))))
out <- nuevo_directorio("06_psicolinguistica_sesion")
for (nm in c("medias", "descriptivos", "incompletos", "resultados", "sensibilidades", "comparacion", "cuentas")) guardar_csv(get(nm), file.path(out, paste0(nm, ".csv")))
guardar_csv(main$pares, file.path(out, "pares.csv")); saveRDS(main$modelos, file.path(out, "modelos.rds"))
dibujar <- function() {
  par(mfrow = c(1, 3), mar = c(4, 7.5, 3, 1), oma = c(2, 0, 0, 0))
  for (tarea in c("F", "actions", "animals")) {
    d <- resultados |> filter(task == tarea, effect == "ON_minus_OFF", variable != "AoA")
    plot(range(c(0, d$conf_low, d$conf_high)), c(.5, nrow(d) + .5), type = "n", yaxt = "n",
      xlab = "ON − OFF ajustado", ylab = "", main = tarea, bty = "l")
    axis(2, seq_len(nrow(d)), d$variable, las = 1, cex.axis = .8)
    abline(v = 0, lty = 3, col = "grey60")
    segments(d$conf_low, seq_len(nrow(d)), d$conf_high, seq_len(nrow(d)), col = "#244B5A", lwd = 2)
    points(d$estimate, seq_len(nrow(d)), pch = 19, col = "#244B5A")
  }
  mtext("Medias por persona · IC95 % marginales · escalas originales de cada variable · AoA aparte", side = 1, outer = TRUE, cex = .8)
}
png(file.path(out, "psicolinguistica_ON_OFF.png"), width = 1800, height = 650, res = 140); dibujar(); dev.off()
svg(file.path(out, "psicolinguistica_ON_OFF.svg"), width = 12.85, height = 4.65); dibujar(); dev.off()
registrar_entorno(out, c(entrada, "analisis/R/06_psicolinguistica_sesion.R", "analisis/R/00_sensibilidad_datos.R",
  "analisis/R/00_utilidades.R", list.files("../parkinson/data", pattern = "\\.xlsx$", full.names = TRUE)))
writeLines(c("# Medias psicolingüísticas — 22-10-2026",
  "Una media de valores disponibles por persona/condición/tarea; una diferencia ON−OFF por paciente.",
  "D~z: intercepto condición y pendiente periodo; pares completos por variable, sin imputación, informar N/cobertura.",
  "Familia principal BH26 (13 resultados × condición/periodo); AoA provisional en BH6 exploratoria separada.",
  "HC descriptivo. Estos modelos no contienen partición y NO sustituyen la futura Tabla 2.",
  "Sensibilidades: HC3, retirar repeticiones dentro de lista y promediar AoA entre diccionarios.",
  "Puntuaciones solo de palabras cubiertas: cobertura incompleta puede sesgar la comparación; medias faltantes no son ceros.",
  "P3/P10 pendientes; datos y manuscrito intactos."), file.path(out, "LEEME.md"), useBytes = TRUE)
print(resultados |> filter(effect == "ON_minus_OFF"), n = 16, width = 140)
print(incompletos); cat("Salida:", out, "\n")
