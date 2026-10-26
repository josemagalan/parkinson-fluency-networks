# Propósito: contrastes pareados de producción ON/OFF y periodo por tarea.
# Entradas: output/03_orden_sesiones/tokens.rds (2.965 tokens; 21 pares PD).
# Salidas: output/04_produccion_sesion[/_vN]/, CSV/RDS, diagnósticos y PNG/SVG.
# E1/E2/E4/E5: participante como unidad, control de periodo, BH de seis contrastes.
# Paquetes: los de 00_utilidades.R; gráficos con graphics/grDevices base.
# Ejecutar desde revision_R3/: Rscript --vanilla analisis/R/04_produccion_sesion.R
# --- English ---
# Purpose: paired ON-OFF word-production contrasts per task adjusted for period; BH over 6 contrasts,
#   HC3 and leave-one-out sensitivities (E1, E2, E4, E5).
# Inputs: 03 tokens.rds. Outputs: analisis/output/04_produccion_sesion[_vN]/ (CSV, RDS, PNG/SVG).
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
set.seed(20261022)
entrada <- "analisis/output/03_orden_sesiones/tokens.rds"
tokens <- readRDS(entrada)
stopifnot(nrow(tokens) == 2965L, !anyDuplicated(tokens$token_id))
conteos <- tokens |> count(participant, cohort, condition, session, order, task, name = "n_words")
stopifnot(nrow(conteos) == 198L, sum(conteos$n_words) == 2965L)
pacientes <- conteos |> filter(cohort == "PD") |>
  mutate(participant = factor(participant), condition = factor(condition, levels = c("OFF", "ON")),
    period2 = as.integer(session == 2L))
diferencias <- pacientes |> select(participant, task, order, condition, n_words) |>
  pivot_wider(names_from = condition, values_from = n_words) |>
  mutate(delta_ON_OFF = ON - OFF, z = if_else(order == "OFF_ON", 1, -1))
stopifnot(nrow(diferencias) == 63L, !anyNA(diferencias),
  all((diferencias |> count(task))$n == 21L),
  all((pacientes |> count(participant, task))$n == 2L))
descriptivos <- conteos |> group_by(cohort, condition, task) |>
  summarise(n = n(), mean = mean(n_words), sd = sd(n_words), median = median(n_words),
    min = min(n_words), max = max(n_words), .groups = "drop")
por_secuencia <- pacientes |> group_by(task, order, condition, session) |>
  summarise(n = n(), mean = mean(n_words), sd = sd(n_words), .groups = "drop")
delta_secuencia <- diferencias |> group_by(task, order) |>
  summarise(n = n(), mean_delta = mean(delta_ON_OFF), sd_delta = sd(delta_ON_OFF), .groups = "drop")
extraer_lm <- function(m, covariance = vcov(m)) {
  b <- coef(m); se <- sqrt(diag(covariance)); df <- df.residual(m)
  tibble(term = names(b), estimate = unname(b), SE = unname(se), df = df,
    conf_low = b - qt(.975, df) * se, conf_high = b + qt(.975, df) * se,
    p = 2 * pt(-abs(b / se), df))
}
resultados <- hc3 <- sin_ajuste <- diagnosticos <- residuos <- loo <- rango <- list()
modelos <- list()
for (tarea in c("F", "actions", "animals")) {
  d <- diferencias |> filter(task == tarea)
  p <- pacientes |> filter(task == tarea)
  m <- lm(delta_ON_OFF ~ z, data = d)
  completo <- lm(n_words ~ participant + condition + period2, data = p)
  esperado <- coef(m)
  observado <- coef(completo)[c("conditionON", "period2")]
  stopifnot(df.residual(m) == 19L, df.residual(completo) == 19L,
    max(abs(esperado - observado)) < 1e-10,
    max(abs(vcov(m) - vcov(completo)[c("conditionON", "period2"), c("conditionON", "period2")])) < 1e-10)
  # Comprobación independiente: media y semidiferencia de los dos grupos de secuencia.
  a <- mean(d$delta_ON_OFF[d$z == 1]); b <- mean(d$delta_ON_OFF[d$z == -1])
  stopifnot(abs(coef(m)[1] - (a + b) / 2) < 1e-10,
    abs(coef(m)[2] - (a - b) / 2) < 1e-10)
  resultados[[tarea]] <- extraer_lm(m) |> mutate(task = tarea)
  X <- model.matrix(m); bread <- solve(crossprod(X))
  # HC3: sensibilidad de varianza a heterocedasticidad, t con 19 gl aproximada.
  scale <- resid(m)^2 / (1 - hatvalues(m))^2
  V3 <- bread %*% crossprod(X, X * scale) %*% bread
  hc3[[tarea]] <- extraer_lm(m, V3) |> mutate(task = tarea)
  unadj <- lm(delta_ON_OFF ~ 1, data = d)
  sin_ajuste[[tarea]] <- extraer_lm(unadj) |> mutate(task = tarea)
  diagnosticos[[tarea]] <- tibble(task = tarea, n_pairs = nrow(d), df = df.residual(m),
    residual_sd = sigma(m), shapiro_W = unname(shapiro.test(resid(m))$statistic),
    shapiro_p = shapiro.test(resid(m))$p.value, max_abs_student = max(abs(rstudent(m))),
    max_cook = max(cooks.distance(m)), n_cook_over_4_n = sum(cooks.distance(m) > 4 / nrow(d)))
  residuos[[tarea]] <- d |> mutate(fitted = fitted(m), residual = resid(m),
    student = rstudent(m), cook = cooks.distance(m))
  loo[[tarea]] <- bind_rows(lapply(seq_len(nrow(d)), function(i) {
    extraer_lm(lm(delta_ON_OFF ~ z, data = d[-i, ])) |>
      mutate(task = tarea, excluded_participant = as.character(d$participant[i]))
  }))
  # El término condition:order no añade un contraste independiente del periodo
  # en presencia de efectos fijos de paciente: comprobar rango, no interpretarlo.
  xx <- model.matrix(~ participant + condition + period2 + condition:order, data = p)
  rango[[tarea]] <- tibble(task = tarea, n_columns = ncol(xx), rank = qr(xx)$rank,
    deficient = qr(xx)$rank < ncol(xx))
  modelos[[tarea]] <- list(differences = m, participant_fixed = completo, unadjusted = unadj)
}
etiquetar <- function(x) x |> mutate(effect = recode(term, `(Intercept)` = "ON_minus_OFF", z = "session2_minus_session1"))
resultados <- bind_rows(resultados) |> etiquetar() |> mutate(p_BH6 = p.adjust(p, "BH"), family = "6 production contrasts")
hc3 <- bind_rows(hc3) |> etiquetar() |> mutate(p_BH6_sensitivity = p.adjust(p, "BH"))
sin_ajuste <- bind_rows(sin_ajuste) |> mutate(p_BH3_secondary = p.adjust(p, "BH"))
loo <- bind_rows(loo) |> etiquetar()
diagnosticos <- bind_rows(diagnosticos); residuos <- bind_rows(residuos); rango <- bind_rows(rango)
stopifnot(nrow(resultados) == 6L, all(rango$deficient))
resumen_loo <- loo |> group_by(task, effect) |>
  summarise(min_estimate = min(estimate), max_estimate = max(estimate), min_p = min(p), max_p = max(p), .groups = "drop")
comparacion <- resultados |> select(task, effect, estimate, conf_low, conf_high, p, p_BH6) |>
  left_join(hc3 |> select(task, effect, HC3_SE = SE, HC3_p = p, HC3_p_BH6 = p_BH6_sensitivity), by = c("task", "effect")) |>
  left_join(resumen_loo, by = c("task", "effect"))
out <- nuevo_directorio("04_produccion_sesion")
for (nm in c("conteos", "diferencias", "descriptivos", "por_secuencia", "delta_secuencia", "resultados", "hc3", "sin_ajuste", "diagnosticos", "residuos", "loo", "rango", "comparacion")) {
  guardar_csv(get(nm), file.path(out, paste0(nm, ".csv")))
}
saveRDS(modelos, file.path(out, "modelos.rds"))
dibujar <- function() {
  par(mfrow = c(2, 3), mar = c(4.5, 6.3, 3.2, 1), oma = c(2.4, 0, 1, 0), family = "sans")
  titulos <- c(F = "Fonológica (F)", actions = "Acciones", animals = "Animales")
  for (tarea in names(titulos)) {
    d <- diferencias |> filter(task == tarea)
    plot(c(.8, 2.2), range(c(d$OFF, d$ON)), type = "n", xaxt = "n", xlab = "Condición",
      ylab = "Palabras por paciente", main = titulos[[tarea]], bty = "l")
    axis(1, 1:2, c("OFF", "ON"))
    cols <- ifelse(d$order == "OFF_ON", "#247BA0AA", "#E07A5FAA")
    for (i in seq_len(nrow(d))) lines(1:2, c(d$OFF[i], d$ON[i]), col = cols[i], type = "o", pch = 16, cex = .5)
    if (tarea == "F") legend("topleft", c("OFF primero", "ON primero"), col = c("#247BA0", "#E07A5F"), lty = 1, bty = "n", cex = .8)
  }
  for (tarea in names(titulos)) {
    d <- resultados |> filter(task == tarea)
    lim <- range(c(0, d$conf_low, d$conf_high)) + c(-.4, .4)
    plot(lim, c(.5, 2.5), type = "n", yaxt = "n", xlab = "Diferencia ajustada (palabras)", ylab = "", bty = "l")
    axis(2, 1:2, c("ON − OFF", "Sesión 2 − 1"), las = 1, cex.axis = .78)
    abline(v = 0, lty = 3, col = "grey60")
    segments(d$conf_low, 1:2, d$conf_high, 1:2, lwd = 2, col = "#244B5A")
    points(d$estimate, 1:2, pch = 19, col = "#244B5A")
    title(main = paste0("IC95 %; p BH: ", paste(format.pval(d$p_BH6, digits = 2, eps = .001), collapse = " / ")), cex.main = .9)
  }
  mtext("21 pacientes por tarea · IC marginales · BH sobre los seis contrastes", outer = TRUE, side = 1, line = .4, cex = .8)
}
png(file.path(out, "produccion_y_sesion.png"), width = 1800, height = 1100, res = 150)
dibujar(); dev.off()
svg(file.path(out, "produccion_y_sesion.svg"), width = 12, height = 7.33)
dibujar(); dev.off()
png(file.path(out, "diagnosticos.png"), width = 1500, height = 950, res = 140)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (tarea in names(modelos)) { qqnorm(resid(modelos[[tarea]]$differences), main = paste(tarea, "Q-Q")); qqline(resid(modelos[[tarea]]$differences)) }
for (tarea in names(modelos)) { d <- residuos |> filter(task == tarea); plot(d$z, d$residual, xlab = "z: ON segunda (+1) / primera (-1)", ylab = "Residuo", main = tarea, pch = 19); abline(h = 0, lty = 3) }
dev.off()
registrar_entorno(out, c(entrada, "analisis/R/04_produccion_sesion.R", "analisis/R/00_utilidades.R"))
writeLines(c(
  "# Producción ON/OFF con ajuste de sesión — 22-10-2026",
  "Unidad: una diferencia ON−OFF por paciente y tarea (21 pares; tres tareas analizadas por separado).",
  "Modelo principal D_i = beta_condicion + beta_periodo*z_i + error_i; z=+1 OFF_ON, -1 ON_OFF.",
  "Intercepto: ON−OFF ajustado; pendiente: sesión 2−1 ajustada. Pruebas t con 19 gl; IC95 % marginales.",
  "Equivalencia numérica comprobada con n ~ participant + condition + period2 y con medias por secuencia.",
  "BH aplicado a una familia de seis p principales (3 tareas × 2 contrastes). No fija familias de futura Tabla 2.",
  "Se cuentan tokens, sin eliminar repeticiones; la celda ambigua de control no entra en los modelos de pacientes.",
  "Los controles se presentan solo como descriptivos. No se mezclan con los pacientes en pruebas pareadas.",
  "Sensibilidades HC3 (p aproximada t19) y exclusión de un paciente por vez; no se elige modelo según p.",
  "El modelo estima diferencias de medias en palabras, no asume que conteos completos sean Poisson. Revisar diagnósticos de diferencias.",
  "condition×order y periodo no se pueden separar con estos efectos fijos de paciente: matriz deficiente verificada.",
  "Sesión puede incluir práctica, fatiga, tiempo y arrastre; no se identifica su mecanismo ni se descarta carryover.",
  "No se ha confirmado aleatorización de secuencia ni fechas/intervalo (P13); no presentar efectos como prueba causal de dopamina.",
  "Este bloque no sustituye E4 para resultados core-periphery ni iguala todavía producción de redes (E1).",
  "Fuentes metodológicas: https://www.bmj.com/content/366/bmj.l4378.full ; https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html",
  "No se han modificado fuentes, resultados previos ni manuscrito."
), file.path(out, "LEEME.md"), useBytes = TRUE)
print(resultados, width = 160); print(diagnosticos); cat("Salida:", out, "\n")
