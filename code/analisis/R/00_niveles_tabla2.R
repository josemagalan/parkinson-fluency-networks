# Propósito: funciones comunes de la Tabla 2 por persona con la estructura de la R2 (D-024), usadas por 21, 22, 45, 47 y 48.
# Entradas: medias por persona/lista/partición (columnas participant, cohort, condition, session, partition, mean).
# Salidas: funciones; no escribe ni consume RNG. Paquetes: nlme, dplyr/tibble (00_utilidades.R).
# --- English ---
# Purpose: shared Table 2 helpers (R2 structure, D-024): nlme model mean ~ group*partition + session,
#   random intercepts participant/list (fallback: participant), and derived linear contrasts with CIs.
# Inputs: per-person/list/partition means. Outputs: functions only (no files, no RNG). Used by 21, 22, 45, 47, 48.
suppressPackageStartupMessages(library(nlme))
# ---- Modelo mixto con la estructura de la R2 ----------------------------------------------------------
terminos <- c(OFF = "groupOFF", ON = "groupON", Periphery = "partitionperiphery",
  OFF_x_Periphery = "groupOFF:partitionperiphery", ON_x_Periphery = "groupON:partitionperiphery", Session = "ses")
contrastes <- list(OFF_minus_HC_periphery = c(groupOFF = 1, `groupOFF:partitionperiphery` = 1),
  ON_minus_HC_periphery = c(groupON = 1, `groupON:partitionperiphery` = 1),
  ON_minus_OFF_core = c(groupON = 1, groupOFF = -1),
  ON_minus_OFF_periphery = c(groupON = 1, `groupON:partitionperiphery` = 1, groupOFF = -1, `groupOFF:partitionperiphery` = -1))
ajustar_uno <- function(d, aleatorio = "anidado") {
  d <- d |> mutate(group = factor(condition, levels = c("HC", "OFF", "ON")), partition = factor(partition, c("core", "periphery")),
    ses = if_else(cohort == "HC", 0, session - 1), lista = paste(participant, condition, sep = ":"))
  f_ran <- if (aleatorio == "anidado") ~ 1 | participant / lista else ~ 1 | participant
  ctrl <- lmeControl(maxIter = 200, msMaxIter = 200, niterEM = 50, returnObject = FALSE)
  ajuste <- function(f) tryCatch(lme(mean ~ group * partition + ses, random = f, data = d, method = "REML", control = ctrl),
    error = function(e) NULL)
  m <- ajuste(f_ran); usado <- aleatorio
  if (is.null(m) && aleatorio == "anidado") { m <- ajuste(~ 1 | participant); usado <- "participante (anidado no convergió)" }
  # Sin ajuste posible (p. ej., una celda grupo × partición vacía en una sensibilidad): fila «no estimable».
  if (is.null(m) || !all(terminos %in% names(fixef(m)))) {
    return(tibble(term = c("Intercept_HC_core", names(terminos), names(contrastes)), estimate = NA_real_, SE = NA_real_,
      df = NA_real_, p = NA_real_, conf_low = NA_real_, conf_high = NA_real_, random_effects = "no estimable",
      n_obs = nrow(d), n_people = n_distinct(d$participant), n_HC = NA_integer_, n_OFF = NA_integer_, n_ON = NA_integer_))
  }
  tt <- summary(m)$tTable; b <- fixef(m); V <- vcov(m)
  filas <- lapply(names(terminos), function(k) {
    j <- terminos[[k]]
    tibble(term = k, estimate = tt[j, "Value"], SE = tt[j, "Std.Error"], df = tt[j, "DF"], p = tt[j, "p-value"])
  })
  for (k in names(contrastes)) {
    w <- contrastes[[k]]; L <- setNames(numeric(length(b)), names(b)); L[names(w)] <- w
    est <- sum(L * b); se <- sqrt(drop(t(L) %*% V %*% L)); df <- min(tt[names(w), "DF"])
    filas[[length(filas) + 1L]] <- tibble(term = k, estimate = est, SE = se, df = df, p = 2 * pt(-abs(est / se), df))
  }
  bind_rows(tibble(term = "Intercept_HC_core", estimate = tt["(Intercept)", "Value"], SE = tt["(Intercept)", "Std.Error"],
      df = tt["(Intercept)", "DF"], p = NA_real_), bind_rows(filas)) |>
    mutate(conf_low = estimate - qt(.975, df) * SE, conf_high = estimate + qt(.975, df) * SE,
      random_effects = usado, n_obs = nrow(d), n_people = n_distinct(d$participant),
      n_HC = n_distinct(d$participant[d$condition == "HC"]), n_OFF = n_distinct(d$participant[d$condition == "OFF"]),
      n_ON = n_distinct(d$participant[d$condition == "ON"]))
}
# Nota: los grados de libertad son los de nlme (regla de contención). Grupo y sesión varían entre listas, por lo que
# reciben los gl del nivel lista (pocos, conservador); partición e interacciones, los del nivel de observación.
