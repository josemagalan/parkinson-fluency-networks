# Propósito: gradiente periferia−core estimado en cada grupo (HC, OFF, ON) con los mismos modelos de la Tabla 2 (21).
#   En 21 el gradiente se estima en HC (columna Periphery) y en los pacientes solo se contrasta su diferencia con HC
#   (interacciones). Para afirmar que el patrón está «presente en todos los grupos» hace falta el gradiente de cada grupo:
#   Periphery + OFF_x_Periphery y Periphery + ON_x_Periphery, con la matriz de covarianzas completa y los gl mínimos de
#   los términos implicados (misma regla que los contrastes derivados de 21). BH por columna en los 22 modelos estimables.
#   Comprobación añadida el 24-10-2026 (D-031) para sustentar la carta/tabla de estado; no cambia la Tabla 2.
# Entradas: output/21_tabla2_niveles_R2_v2/{niveles_persona,resultados}.csv.
# Salidas: output/48_gradiente_por_grupo[_vN]/ (gradientes, resumen, verificación frente a 21, LEEME).
# E2/E5/E6. Paquetes: 00_utilidades.R + nlme (vía 00_niveles_tabla2.R).
# Ejecutar desde revision_R3/: Rscript --vanilla analisis/R/48_gradiente_por_grupo.R
# --- English ---
# Purpose: periphery-minus-core gradient within HC, OFF and ON from the script 21 models, full covariance,
#   BH per column over the 22 estimable models (D-031) (E2, E5, E6).
# Inputs: 21_tabla2_niveles_R2_v2/ CSVs. Outputs: analisis/output/48_gradiente_por_grupo[_vN]/.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
source("analisis/R/00_niveles_tabla2.R", encoding = "UTF-8")
set.seed(20261024)
fnv <- "analisis/output/21_tabla2_niveles_R2_v2/niveles_persona.csv"
f21 <- "analisis/output/21_tabla2_niveles_R2_v2/resultados.csv"
nv <- read_csv(fnv, show_col_types = FALSE)
r21 <- read_csv(f21, show_col_types = FALSE)
# Contrastes añadidos a los de 00_niveles_tabla2.R (ajustar_uno usa la lista global `contrastes`).
contrastes <- c(contrastes, list(
  Gradient_OFF = c(partitionperiphery = 1, `groupOFF:partitionperiphery` = 1),
  Gradient_ON = c(partitionperiphery = 1, `groupON:partitionperiphery` = 1)))
res <- nv |> group_by(task, network, variable) |> group_modify(~ ajustar_uno(.x)) |> ungroup()
no_id <- r21 |> distinct(task, network, variable, no_identificable)
ver <- res |> inner_join(r21 |> select(task, network, variable, term, estimate_21 = estimate), by = c("task", "network", "variable", "term")) |>
  mutate(dif = abs(estimate - estimate_21))
stopifnot(max(ver$dif, na.rm = TRUE) < 1e-6)
grad <- res |> filter(term %in% c("Periphery", "Gradient_OFF", "Gradient_ON")) |>
  mutate(term = recode(term, Periphery = "Gradient_HC")) |>
  left_join(no_id, by = c("task", "network", "variable")) |>
  mutate(family = paste0(if_else(variable == "AoA", "AoA_", ""), term)) |>
  group_by(family) |> mutate(p_BH = if_else(no_identificable, NA_real_, p.adjust(if_else(no_identificable, NA_real_, p), "BH"))) |> ungroup()
resumen <- grad |> filter(!no_identificable) |> group_by(family) |>
  summarise(m = n(), negativos = sum(estimate < 0), positivos = sum(estimate > 0), sig_BH = sum(p_BH < .05), .groups = "drop")
out <- nuevo_directorio("48_gradiente_por_grupo")
guardar_csv(grad, file.path(out, "gradientes.csv"))
guardar_csv(resumen, file.path(out, "resumen.csv"))
guardar_csv(ver, file.path(out, "verificacion_21.csv"))
registrar_entorno(out, c(fnv, f21, "analisis/R/48_gradiente_por_grupo.R", "analisis/R/00_niveles_tabla2.R", "analisis/R/00_utilidades.R"))
writeLines(c("# 48: gradiente periferia−core por grupo (D-031)",
  "Mismos modelos que 21 (media ~ grupo × partición + sesión; 1 | participante/lista). Gradient_HC = Periphery de 21;",
  "Gradient_OFF/ON = Periphery + interacción, con covarianza completa y gl mínimos. BH por columna (22 modelos; AoA aparte, 5).",
  sprintf("Verificación: los términos de 21 se reproducen (diferencia máxima %.1e).", max(ver$dif, na.rm = TRUE)),
  "F secuencia queda fuera (particiones de pacientes no estables). Para AoA el signo esperado es positivo (periferia más tardía).",
  paste(capture.output(print(as.data.frame(resumen))), collapse = "\n")), file.path(out, "LEEME.md"), useBytes = TRUE)
print(resumen); cat("Salida:", out, "\n")
