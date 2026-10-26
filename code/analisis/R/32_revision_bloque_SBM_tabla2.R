# Propósito: revisar resultados 18–20, reproducción Windows y N efectivos de Tabla 2.
# Entradas: outputs 19, 20 (Linux) y 20_v2 (repetición Windows), sin modificarlos.
# Salidas: output/32_revision_bloque_SBM_tabla2[_vN]: comparación, N y tabla rotulada.
# E2/E4/E5/E6. No cambia particiones, tests ni familias de D-023; revisión condicional.
# Paquetes: 00_utilidades.R; sin paquetes nuevos. Ejecutar desde revision_R3/.
# --- English ---
# Purpose (optional): audit scripts 19-20: Linux run vs Windows re-run (20_v2), BH recomputation and
#   effective N per Table 2 cell (E2, E4, E5, E6).
# Inputs: 19 summary, 20 and 20_v2 outputs. Outputs: analisis/output/32_revision_bloque_SBM_tabla2[_vN]/.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
set.seed(20261022)
d20 <- "analisis/output/20_tabla2_core_periferia"
dwin <- "analisis/output/20_tabla2_core_periferia_v2"
d19 <- "analisis/output/19_SBM_redes_limpias"
leer <- function(path) read_csv(path, show_col_types = FALSE, trim_ws = FALSE)
ficheros <- c("resultados.csv", "sensibilidades.csv", "comparacion_sensibilidades.csv",
  "brechas_persona.csv", "descriptivos_particion.csv", "cobertura.csv", "tabla2_borrador.csv")
comparacion <- bind_rows(lapply(ficheros, function(f) {
  a <- leer(file.path(d20, f)); b <- leer(file.path(dwin, f))
  stopifnot(identical(names(a), names(b)), nrow(a) == nrow(b))
  num <- names(a)[vapply(a, is.numeric, logical(1))]
  no_num <- setdiff(names(a), num)
  numeric_ok <- all(vapply(num, function(nm)
    isTRUE(all.equal(a[[nm]], b[[nm]], tolerance = 1e-10)), logical(1)))
  text_ok <- all(vapply(no_num, function(nm) identical(a[[nm]], b[[nm]]), logical(1)))
  delta <- unlist(lapply(num, function(nm) abs(a[[nm]] - b[[nm]])))
  delta <- delta[is.finite(delta)]
  tibble(file = f, rows = nrow(a), numeric_equal_tol_1e10 = numeric_ok,
    nonnumeric_identical = text_ok, max_absolute_difference = if (length(delta)) max(delta) else 0)
}))
stopifnot(all(comparacion$numeric_equal_tol_1e10), all(comparacion$nonnumeric_identical))
r <- leer(file.path(d20, "resultados.csv"))
b <- leer(file.path(d20, "brechas_persona.csv"))
bh <- r |> group_by(family) |> summarise(m = n(),
  max_BH_difference = max(abs(p_BH - p.adjust(p, "BH"))), .groups = "drop")
stopifnot(all(bh$max_BH_difference < 1e-12))
ns <- r |> group_by(family) |> summarise(tests = n(), n_min = min(n), n_max = max(n),
  significant_BH05 = sum(p_BH < .05), positive_estimates = sum(estimate > 0), .groups = "drop")
sesiones <- b |> filter(reference == "conjunta_3", cohort == "PD", is.finite(gap)) |>
  group_by(task, variable, network, participant) |>
  summarise(n_available_conditions = n_distinct(condition), .groups = "drop") |>
  group_by(task, variable, network) |>
  summarise(n_PD_both_conditions = sum(n_available_conditions == 2),
    n_PD_single_condition = sum(n_available_conditions == 1), .groups = "drop")
# Corregir solo rotulación: mismas estimaciones, IC y p; N explícito por contraste.
tabla <- r |> filter(effect != "session2_minus_session1") |>
  mutate(column = recode(effect, core_minus_periphery = "Core_minus_periphery",
    ON_minus_OFF = "ON_minus_OFF_gap", HC_minus_OFF = "HC_minus_OFF_gap", HC_minus_ON = "HC_minus_ON_gap"),
    N_description = if_else(effect == "ON_minus_OFF", paste0("paired N=", n),
      paste0("N=", n, " (HC=", n_HC, "; PD=", n_PD, ")")),
    cell = sprintf("%.2f [%.2f, %.2f]; p_BH=%s; %s", estimate, conf_low, conf_high,
      formatC(p_BH, format = "g", digits = 2), N_description)) |>
  select(task, variable, network, column, cell) |>
  pivot_wider(names_from = column, values_from = cell) |>
  mutate(exploratory = variable == "AoA")
rr <- leer(file.path(d19, "resumen_redes.csv"))
ganadores <- rr |> filter(family != "igualada", modelo_ganador != "hubspoke") |>
  select(graph, family, modelo_ganador, mdl_hubspoke_min, mdl_capas_min, mdl_sin_core)
out <- nuevo_directorio("32_revision_bloque_SBM_tabla2")
guardar_csv(comparacion, file.path(out, "replica_Windows_vs_Linux.csv"))
guardar_csv(ns, file.path(out, "N_y_resultados_por_familia.csv"))
guardar_csv(bh, file.path(out, "verificacion_BH.csv"))
guardar_csv(sesiones, file.path(out, "condiciones_disponibles_gradiente_PD.csv"))
guardar_csv(tabla, file.path(out, "tabla2_rotulada_con_N.csv"))
guardar_csv(ganadores, file.path(out, "modelos_ganadores_distintos_HS.csv"))
registrar_entorno(out, c(file.path(d20, ficheros), file.path(dwin, ficheros),
  file.path(d19, "resumen_redes.csv"), "analisis/R/32_revision_bloque_SBM_tabla2.R"))
writeLines(c("# Revisión numérica/documental del bloque SBM y Tabla2",
  "Siete CSV reproducidos en Windows20_v2 frente a Linux20, tolerancia1e-10; BH verificado.",
  "N efectivos principales: G39–45; C/S14–21 pares; H35–45 personas en total.",
  "Tabla rotulada con N usa los mismos resultados20; no cambia el método ni acepta D-023.",
  "La partición se estima en los mismos participantes: inferencia condicional a ella,",
  "sin incorporar incertidumbre de reestimar redes/particiones en otras muestras.",
  "No detectar cambios en brechas no demuestra equivalencia de grupos ni de redes.",
  "Igualadas19:10 extracciones y3 cadenas; no prueba de igualdad de proporciones.",
  "Pendientes de la revisión metodológica: ver notas/revision_D023_SBM_tabla2_20261022.md."
), file.path(out, "LEEME.md"), useBytes = TRUE)
print(comparacion)
print(ns)
cat("Salida:", out, "\n")
