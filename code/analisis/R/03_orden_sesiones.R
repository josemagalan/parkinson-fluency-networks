# Propósito: incorporar el orden ON/OFF aportado por Elena y derivar la sesión.
# Entradas: CSV de Elena (raíz), copia extraída del MSG, demografía histórica,
# y tokens.rds de 01_datos_tidy (congelado, no se sobrescribe).
# Salidas: analisis/output/03_orden_sesiones[/_vN]/, tablas CSV/RDS y controles.
# E4 (también prepara E2): integración de metadatos, sin ajustar modelos nuevos.
# Paquetes: readr, dplyr, tidyr, readxl, purrr, broom (00_utilidades.R).
# Ejecutar desde revision_R3/: Rscript --vanilla analisis/R/03_orden_sesiones.R
# --- English ---
# Purpose: add the individual ON/OFF order (11 OFF->ON, 10 ON->OFF) and derive session 1/2 (E4; prepares E2).
# Inputs: order CSV, historical demographics, 01 tokens.rds. Outputs: analisis/output/03_orden_sesiones[_vN]/
#   (tokens.csv/.rds = tidy dataset used by later scripts). No models are fitted.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
set.seed(20261022)
fuente <- "../datos sexo_edad_parkinson_control.csv"
adjunto <- "historial/elena_20261022_v3/datos sexo_edad_parkinson_control.csv"
mensaje <- "../RE_ Major revision requested_ TLCD-2025-0015_elena.msg"
historico <- "../parkinson/data/datos sexo_edad_parkinson_control.csv"
base_tokens <- "analisis/output/01_datos_tidy/tokens.rds"
stopifnot(unname(tools::md5sum(fuente)) == unname(tools::md5sum(adjunto)))
clinica <- read_delim(fuente, delim = ";", locale = locale(encoding = "UTF-8"), show_col_types = FALSE) |>
  transmute(participant_raw = as.integer(`ID parkinson`), sexo, edad,
    escolaridad = `años de escolarización`, MMSE, first_condition = normalizar(`primera evaluacion`),
    source_row = row_number() + 1L)
antigua <- read_delim(historico, delim = ";", locale = locale(encoding = "UTF-8"), show_col_types = FALSE) |>
  rename(participant_raw = participante)
stopifnot(nrow(clinica) == 21L, !anyDuplicated(clinica$participant_raw),
  setequal(clinica$participant_raw, 1:21), setequal(antigua$participant_raw, 1:21),
  all(clinica$first_condition %in% c("on", "off")))
# Demografía: los 21 ID y los cuatro campos coinciden con la tabla histórica.
validacion_ids <- full_join(clinica, antigua, by = "participant_raw", suffix = c("_Elena", "_R2")) |>
  transmute(participant_raw, sexo_igual = sexo_Elena == sexo_R2, edad_igual = edad_Elena == edad_R2,
    escolaridad_igual = escolaridad_Elena == escolaridad_R2, MMSE_igual = MMSE_Elena == MMSE_R2)
stopifnot(nrow(validacion_ids) == 21L, all(as.matrix(validacion_ids[-1])))
orden_pacientes <- clinica |> transmute(participant_raw = as.character(participant_raw),
  participant = paste0("PD_", sprintf("%02d", as.integer(participant_raw))),
  first_condition = toupper(first_condition),
  order = if_else(first_condition == "OFF", "OFF_ON", "ON_OFF"),
  source_file = fuente, source_field = "primera evaluacion", source_row)
stopifnot(sum(orden_pacientes$order == "OFF_ON") == 11L,
  sum(orden_pacientes$order == "ON_OFF") == 10L)
previos <- readRDS(base_tokens)
stopifnot(nrow(previos) == 2965L, !anyDuplicated(previos$token_id),
  all(is.na(previos$session)), all(is.na(previos$order)),
  setequal(unique(previos$participant[previos$cohort == "PD"]), orden_pacientes$participant))
tokens <- previos |> select(-session, -order) |>
  left_join(orden_pacientes |> select(participant, first_condition, order),
    by = "participant", relationship = "many-to-one") |>
  mutate(session = case_when(cohort == "HC" ~ 1L,
    condition == first_condition ~ 1L, cohort == "PD" ~ 2L),
    order_source = if_else(cohort == "PD", fuente, NA_character_),
    session_source = if_else(cohort == "PD", "condition + primera evaluacion (Elena)", "diseño: evaluación única HC")) |>
  relocate(session, order, first_condition, .after = condition)
# Invariantes: el enriquecimiento solo cambia metadatos de sesión/orden.
stopifnot(nrow(tokens) == nrow(previos), identical(tokens$token_id, previos$token_id),
  identical(tokens[names(previos)[!names(previos) %in% c("session", "order")]],
    previos[names(previos)[!names(previos) %in% c("session", "order")]]),
  all(!is.na(tokens$session)), all(!is.na(tokens$order[tokens$cohort == "PD"])),
  all(is.na(tokens$order[tokens$cohort == "HC"])))
sesiones <- tokens |> distinct(participant, cohort, group, condition, session, order, first_condition)
listas <- tokens |> distinct(participant, task, condition, session, order)
pares <- sesiones |> filter(cohort == "PD") |> group_by(participant) |>
  summarise(n = n(), n_sesion = n_distinct(session), n_condicion = n_distinct(condition), .groups = "drop")
stopifnot(nrow(sesiones) == 66L, nrow(listas) == 198L, nrow(pares) == 21L,
  all(pares$n == 2L & pares$n_sesion == 2L & pares$n_condicion == 2L),
  all(sesiones$session[sesiones$condition == sesiones$first_condition & sesiones$cohort == "PD"] == 1L),
  all((sesiones |> filter(cohort == "HC"))$session == 1L))
resumen_orden <- orden_pacientes |> count(order, name = "n_pacientes")
cruce_sesion <- sesiones |> count(cohort, order, condition, session, name = "n_evaluaciones")
out <- nuevo_directorio("03_orden_sesiones")
for (nm in c("tokens", "orden_pacientes", "sesiones", "listas", "validacion_ids", "resumen_orden", "cruce_sesion")) {
  guardar_csv(get(nm), file.path(out, paste0(nm, ".csv")))
}
saveRDS(tokens, file.path(out, "tokens.rds"))
# Relectura exacta de las palabras y del orden de filas del CSV exportado.
releido <- read_csv(file.path(out, "tokens.csv"), trim_ws = FALSE, show_col_types = FALSE)
stopifnot(identical(releido$word_raw, tokens$word_raw), identical(releido$word_norm, tokens$word_norm),
  identical(releido$token_id, tokens$token_id), identical(as.integer(releido$session), tokens$session))
writeLines(c(
  "# Orden de sesiones incorporado — 22-10-2026",
  "Fuente: CSV adjunto al correo de Elena, idéntico por MD5 al CSV recibido en raíz.",
  "21 pacientes: PD_01–PD_11 OFF→ON (11); PD_12–PD_21 ON→OFF (10).",
  "La R2 describe el reparto inverso (11 ON primero/10 OFF primero): corregirlo al preparar R3, sin editar ahora el manuscrito.",
  "session=1 si la condición corresponde a primera evaluacion; session=2 para la otra condición.",
  "order=OFF_ON u ON_OFF. first_condition=OFF u ON. Los tres tipos de tarea heredan la misma sesión de cada evaluación.",
  "HC: session=1 por diseño de evaluación única; order y first_condition=NA porque no hay cruce ON/OFF.",
  "Los 21 ID y sexo/edad/escolaridad/MMSE coinciden con la tabla demográfica histórica.",
  "Se conservan los 2.965 tokens, sus ID, palabras y variables psicolingüísticas: solo se añaden/completan metadatos.",
  "La celda partir<TAB>llegar continúa pendiente; no se cambia la limpieza ni la asignación AoA provisional.",
  "No hay fechas, intervalo ni confirmación de si las dos evaluaciones fueron el mismo día; no se infieren del orden.",
  "CSV UTF-8: importar con trim_ws=FALSE para preservar word_raw; tokens.rds conserva los tipos de R.",
  "El fichero de tokens de 01_datos_tidy y los resultados históricos permanecen intactos.",
  "No se han ajustado modelos ni modificado el manuscrito."
), file.path(out, "LEEME.md"), useBytes = TRUE)
registrar_entorno(out, c(fuente, adjunto, mensaje, historico, base_tokens,
  "analisis/R/03_orden_sesiones.R", "analisis/R/00_utilidades.R"))
print(resumen_orden); print(cruce_sesion); cat("Salida:", out, "\n")
