# Propósito: datos largos reproducibles de palabras-token, limpieza mínima y cobertura.
# Entradas: ../parkinson/data/FLUIDEZ PARKINSON.xlsx y diccionarios originales.
# Salidas: analisis/output/01_datos_tidy[/_vN] (tokens.csv/RDS, limpieza, cobertura, conteos, incidencias, dic_duplicados).
# E1, E2, E4, E5. Paquetes: los de 00_utilidades.R. No calcula modelos nuevos.
# Ejecutar desde revision_R3/: Rscript --vanilla analisis/R/01_datos_tidy.R
# --- English ---
# Purpose: long token table (one row per produced word) with minimal cleaning, psycholinguistic annotation
#   and dictionary coverage (E1, E2, E4, E5).
# Inputs: ../parkinson/data/FLUIDEZ PARKINSON.xlsx and dictionaries. Outputs: analisis/output/01_datos_tidy[_vN]/.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
set.seed(20261022)
config <- tibble(
  sheet = c("acciones park", "animales_park", "palab F_park", "acciones cont", "animales_cont", "palab F_Cont"),
  test = rep(c("Acciones", "Animales", "Fpalabra"), 2),
  task = rep(c("actions", "animals", "F"), 2), cohort = rep(c("PD", "HC"), each = 3))
entrada <- "../parkinson/data/FLUIDEZ PARKINSON.xlsx"
tokens <- bind_rows(lapply(seq_len(nrow(config)), function(i) {
  cfg <- config[i, ]
  d <- suppressMessages(read_excel(entrada, sheet = cfg$sheet, col_types = "text", trim_ws = FALSE))
  stopifnot(names(d)[1] == "participante", names(d)[2] == "grupo")
  d |> mutate(source_row = row_number() + 1L) |>
    pivot_longer(cols = all_of(names(d)[3:ncol(d)]), names_to = "source_column", values_to = "word_raw") |>
    filter(!is.na(word_raw), nzchar(word_raw)) |>
    transmute(participant_raw = participante,
      participant = paste0(cfg$cohort, "_", sprintf("%02d", as.integer(participante))),
      cohort = cfg$cohort, group = normalizar(grupo),
      condition = ifelse(group == "control", "HC", toupper(group)),
      session = NA_integer_, order = NA_character_, task = cfg$task, test = cfg$test,
      position = match(source_column, names(d)) - 2L,
      word_raw, word_norm = normalizar(word_raw), source_sheet = cfg$sheet, source_row, source_column)
}))
stopifnot(nrow(tokens) == 2965L, all(nzchar(tokens$word_norm)),
  n_distinct(tokens$participant) == 45L,
  !anyDuplicated(tokens[c("participant", "group", "task", "position")]))
sesiones <- tokens |> distinct(participant, cohort, group, task)
stopifnot(nrow(sesiones) == 198L,
  all((sesiones |> count(task, group))$n %in% c(21L, 24L)),
  all((sesiones |> filter(cohort == "PD") |> count(participant, task))$n == 2L))
dic <- leer_diccionarios()
tokens <- anotar_limpio(tokens, dic)
tokens <- tokens |> mutate(token_id = paste(participant, group, task, position, sep = ":"), .before = 1)
limpieza <- tokens |> filter(word_raw != word_norm) |>
  select(token_id, participant, group, task, position, word_raw, word_norm)
cobertura <- tokens |> pivot_longer(all_of(variables_psico), names_to = "variable", values_to = "valor") |>
  group_by(task, group, variable) |> summarise(n_tokens = n(), n_con_dato = sum(!is.na(valor)),
    n_NA = sum(is.na(valor)), porcentaje_con_dato = 100 * mean(!is.na(valor)), .groups = "drop")
conteos <- tokens |> count(participant, cohort, group, task, name = "n_tokens")
incidencias <- tokens |> filter(grepl("[\\h\\v]", word_norm, perl = TRUE)) |>
  select(token_id, word_raw, word_norm, source_sheet, source_row, source_column)
dic_duplicados <- bind_rows(lapply(names(dic), function(nm) {
  dic[[nm]] |> mutate(palabra = normalizar(palabra)) |> count(palabra) |> filter(n > 1) |> mutate(diccionario = nm)
}))
out <- nuevo_directorio("01_datos_tidy")
for (nm in c("tokens", "limpieza", "cobertura", "conteos", "incidencias", "dic_duplicados")) {
  guardar_csv(get(nm), file.path(out, paste0(nm, ".csv")))
}
saveRDS(tokens, file.path(out, "tokens.rds"))
writeLines(c(
  "# Datos tidy — 22-10-2026",
  "Una fila por celda no vacía de producción: 2.965 tokens; 45 personas; 66 evaluaciones; 198 listas.",
  "participant: PD/HC + identificador original. El prefijo evita colisiones paciente/control.",
  "group: control/off/on; cohort: HC/PD; condition: HC/OFF/ON.",
  "position conserva la columna de producción original (columna Excel menos dos).",
  "word_raw conserva mayúsculas/tabuladores; word_norm aplica trimws Unicode (incluye U+00A0) + tolower (D-006/D-008).",
  "No se eliminan repeticiones ni se lematiza. No se separan palabras con tabulador interno.",
  "session y order permanecen NA pendientes de P1; no se infieren del orden de filas.",
  "AoA: verbos en acciones; sustantivos en animales; sustantivos con respaldo de verbos en F (D-005 provisional).",
  "AoA_source indica el diccionario. Las claves duplicadas de Lexesp se promedian por variable.",
  "Motor se anota donde existe, aunque solo sea pertinente su interpretación en acciones.",
  "NA = dato no disponible; no se imputa. tokens.csv es UTF-8; tokens.rds conserva tipos y tabuladores.",
  "Importar el CSV con readr::read_csv(..., trim_ws = FALSE) para preservar word_raw; el valor por defecto recorta los tabuladores.",
  "Las incidencias con espacios internos requieren consulta clínica (no son errores corregidos automáticamente)."
), file.path(out, "LEEME.md"), useBytes = TRUE)
registrar_entorno(out, c(entrada, list.files("../parkinson/data", pattern = "\\.xlsx$", full.names = TRUE)))
print(tokens |> count(task, group)); cat("Tokens:", nrow(tokens), "Limpiezas:", nrow(limpieza), "\nSalida:", out, "\n")
