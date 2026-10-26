# Propósito: verificar si Sequences.R conserva el orden de producción de Figura 4.
# Entradas: Sequences.R (solo se evalúa su función pura), Excel original y diccionarios.
# Salidas: output/09_auditoria_secuencia[/_vN]/, correspondencias de posición y resumen.
# E2/E5/E6: auditoría documental/computacional; no modifica ni ejecuta todo el script histórico.
# Paquetes: 00_utilidades.R. set.seed explícito; ninguna simulación.
# --- English ---
# Purpose: check whether R2 Sequences.R kept the production order used in Figure 4 (it did not)
#   by evaluating only its word-mapping function (E2, E5, E6).
# Inputs: ../parkinson/Sequences.R, raw Excel, AoA dictionaries. Outputs: analisis/output/09_auditoria_secuencia[_vN]/.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
set.seed(20261022)
ast <- parse("../parkinson/Sequences.R", encoding = "UTF-8")
found <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
  identical(e[[2]], as.name("transformar_palabras")), as.list(ast))
stopifnot(length(found) == 1L)
env <- new.env(parent = baseenv()); eval(found[[1]], envir = env)
cfg <- tibble(sheet = c("acciones park", "animales_park", "palab F_park", "acciones cont", "animales_cont", "palab F_Cont"),
  task = rep(c("actions", "animals", "F"), 2))
dic <- leer_diccionarios()
mapping <- list(); cursor <- 0L
for (i in seq_len(nrow(cfg))) {
  d <- suppressMessages(read_excel("../parkinson/data/FLUIDEZ PARKINSON.xlsx", sheet = cfg$sheet[i]))
  for (j in seq_len(nrow(d))) {
    w <- unlist(d[j, 3:ncol(d)], use.names = FALSE); w <- w[!is.na(w)]
    dictionary <- if (cfg$task[i] == "actions") dic$aoa_verbs else dic$aoa_nouns
    literal <- env$transformar_palabras(w, dictionary)
    m <- merge(data.frame(palabra = w, original_rank = seq_along(w)), dictionary,
      by = "palabra", all.x = TRUE)
    stopifnot(isTRUE(all.equal(unname(literal), unname(m$AoA))))
    cursor <- cursor + 1L
    mapping[[cursor]] <- as_tibble(m) |> mutate(task = cfg$task[i], group = d$grupo[j],
      participant = paste0(if_else(d$grupo[j] == "control", "HC_", "PD_"), sprintf("%02d", d$participante[j])),
      assigned_position_by_Sequences = row_number(),
      word_in_original_position = w[assigned_position_by_Sequences],
      wrong_word_at_position = palabra != word_in_original_position)
  }
}
mapping <- bind_rows(mapping)
por_lista <- mapping |> group_by(participant, group, task) |>
  summarise(n_positions = n(), n_wrong_words_at_position = sum(wrong_word_at_position),
    any_reordered = any(wrong_word_at_position), .groups = "drop")
resumen <- por_lista |> group_by(task) |>
  summarise(n_lists = n(), n_reordered = sum(any_reordered), n_positions = sum(n_positions),
    n_wrong_positions = sum(n_wrong_words_at_position), wrong_pct = 100 * n_wrong_positions / n_positions,
    .groups = "drop")
stopifnot(nrow(por_lista) == 198L)
out <- nuevo_directorio("09_auditoria_secuencia")
for (nm in c("mapping", "por_lista", "resumen")) guardar_csv(get(nm), file.path(out, paste0(nm, ".csv")))
registrar_entorno(out, c("../parkinson/Sequences.R", "../parkinson/data/FLUIDEZ PARKINSON.xlsx",
  "analisis/R/09_auditoria_secuencia.R", list.files("../parkinson/data", pattern = "\\.xlsx$", full.names = TRUE)))
writeLines(c("# Auditoría de posición en Sequences.R — 22-10-2026",
  "La función histórica transformar_palabras usa merge(...,all.x=TRUE) con sort=TRUE por defecto.",
  "Se evalúa SOLO esa función literal, nunca el script con escrituras; se comprueba su salida con el mapa de palabras/posición original.",
  "Después, obtener_datos_tidy asigna row_number sobre el vector ya ordenado por palabra: no conserva la posición de producción.",
  "El hallazgo se cuantifica con AoA (diccionario sin claves duplicadas por tarea); afecta al procedimiento de todas las variables de la figura.",
  "La Figura 4 no debe interpretarse como evolución temporal sin regenerarla a partir de position del tidy.",
  "Esto no implica que mainParkinson_Dirigidas.R ordene alfabéticamente: ese script construye adyacencias directamente desde words.",
  "No se han modificado fuente, figura ni manuscrito."), file.path(out, "LEEME.md"), useBytes = TRUE)
print(resumen); cat("Salida:", out, "\n")
