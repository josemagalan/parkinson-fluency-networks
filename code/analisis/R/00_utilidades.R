# Propósito: utilidades comunes de lectura, diccionarios y trazabilidad R3.
# Entradas: ../parkinson/data (solo lectura). Salidas: funciones; no escribe al cargar.
# Puntos editoriales: E2, E5. Paquetes: readxl, readr, dplyr, tidyr, purrr, broom.
# Ejecutar los scripts llamantes desde revision_R3/. Semilla común, sin simulación.
# --- English ---
# Purpose: shared helpers (text normalisation, versioned output folders, CSV export, dictionary reading and
#   annotation, environment/MD5 logging). Reads ../parkinson/data (read-only); writes nothing when sourced. Seed 20261022.
set.seed(20261022)
user_library <- "C:/Users/josem/AppData/Local/R/win-library/4.5"
if (dir.exists(user_library)) .libPaths(c(user_library, .libPaths()))
suppressPackageStartupMessages({
  library(readxl); library(readr); library(dplyr); library(tidyr)
  library(purrr); library(broom)
})
stopifnot(dir.exists("../parkinson/data"))
# Incluye los tres espacios no separables (U+00A0) del Excel; trimws por
# defecto solo elimina espacio ASCII/tabulador/CR/LF.
normalizar <- function(x) tolower(trimws(x, whitespace = "[\\h\\v]"))
nuevo_directorio <- function(nombre) {
  base <- file.path("analisis/output", nombre)
  dest <- base; v <- 1L
  while (dir.exists(dest)) { v <- v + 1L; dest <- paste0(base, "_v", v) }
  dir.create(dest, recursive = TRUE)
  dest
}
# Exportación explícitamente entrecomillada. Al importar word_raw se necesita
# trim_ws=FALSE: readr recorta tabuladores incluso dentro de campos citados.
guardar_csv <- function(x, path) readr::write_csv(x, path, na = "NA", quote = "all")
leer_diccionarios <- function() {
  leer <- function(file, sheet) read_excel(file.path("../parkinson/data", file), sheet = sheet)
  lex <- leer("Lexesp data.xlsx", "Hoja1")
  list(
    logFreq = leer("Buscador frecuencia oral Subtlex-Esp.xlsx", "Subtlex-Esp") |>
      select(palabra = Word, logFreq = `Log freq.`),
    motor = leer("buscador motor content.xlsx", "motor_content") |>
      rename(palabra = Verbs, motor = `Average motor content`) |>
      mutate(palabra = tolower(palabra)),
    imageability = lex |> transmute(palabra = .data[[names(lex)[1]]], imageability = as.numeric(Imageability)),
    concreteness = lex |> transmute(palabra = .data[[names(lex)[1]]], concreteness = as.numeric(Concreteness)),
    familiarity = lex |> transmute(palabra = .data[[names(lex)[1]]], familiarity = as.numeric(Familiarity)),
    aoa_verbs = leer("AoA_verbos.xlsx", "Hoja1") |> select(palabra = Word, AoA = AoA_average),
    aoa_nouns = leer("AoA_sustantivos.xlsx", "Hoja1") |> select(palabra = word, AoA = averageAoA)
  )
}
diccionarios_unicos <- function(dic, normalizar_claves = TRUE) {
  lapply(dic, function(d) {
    if (normalizar_claves) d$palabra <- normalizar(d$palabra)
    d |> group_by(palabra) |> summarise(across(where(is.numeric),
      ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)), .groups = "drop")
  })
}
anotar_limpio <- function(datos, dic, normalizar_claves = TRUE) {
  dic <- diccionarios_unicos(dic, normalizar_claves)
  n <- nrow(datos)
  for (nm in c("logFreq", "motor", "imageability", "concreteness", "familiarity")) {
    datos <- left_join(datos, dic[[nm]], by = c("word_norm" = "palabra"), relationship = "many-to-one")
  }
  vi <- match(datos$word_norm, dic$aoa_verbs$palabra)
  ni <- match(datos$word_norm, dic$aoa_nouns$palabra)
  av <- dic$aoa_verbs$AoA[vi]; an <- dic$aoa_nouns$AoA[ni]
  datos$AoA <- ifelse(datos$test == "Acciones", av,
    ifelse(datos$test == "Fpalabra" & is.na(an), av, an))
  datos$AoA_source <- ifelse(datos$test == "Acciones", "verbs",
    ifelse(datos$test == "Fpalabra" & is.na(an), "verbs_fallback", "nouns"))
  datos$AoA_source[is.na(datos$AoA)] <- NA_character_
  stopifnot(nrow(datos) == n)
  datos
}
variables_psico <- c("logFreq", "motor", "imageability", "concreteness", "familiarity", "AoA")
registrar_entorno <- function(out, entradas) {
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  guardar_csv(tibble(path = entradas, md5 = unname(tools::md5sum(entradas))), file.path(out, "entradas_md5.csv"))
}
