# Propósito: auditar exportaciones ARFF y preparar folds por persona sin elegir nueva vía ML.
# Entradas: ARFF históricos, fluidez original, tokens y medias nuevas.
# Salidas: output/08_auditoria_ML[/_vN]/, inventario, datos preparados y 45 folds sin fuga.
# E3/E2: no replica resultados Weka ni entrena clasificadores sin decisión de coautores.
# Paquetes: 00_utilidades.R + foreign (R recomendado).
# --- English ---
# Purpose: audit the R2 Weka ARFF exports and prepare 45 leave-one-person-out folds without leakage;
#   no classifiers are trained (E3, E2).
# Inputs: ../parkinson/results/*.arff, raw Excel, 03 and 06 outputs. Outputs: analisis/output/08_auditoria_ML[_vN]/.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
set.seed(20261022)
files <- sort(list.files("../parkinson/results", pattern = "\\.arff$", full.names = TRUE))
inventario <- faltantes <- list()
for (f in files) {
  d <- foreign::read.arff(f)
  gcol <- which(names(d) %in% c("grupo", "y"))
  if (!length(gcol)) gcol <- which(vapply(d, function(x) {
    y <- unique(as.character(x[!is.na(x)])); length(y) >= 2L && all(y %in% c("control", "off", "on"))
  }, logical(1)))
  stopifnot(length(gcol) == 1L)
  predictors <- setdiff(seq_len(ncol(d)), gcol)
  numerical <- predictors[vapply(d[predictors], is.numeric, logical(1))]
  mat <- as.matrix(d[numerical])
  grupos <- table(d[[gcol]])
  inventario[[basename(f)]] <- tibble(file = basename(f), n_rows = nrow(d), n_predictors = length(predictors),
    class_name = names(d)[gcol], class_position = gcol, n_columns = ncol(d), class_is_last = gcol == ncol(d),
    has_participant = any(grepl("participante|participant|subject", names(d), ignore.case = TRUE)),
    n_HC = sum(d[[gcol]] == "control"), n_OFF = sum(d[[gcol]] == "off"), n_ON = sum(d[[gcol]] == "on"),
    NA_numeric = sum(is.na(mat)), infinite_numeric = sum(is.infinite(mat)),
    rows_complete_numeric = sum(rowSums(is.na(mat) | !is.finite(mat)) == 0),
    rows_repeated_predictors = sum(duplicated(d[predictors])))
  faltantes[[basename(f)]] <- bind_rows(lapply(predictors, function(j) tibble(file = basename(f),
    feature = names(d)[j], n_missing = sum(is.na(d[[j]])),
    n_infinite = if (is.numeric(d[[j]])) sum(is.infinite(d[[j]])) else 0L)))
}
inventario <- bind_rows(inventario); faltantes <- bind_rows(faltantes)
# Validar el mapa de las 66 filas de cada ARFF de características contra ID/condición/conteo original.
configs <- tibble(task = c("actions", "animals", "F"), park = c("acciones park", "animales_park", "palab F_park"),
  hc = c("acciones cont", "animales_cont", "palab F_Cont"), prefix = c("ac", "an", "f"),
  file = c("accionesWeka.arff", "animalesWeka.arff", "FpalabraWeka.arff"))
mapas <- list()
for (i in seq_len(nrow(configs))) {
  cfg <- configs[i, ]
  raw <- bind_rows(lapply(c(cfg$park, cfg$hc), function(sheet) {
    d <- suppressMessages(read_excel("../parkinson/data/FLUIDEZ PARKINSON.xlsx", sheet = sheet))
    tibble(participant = paste0(if_else(d$grupo == "control", "HC_", "PD_"), sprintf("%02d", as.integer(d$participante))),
      group = d$grupo, n_original = rowSums(!is.na(as.data.frame(d[3:ncol(d)]))))
  }))
  d <- foreign::read.arff(file.path("../parkinson/results", cfg$file))
  stopifnot(nrow(d) == nrow(raw), all(as.character(d$grupo) == raw$group),
    all(d[[paste0(cfg$prefix, "_longitud")]] == raw$n_original))
  mapas[[cfg$task]] <- raw |> mutate(file = cfg$file, row = row_number(), task = cfg$task)
}
mapas <- bind_rows(mapas)
personas <- mapas |> distinct(participant) |> arrange(participant) |> mutate(fold = row_number())
folds <- mapas |> left_join(personas, by = "participant", relationship = "many-to-one")
validacion <- bind_rows(lapply(personas$fold, function(f) {
  test <- folds |> filter(fold == f); train <- folds |> filter(fold != f)
  tibble(fold = f, test_person = unique(test$participant), train_test_overlap = length(intersect(train$participant, test$participant)),
    n_test_rows_per_task = nrow(test) / 3L)
}))
stopifnot(nrow(personas) == 45L, all(validacion$train_test_overlap == 0),
  all(validacion$n_test_rows_per_task %in% c(1, 2)))
# Demostración determinista de la fuga de un LOO por filas en los 42 registros de pacientes.
fuga <- mapas |> filter(task == "actions") |> group_by(participant) |>
  mutate(same_person_remains_in_train_if_row_LOO = n() > 1L) |> ungroup()
# Matriz nueva, solo preparación: no edad/sexo/MMSE, orden o sesión como predictores;
# AoA fuera hasta decisión, medias de cobertura explícita para imputar DENTRO de cada fold si se usa.
medias <- read_csv("analisis/output/06_psicolinguistica_sesion/medias.csv", show_col_types = FALSE)
caracteristicas <- medias |> filter(variable != "AoA") |>
  select(participant, condition, task, variable, mean, coverage) |>
  pivot_wider(names_from = c(task, variable), values_from = c(mean, coverage))
produccion <- readRDS("analisis/output/03_orden_sesiones/tokens.rds") |>
  count(participant, condition, task, name = "n_tokens") |> pivot_wider(names_from = task, values_from = n_tokens, names_prefix = "count_")
caracteristicas <- left_join(produccion, caracteristicas, by = c("participant", "condition"), relationship = "one-to-one") |>
  left_join(personas, by = "participant", relationship = "many-to-one")
stopifnot(nrow(caracteristicas) == 66L, !anyNA(caracteristicas$fold))
out <- nuevo_directorio("08_auditoria_ML")
for (nm in c("inventario", "faltantes", "mapas", "folds", "validacion", "fuga", "caracteristicas")) guardar_csv(get(nm), file.path(out, paste0(nm, ".csv")))
registrar_entorno(out, c(files, "../parkinson/processing.R", "../parkinson/naiveBayes.R",
  "../parkinson/data/FLUIDEZ PARKINSON.xlsx", "analisis/output/06_psicolinguistica_sesion/medias.csv", "analisis/R/08_auditoria_ML.R"))
writeLines(c("# Auditoría ML — 22-10-2026",
  "Se leen los ocho ARFF guardados; su existencia no demuestra qué archivo/configuración se usó en Weka.",
  "La columna participante se eliminó en processing.R. Los tres ARFF por tarea se vinculan al orden original por grupo+longitud exactos.",
  "Los 42 registros PD de un LOO por FILAS dejan la otra condición del paciente en entrenamiento. El comentario del código naiveBayes.R no equivale a LOO por persona.",
  "Class_is_last=FALSE exige seleccionar la clase explícitamente en Weka; no se infiere que se seleccionara mal.",
  "Los ARFF de características contienen demografía; faltantes y no finitos se cuantifican por columna.",
  "Folds preparados: dejar fuera una persona completa (45 folds), nunca una sola fila; validado solapamiento cero.",
  "Si se rehace ML: fijar objetivo/modelo antes de evaluar, imputar/escalar/seleccionar solo en entrenamiento, tuning anidado por persona, intervalos por persona y prueba nula compatible con las etiquetas repetidas.",
  "No se han entrenado modelos ni decidido eliminar ML. P4/P7 siguen pendientes; no añadir nuevos clasificadores para elegir el mejor tras ver accuracy.",
  "La matriz nueva excluye demografía, AoA, orden/session como predictores; sus NA requieren imputación dentro del fold, no global.",
  "Sin cambios en ARFF, datos ni manuscrito."), file.path(out, "LEEME.md"), useBytes = TRUE)
print(inventario, width = 150); print(validacion |> count(n_test_rows_per_task)); cat("Salida:", out, "\n")
