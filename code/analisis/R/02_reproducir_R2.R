# Propósito: reproducir Tabla 2/3 R2 y cuantificar duplicación/codificación/limpieza.
# Entradas: corePeriphery.R (lógica reimplementada; solo MD5), particiones CSV y diccionarios del repo histórico;
# resultsCorePeriphery/results_table.xlsx, tabla_resumen_GLM.xlsx y descriptivos.
# Salidas: analisis/output/02_reproducir_R2[/_vN], CSV, RDS, informe y sessionInfo.
# E2/E5: auditoría del análisis histórico; NO inferencia nueva ni reajuste del SBM.
# Paquetes: los de 00_utilidades.R. Ejecutar: Rscript --vanilla analisis/R/02_reproducir_R2.R
# --- English ---
# Purpose: reproduce R2 Tables 2/3 (GLMs on the R2 SBM partitions) and quantify merge-duplication,
#   encoding and cleaning errors (E2, E5 audit; no new inference, no SBM refit).
# Inputs: R2 partition CSVs, dictionaries and R2 Excel tables. Outputs: analisis/output/02_reproducir_R2[_vN]/.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
set.seed(20261022)
dic <- leer_diccionarios()
leer_nodos <- function(encoding, reparar = FALSE) {
  bind_rows(lapply(c("directed", "undirected"), function(net) {
    files <- sort(list.files(file.path("../parkinson/hub and spoke model", net, "csvs"),
      pattern = "\\.csv$", full.names = TRUE))
    bind_rows(lapply(files, function(f) {
      enc <- encoding
      if (encoding == "auto") enc <- if (is.na(iconv(rawToChar(read_file_raw(f)), from = "UTF-8", to = "UTF-8"))) "latin1" else "UTF-8"
      if (reparar) {
        lineas <- read_lines(f, locale = locale(encoding = enc))
        # Una fila completa está entrecomillada: recuperar sus cuatro campos,
        # conservando literalmente el tabulador interno de partir\tllegar.
        lineas <- sub('^"([0-9]+,n[0-9]+,.*,[01])"$', "\\1", lineas)
        d <- suppressMessages(read_csv(I(paste(lineas, collapse = "\n")), show_col_types = FALSE))
      } else d <- suppressMessages(read_csv(f, locale = locale(encoding = enc), show_col_types = FALSE))
      # Índice auxiliar no usado: readr actual infiere character en el CSV mal formado.
      # Coerción explícita para que bind_rows no falle; no altera nodos/modelos.
      d[[1]] <- as.character(d[[1]])
      d$test <- regmatches(basename(f), regexpr("Acciones|Animales|Fpalabra", basename(f)))
      d$group <- sub(".*(control|off|on)\\.graphml\\.csv$", "\\1", basename(f))
      d$network <- net
      d
    }))
  }))
}
merge_historico <- function(nodos) {
  dic_legacy <- c(dic[c("logFreq", "motor", "imageability", "concreteness", "familiarity")],
    list(AoA = bind_rows(dic$aoa_verbs, dic$aoa_nouns)))
  etapas <- tibble(etapa = "nodos", n = nrow(nodos))
  for (nm in names(dic_legacy)) {
    nodos <- merge(nodos, dic_legacy[[nm]], by.x = "name", by.y = "palabra", all.x = TRUE)
    etapas <- bind_rows(etapas, tibble(etapa = nm, n = nrow(nodos)))
  }
  list(datos = nodos, etapas = etapas)
}
ajustar <- function(d) {
  d$group <- factor(d$group, levels = c("control", "off", "on"))
  d$hs_partition <- factor(d$hs_partition, levels = c(0, 1))
  bind_rows(lapply(variables_psico, function(var) {
    # Misma fórmula, family gaussian, filtrado de NA y summary/glm que la R2.
    # No hay modelo motor de animales en los Excel históricos (34 modelos).
    d |> filter(!is.na(.data[[var]]), !(var == "motor" & test == "Animales")) |> group_by(network, test) |>
      group_modify(~ {
        m <- tryCatch(glm(as.formula(paste(var, "~ group + hs_partition + group:hs_partition")),
          data = .x, family = "gaussian"), error = function(e) stop(paste(var, paste(.y, collapse = "/"),
            paste(unique(.x$group), collapse = "/"), paste(unique(.x$hs_partition), collapse = "/"), conditionMessage(e))))
        tidy(m) |> mutate(n = nobs(m), df_residual = df.residual(m))
      }) |> ungroup() |> mutate(DependentVar = var)
  }))
}
estrellas <- function(p) ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "")))
presentar <- function(r) {
  r |> mutate(
    network = recode(network, directed = "Sequence", undirected = "Co-occurrence"),
    test = recode(test, Acciones = "Actions fluency", Animales = "Semantic fluency", Fpalabra = "Phonological fluency"),
    DependentVar = recode(DependentVar, logFreq = "Log Frequency", motor = "Motor Content",
      imageability = "Imageability", concreteness = "Concreteness", familiarity = "Familiarity", AoA = "Age of acquisition"),
    term = recode(term, `(Intercept)` = "Intercept", groupoff = "Off", groupon = "On", hs_partition1 = "Periphery",
      `groupoff:hs_partition1` = "Off:Periphery", `groupon:hs_partition1` = "On:Periphery"),
    estimate_signif = paste0(round(estimate, 3), " ", estrellas(p.value)))
}
describir <- function(d) {
  d |> pivot_longer(all_of(variables_psico), names_to = "variable", values_to = "valor") |>
    group_by(network, test, group, hs_partition, variable) |>
    summarise(n = sum(!is.na(valor)), mean = mean(valor, na.rm = TRUE), sd = sd(valor, na.rm = TRUE),
      se = sd / sqrt(n), .groups = "drop")
}
utf <- leer_nodos("UTF-8"); lat <- leer_nodos("latin1")
mixto <- leer_nodos("auto"); reparado <- leer_nodos("auto", reparar = TRUE)
sin_shy <- reparado |> mutate(name_csv = name, name = gsub("\u00ad", "", name, fixed = TRUE))
reparaciones <- sin_shy |> filter(name != name_csv) |> select(network, test, group, node_id, name_csv, name)
graphml <- "../parkinson/networks/graphAnimales_directedcontrol.graphml"
xml <- paste(read_lines(graphml), collapse = "\n")
stopifnot(nrow(reparaciones) == 5L,
  all(vapply(reparaciones$name, function(w) grepl(paste0(">", w, "</data>"), xml, fixed = TRUE), logical(1))),
  !anyNA(reparado$node_id), !anyNA(reparado$name), !anyNA(reparado$hs_partition))
stopifnot(nrow(lat) == 2856L)
u <- merge_historico(utf); l <- merge_historico(lat); m <- merge_historico(mixto); r <- merge_historico(reparado)
h <- merge_historico(sin_shy)
# Mantener las particiones originales: limpiar etiquetas NO equivale a volver a ajustar SBM.
fijo <- anotar_limpio(mutate(sin_shy, word_norm = name), dic, FALSE)
limpio <- anotar_limpio(mutate(sin_shy, word_norm = normalizar(name)), dic, TRUE)
stopifnot(nrow(fijo) == 2856L, nrow(limpio) == 2856L)
colisiones <- limpio |> group_by(network, test, group, word_norm) |>
  summarise(n = n(), partitions = n_distinct(hs_partition), .groups = "drop") |> filter(n > 1)
# No se colapsan las ocho parejas: fuego/fuego+U+00A0 tiene particiones
# distintas en F-OFF-dirigida. Resolverla exige volver a ajustar el SBM (P12).
stopifnot(nrow(colisiones) == 8L, sum(colisiones$partitions > 1L) == 1L)
escenarios <- list(original_utf8 = u$datos, original_latin1 = l$datos, original_mixto = m$datos,
  csv_reparado = r$datos,
  csv_sin_soft_hyphen = h$datos,
  union_corregida = fijo, union_y_limpieza = limpio)
modelos <- bind_rows(lapply(names(escenarios), function(s) ajustar(escenarios[[s]]) |> mutate(escenario = s)))
descriptivos <- bind_rows(lapply(names(escenarios), function(s) describir(escenarios[[s]]) |> mutate(escenario = s)))
ref <- read_excel("../parkinson/resultsCorePeriphery/results_table.xlsx") |>
  pivot_longer(-c(network, test, DependentVar), names_to = "term", values_to = "publicado")
comparacion <- presentar(modelos) |> left_join(ref, by = c("network", "test", "DependentVar", "term")) |>
  mutate(coincide = trimws(estimate_signif) == trimws(publicado),
    publicado_num = as.numeric(sub(" .*", "", trimws(publicado))), delta_publicado = round(estimate, 3) - publicado_num)
resumen <- comparacion |> group_by(escenario) |> summarise(celdas = n(), coinciden = sum(coincide, na.rm = TRUE),
  max_diferencia = max(abs(delta_publicado), na.rm = TRUE), .groups = "drop")
ref_glm <- read_excel("../parkinson/resultsCorePeriphery/tabla_resumen_GLM.xlsx")
precision <- modelos |> left_join(ref_glm |> rename_with(~paste0(.x, "_guardado"),
    -all_of(c("network", "test", "DependentVar", "term"))), by = c("network", "test", "DependentVar", "term")) |>
  mutate(delta_estimate = estimate - estimate_guardado, delta_p = p.value - p.value_guardado)
impacto <- modelos |> select(network, test, DependentVar, term, escenario, estimate, std.error, p.value, n) |>
  pivot_wider(names_from = escenario, values_from = c(estimate, std.error, p.value, n)) |>
  mutate(delta_coef_csv = estimate_csv_sin_soft_hyphen - estimate_original_mixto,
    delta_coef_union = estimate_union_corregida - estimate_csv_sin_soft_hyphen,
    delta_coef_limpieza = estimate_union_y_limpieza - estimate_union_corregida,
    cambio_p05_union = (p.value_csv_sin_soft_hyphen < .05) != (p.value_union_corregida < .05),
    cambio_p05_total = (p.value_csv_sin_soft_hyphen < .05) != (p.value_union_y_limpieza < .05),
    delta_coef_total = estimate_union_y_limpieza - estimate_csv_sin_soft_hyphen)
ref_desc <- read_excel("../parkinson/resultsCorePeriphery/tabla_resumen_corePeriphery.xlsx") |>
  pivot_longer(matches("^(mean|sd|se)_"), names_to = c("estadistico", "variable"),
    names_pattern = "(mean|sd|se)_(.*)", values_to = "guardado")
precision_desc <- descriptivos |> pivot_longer(c(mean, sd, se), names_to = "estadistico", values_to = "calculado") |>
  left_join(ref_desc, by = c("network", "test", "group", "hs_partition", "variable", "estadistico")) |>
  mutate(delta = calculado - guardado)
validacion <- precision |> filter(escenario == "csv_sin_soft_hyphen")
validacion_desc <- precision_desc |> filter(escenario == "csv_sin_soft_hyphen", !is.na(delta))
stopifnot(nrow(validacion) == 204L,
  all(abs(validacion$delta_estimate) < 1e-10), all(abs(validacion$delta_p) < 1e-10),
  all(comparacion$coincide[comparacion$escenario == "csv_sin_soft_hyphen"]),
  nrow(validacion_desc) == 612L, all(abs(validacion_desc$delta) < 1e-10))
impacto_desc <- descriptivos |> select(-se) |>
  pivot_wider(names_from = escenario, values_from = c(n, mean, sd)) |>
  mutate(delta_media = mean_union_y_limpieza - mean_csv_sin_soft_hyphen,
    delta_sd = sd_union_y_limpieza - sd_csv_sin_soft_hyphen)
tabla2_antes_despues <- impacto |> transmute(network, test, DependentVar, term,
  publicado_manuscrito = !(DependentVar == "motor" & test == "Fpalabra"),
  coef_R2 = estimate_csv_sin_soft_hyphen, coef_corregido = estimate_union_y_limpieza,
  delta_coef = delta_coef_total, SE_R2 = std.error_csv_sin_soft_hyphen, SE_corregido = std.error_union_y_limpieza,
  p_R2 = p.value_csv_sin_soft_hyphen, p_corregido = p.value_union_y_limpieza,
  n_R2 = n_csv_sin_soft_hyphen, n_corregido = n_union_y_limpieza, cambio_p05 = cambio_p05_total)
tabla3_antes_despues <- impacto_desc |> filter(!is.na(hs_partition)) |> transmute(network, test, group, hs_partition, variable,
  n_R2 = n_csv_sin_soft_hyphen, n_corregido = n_union_y_limpieza,
  media_R2 = mean_csv_sin_soft_hyphen, media_corregida = mean_union_y_limpieza,
  DE_R2 = sd_csv_sin_soft_hyphen, DE_corregida = sd_union_y_limpieza, delta_media)
filas <- bind_rows(lapply(names(escenarios), function(s) escenarios[[s]] |>
  count(network, test, group, name = "n_filas") |> mutate(escenario = s)))
inflacion <- filas |> filter(escenario %in% c("csv_sin_soft_hyphen", "union_corregida", "union_y_limpieza")) |>
  pivot_wider(names_from = escenario, values_from = n_filas) |>
  mutate(filas_extra = csv_sin_soft_hyphen - union_corregida,
    factor_inflacion = csv_sin_soft_hyphen / union_corregida,
    porcentaje_filas_extra = 100 * filas_extra / csv_sin_soft_hyphen)
encoding_csv <- tibble(path = list.files("../parkinson/hub and spoke model", pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)) |>
  mutate(encoding = vapply(path, function(f) if(is.na(iconv(rawToChar(read_file_raw(f)), from = "UTF-8", to = "UTF-8"))) "latin1" else "UTF-8", character(1)))
out <- nuevo_directorio("02_reproducir_R2")
for (nm in c("modelos", "descriptivos", "comparacion", "resumen", "precision", "impacto", "filas", "colisiones", "precision_desc", "impacto_desc", "reparaciones", "inflacion", "encoding_csv", "tabla2_antes_despues", "tabla3_antes_despues")) {
  guardar_csv(get(nm), file.path(out, paste0(nm, ".csv")))
}
guardar_csv(bind_rows(mutate(u$etapas, encoding = "UTF-8"), mutate(l$etapas, encoding = "latin1")), file.path(out, "etapas_merge.csv"))
guardar_csv(bind_rows(mutate(m$etapas, escenario = "original_mixto"), mutate(r$etapas, escenario = "csv_reparado")), file.path(out, "etapas_merge_mixto.csv"))
saveRDS(escenarios, file.path(out, "nodos_escenarios.rds"))
entradas <- c("../parkinson/corePeriphery.R", list.files("../parkinson/data", pattern = "\\.xlsx$", full.names = TRUE),
  list.files("../parkinson/hub and spoke model", pattern = "\\.csv$", recursive = TRUE, full.names = TRUE),
  "../parkinson/resultsCorePeriphery/results_table.xlsx", "../parkinson/resultsCorePeriphery/tabla_resumen_GLM.xlsx")
entradas <- c(entradas, graphml, "../parkinson/resultsCorePeriphery/tabla_resumen_corePeriphery.xlsx")
registrar_entorno(out, entradas)
writeLines(c(
  "# Reproducción R2 — 22-10-2026",
  "La lógica de merge/glm histórica se mantiene; solo se omiten gráficos y se tipa como texto el índice auxiliar del CSV.",
  "Hay 16 CSV UTF-8 y 2 Latin-1. No se debe forzar Latin-1 a todos.",
  "Una fila completamente entrecomillada en acciones-control-dirigida se recupera en memoria.",
  "Cinco nombres de animales-control-dirigida contienen U+00AD; se verifica su forma sin ese carácter en el GraphML original.",
  "Escenario csv_sin_soft_hyphen: 204/204 celdas de results_table.xlsx, con coeficientes y asteriscos idénticos.",
  sprintf("Precisión frente a tabla_resumen_GLM.xlsx: error máximo de coeficientes %.3g; p %.3g.", max(abs(validacion$delta_estimate)), max(abs(validacion$delta_p))),
  sprintf("Tabla 3: 612 medias/DE/EE no ausentes coinciden; error máximo %.3g.", max(abs(validacion_desc$delta))),
  "El Excel histórico contiene 34 modelos (204 términos): 32 del manuscrito y dos de motor en F. No contiene motor en animales.",
  "El bloque GLM original intentaría motor en animales y falla por insuficiencia de niveles; aquí se limita a los modelos guardados.",
  "Las uniones históricas pasan de 2.856 a 3.500 filas (14 adicionales tras Lexesp y otras 630 tras AoA).",
  "union_corregida: D-005 provisional, una fila por node_id (2.856); Lexesp promediado por palabra/variable.",
  "union_y_limpieza: D-006/D-008, 2.856 node_id conservados. Hay 2.848 claves normalizadas distintas por red/tarea/grupo.",
  "Ocho parejas colisionan al normalizar: pensar, figura, fuera y fuego en ambas redes. En F-OFF-dirigida, fuego tiene dos particiones distintas.",
  "No se colapsan los node_id ni se asigna una partición inventada: el resultado es una sensibilidad con particiones R2 fijas, pendiente del SBM limpio.",
  sprintf("Cruces de p=0,05 tras corregir: %d/204 términos. Son p sin corrección múltiple del modelo histórico.", sum(impacto$cambio_p05_total)),
  "impacto.csv compara coeficientes, errores, p y n de todos los escenarios; impacto_desc.csv compara medias/DE/n.",
  "Las correcciones son una auditoría/sensibilidad y no resuelven la dependencia de palabras/personas ni el pareado ON/OFF.",
  "El token partir<TAB>llegar se conserva unido hasta respuesta clínica; no se inventan dos tokens.",
  "No se ha reconstruido el SBM, ni calculado modelos nuevos ni modificado el manuscrito."
), file.path(out, "LEEME.md"), useBytes = TRUE)
print(resumen); print(bind_rows(mutate(u$etapas, encoding = "UTF-8"), mutate(l$etapas, encoding = "latin1")))
cat("Salida:", out, "\n")
