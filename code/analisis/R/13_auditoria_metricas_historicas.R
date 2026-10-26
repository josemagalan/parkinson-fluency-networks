# Propósito: comprobar métricas de GraphML históricos frente a sus Excel y aclarar pesos/distancias.
# Entradas: funciones puras de mainParkinson*.R, GraphML y tablas Excel (solo lectura).
# Salidas: output/13_auditoria_metricas_historicas[_vN]/, comparación y diagnóstico.
# E1/E5: sin SBM ni escrituras en repositorio histórico; paquetes: 00_utilidades.R + igraph.
# --- English ---
# Purpose: recompute R2 network metrics (function taken from mainParkinson*.R) on the R2 GraphML files
#   and compare with the R2 Excel tables; path lengths under three weightings (E1, E5).
# Inputs: ../parkinson GraphML, Excel and scripts. Outputs: analisis/output/13_auditoria_metricas_historicas[_vN]/.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
suppressPackageStartupMessages(library(igraph))
set.seed(20261022)
extraer_funcion <- function(path, nombre) {
  a <- as.list(parse(path, encoding = "UTF-8"))
  f <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name(nombre)), a)
  stopifnot(length(f) == 1L)
  e <- new.env(parent = globalenv()); eval(f[[1]], e); get(nombre, envir = e)
}
res <- distancias <- list()
fuentes <- c("../parkinson/mainParkinson.R", "../parkinson/mainParkinson_Dirigidas.R")
for (net in c("undirected", "directed")) {
  f <- fuentes[if (net == "directed") 2L else 1L]
  funcion <- extraer_funcion(f, "red_caracterizacion")
  xls <- if (net == "directed") "../parkinson/results/tabla redes dirigidas.xlsx" else "../parkinson/results/tabla redes no dirigidas.xlsx"
  tabla <- read_excel(xls)
  for (i in seq_len(nrow(tabla))) {
    path <- paste0("../parkinson/", tabla$nom_grafo[i])
    g <- read_graph(path, format = "graphml")
    r <- suppressWarnings(funcion(g))
    for (nm in setdiff(names(tabla), "nom_grafo")) {
      saved <- as.numeric(tabla[[nm]][i]); calc <- as.numeric(r[[nm]])
      res[[paste(net, i, nm)]] <- tibble(network = net, graph = tabla$nom_grafo[i], metric = nm,
        saved = saved, recalculated = calc, difference = calc - saved,
        matches = (is.na(saved) && is.na(calc)) || isTRUE(all.equal(saved, calc, tolerance = 1e-10)))
    }
    for (strategy in c("frequency_as_length_legacy", "unweighted", "inverse_frequency")) {
      w <- switch(strategy, frequency_as_length_legacy = E(g)$weight, unweighted = NA,
        inverse_frequency = 1 / E(g)$weight)
      distancias[[paste(net, i, strategy)]] <- tibble(network = net, graph = tabla$nom_grafo[i], strategy = strategy,
        mean_distance = mean_distance(g, directed = is_directed(g), weights = w, unconnected = TRUE),
        diameter = diameter(g, directed = is_directed(g), weights = w, unconnected = TRUE),
        mean_edge_weight = mean(E(g)$weight), mean_weighted_strength = mean(strength(g, mode = "all")))
    }
  }
}
res <- bind_rows(res); distancias <- bind_rows(distancias)
out <- nuevo_directorio("13_auditoria_metricas_historicas")
guardar_csv(res, file.path(out, "comparacion.csv")); guardar_csv(distancias, file.path(out, "distancias_pesos.csv"))
files <- c(fuentes, "../parkinson/results/tabla redes dirigidas.xlsx", "../parkinson/results/tabla redes no dirigidas.xlsx",
  paste0("../parkinson/", unique(res$graph)), "analisis/R/13_auditoria_metricas_historicas.R")
registrar_entorno(out, files)
writeLines(c("# Auditoría de métricas históricas",
  "Solo se evalúa la función pura red_caracterizacion, nunca los scripts históricos completos con escrituras.",
  "Comparación contra los GraphML guardados en results y sus Excel; no demuestra por sí sola qué versión usó cada figura/partición.",
  "avg_degree_weighted del código es mean(E(g)$weight): peso medio de arista, no grado/strength medio de nodo.",
  "average.path.length y diameter usan por defecto el atributo weight como longitud; los pesos de frecuencia mayores se interpretan como distancias mayores.",
  "Se guardan variantes explícitas sin pesos e inverso de frecuencia; ninguna se incorpora al manuscrito sin explicar el estimando.",
  "Comprobar diferencias de exportación/versiones; no sustituir resultados silenciosamente. Sin SBM ni cambios en originales."), file.path(out, "LEEME.md"))
print(res |> group_by(network) |> summarise(n_cells = n(), n_matches = sum(matches), .groups = "drop"))
print(res |> filter(!matches), n = 40); cat("Salida:", out, "\n")
