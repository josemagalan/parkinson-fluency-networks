# Propósito: paquete local de entrada para integrar el código de Virginia cuando llegue.
# Entradas: 18 GraphML normalizados de 05 y mapa de particiones históricas de 02.
# Salidas: output/14_preparar_entrada_SBM[_vN]/, GraphML, nodos/aristas y mapa auditado.
# E1/E2/E5: sin ejecutar SBM ni trasladar automáticamente particiones a nodos fusionados.
# Paquetes: 00_utilidades.R + igraph. Sin RNG salvo semilla documental.
# --- English ---
# Purpose: package the 18 normalised networks as SBM input and map R2 partitions for audit only
#   (never transferred) (E1, E2, E5). No SBM is fitted.
# Inputs: 05 GraphML, 02 node map. Outputs: analisis/output/14_preparar_entrada_SBM[_vN]/ (GraphML, nodes, edges).
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
suppressPackageStartupMessages(library(igraph))
set.seed(20261022)
directory <- "analisis/output/05_redes_igualadas_preparacion/graphml_completos"
files <- sort(list.files(directory, pattern = "graphml$", full.names = TRUE))
stopifnot(length(files) == 18L)
old <- readRDS("analisis/output/02_reproducir_R2/nodos_escenarios.rds")$union_y_limpieza |>
  transmute(task = recode(test, Acciones = "actions", Animales = "animals", Fpalabra = "F"),
    condition = recode(group, control = "HC", off = "OFF", on = "ON"), network,
    old_node_id = node_id, old_name_csv = name_csv, repaired_name = name, word_norm,
    historical_partition = hs_partition)
nodes <- edges <- index <- list()
for (f in files) {
  name <- sub("\\.graphml$", "", basename(f)); info <- strsplit(name, "_", fixed = TRUE)[[1]]
  stopifnot(length(info) == 3L)
  g <- read_graph(f, format = "graphml")
  nodes[[name]] <- tibble(graph = name, task = info[1], condition = info[2], network = info[3],
    new_node_id = paste0("v", seq_len(vcount(g))), word_norm = V(g)$name)
  e <- as_data_frame(g, what = "edges")
  edges[[name]] <- as_tibble(e) |> mutate(graph = name, is_loop = from == to, .before = 1)
  index[[name]] <- tibble(graph = name, task = info[1], condition = info[2], network = info[3],
    nodes = vcount(g), edges_including_loops = ecount(g), loops = sum(which_loop(g)),
    HC_provisional_P11 = info[2] == "HC" && info[1] == "actions")
}
nodes <- bind_rows(nodes); edges <- bind_rows(edges); index <- bind_rows(index)
mapa_historico <- left_join(old, nodes, by = c("task", "condition", "network", "word_norm"), relationship = "many-to-one")
stopifnot(nrow(nodes) == 2848L, nrow(mapa_historico) == 2856L, !anyNA(mapa_historico$new_node_id))
colisiones <- mapa_historico |> group_by(graph, new_node_id, word_norm) |>
  summarise(n_old_nodes = n(), n_historical_partitions = n_distinct(historical_partition), .groups = "drop") |>
  filter(n_old_nodes > 1L)
stopifnot(nrow(colisiones) == 8L, sum(colisiones$n_historical_partitions > 1) == 1L)
out <- nuevo_directorio("14_preparar_entrada_SBM"); dir.create(file.path(out, "graphml"))
ok <- file.copy(files, file.path(out, "graphml", basename(files)))
stopifnot(all(ok), identical(unname(tools::md5sum(files)), unname(tools::md5sum(file.path(out, "graphml", basename(files))))))
for (nm in c("nodes", "edges", "index", "mapa_historico", "colisiones")) guardar_csv(get(nm), file.path(out, paste0(nm, ".csv")))
writeLines(c("# Entrada para integrar SBM — paquete local, no enviado",
  "18 redes completas normalizadas: tres tareas × tres condiciones × dos tipos. Revisar index.csv.",
  "GraphML idénticos por MD5 a salida 05; nodes.csv/edges.csv conservan etiquetas y pesos. new_node_id solo identifica esta exportación.",
  "Los GraphML contienen frecuencia como peso y algunos bucles; NO asumir que el SBM admite esos atributos sin revisar su implementación.",
  "Confirmar: versión de librería/dependencias, semilla, grafos dirigidos/no dirigidos, tratamiento de peso/bucle, número de capas y criterion de comparación.",
  "Guardar por red el vector de partición, longitud de descripción de cada modelo y configuración; no mezclar longitudes calculadas sobre entradas distintas.",
  "mapa_historico.csv sirve para auditoría, no para asignar particiones nuevas. Ocho colisiones de etiquetas; fuego F-OFF-dirigida tiene particiones históricas distintas.",
  "HC-acciones contiene una celda ambigua (P11); validar antes de cerrar resultados. Las redes submuestreadas se regeneran con scripts 05/07.",
  "Pendiente de integrar código de Virginia; no se ha ajustado ni seleccionado ningún modelo SBM.",
  "No se ha enviado este paquete ni modificado el manuscrito."), file.path(out, "LEEME_Virginia.md"), useBytes = TRUE)
registrar_entorno(out, c(files, "analisis/output/02_reproducir_R2/nodos_escenarios.rds", "analisis/R/14_preparar_entrada_SBM.R"))
print(index |> count(network)); print(colisiones); cat("Salida:", out, "\n")
