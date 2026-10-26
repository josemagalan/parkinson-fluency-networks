# Propósito: cálculo rápido de métricas no ponderadas en redes de listas seleccionadas.
# Entradas: listas de índices de palabra, selección de posiciones, tamaño del vocabulario.
# Salida: seis métricas; E1. Paquetes: igraph. No escribe; seleccionar_posiciones consume RNG (sample.int).
# --- English ---
# Purpose: fast unweighted network metrics (nodes, edges, density, transitivity, mean degree, largest
#   component) for networks built from selected list positions (E1). Functions only; seleccionar_posiciones() uses RNG.
red_metricas <- function(ws, ix, V, modo) {
  selected <- unique(unlist(Map(function(w, i) w[i], ws, ix), use.names = FALSE))
  ed <- Map(function(w, i) {
    if (length(i) < 2L) return(matrix(integer(), ncol = 2))
    if (modo == "cooccurrence") return(t(combn(w[i], 2)))
    e <- cbind(w[head(i, -1L)], w[tail(i, -1L)])
    if (modo == "sequence_original") e <- e[diff(i) == 1L, , drop = FALSE]
    e
  }, ws, ix)
  ed <- do.call(rbind, ed)
  g <- igraph::make_empty_graph(V, directed = modo != "cooccurrence")
  if (nrow(ed)) g <- igraph::add_edges(g, as.vector(t(ed)))
  g <- igraph::induced_subgraph(g, selected)
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  c(nodes = igraph::vcount(g), edges = igraph::ecount(g),
    density = igraph::edge_density(g, loops = FALSE),
    transitivity = igraph::transitivity(igraph::as_undirected(g, mode = "collapse"), type = "global", isolates = "zero"),
    mean_degree_total = mean(igraph::degree(g, mode = "all")),
    largest_weak_component_pct = 100 * max(igraph::components(g, mode = "weak")$csize) / igraph::vcount(g))
}
seleccionar_posiciones <- function(ws, k) {
  stopifnot(length(ws) == length(k), all(lengths(ws) >= k), all(k >= 2L))
  lapply(seq_along(ws), function(i) sort(sample.int(length(ws[[i]]), k[i])))
}
