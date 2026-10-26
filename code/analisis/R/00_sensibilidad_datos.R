# Propósito: escenarios explícitos para P3/P11; no escribe ni cambia datos base.
# Entradas: tokens tidy y diccionarios (solo al separar hipotéticamente la celda).
# Salidas: función escenario_datos. E1/E2/E4; paquetes de 00_utilidades.R.
# --- English ---
# Purpose: data scenarios for sensitivity checks P3/P11 (as recorded, within-list deduplicated, ambiguous
#   cell HC_22/actions/pos. 16 dropped or split). Function escenario_datos(); writes nothing. Source after 00_utilidades.R.
escenario_datos <- function(tokens, escenario = "recorded") {
  d <- tokens |> mutate(subposition = 0L)
  if (escenario == "deduplicated") {
    d <- d |> arrange(position) |> distinct(participant, task, condition, word_norm, .keep_all = TRUE)
  } else if (escenario %in% c("ambiguous_drop", "ambiguous_split")) {
    sel <- d$participant == "HC_22" & d$task == "actions" & d$position == 16L
    stopifnot(sum(sel) == 1L, d$word_norm[sel] == "partir\tllegar")
    if (escenario == "ambiguous_split") {
      nuevos <- bind_rows(d[sel, ], d[sel, ]) |>
        select(-all_of(variables_psico), -AoA_source) |>
        mutate(word_norm = c("partir", "llegar"), subposition = 1:2,
          token_id = paste0(token_id, ":hypothesis_split:", subposition))
      nuevos <- anotar_limpio(nuevos, leer_diccionarios())
      d <- bind_rows(d[!sel, ], nuevos)
    } else d <- d[!sel, ]
  } else stopifnot(escenario == "recorded")
  d |> arrange(participant, task, condition, position, subposition)
}
