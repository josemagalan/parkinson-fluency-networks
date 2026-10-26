# Propósito: robustez de la inferencia de la Tabla 2 (script 21) frente a la no normalidad de residuos (Shapiro p<0,05
#   en 21/32 modelos principales; 21_v2/diagnosticos_modelos.csv). Plan fijado antes de calcular (DECISIONES D-029):
#   (1) principal: errores típicos robustos por conglomerado (participante), CR2 con gl de Satterthwaite (clubSandwich;
#       Pustejovsky y Tipton, 2018). Las estimaciones puntuales no cambian; solo SE, gl, IC y p.
#   (2) secundario: bootstrap de participantes estratificado por cohorte (24 HC, 21 EP; cada paciente conserva sus dos
#       listas), reajustando el mismo modelo mixto; IC percentil 95 % y p bootstrap bilateral.
#   Mismas familias BH que 21 (22 modelos estimables por término; AoA aparte, 5; F secuencia n.i. fuera; intercepto sin test).
#   Se compara signo y decisión BH 0,05 con el modelo principal de 21.
# Entradas: analisis/output/21_tabla2_niveles_R2_v2/{niveles_persona,soporte_particiones,resultados}.csv.
# Salidas: analisis/output/45_tabla2_robustez_residuos[_vN]/ (robustez_CR2, robustez_bootstrap, comparacion_21, resumen, cambios_decision, LEEME).
# Puntos editoriales: E2/E5 (supuestos de los modelos de la Tabla 2). No modifica 21 ni sus salidas.
# Paquetes: 00_utilidades.R (readr, dplyr, tidyr, purrr), nlme, clubSandwich (>= 0.5.6; linear_contrast).
#   En Windows, si falta clubSandwich, instalarlo en analisis/R/library (como randomForest), no en la biblioteca global.
# Ejecutar desde revision_R3/: Rscript --vanilla analisis/R/45_tabla2_robustez_residuos.R
#   Variable de entorno opcional R3_B45 = número de réplicas bootstrap (por defecto 1999; solo para pruebas).
# --- English ---
# Purpose: robustness of Table 2 (script 21) to non-normal residuals (D-029): CR2 cluster-robust SEs
#   (clubSandwich) and stratified participant bootstrap (B = 1999), same BH families (E2, E5).
# Inputs: 21_tabla2_niveles_R2_v2/ CSVs. Outputs: analisis/output/45_tabla2_robustez_residuos[_vN]/.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
lib_r3 <- "analisis/R/library"
if (dir.exists(lib_r3)) .libPaths(c(normalizePath(lib_r3), .libPaths()))
suppressPackageStartupMessages(library(clubSandwich))
source("analisis/R/00_niveles_tabla2.R", encoding = "UTF-8")
set.seed(20261024)
B <- as.integer(Sys.getenv("R3_B45", "1999"))
d21 <- "analisis/output/21_tabla2_niveles_R2_v2"
f_nv <- file.path(d21, "niveles_persona.csv"); f_sop <- file.path(d21, "soporte_particiones.csv"); f_res <- file.path(d21, "resultados.csv")
nv <- read_csv(f_nv, show_col_types = FALSE)
soporte <- read_csv(f_sop, show_col_types = FALSE)
res21 <- read_csv(f_res, show_col_types = FALSE)
stopifnot(nrow(nv) == 4060L)

preparar <- function(d) d |> mutate(group = factor(condition, levels = c("HC", "OFF", "ON")),
  partition = factor(partition, c("core", "periphery")), ses = if_else(cohort == "HC", 0, session - 1),
  lista = paste(participant, condition, sep = ":"))
ctrl <- lmeControl(maxIter = 200, msMaxIter = 200, niterEM = 50, returnObject = FALSE)
# Mismo ajuste y misma regla de reserva que ajustar_uno() (00_niveles_tabla2.R).
ajustar_lme <- function(d) {
  f <- function(r) tryCatch(lme(mean ~ group * partition + ses, random = r, data = d, method = "REML", control = ctrl),
    error = function(e) NULL)
  m <- f(~ 1 | participant / lista); usado <- "anidado"
  if (is.null(m)) { m <- f(~ 1 | participant); usado <- "participante" }
  if (!is.null(m) && !all(terminos %in% names(fixef(m)))) m <- NULL
  list(m = m, usado = usado)
}
nombres <- c(names(terminos), names(contrastes))
matriz_L <- function(b) {
  L <- matrix(0, length(nombres), length(b), dimnames = list(nombres, names(b)))
  for (k in names(terminos)) L[k, terminos[[k]]] <- 1
  for (k in names(contrastes)) { w <- contrastes[[k]]; L[k, names(w)] <- w }
  L
}

# ---- 0. Comprobación: el reajuste desde niveles_persona.csv reproduce 21 -------------------------------------
reproduccion <- nv |> group_by(task, network, variable) |> group_modify(~ ajustar_uno(.x)) |> ungroup() |>
  filter(term != "Intercept_HC_core") |>
  inner_join(res21 |> select(task, network, variable, term, est21 = estimate, p21 = p), by = c("task", "network", "variable", "term"))
max_dif <- max(abs(reproduccion$estimate - reproduccion$est21), na.rm = TRUE)
max_dif_p <- max(abs(reproduccion$p - reproduccion$p21), na.rm = TRUE)
stopifnot(nrow(reproduccion) == 320L, max_dif < 1e-6, max_dif_p < 1e-6)

# ---- 1. CR2 + Satterthwaite -----------------------------------------------------------------------------------
cr2_uno <- function(d) {
  a <- ajustar_lme(preparar(d))
  if (is.null(a$m)) return(tibble(term = nombres, estimate = NA_real_, SE_CR2 = NA_real_, df_CR2 = NA_real_,
    p_CR2 = NA_real_, conf_low_CR2 = NA_real_, conf_high_CR2 = NA_real_, random_effects = "no estimable"))
  L <- matriz_L(fixef(a$m))
  lc <- linear_contrast(a$m, vcov = "CR2", contrasts = L, cluster = preparar(d)$participant, test = "Satterthwaite", p_values = TRUE)
  tibble(term = nombres, estimate = lc$Est, SE_CR2 = lc$SE, df_CR2 = lc$df, p_CR2 = lc$p_val,
    conf_low_CR2 = lc$CI_L, conf_high_CR2 = lc$CI_U, random_effects = a$usado)
}
cr2 <- nv |> group_by(task, network, variable) |> group_modify(~ cr2_uno(.x)) |> ungroup()

# ---- 2. Bootstrap de participantes estratificado por cohorte ---------------------------------------------------
est_uno <- function(d) {
  a <- ajustar_lme(d); if (is.null(a$m)) return(rep(NA_real_, length(nombres)))
  b <- fixef(a$m); drop(matriz_L(b) %*% b)
}
claves <- nv |> distinct(task, network, variable) |> arrange(task, network, variable)
ids <- nv |> distinct(participant, cohort)
ids_HC <- ids$participant[ids$cohort == "HC"]; ids_PD <- ids$participant[ids$cohort == "PD"]
stopifnot(length(ids_HC) == 24L, length(ids_PD) == 21L)
nv_p <- preparar(nv)
por_modelo <- split(nv_p, interaction(nv_p$task, nv_p$network, nv_p$variable, drop = TRUE, lex.order = TRUE))
# Todas las extracciones se generan antes de ajustar (el ajuste no consume RNG): el resultado no depende del número de
# procesos. En Unix se reparten con parallel::mclapply (R3_CORES45, por defecto 2); en Windows, en serie.
draws <- lapply(seq_len(B), function(r) c(sample(ids_HC, replace = TRUE), sample(ids_PD, replace = TRUE)))
una_replica <- function(r) {
  draw <- tibble(participant = draws[[r]]) |> mutate(nuevo = paste0(participant, "_b", row_number()))
  t(vapply(seq_len(nrow(claves)), function(j) {
    k <- paste(claves$task[j], claves$network[j], claves$variable[j], sep = ".")
    d <- inner_join(draw, por_modelo[[k]], by = "participant", relationship = "many-to-many") |>
      mutate(participant = nuevo, lista = paste(participant, condition, sep = ":"))
    est_uno(d)
  }, numeric(length(nombres))))
}
nucleos <- if (.Platform$OS.type == "windows") 1L else as.integer(Sys.getenv("R3_CORES45", "2"))
t0 <- Sys.time()
lista_boot <- if (nucleos > 1L) parallel::mclapply(seq_len(B), una_replica, mc.cores = nucleos) else lapply(seq_len(B), una_replica)
stopifnot(all(vapply(lista_boot, is.matrix, logical(1))))
boot <- array(NA_real_, c(B, nrow(claves), length(nombres)))
for (r in seq_len(B)) boot[r, , ] <- lista_boot[[r]]
resumen_boot <- bind_rows(lapply(seq_len(nrow(claves)), function(j) {
  bind_rows(lapply(seq_along(nombres), function(i) {
    x <- boot[, j, i]; ok <- x[!is.na(x)]; n_ok <- length(ok)
    tibble(task = claves$task[j], network = claves$network[j], variable = claves$variable[j], term = nombres[i],
      B_ok = n_ok, B_fail = B - n_ok,
      conf_low_boot = if (n_ok > 0) unname(quantile(ok, .025)) else NA_real_,
      conf_high_boot = if (n_ok > 0) unname(quantile(ok, .975)) else NA_real_,
      p_boot = if (n_ok > 0) min(1, 2 * min((sum(ok <= 0) + 1) / (n_ok + 1), (sum(ok >= 0) + 1) / (n_ok + 1))) else NA_real_)
  }))
}))
minutos <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)

# ---- 3. Familias BH (como 21) y comparación ---------------------------------------------------------------------
base21 <- res21 |> filter(term != "Intercept_HC_core") |>
  select(task, network, variable, term, estimate, p, p_BH, conf_low, conf_high, no_identificable, family)
bh <- function(x, p) x |> group_by(family) |> mutate(out = if_else(no_identificable, NA_real_, p.adjust(if_else(no_identificable, NA_real_, .data[[p]]), "BH"))) |>
  ungroup() |> pull(out)
comparacion <- base21 |>
  left_join(cr2 |> select(-estimate), by = c("task", "network", "variable", "term")) |>
  left_join(resumen_boot, by = c("task", "network", "variable", "term"))
comparacion$p_BH_CR2 <- bh(comparacion, "p_CR2")
comparacion$p_BH_boot <- bh(comparacion, "p_boot")
comparacion <- comparacion |> mutate(
  sig_main = !no_identificable & p_BH < .05, sig_CR2 = !no_identificable & p_BH_CR2 < .05, sig_boot = !no_identificable & p_BH_boot < .05,
  same_decision_CR2 = sig_main == sig_CR2, same_decision_boot = sig_main == sig_boot,
  boot_CI_excludes_0 = conf_low_boot > 0 | conf_high_boot < 0)
estim <- comparacion |> filter(!no_identificable)
resumen <- estim |> mutate(primary = !grepl("^AoA_", family)) |> group_by(primary, term) |>
  summarise(m = n(), sig_main = sum(sig_main), sig_CR2 = sum(sig_CR2), sig_boot = sum(sig_boot),
    cambios_CR2 = sum(!same_decision_CR2), cambios_boot = sum(!same_decision_boot),
    min_B_ok = min(B_ok), .groups = "drop") |> arrange(desc(primary), match(term, nombres))
cambios <- estim |> filter(!same_decision_CR2 | !same_decision_boot) |>
  select(task, network, variable, term, estimate, p_BH, p_BH_CR2, p_BH_boot, conf_low_boot, conf_high_boot)

out <- nuevo_directorio("45_tabla2_robustez_residuos")
guardar_csv(cr2, file.path(out, "robustez_CR2.csv"))
guardar_csv(resumen_boot, file.path(out, "robustez_bootstrap.csv"))
guardar_csv(comparacion, file.path(out, "comparacion_21.csv"))
guardar_csv(resumen, file.path(out, "resumen.csv"))
guardar_csv(cambios, file.path(out, "cambios_decision.csv"))
registrar_entorno(out, c(f_nv, f_sop, f_res, "analisis/R/45_tabla2_robustez_residuos.R", "analisis/R/00_niveles_tabla2.R",
  "analisis/R/00_utilidades.R"))
writeLines(c("# Robustez de la Tabla 2 frente a la no normalidad de residuos (D-029)",
  "Datos: 21_tabla2_niveles_R2_v2/niveles_persona.csv. Reajuste comprobado frente a 21: máx. dif. estimación " ,
  sprintf("  %.2e; máx. dif. p %.2e (320 términos).", max_dif, max_dif_p),
  "CR2: clubSandwich::linear_contrast, conglomerado = participante, gl Satterthwaite; estimaciones idénticas a 21.",
  sprintf("Bootstrap: %d réplicas de participantes estratificadas por cohorte (24 HC, 21 EP), reajuste del mismo lme y misma", B),
  "  regla de reserva; IC percentil 95 %; p bilateral = min(1, 2·min((#≤0+1)/(B+1), (#≥0+1)/(B+1))). Semilla 20261024.",
  sprintf("  Tiempo: %.1f min. Réplicas no estimables por modelo en robustez_bootstrap.csv (B_fail).", minutos),
  "Familias BH idénticas a 21 (22 estimables; AoA 5; F secuencia n.i. fuera; intercepto sin test).",
  "resumen.csv: significativos tras BH por término con el modelo principal, CR2 y bootstrap; cambios_decision.csv: discrepancias.",
  "Los p bootstrap tienen resolución 1/(B+1); con BH pueden ser conservadores en las familias grandes.",
  "Inferencia condicional a las particiones estimadas (como 21); la robustez no propaga la incertidumbre de la partición."),
  file.path(out, "LEEME.md"), useBytes = TRUE)
print(resumen, n = 50); print(cambios, n = 50); cat("Salida:", out, "| bootstrap", minutos, "min\n")
