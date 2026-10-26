# Propósito: exportar sin regenerar los datos y planes de 15 para TabPFN v2 local.
# Entradas: output/15_ML_agrupado/{features.csv,planes.rds,predictions.csv}.
# Salidas: output/40_preparar_TabPFN[_vN]: planes CSV y referencias de escala R.
# E3/E2: igualdad de participantes, predictores, folds y remuestreo entre modelos.
# Paquetes: los de 00_utilidades.R; no se ajustan modelos ni se instala nada.
# --- English ---
# Purpose (optional TabPFN add-on, D-026): export script 15 data, folds and resampling plans, plus R
#   preprocessing references to check the Python pipeline (E3). No model is fitted.
# Inputs: 15 features/plans/predictions. Outputs: analisis/output/40_preparar_TabPFN[_vN]/.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding="UTF-8")
set.seed(20261022)
base <- "analisis/output/15_ML_agrupado"
f <- read_csv(file.path(base,"features.csv"), show_col_types=FALSE)
p <- readRDS(file.path(base,"planes.rds"))
old <- read_csv(file.path(base,"predictions.csv"), show_col_types=FALSE)
stopifnot(nrow(f)==66L, n_distinct(f$participant)==45L,
          identical(lengths(p$sets), c(counts=3L,distributions=107L,distributions_AoA=131L)))
out <- nuevo_directorio("40_preparar_TabPFN")
guardar_csv(bind_rows(lapply(names(p$sets),function(k)
  tibble(panel=k,position=seq_along(p$sets[[k]]),feature=p$sets[[k]]))),file.path(out,"panels.csv"))
folds <- bind_rows(lapply(c("three_class","ON_OFF"),function(e) {
  dat <- if(e=="ON_OFF") filter(f,condition!="HC") else f
  tibble(endpoint=e,row=seq_len(nrow(dat)),participant=dat$participant,
         truth=dat$condition,fold=match(dat$participant,unique(dat$participant)))
}))
check <- old |> filter(panel=="distributions",algorithm=="random_forest") |>
  arrange(endpoint,row) |> select(endpoint,row,participant,truth,fold)
stopifnot(isTRUE(all.equal(as.data.frame(arrange(folds,endpoint,row)),as.data.frame(check),check.attributes=FALSE)))
guardar_csv(folds,file.path(out,"folds.csv"))
guardar_csv(bind_rows(lapply(seq_along(p$bootstrap_plans),function(b)
  tibble(replicate=b,position=seq_along(p$bootstrap_plans[[b]]),participant=p$bootstrap_plans[[b]]))),
  file.path(out,"bootstrap_plans.csv"))
guardar_csv(bind_rows(lapply(seq_along(p$grouped_CV_plans),function(b)
  tibble(repeat_id=b,participant=names(p$grouped_CV_plans[[b]]),fold=unname(p$grouped_CV_plans[[b]])))),
  file.path(out,"repeated_plans.csv"))
guardar_csv(bind_rows(lapply(seq_along(p$sequence_permutations),function(b)
  tibble(permutation=b,participant=p$sequence_permutations[[b]]))),file.path(out,"sequence_permutations.csv"))
# Referencia numérica independiente: primer paciente excluido de cada objetivo,
# los tres paneles; Python debe reproducir la escala R antes de ajustar TabPFN.
for(e in c("three_class","ON_OFF")) for(panel in names(p$sets)) {
  dat <- if(e=="ON_OFF") filter(f,condition!="HC") else f
  test <- which(dat$participant==dat$participant[grepl("^PD_",dat$participant)][1])
  train <- setdiff(seq_len(nrow(dat)),test)
  a <- as.matrix(dat[train,p$sets[[panel]],drop=FALSE])
  z <- as.matrix(dat[test,p$sets[[panel]],drop=FALSE])
  keep <- vapply(seq_len(ncol(a)),function(j) length(unique(a[is.finite(a[,j]),j]))>1L,logical(1))
  a <- a[,keep,drop=FALSE]; z <- z[,keep,drop=FALSE]
  for(j in seq_len(ncol(a))) {
    med <- median(a[is.finite(a[,j]),j]); a[!is.finite(a[,j]),j] <- med; z[!is.finite(z[,j]),j] <- med
  }
  mu <- colMeans(a); s <- apply(a,2,sd); s[!is.finite(s)|s==0] <- 1
  a <- sweep(sweep(a,2,mu),2,s,"/"); z <- sweep(sweep(z,2,mu),2,s,"/")
  xy <- bind_rows(bind_cols(tibble(role="train",row=train),as_tibble(a)),
                  bind_cols(tibble(role="test",row=test),as_tibble(z)))
  guardar_csv(xy,file.path(out,paste0("reference_",e,"_",panel,".csv")))
}
registrar_entorno(out,c(file.path(base,c("features.csv","planes.rds","predictions.csv")),
  "analisis/R/40_preparar_TabPFN.R","analisis/R/00_utilidades.R"))
writeLines(c("# Preparación TabPFN v2 — D-026",
  "Exportación de los planes guardados por15, no una nueva generación aleatoria.",
  "Folds verificados contra predictions.csv15. Datos66/45, ON/OFF42/21; paneles3/107/131.",
  "2000 bootstrap por persona estratificados PD/HC,20 CV agrupadas de5 folds,199 secuencias11/10.",
  "Las seis referencias numéricas permiten verificar la igualdad del preprocesamiento R/Python.",
  "No contiene resultados del nuevo modelo. AoA provisional; no implica adopción clínica."),
  file.path(out,"LEEME.md"),useBytes=TRUE)
cat("Preparación verificada:",out,"\n")
