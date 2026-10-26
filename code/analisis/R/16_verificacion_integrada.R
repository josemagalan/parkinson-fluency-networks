# Propósito: comprobaciones cruzadas de los bloques independientes del SBM.
# Entradas: scripts R numerados y outputs 03–15; fuentes históricas solo lectura.
# Salidas: output/16_verificacion_integrada[_vN], aserciones y huellas contrastadas.
# E1–E5: verifica parejas, planes, métricas, multiplicidad y ausencia de fuga.
# Paquetes: 00_utilidades.R, igraph. No vuelve a entrenar ni elige resultados.
# --- English ---
# Purpose: deterministic cross-checks of outputs 03-15 and of every recorded input MD5; stops on failure
#   (E1-E5). Note: it compares the MD5 of the scripts recorded when each output was produced.
# Inputs: outputs 03-15 and analisis/output/*/entradas_md5.csv. Outputs: analisis/output/16_verificacion_integrada[_vN]/.
if (.Platform$OS.type=="windows") invisible(Sys.setlocale("LC_CTYPE","English_United States.utf8"))
source("analisis/R/00_utilidades.R",encoding="UTF-8")
source("analisis/R/00_redes_metricas.R",encoding="UTF-8")
set.seed(20261022)
checks <- list()
check <- function(name, value, detail="") {stopifnot(isTRUE(value)); checks[[length(checks)+1L]] <<- tibble(check=name,passed=value,detail=detail)}
scripts <- list.files("analisis/R",pattern="\\.R$",full.names=TRUE)
for(f in scripts) parse(f,encoding="UTF-8")
check("Sintaxis de todos los scripts R",TRUE,as.character(length(scripts)))
t <- readRDS("analisis/output/03_orden_sesiones/tokens.rds")
check("Dataset y listas",nrow(t)==2965L && n_distinct(t$participant)==45L && nrow(distinct(t,participant,condition,task))==198L)
pd <- t |> filter(cohort=="PD") |> distinct(participant,condition,session,order)
check("21 pares completos con periodo",nrow(pd)==42L && all(table(pd$participant)==2L) && !anyNA(pd))
q <- pd |> distinct(participant,order)
check("11 OFF_ON y 10 ON_OFF",sum(q$order=="OFF_ON")==11L && sum(q$order=="ON_OFF")==10L)
# Recalcular las seis métricas del helper rápido y compararlas a las exportadas
# por el constructor independiente de GraphML de 05 (18 grafos completos).
m <- read_csv("analisis/output/05_redes_igualadas_preparacion/full_metrics.csv",show_col_types=FALSE)
for(i in seq_len(nrow(m))) {
  s <- t |> filter(task==m$task[i],condition==m$condition[i]) |> arrange(participant,position)
  vocab <- unique(s$word_norm); ws <- split(match(s$word_norm,vocab),s$participant)
  mode <- if(m$network[i]=="directed") "sequence_rebuilt" else "cooccurrence"
  v <- red_metricas(ws,lapply(ws,seq_along),length(vocab),mode)
  expected <- unlist(m[i,c("nodes","edges","density","transitivity_undirected","mean_degree_total","largest_weak_component_pct")])
  check(paste("Métricas independientes",m$task[i],m$condition[i],m$network[i]),max(abs(v-expected))<1e-10)
}
# Portabilidad (22-10-2026): la salida vigente de 07 es _v2 en la máquina de
# referencia; en una ejecución limpia (repositorio) es la carpeta sin sufijo.
d07 <- if (dir.exists("analisis/output/07_redes_HC_bootstrap_v2")) "analisis/output/07_redes_HC_bootstrap_v2" else "analisis/output/07_redes_HC_bootstrap"
r <- readRDS(file.path(d07,"simulaciones.rds"))
b <- read_csv(file.path(d07,"presupuestos.csv"),show_col_types=FALSE)
check("11.000 extracciones y 594.000 métricas",nrow(b)==11000L && nrow(r$metrics)==594000L && !anyNA(r$metrics$value))
plans <- readRDS(file.path(d07,"plan_personas.rds"))
check("Bootstrap estratificado por secuencia",all(apply(plans$PD_boot,1,function(x)sum(plans$PD$order[x]=="OFF_ON")==11L)))
check("HC condicional sin reemplazo",all(apply(plans$HC_conditional,1,function(x)length(unique(x))==21L)))
nb <- b |> filter(phase=="participant_bootstrap") |> select(task,draw,unique_PD,unique_HC) |> pivot_wider(names_from=task,values_from=c(unique_PD,unique_HC))
check("Mismas personas bootstrap en las tres tareas",all(nb$unique_PD_F==nb$unique_PD_actions & nb$unique_PD_F==nb$unique_PD_animals & nb$unique_HC_F==nb$unique_HC_actions & nb$unique_HC_F==nb$unique_HC_animals))
p6 <- read_csv("analisis/output/06_psicolinguistica_sesion/resultados.csv",show_col_types=FALSE)
check("BH26 y BH6 de medias",sum(p6$family=="psycholinguistic_26")==26L && nrow(p6)==32L)
for(f in unique(p6$family)) {v<-filter(p6,family==f);check(paste("BH recalculado",f),max(abs(p.adjust(v$p,"BH")-v$p_BH))<1e-12)}
p12 <- readRDS("analisis/output/12_permutacion_secuencia_sensibilidad/permutaciones.rds")
check("999 permutaciones de secuencia balanceadas",nrow(p12$assignments)==1000L && all(rowSums(p12$assignments)==11L))
p13 <- read_csv("analisis/output/13_auditoria_metricas_historicas/comparacion.csv",show_col_types=FALSE)
check("Auditoría histórica 233/234",nrow(p13)==234L && sum(p13$matches)==233L,"Radio acciones OFF dirigida: 6 guardado, 7 recalculado")
f15 <- read_csv("analisis/output/15_ML_agrupado/features.csv",show_col_types=FALSE)
pr <- read_csv("analisis/output/15_ML_agrupado/predictions.csv",show_col_types=FALSE)
check("ML doce evaluaciones completas LOPO",nrow(pr)==6L*66L+6L*42L && !anyNA(pr$predicted))
for(g in group_split(pr,endpoint,panel,algorithm)) {
  check(paste("ML folds indivisibles",g$endpoint[1],g$panel[1],g$algorithm[1]),
    all(vapply(split(g$fold,g$participant),function(v)length(unique(v))==1L,logical(1))))
  prob <- as.matrix(g[c("p_HC","p_OFF","p_ON")])
  check(paste("ML probabilidades",g$endpoint[1],g$panel[1],g$algorithm[1]),max(abs(rowSums(prob)-1))<1e-7 && min(prob)>=0)
}
null <- read_csv("analisis/output/15_ML_agrupado/null.csv",show_col_types=FALSE)
check("ML 199 nulos por algoritmo",nrow(null)==398L && all(table(null$algorithm)==199L))
pl15 <- readRDS("analisis/output/15_ML_agrupado/planes.rds")
check("ML sin metadatos como predictores",!any(c("participant","condition","session","order","fold") %in% unlist(pl15$sets)))
check("ML paneles fijos 3/107/131",identical(as.integer(lengths(pl15$sets)),c(3L,107L,131L)))
# Contrastar las entradas ya registradas por cada bloque. No inferir que un
# git dirty histórico pertenece a esta sesión; comparar huellas de archivos.
dirs <- list.dirs("analisis/output",recursive=FALSE,full.names=TRUE)
# Portabilidad: no contrastar manifiestos de ejecuciones previas de este mismo script.
dirs <- dirs[!grepl("^16_verificacion_integrada", basename(dirs))]
man <- bind_rows(lapply(dirs,function(dir) {
  f <- file.path(dir,"entradas_md5.csv"); if(!file.exists(f)) return(NULL)
  x <- read_csv(f,show_col_types=FALSE); if(!all(c("path","md5")%in%names(x)))return(NULL)
  x |> mutate(manifest=f,current_md5=unname(tools::md5sum(path)),
    # Un fichero opcional ausente en ambos momentos (p. ej., el correo privado) cuenta como sin cambios.
    unchanged=(!is.na(md5) & !is.na(current_md5) & md5==current_md5) | (is.na(md5) & is.na(current_md5)))
}))
check("Fuentes de manifiestos preservadas",all(man$unchanged),paste(nrow(man),"huellas"))
out <- nuevo_directorio("16_verificacion_integrada")
guardar_csv(bind_rows(checks),file.path(out,"comprobaciones.csv"));guardar_csv(man,file.path(out,"huellas.csv"))
registrar_entorno(out,c("analisis/R/16_verificacion_integrada.R",scripts))
writeLines(c("# Verificación integrada","Comprobaciones deterministas completadas; no equivalen a validación de supuestos clínicos ni de calibración estadística.",
  "No resuelve SBM, P3/P10/P11/P13/P14, ni convierte intervalos de grafos/OOF en IC exactos.","Figuras 06 y10 inspeccionadas visualmente en la sesión."),file.path(out,"LEEME.md"))
print(bind_rows(checks),n=100); cat("Salida:",out,"\n")
