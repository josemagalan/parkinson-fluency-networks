# Propósito: auditar las repeticiones Windows de 21/22 frente a Linux, sin cambiar modelos.
# Entradas: outputs 21/22 originales (Linux) y sus _v2 (Windows), plan D-024 y scripts.
# Salidas: output/44_verificar_21_22_Windows[_vN]: diferencias, BH, diagnósticos y resumen 22
#   enlazado al 21 Windows (22 sin modificar conserva como referencia el 21 original).
# E1/E2/E4/E5/E6; D-024/D-027. Paquetes: 00_utilidades.R (nlme no se usa). No se reestima el SBM.
# --- English ---
# Purpose (optional): compare Linux runs of 21/22 with their Windows re-runs (_v2): CSVs, BH decisions,
#   diagnostics and input MD5s; models are not refitted (E1, E2, E4-E6; D-024/D-027).
# Inputs: 21, 21_v2, 22, 22_v2 outputs. Outputs: analisis/output/44_verificar_21_22_Windows[_vN]/.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R",encoding="UTF-8")
set.seed(20261023)
dirs <- c(`21`="analisis/output/21_tabla2_niveles_R2",`22`="analisis/output/22_tabla2_niveles_igualadas")
leer <- function(f) read_csv(f,show_col_types=FALSE,trim_ws=FALSE)
comparaciones <- list(); columnas <- list(); usados <- character()
tol <- 1e-10 # umbral estricto de comparación; los cambios no se ocultan aumentando tolerancia.
for (bloque in names(dirs)) {
  ref <- dirs[[bloque]]; win <- paste0(ref,"_v2")
  archivos <- setdiff(list.files(ref,pattern="[.]csv$"),"entradas_md5.csv")
  for (f in archivos) {
    pa <- file.path(ref,f); pb <- file.path(win,f); usados <- c(usados,pa,pb)
    a <- leer(pa); b <- leer(pb)
    stopifnot(identical(names(a),names(b)),nrow(a)==nrow(b))
    nums <- names(a)[vapply(a,is.numeric,logical(1))]
    resto <- setdiff(names(a),nums)
    texto_ok <- all(vapply(resto,function(k) identical(a[[k]],b[[k]]),logical(1)))
    na_ok <- identical(is.na(as.data.frame(a)),is.na(as.data.frame(b)))
    cc <- bind_rows(lapply(nums,function(k) {
      delta <- abs(a[[k]]-b[[k]]); dd <- delta[is.finite(delta)]
      tibble(block=bloque,file=f,column=k,max_abs_diff=max(c(0,dd)),n_over_1e10=sum(dd>tol),
        same_NA=identical(is.na(a[[k]]),is.na(b[[k]])))
    }))
    if(!length(nums)) cc <- tibble(block=character(),file=character(),column=character(),
      max_abs_diff=numeric(),n_over_1e10=integer(),same_NA=logical())
    columnas[[length(columnas)+1L]] <- cc
    comparaciones[[length(comparaciones)+1L]] <- tibble(block=bloque,file=f,rows=nrow(a),
      same_nonnumeric=texto_ok,same_NA=na_ok,max_abs_diff=max(c(0,cc$max_abs_diff)),
      all_numeric_within_1e10=all(cc$n_over_1e10==0),
      all_equal=identical(as.data.frame(a),as.data.frame(b)))
    stopifnot(texto_ok,na_ok)
  }
}
comparacion <- bind_rows(comparaciones); por_columna <- bind_rows(columnas)
entradas <- bind_rows(lapply(names(dirs),function(b) bind_rows(lapply(c("Linux","Windows"),function(s) {
  dir <- paste0(dirs[[b]],if(s=="Windows")"_v2"else"")
  p <- file.path(dir,"entradas_md5.csv"); usados <<- c(usados,p)
  leer(p) |> mutate(block=b,platform=s,current_md5=unname(tools::md5sum(path)),matches=current_md5==md5)
}))))
stopifnot(all(entradas$matches))
l21 <- leer(file.path(dirs[["21"]],"resultados.csv"))
w21 <- leer(file.path(paste0(dirs[["21"]],"_v2"),"resultados.csv"))
l22 <- leer(file.path(dirs[["22"]],"resultados_por_extraccion.csv"))
w22 <- leer(file.path(paste0(dirs[["22"]],"_v2"),"resultados_por_extraccion.csv"))
decidir <- function(a,b,nombre,ajustada=FALSE) {
  keep <- !a$no_identificable & a$term!="Intercept_HC_core"
  pcol <- if(ajustada)"p_BH"else"p"
  tibble(block=nombre,n=sum(keep),changed_sign=sum(sign(a$estimate[keep])!=sign(b$estimate[keep]),na.rm=TRUE),
    changed_decision=sum((a[[pcol]][keep]<.05)!=(b[[pcol]][keep]<.05),na.rm=TRUE),
    changed_random_effects=sum(a$random_effects[keep]!=b$random_effects[keep]),
    max_estimate_diff=max(abs(a$estimate[keep]-b$estimate[keep]),na.rm=TRUE),
    max_p_diff=max(abs(a[[pcol]][keep]-b[[pcol]][keep]),na.rm=TRUE))
}
decisiones <- bind_rows(decidir(l21,w21,"21 principal BH",TRUE),decidir(l22,w22,"22 nominal descriptivo"))
ls <- leer(file.path(dirs[["21"]],"sensibilidades.csv"))
ws <- leer(file.path(paste0(dirs[["21"]],"_v2"),"sensibilidades.csv"))
sens <- tibble(sensitivity=ws$sensitivity,valid=!ws$no_identificable,
  decision_changed=(!is.na(ls$p_BH)&ls$p_BH<.05)!=(!is.na(ws$p_BH)&ws$p_BH<.05),
  sign_changed=sign(ls$estimate)!=sign(ws$estimate)) |>
  filter(valid) |> group_by(sensitivity) |>
  summarise(n=n(),changed_decision=sum(decision_changed),changed_sign=sum(sign_changed,na.rm=TRUE),.groups="drop")
stopifnot(all(decisiones$changed_decision==0),all(decisiones$changed_sign==0),
  all(decisiones$changed_random_effects==0),all(sens$changed_decision==0),all(sens$changed_sign==0))
bh <- w21 |> filter(!no_identificable,term!="Intercept_HC_core") |> group_by(family) |>
  summarise(m=n(),significant_BH=sum(p_BH<.05),positive=sum(estimate>0),
    max_BH_diff=max(abs(p_BH-p.adjust(p,"BH"))),.groups="drop")
stopifnot(all(bh$m==ifelse(grepl("^AoA_",bh$family),5L,22L)),all(bh$max_BH_diff<1e-12),
  all(is.na(w21$p_BH[w21$no_identificable])))
diag <- bind_rows(lapply(c("Linux","Windows"),function(s) {
  x <- leer(file.path(paste0(dirs[["21"]],if(s=="Windows")"_v2"else""),"diagnosticos_modelos.csv"))
  tibble(platform=s,models=nrow(x),converged=sum(x$converged_nested),shapiro_p05=sum(x$shapiro_p_resid<.05),
    max_normalised_residual=max(x$max_abs_norm_resid),n_abs_residual_gt3=sum(x$n_abs_resid_gt3))
}))
stopifnot(all(diag$models==32L),all(diag$converged==32L))
cero_gl <- ws |> filter(!is.na(df),df<=0) |>
  select(sensitivity,task,network,variable,term,df,p,p_BH,no_identificable)
stopifnot(all(cero_gl$no_identificable))
# Reenlazar el resumen de extracciones Windows con el principal Windows, sin alterar22/22_v2.
principal <- w21 |> select(task,network,variable,term,estimate_main=estimate,p_BH_main=p_BH,no_identificable_main=no_identificable)
resumen_windows <- w22 |> filter(term!="Intercept_HC_core") |> group_by(task,network,variable,term) |>
  summarise(draws_estimables=sum(!no_identificable),
    mean_estimate=if(any(!no_identificable))mean(estimate[!no_identificable])else NA_real_,
    min_estimate=if(any(!no_identificable))min(estimate[!no_identificable])else NA_real_,
    max_estimate=if(any(!no_identificable))max(estimate[!no_identificable])else NA_real_,
    draws_p05=sum(!no_identificable&p<.05),.groups="drop") |>
  left_join(principal,by=c("task","network","variable","term")) |>
  mutate(same_sign_mean=sign(mean_estimate)==sign(estimate_main),
    main_sig_BH=!no_identificable_main&!is.na(p_BH_main)&p_BH_main<.05)
resumen_anterior <- leer(file.path(paste0(dirs[["22"]],"_v2"),"resumen.csv"))
stopifnot(identical(resumen_windows$same_sign_mean,resumen_anterior$same_sign_mean),
  identical(resumen_windows$main_sig_BH,resumen_anterior$main_sig_BH))
res22 <- resumen_windows |> filter(variable!="AoA") |> group_by(term) |>
  summarise(models=n(),estimable=sum(draws_estimables>0),same_sign=sum(same_sign_mean,na.rm=TRUE),
    majority_p05=sum(draws_p05>=6),.groups="drop")
d027 <- w21 |> filter(!no_identificable,variable!="AoA",term%in%c("OFF","ON","OFF_minus_HC_periphery","ON_minus_HC_periphery")) |>
  summarise(n=n(),nominal_p05=sum(p<.05),BH_p05=sum(p_BH<.05))
motor <- ws |> filter(sensitivity=="contrastes_simples",task=="actions",network=="sequence",variable=="motor",term=="ON_minus_HC_core")
stopifnot(d027$n==88L,d027$nominal_p05==7L,d027$BH_p05==0L)
test_grupo <- w22 |> filter(!no_identificable,variable!="AoA",term%in%c("OFF","ON","OFF_x_Periphery","ON_x_Periphery",
  "OFF_minus_HC_periphery","ON_minus_HC_periphery","ON_minus_OFF_core","ON_minus_OFF_periphery")) |>
  summarise(n=n(),nominal_p05=sum(p<.05))
out <- nuevo_directorio("44_verificar_21_22_Windows")
for(nm in c("comparacion","por_columna","entradas","decisiones","sens","bh","diag","cero_gl","resumen_windows","res22","d027","motor","test_grupo"))
  guardar_csv(get(nm),file.path(out,paste0(nm,".csv")))
usados <- unique(c(usados,"analisis/R/21_tabla2_niveles_R2.R","analisis/R/22_tabla2_niveles_igualadas.R",
  "analisis/R/00_niveles_tabla2.R","analisis/R/44_verificar_21_22_Windows.R"))
registrar_entorno(out,usados)
max_por <- function(block,col) max(por_columna$max_abs_diff[por_columna$block==block & por_columna$column==col])
writeLines(c("# Repetición Windows21/22 frente a Linux — 23-10-2026",
  "Scripts21/22 y entradas sin modificar, comprobados con las20 huellas de las ejecuciones Linux.",
  "Windows R4.5.2/nlme3.1-168; Linux R4.3.3/nlme3.1-164. Se conservan las salidas originales; nuevas21_v2 y22_v2.",
  "Comparados13 CSV: misma estructura, textos, datos ausentes y decisiones. No identidad numérica estricta a1e-10.",
  sprintf("Máxima diferencia en estimate:21 %.3g;22 %.3g. Máximo global (DE del efecto de participante):%.3g.",
    max_por("21","estimate"),max_por("22","estimate"),max(comparacion$max_abs_diff)),
  "Tabla21 en formato R2: idéntica, incluidos coeficientes/IC redondeados y asteriscos. Cero cambios de signo/decisión BH en principal/sensibilidades; cero cambios p<0,05 en22.",
  "Principal:Periphery17/22;HC0/88 tras BH (7/88 nominales);ON−OFF0/22 por partición. OFF×Periphery concreción acciones secuencia y sesión familiaridad F co-ocurrencia se mantienen.",
  "32/32 modelos principales convergen;Shapiro rechaza normalidad en21/32. Los avisos NaN corresponden a8 términos df=0 de AoA F secuencia con solo nodos estables:ya no identificable y fuera de BH.",
  "22:idénticos tokens retenidos, N, soporte, estimabilidad y uso del aleatorio. Gradiente:signo21/22 y mayoría p nominal<0,05 en13/22. Las frecuencias por extracción siguen siendo descriptivas.",
  "22 sin modificar lee el21 original para su referencia del principal. resumen_windows.csv reenlaza los ajustes22 Windows con21 Windows, sin cambios de signos ni decisiones;los otros CSV22_v2 se conservan.",
  "Reproducibilidad numérica de la inferencia/tabla, no validación de supuestos ni equivalencia entre grupos. No se cambian modelos, familias o reglas D-024/D-027."),
  file.path(out,"RESULTADOS.md"),useBytes=TRUE)
print(comparacion,n=20);print(decisiones);print(bh,n=20);print(diag);print(d027);print(test_grupo)
cat("Salida:",out,"\n")
