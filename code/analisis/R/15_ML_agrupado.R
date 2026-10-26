# Propósito: reevaluar la afirmación de clasificación modesta con validación por persona.
# Entradas: tokens con sesión de 03, diccionarios ya anotados; sin particiones SBM.
# Salidas: output/15_ML_agrupado[_vN], predicciones externas, matrices, métricas e incertidumbre.
# E2/E3/E4: LOPO; imputación/escala solo train; modelos fijos, sin buscar el mejor.
# Paquetes: utilidades; nnet 7.3-20; randomForest 4.7-1.2 (biblioteca local de revisión).
# --- English ---
# Purpose: re-evaluate classification with leave-one-person-out CV (fixed ridge multinomial and random
#   forest), person-bootstrap intervals, repeated grouped CV and permutation nulls (E2, E3, E4).
# Inputs: 03 tokens.rds. Outputs: analisis/output/15_ML_agrupado[_vN]/. Packages: nnet, randomForest.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R", encoding = "UTF-8")
.libPaths(c(normalizePath("analisis/R/library"), .libPaths()))
stopifnot(requireNamespace("nnet", quietly = TRUE), requireNamespace("randomForest", quietly = TRUE))
set.seed(20261022)
entrada <- "analisis/output/03_orden_sesiones/tokens.rds"
d <- readRDS(entrada)
# Características comparables al esquema de distribuciones del manuscrito,
# reparando uniones y no finitos. AoA solo sensibilidad; motor solo acciones.
resumir <- function(x) {
  v <- x[is.finite(x)]; n <- length(v)
  if (!n) return(c(mean=NA, sd=NA, min=NA, max=NA, IQR=NA, skew=NA, kurtosis=NA, coverage=0))
  s <- if (n > 1) sd(v) else NA_real_; m2 <- mean((v - mean(v))^2)
  c(mean=mean(v), sd=s, min=min(v), max=max(v), IQR=IQR(v),
    skew=if (n>2 && s>0) mean((v-mean(v))^3)/s^3 else NA_real_,
    kurtosis=if (n>3 && m2>0) mean((v-mean(v))^4)/m2^2 else NA_real_, coverage=n/length(x))
}
long <- d |> pivot_longer(all_of(variables_psico), names_to="variable", values_to="score") |>
  filter(variable != "motor" | task == "actions") |>
  group_by(participant, condition, task, variable) |>
  summarise(stats=list(resumir(score)), .groups="drop") |>
  unnest_wider(stats) |> pivot_longer(mean:coverage, names_to="stat", values_to="value")
features <- long |> pivot_wider(names_from=c(task, variable, stat), values_from=value)
counts <- d |> count(participant, condition, task, name="n") |>
  pivot_wider(names_from=task, values_from=n, names_prefix="count_")
meta <- d |> distinct(participant, condition, session, order)
features <- meta |> left_join(counts, by=c("participant","condition"), relationship="one-to-one") |>
  left_join(features, by=c("participant","condition"), relationship="one-to-one") |> arrange(participant, condition)
stopifnot(nrow(features)==66L, n_distinct(features$participant)==45L)
count_cols <- grep("^count_", names(features), value=TRUE)
all_cols <- setdiff(names(features), c("participant","condition","session","order"))
sets <- list(counts=count_cols, distributions=setdiff(all_cols, grep("_AoA_",all_cols,value=TRUE)),
  distributions_AoA=all_cols)
algorithms <- c("ridge_multinomial", "random_forest")
# No selección global de variables ni imputación global. Se eliminan variables
# sin variación/finitos solo en train; el test usa la mediana y escala del train.
preparar <- function(X, train, test) {
  a <- X[train,,drop=FALSE]; b <- X[test,,drop=FALSE]
  keep <- vapply(seq_len(ncol(a)), function(j) length(unique(a[is.finite(a[,j]),j]))>1L, logical(1))
  a <- a[,keep,drop=FALSE]; b <- b[,keep,drop=FALSE]
  for (j in seq_len(ncol(a))) {
    med <- median(a[is.finite(a[,j]),j]); a[!is.finite(a[,j]),j] <- med; b[!is.finite(b[,j]),j] <- med
  }
  mu <- colMeans(a); s <- apply(a,2,sd); s[!is.finite(s)|s==0] <- 1
  list(a=sweep(sweep(a,2,mu),2,s,"/"), b=sweep(sweep(b,2,mu),2,s,"/"), n_features=ncol(a))
}
fit_predict <- function(x,y,z,algorithm,seed) {
  set.seed(seed)
  if (algorithm=="ridge_multinomial") {
    train <- data.frame(y=y,x,check.names=FALSE); test <- data.frame(z,check.names=FALSE)
    fit <- nnet::multinom(y~., data=train, decay=1, maxit=1000, trace=FALSE, MaxNWts=10000)
    stopifnot(fit$convergence==0L)
    prob <- predict(fit,newdata=test,type="probs")
    if (is.null(dim(prob))) {
      if (nlevels(y)==2L) prob <- cbind(1-prob,prob) else prob <- matrix(prob,nrow=1L)
    }
    colnames(prob) <- levels(y)
  } else {
    fit <- randomForest::randomForest(x=x,y=y,ntree=200L,mtry=max(1L,floor(sqrt(ncol(x)))),nodesize=3L)
    prob <- predict(fit,z,type="prob")
  }
  prob <- as.matrix(prob); colnames(prob) <- levels(y)
  list(pred=levels(y)[max.col(prob,ties.method="first")],prob=prob)
}
score <- function(y,pred,w=rep(1,length(y))) {
  recall <- vapply(sort(unique(y)),function(g) sum(w[y==g]*(pred[y==g]==g))/sum(w[y==g]),numeric(1))
  c(accuracy=sum(w*(y==pred))/sum(w), balanced_accuracy=mean(recall), setNames(recall,paste0("recall_",names(recall))))
}
cv <- function(dat, cols, target, algorithm, fold, seed=7321L) {
  y <- factor(target); X <- as.matrix(dat[cols]); results <- list()
  for (f in sort(unique(fold))) {
    test <- which(fold==f); train <- which(fold!=f)
    stopifnot(!length(intersect(dat$participant[train],dat$participant[test])), all(table(y[train])>0))
    xy <- preparar(X,train,test); pr <- fit_predict(xy$a,y[train],xy$b,algorithm,seed+f)
    prob <- matrix(0,nrow=length(test),ncol=3,dimnames=list(NULL,c("HC","OFF","ON")))
    prob[,levels(y)] <- pr$prob
    results[[length(results)+1L]] <- tibble(row=test,participant=dat$participant[test],truth=as.character(y[test]),
      predicted=pr$pred, fold=f,n_train=length(train),n_features=xy$n_features,
      p_HC=prob[,"HC"],p_OFF=prob[,"OFF"],p_ON=prob[,"ON"])
  }
  bind_rows(results) |> arrange(row)
}
predictions <- list(); counter <- 0L
for (endpoint in c("three_class","ON_OFF")) {
  dat <- if(endpoint=="ON_OFF") filter(features,condition!="HC") else features
  fold <- match(dat$participant,unique(dat$participant))
  for (panel in names(sets)) for (algorithm in algorithms) {
    p <- cv(dat,sets[[panel]],dat$condition,algorithm,fold)
    counter <- counter+1L; predictions[[counter]] <- mutate(p,endpoint=endpoint,panel=panel,algorithm=algorithm)
    cat("LOPO",endpoint,panel,algorithm, paste(round(score(p$truth,p$predicted),3),collapse=" "),"\n"); flush.console()
  }
}
predictions <- bind_rows(predictions)
# Referencia sin características: predice la clase mayoritaria de train.
# BA = 1/C exactamente para predicción constante; accuracy de HC = 24/66.
baseline <- tibble(endpoint=c("three_class","ON_OFF"),balanced_accuracy=c(1/3,.5),
  constant_accuracy=c(24/66,.5),description=c("Predicción constante HC","Predicción constante OFF"))
# Bootstrap de las predicciones OOF por persona: condicionado a modelos ya
# ajustados. No equivale a bootstrap anidado de TODO el entrenamiento.
B <- 2000L
people <- unique(features$participant); PD <- people[grepl("PD_",people)]; HC <- people[grepl("HC_",people)]
set.seed(9312)
boot_plans <- lapply(seq_len(B),function(i)c(sample(PD,length(PD),replace=TRUE),sample(HC,length(HC),replace=TRUE)))
summary <- list(); boot <- list()
groups <- predictions |> group_split(endpoint,panel,algorithm)
for (j in seq_along(groups)) {
  p <- groups[[j]]; point <- score(p$truth,p$predicted)
  vals <- t(vapply(boot_plans,function(ids) score(p$truth,p$predicted,as.numeric(table(factor(ids,levels=people)))[match(p$participant,people)]),point))
  summary[[j]] <- tibble(endpoint=p$endpoint[1],panel=p$panel[1],algorithm=p$algorithm[1],metric=names(point),estimate=unname(point),
    lo=apply(vals,2,quantile,.025),hi=apply(vals,2,quantile,.975))
  boot[[j]] <- vals
}
summary <- bind_rows(summary)
# Variación de partición/entrenamiento: 20 CV de cinco folds agrupados,
# cada fold recibe HC y PD; mismo plan para paneles/algoritmos.
set.seed(3151); plans <- lapply(seq_len(20),function(i) {
  m <- setNames(integer(length(people)),people)
  m[PD] <- sample(rep(1:5,length.out=length(PD))); m[HC] <- sample(rep(1:5,length.out=length(HC))); m
})
repeated <- list(); k <- 0L
for (r in seq_along(plans)) for(endpoint in c("three_class","ON_OFF")) {
  dat <- if(endpoint=="ON_OFF") filter(features,condition!="HC") else features
  for (panel in c("counts","distributions")) for (algorithm in algorithms) {
    p <- cv(dat,sets[[panel]],dat$condition,algorithm,unname(plans[[r]][dat$participant]),seed=5000+r*10)
    k <- k+1L; repeated[[k]] <- tibble(repeat_id=r,endpoint=endpoint,panel=panel,algorithm=algorithm,
      balanced_accuracy=unname(score(p$truth,p$predicted)["balanced_accuracy"]))
  }
  if(endpoint=="ON_OFF") {cat("CV agrupada repetida",r,"/20\n");flush.console()}
}
repeated <- bind_rows(repeated)
# Nulo exploratorio ON/OFF que conserva periodo: reasignar las 11 secuencias
# OFF_ON entre los 21 pacientes. Reentrenar TODO el LOPO para cada permutación.
# No es un test aleatorizado documentado; requiere intercambiabilidad de secuencia.
# No se fabrica un nulo global HC/PD intercambiando filas de tamaños de cluster distintos.
M <- 199L; dat <- filter(features,condition!="HC"); ids <- unique(dat$participant)
fold <- match(dat$participant,ids); set.seed(82342)
perms <- replicate(M,sample(ids,11L),simplify=FALSE)
null <- list(); k <- 0L
for (b in seq_len(M)) {
  off_first <- dat$participant %in% perms[[b]]
  label <- ifelse((dat$session==1L & off_first)|(dat$session==2L & !off_first),"OFF","ON")
  for(algorithm in algorithms) {
    p <- cv(dat,sets$distributions,label,algorithm,fold,seed=7321L)
    k <- k+1L; null[[k]] <- tibble(permutation=b,algorithm=algorithm,balanced_accuracy=unname(score(p$truth,p$predicted)["balanced_accuracy"]))
  }
  if(b%%20L==0L) {cat("Permutación ML",b,"/",M,"\n");flush.console()}
}
null <- bind_rows(null)
permutation_tests <- bind_rows(lapply(algorithms,function(a) {
  observed <- summary |> filter(endpoint=="ON_OFF",panel=="distributions",algorithm==a,metric=="balanced_accuracy") |> pull(estimate)
  vals <- null$balanced_accuracy[null$algorithm==a]
  tibble(algorithm=a,observed=observed,B=M,p_one_sided=(1+sum(vals>=observed-1e-12))/(M+1),null_mean=mean(vals))
})) |> mutate(p_BH2=p.adjust(p_one_sided,"BH"))
confusion <- predictions |> count(endpoint,panel,algorithm,truth,predicted,name="n")
out <- nuevo_directorio("15_ML_agrupado")
for(nm in c("features","predictions","summary","baseline","repeated","null","permutation_tests","confusion")) guardar_csv(get(nm),file.path(out,paste0(nm,".csv")))
saveRDS(list(sets=sets,bootstrap_plans=boot_plans,bootstrap_metrics=boot,grouped_CV_plans=plans,sequence_permutations=perms),file.path(out,"planes.rds"))
registrar_entorno(out,c(entrada,"analisis/R/15_ML_agrupado.R","analisis/R/00_utilidades.R"))
writeLines(c("# Reevaluación ML agrupada — 22-10-2026",
  "Autorización: Josema solicita rehacer ML para comprobar la frase de resultados modestos.",
  "No es réplica exacta de Weka: configuración/resultados históricos no recuperados. No se reevalúan XGBoost/superlearner/rotation forest.",
  "Objetivo tres clases HC/OFF/ON (66 evaluaciones,45 personas); secundario ON/OFF (42,21). Sin demografía, ID, orden ni sesión como predictores.",
  "Panel principal: 3 conteos + 13 distribuciones x 8 características (media,sd,min,max,IQR,asimetría,kurtosis,cobertura)=107. AoA añade24 solo sensibilidad.",
  "Referencia con solo3 conteos y referencia constante. No se selecciona panel/modelo según rendimiento.",
  "Regresión multinomial con penalización L2 decay=1; RF 200 árboles,mtry=floor(sqrt(p)),nodesize=3. Parámetros fijos: no hay tuning.",
  "LOPO: se excluyen ambas evaluaciones del paciente. Imputación por mediana, retirada de constantes y estandarización solo en train.",
  "2000 bootstrap de personas de predicciones externas, PD/HC estratificados: intervalos percentiles aproximados, condicionales a modelos ajustados, no incluyen toda variación de entrenamiento.",
  "20 CV adicionales de5 folds por persona describen sensibilidad a particiones y menor train. Su rango tampoco es IC.",
  "Permutación199 ON/OFF: conserva periodo y reasigna secuencias11/10; reentrena LOPO completo para2 algoritmos, BH2. Depende de intercambiabilidad no confirmada (P14), sin interpretación causal.",
  "No hay p de azar global de tres clases: distintos tamaños de cluster no permiten una permutación ingenua de filas/etiquetas HC/OFF/ON.",
  "Los modelos predicen condición/grupo; no prueban diagnóstico clínico ni separación del efecto de medicación y práctica. Incidencias P3/P11 y cobertura siguen pendientes.",
  "Todo se informa sin afirmar equivalencia al azar por un p no significativo; revisar métricas antes de conservar la frase modestos."),file.path(out,"LEEME.md"),useBytes=TRUE)
print(summary |> filter(metric=="balanced_accuracy"),n=20); print(permutation_tests); cat("Salida:",out,"\n")
