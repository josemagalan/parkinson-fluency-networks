# Propósito: resumen gráfico/tabular del complemento TabPFN v2, sin nuevos tests.
# Entradas: output/41_TabPFN_v2 (métricas de los tres algoritmos y diferencias).
# Salidas: output/43_resumen_TabPFN[_vN]: figura PNG/PDF/SVG y tabla legible.
# E3: informe completo de paneles/objetivos, conservando la incertidumbre por persona.
# Paquetes: 00_utilidades.R y ggplot2; inferencia ya terminada en 41.
# --- English ---
# Purpose (optional TabPFN add-on): summary table and figure of balanced accuracy for ridge, random forest
#   and TabPFN v2; no new tests (E3).
# Inputs: 41_TabPFN_v2 outputs. Outputs: analisis/output/43_resumen_TabPFN[_vN]/ (CSV, PNG/PDF/SVG).
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
source("analisis/R/00_utilidades.R",encoding="UTF-8")
library(ggplot2)
set.seed(20261022)
base <- "analisis/output/41_TabPFN_v2"
inputs <- file.path(base,c("summary_all_algorithms.csv","paired_differences.csv","repeated_summary.csv","verification.csv"))
stopifnot(all(file.exists(inputs)))
s <- read_csv(inputs[1],show_col_types=FALSE) |> filter(metric=="balanced_accuracy")
v <- read_csv(inputs[4],show_col_types=FALSE)
stopifnot(nrow(s)==18L,all(v$passed))
out <- nuevo_directorio("43_resumen_TabPFN")
labels <- c(ridge_multinomial="Ridge multinomial",random_forest="Random forest",tabpfn_v2="TabPFN v2")
panel_labels <- c(counts="3 word counts",distributions="107 variables",distributions_AoA="131 incl. AoA (exploratory)")
s <- s |> mutate(algorithm_label=unname(labels[algorithm]),
  algorithm_label=factor(algorithm_label,levels=rev(unname(labels))),
  panel_label=factor(unname(panel_labels[panel]),levels=unname(panel_labels)),
  endpoint_label=factor(if_else(endpoint=="three_class","HC / OFF / ON","ON / OFF"),
                        levels=c("HC / OFF / ON","ON / OFF")),
  baseline=if_else(endpoint=="three_class",1/3,.5))
tabla <- s |> mutate(BA_IC95=sprintf("%.1f%% [%.1f, %.1f]",100*estimate,100*lo,100*hi)) |>
  select(endpoint,panel,algorithm,BA_IC95) |> pivot_wider(names_from=algorithm,values_from=BA_IC95)
guardar_csv(tabla,file.path(out,"tabla_ML.csv"))
guardar_csv(s |> select(endpoint,panel,algorithm,estimate,lo,hi),file.path(out,"balanced_accuracy.csv"))
p <- ggplot(s,aes(x=estimate,y=algorithm_label,colour=algorithm_label))+
  geom_vline(data=distinct(s,endpoint_label,panel_label,baseline),aes(xintercept=baseline),
             linetype="dotted",colour="grey45",inherit.aes=FALSE)+
  geom_segment(aes(x=lo,xend=hi,yend=algorithm_label),linewidth=.8)+
  geom_point(size=2.4)+
  facet_grid(endpoint_label~panel_label)+
  scale_x_continuous(limits=c(0,1),breaks=c(0,.25,.5,.75,1),labels=function(z)paste0(100*z,"%"))+
  scale_colour_manual(values=c("Ridge multinomial"="#1B6CA8","Random forest"="#C65B27","TabPFN v2"="#3A7D44"))+
  labs(x="Balanced accuracy (95% conditional bootstrap interval)",y=NULL,
       caption="LOPO by participant; n=45 people (three classes), n=21 (ON/OFF). Dotted: constant-class BA.\n2,000 participant bootstraps of fixed OOF predictions; no external validation. AoA provisional.")+
  theme_bw(base_size=10)+
  theme(legend.position="none",strip.background=element_rect(fill="grey93"),
        panel.grid.minor=element_blank(),panel.grid.major.y=element_blank(),
        plot.caption=element_text(hjust=0,size=8),axis.text.y=element_text(size=9))
ggsave(file.path(out,"comparacion_ML.png"),p,width=10,height=5,dpi=300,bg="white")
ggsave(file.path(out,"comparacion_ML.pdf"),p,width=10,height=5)
grDevices::svg(file.path(out,"comparacion_ML.svg"),width=10,height=5); print(p); invisible(dev.off())
registrar_entorno(out,c(inputs,"analisis/R/43_resumen_TabPFN.R","analisis/R/00_utilidades.R"))
writeLines(c("# Comparación completa de ML (E3)",
  "Figura descriptiva: BA e IC95% por personas de41, mismos datos/folds de15.",
  "Los intervalos están condicionados a las predicciones OOF ya ajustadas; no son simultáneos.",
  "El solapamiento/no solapamiento de intervalos no es un test entre algoritmos.",
  "Las diferencias emparejadas están en41/paired_differences.csv; variación de20 CV en41/repeated_summary.csv.",
  "No seleccionar el panel/algoritmo por el resultado. AoA provisional. No se estima efecto clínico/causal.",
  "Figura y tabla no sustituyen ni modifican Figura4 o Tabla2 del manuscrito."),
  file.path(out,"LEEME.md"),useBytes=TRUE)
print(tabla,n=6);cat("Salida:",out,"\n")
