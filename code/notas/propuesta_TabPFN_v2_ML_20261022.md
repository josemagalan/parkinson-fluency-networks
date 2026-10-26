# TabPFN v2 como comprobación adicional de ML — 22-10-2026

Propuesta de Josema: completar la reevaluación ML con TabPFN v2 para valorar mejor la descripción de rendimiento modesto. Valoración de Codex: pertinente como sensibilidad E3 del objetivo predictivo existente, compatible con D-025; no convierte ML en el foco central. No se ha instalado ni ejecutado TabPFN en esta consulta.

## Fundamento y precisión de la afirmación

La publicación original presenta TabPFN como modelo preentrenado para tablas pequeñas/medianas, con buenos resultados sin una búsqueda de hiperparámetros por dataset. Nuestros66 registros/45 personas y107 predictores principales están dentro de las dimensiones publicadas. Fuente: Hollmann et al., *Nature*637,319–326 (2025), [artículo](https://www.nature.com/articles/s41586-024-08328-6).

V2 ya no es la versión más reciente según la [documentación de modelos](https://docs.priorlabs.ai/models). Elegirla por su publicación y reproducibilidad es defendible; evitar llamarla universalmente «el mejor modelo actual». El [repositorio oficial](https://github.com/PriorLabs/TabPFN) permite inferencia local y selección explícita de versión; su constructor por defecto actual no selecciona necesariamentev2. Deben fijarse versión de paquete, checkpoint y configuración, aunque no haya tuning ni fine-tuning en esta muestra.

Aquí se comprueba rendimiento predictivo bajo una representación de datos, no un tamaño de efecto clínico/dopaminérgico. Añadir una familia potente reduce la dependencia de la conclusión respecto a los dos algoritmos anteriores, pero no aumenta las45 personas ni elimina su incertidumbre. Un resultado modesto no demuestra ausencia de información; una mejora obligaría a revisar la frase.

## Protocolo propuesto antes de ver los resultados de TabPFN

1. **Objetivo principal:** HC/OFF/ON,66 evaluaciones de45 personas, los mismos paneles que15: distribuciones107 variables y referencia de tres conteos. **Secundario:** ON/OFF,42 evaluaciones/21 personas. AoA añade24 variables solo como sensibilidad provisional separada.
2. **Mismos folds:** LOPO por persona, excluyendo juntas ambas evaluaciones PD. Mantener fuera ID, demografía, orden y sesión como predictores, como15. Construir exactamente las mismas variables y verificar mapa de personas; no escoger subconjuntos según rendimiento.
3. **Configuración fija:** TabPFNClassifier v2 explícito, semilla/checkpoint/configuración registrados; sin búsqueda de hiperparámetros ni fine-tuning. Usar el mismo preprocesamiento comparativo que15, aprendido dentro de cada train, y documentar la configuración interna del modelo. Un manejo nativo de NA alternativo, si se estudia, debe quedar como sensibilidad definida, no elección según el resultado.
4. **Comparación emparejada:** BA como métrica principal, con matrices y probabilidades OOF. Comparar con regresión penalizada/RF y conteos en los mismos folds. Reportar los resultados del nuevo modelo aunque contradigan «modest».
5. **Incertidumbre:** bootstrap de personas sobre las predicciones OOF y diferencias de BA emparejadas con los otros algoritmos; etiquetar los intervalos como condicionales a modelos ya ajustados. Reutilizar las particiones agrupadas repetidas de15 para describir variación de entrenamiento, sin llamar IC al rango de particiones. No afirmar significación porque los IC se solapen o no.
6. **Nulos:** no inventar una permutación por filas para tres clases. Si se amplía el nulo ON/OFF de15, conservar periodo, secuencias11/10 y reentrenamiento completo, y actualizar explícitamente la familia de algoritmos/contrastes para multiplicidad; no reutilizar BH2 como si el tercer algoritmo ya estuviera incluido.
7. **Ejecución local:** biblioteca Python en entorno aislado de la revisión, pesos descargados/identificados y datos locales; no usar una API alojada para enviar registros de participantes. Comprobar dependencias/hardware antes de fijar tiempos. Preservar15 y generar script/output nuevo, con sesión/configuración y sincronización del repositorio al implementarlo.

La propuesta se escribe después de conocer15, pero antes de observar resultados TabPFN. Es una ampliación/sensibilidad posterior, no preespecificación original del estudio. No seleccionar una versión más nueva porque dé el resultado deseado.

## Adopción y concreción antes de resultados — 22-10-2026

Josema autoriza la ejecución («vale. Te puedes encargar de ese análisis?»). Se conserva arriba la propuesta inicial. Se ejecutarán los seis LOPO (dos objetivos × tres paneles), los 2.000 bootstrap emparejados y las 20 CV agrupadas de cinco folds en conteos/distribuciones, reutilizando los planes guardados por15.

Configuración fijada tras inspeccionar el código distribuido de `tabpfn==2.2.1`: checkpoint original `Prior-Labs/TabPFN-v2-clf/tabpfn-v2-classifier.ckpt` (no la variante afinada que2.2.1 elige por defecto), ocho miembros del ensemble, temperatura0,9, sin equilibrado de probabilidades, promedio después de softmax, precisión float32, `fit_preprocessors`, `memory_saving_mode=auto`, un proceso de preprocesamiento, semilla7321+fold en LOPO y5000+10×repetición+fold en CV, como15. GPU local RTX5070Ti; Python3.12.14 aislado y PyTorch2.8.0/CUDA12.8. Sin búsqueda, ajuste de pesos, elección por rendimiento ni predictores nuevos. Se registrarán el hash de pesos y la configuración interna efectiva. Las dependencias exactas se congelarán tras instalación, antes de ajustar datos del estudio. Telemetría desactivada.

Precisión documental tras descargar la ficha de los autores: el nombre por defecto de2.2.1 contiene `finetuned`, pero la ficha actual explica que son alias de checkpoints base; ese nombre no demuestra un ajuste en datos reales. La elección explícita aquí sigue siendo `tabpfn-v2-classifier.ckpt`, sin explorar otros checkpoints. Núcleo numérico fijado antes de ajustar: numpy2.2.6, pandas2.2.3, scipy1.15.3 y scikit-learn1.6.1 (compatible con2.2.1).

El nulo ON/OFF de199 reasignaciones de15 **no se amplía en esta sensibilidad**: el alcance será comparación descriptiva emparejada con incertidumbre, sin p de TabPFN ni afirmación de significación/ausencia de efecto a partir de él. Por tanto, tampoco se presentan los p BH2 anteriores como una familia que incluya tres algoritmos. Esta decisión precede a cualquier resultado TabPFN y evita introducir una prueba nueva que no es necesaria para la pregunta de robustez entre clasificadores. Se exportan las secuencias de15 solo para trazabilidad; no se ejecutarán ni se elegirán permutaciones según los resultados.

Comprobaciones previas: equivalencia numérica R/Python del preprocesamiento en seis referencias independientes (ambos objetivos y tres paneles), identidad de folds con15 y separación completa de personas entre train/test. Outputs nuevos40/41;15 se conserva.

## Interpretación editorial prevista

- Si el rendimiento sigue siendo modesto con incertidumbre amplia: se refuerza su robustez frente a la familia de clasificador, con alcance a estas variables, muestra y validación.
- Si mejora: actualizar la descripción y valorar cuánto aporta frente a solo conteos; conservar el resultado y sus límites.
- En ambos casos: no inferir efecto causal de medicación, diagnóstico clínico ni equivalencia de grupos a partir de clasificación. Ubicación preferente como ampliación breve/suplementaria de E3, manteniendo el enfoque de redes.
