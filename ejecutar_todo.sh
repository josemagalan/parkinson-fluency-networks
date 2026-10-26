#!/usr/bin/env bash
# Ejecuta la cadena completa de la revisión R3 en orden (Linux/macOS).
# Uso, desde la raíz del repositorio:  bash ejecutar_todo.sh [--sin-replica] [--desde NN]
#   --sin-replica  omite la réplica histórica del SBM (17 y su exploración, ~1,5 h); 18–22 no dependen de ella
#   --desde NN     empieza en el paso NN (p. ej., 18) si los anteriores ya están en analisis/output/
# Orden: 01–16 (R) · 17 réplica SBM R2 (Python) · 18 redes para el SBM (R) · 19 SBM R3 (Python, ~2 h) · 20 brechas (R) · 21 Tabla 2 fiel a la R2 (R) · 22 E1 niveles igualados (R)
# Requisitos: Rscript y python en el PATH (ver ENTORNO.md) y datos según DATOS.md.
# 45–48 y los opcionales (30–32, 40–44) se ejecutan a mano después (ver README.md).
# EN: runs the main R3 pipeline (01–22) in order on Linux/macOS; then run 46, 47, 48 and 45 by hand (see README.md).
set -euo pipefail
SIN_REPLICA=0; DESDE=01
while [ $# -gt 0 ]; do
  case "$1" in --sin-replica|--sin-sbm) SIN_REPLICA=1 ;; --desde) DESDE="$2"; shift ;; *) echo "opción desconocida: $1"; exit 2 ;; esac
  shift
done
# Los scripts R leen ficheros UTF-8 con source(..., encoding="UTF-8"): hace falta una locale UTF-8.
# Orden alfabético C, como en la sesión de referencia (Windows, LC_COLLATE=C): con colación ICU
# cambian sort()/merge() y, con ellos, la auditoría 09 de Sequences.R (ver ENTORNO.md).
export LANG="${LANG:-C.UTF-8}"; case "$LANG" in *UTF-8*|*utf8*) ;; *) export LANG=C.UTF-8 ;; esac
unset LC_ALL; export LC_CTYPE="$LANG" LC_COLLATE=C
cd "$(dirname "$0")/code"
mkdir -p analisis/output logs
PASOS="01:R:01_datos_tidy 02:R:02_reproducir_R2 03:R:03_orden_sesiones 04:R:04_produccion_sesion
05:R:05_redes_igualadas_preparacion 06:R:06_psicolinguistica_sesion 07:R:07_redes_HC_bootstrap
08:R:08_auditoria_ML 09:R:09_auditoria_secuencia 10:R:10_perfiles_posicion_corregidos
11:R:11_HC_primera_sesion_influencia 12:R:12_permutacion_secuencia_sensibilidad
13:R:13_auditoria_metricas_historicas 14:R:14_preparar_entrada_SBM 15:R:15_ML_agrupado
16:R:16_verificacion_integrada 17:PY:17_replica_SBM_historico 18:R:18_exportar_redes_SBM
19:PY:19_SBM_redes_limpias 20:R:20_tabla2_core_periferia 21:R:21_tabla2_niveles_R2
22:R:22_tabla2_niveles_igualadas"
for paso in $PASOS; do
  n=${paso%%:*}; resto=${paso#*:}; lang=${resto%%:*}; nombre=${resto#*:}
  [ "$n" \< "$DESDE" ] && continue
  [ "$n" = "17" ] && [ "$SIN_REPLICA" -eq 1 ] && continue
  echo "== $(date +%H:%M:%S) $nombre"
  if [ "$lang" = "R" ]; then
    Rscript --vanilla "analisis/R/$nombre.R" > "logs/$nombre.log" 2>&1 || { echo "FALLO en $nombre (ver logs/)"; exit 1; }
  else
    python -I "analisis/python/$nombre.py" > "logs/$nombre.log" 2>&1 || { echo "FALLO en $nombre (ver logs/)"; exit 1; }
    if [ "$n" = "17" ]; then bash analisis/python/17_explorar_configuracion.sh > logs/17_explorar_configuracion.log 2>&1; fi
  fi
done
echo "== Comparación con la referencia"
python -I analisis/python/98_comparar_outputs.py --informe logs/comparacion_referencia.csv || \
  echo "Hay diferencias con la referencia: revisar logs/comparacion_referencia.csv (ver ENTORNO.md)"
