#!/usr/bin/env bash
# Sensibilidad de la réplica SBM (script 17) a la configuración del muestreador.
# Ejecutar desde revision_R3/. Configuración principal (README de la librería):
# mcmc_init=random, n_mcmc=10*N. Aquí: previous/10N, random/50N, previous/50N
# en tres redes: denso (Animales co-ocurrencia OFF), secuencia típica
# (Acciones secuencia control) y la más dispersa (F secuencia OFF).
# Uso: bash analisis/python/17_explorar_configuracion.sh [REDES] [PARTICIONES]
# --- English ---
# Sampler-sensitivity check for script 17 on three R2 networks (initialisation and chain length).
# Writes analisis/output/17_replica_SBM_historico/exploracion_configuracion/.
set -euo pipefail
REDES="${1:-../parkinson/networks}"
PART="${2:-../parkinson/hub and spoke model}"
BASE="analisis/output/17_replica_SBM_historico/exploracion_configuracion"
for red in graphAnimalesoff Acciones_directedcontrol Fpalabra_directedoff; do
  for cfg in "previous 10" "random 50" "previous 50"; do
    set -- $cfg
    python -I analisis/python/17_replica_SBM_historico.py --redes "$REDES" --particiones "$PART" \
      --salida "$BASE/${red}_${1}_${2}N" --solo "$red" --cadenas 4 --capas "" --n-mdl 20000 \
      --mcmc-init "$1" --mcmc-por-nodo "$2"
  done
done
