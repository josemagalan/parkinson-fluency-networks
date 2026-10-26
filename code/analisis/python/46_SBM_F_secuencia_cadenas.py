# -*- coding: utf-8 -*-
# --- English ---
# F-sequence ON/OFF networks: 50 SBM chains each with script 19's procedure, saving every chain's
# consensus and 5 posterior samples (D-030; feeds script 47) (E2, E5).
# Inputs: 18 and 19 outputs. Outputs: analisis/output/46_SBM_F_secuencia_cadenas/ (overwritten).
"""
46_SBM_F_secuencia_cadenas.py

Propósito
    Cuantificar la incertidumbre de las particiones propias de F secuencia ON y OFF, que el
    script 19 declaró no identificables (cadenas discordantes; MDL hub-and-spoke ≈ «sin core»).
    Se vuelve a ajustar el SBM hub-and-spoke con el MISMO procedimiento que 19 (grafo no
    dirigido y binario, sin bucles; n_gibbs=100, n_mcmc=10·N, last_n_samples=50), con más
    cadenas, y se GUARDAN las particiones de cada cadena:
      - consenso de cada cadena (últimas 50 muestras; misma definición que 19);
      - muestras posteriores individuales (últimas 50, una de cada 10 → 5 por cadena).
    Las cadenas 1–10 usan las semillas de 19 y deben reproducir sus tamaños de core.
    El script 47 estima con estas particiones los contrastes de la Tabla 2 en F secuencia.

Puntos editoriales: E5 (qué pasa con las dos diferencias de grupo de la R2 en F secuencia), E2.
    Sensibilidad añadida el 24-10-2026 (D-030), después de conocer los resultados de 19/21.

Entradas (desde revision_R3/)
    analisis/output/18_exportar_redes_SBM/{index,edges,nodes}.csv
    analisis/output/19_SBM_redes_limpias/cadenas.csv          (verificación de semillas 1–10)
    analisis/python/vendor/core_periphery_sbm/

Salidas  analisis/output/46_SBM_F_secuencia_cadenas/
    particiones.csv     graph, semilla, tipo (consenso_cadena | muestra), muestra, word_norm, hs (0 core, 1 periferia)
    cadenas.csv         por red y cadena: tamaño del core (consenso), log-posterior final, aceptación
    verificacion_19.csv tamaño del core de las semillas 1–10 frente a 19
    config.json, entorno.txt, entradas_md5.csv, LEEME.md

Uso
    python -I analisis/python/46_SBM_F_secuencia_cadenas.py [--cadenas 50] [--procesos 2]
"""
import argparse
import csv
import hashlib
import json
import os
import platform
import sys
import time
import warnings
from multiprocessing import Pool

import numpy as np
import networkx as nx
import scipy

warnings.filterwarnings("ignore")
COMMIT_LIB = "90fb2896172d0489490c3574b4761af7e7f778f9"
REDES = ["F_OFF_directed", "F_ON_directed"]
A = None
G_CACHE = {}


def args_cli():
    p = argparse.ArgumentParser()
    p.add_argument("--entrada", default="analisis/output/18_exportar_redes_SBM")
    p.add_argument("--ref19", default="analisis/output/19_SBM_redes_limpias/cadenas.csv")
    p.add_argument("--lib", default="analisis/python/vendor/core_periphery_sbm")
    p.add_argument("--salida", default="analisis/output/46_SBM_F_secuencia_cadenas")
    p.add_argument("--cadenas", type=int, default=50)
    p.add_argument("--n-gibbs", type=int, default=100)
    p.add_argument("--mcmc-por-nodo", type=int, default=10)
    p.add_argument("--ultimas", type=int, default=50)
    p.add_argument("--cada", type=int, default=10, help="muestras posteriores: una de cada k de las últimas")
    p.add_argument("--procesos", type=int, default=2)
    return p.parse_args()


def leer_csv(ruta):
    with open(ruta, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def escribir(ruta, filas):
    with open(ruta, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0].keys()), lineterminator="\n")
        w.writeheader(); w.writerows(filas)


def md5(ruta):
    with open(ruta, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


def cargar(a):
    grafos = {g: nx.Graph() for g in REDES}
    for r in leer_csv(os.path.join(a.entrada, "nodes.csv")):
        if r["graph"] in grafos:
            grafos[r["graph"]].add_node(r["word_norm"])
    for r in leer_csv(os.path.join(a.entrada, "edges.csv")):
        if r["graph"] in grafos and r["from"] != r["to"]:
            grafos[r["graph"]].add_edge(r["from"], r["to"])  # dirigida → no dirigida, como en 19 (D-021)
    return grafos


def init_hijo(a, grafos):
    global A, cp, G_CACHE
    A = a; G_CACHE = grafos
    sys.path.insert(0, os.path.abspath(a.lib))
    from core_periphery_sbm import core_periphery as cp_  # noqa
    cp = cp_


def t_cadena(job):
    g, semilla = job
    G = G_CACHE[g]
    m = cp.HubSpokeCorePeriphery(n_gibbs=A.n_gibbs, n_mcmc=A.mcmc_por_nodo * len(G), seed=semilla)
    m.infer(G)
    lab = m.get_labels(last_n_samples=A.ultimas, prob=False, return_dict=True)
    # Muestras individuales: filas de node_labels (índice interno → nodo con index2node)
    filas_m = list(range(A.n_gibbs - A.ultimas, A.n_gibbs))[A.cada - 1::A.cada]
    muestras = []
    for k, t in enumerate(filas_m, 1):
        v = m.node_labels[t]
        muestras.append({m.index2node[i]: int(v[i]) for i in range(len(v))})
    return dict(graph=g, semilla=semilla, consenso={n: int(x) for n, x in lab.items()}, muestras=muestras,
                log_post=float(m.log_posteriors[-1]), aceptacion=float(np.mean(m.accepts)))


def main():
    a = args_cli()
    os.makedirs(a.salida, exist_ok=True)
    grafos = cargar(a)
    jobs = [(g, s) for g in REDES for s in range(1, a.cadenas + 1)]
    t0 = time.time()
    with Pool(a.procesos, initializer=init_hijo, initargs=(a, grafos)) as pool:
        res = []
        for k, x in enumerate(pool.imap_unordered(t_cadena, jobs, chunksize=1), 1):
            res.append(x)
            if k % 10 == 0:
                print(f"  {k}/{len(jobs)} {time.time()-t0:.0f}s", flush=True)
    res.sort(key=lambda x: (x["graph"], x["semilla"]))
    filas_p, filas_c = [], []
    for x in res:
        nodos = sorted(grafos[x["graph"]].nodes())
        for n in nodos:
            filas_p.append(dict(graph=x["graph"], semilla=x["semilla"], tipo="consenso_cadena", muestra=0, word_norm=n, hs=x["consenso"][n]))
        for k, mu in enumerate(x["muestras"], 1):
            for n in nodos:
                filas_p.append(dict(graph=x["graph"], semilla=x["semilla"], tipo="muestra", muestra=k, word_norm=n, hs=mu[n]))
        filas_c.append(dict(graph=x["graph"], semilla=x["semilla"], nodos=len(nodos),
                            core=sum(v == 0 for v in x["consenso"].values()),
                            core_muestras="/".join(str(sum(v == 0 for v in mu.values())) for mu in x["muestras"]),
                            log_posterior_final=round(x["log_post"], 3), aceptacion=round(x["aceptacion"], 4)))
    ref = {(r["graph"], int(r["semilla"])): int(r["core"]) for r in leer_csv(a.ref19) if r["graph"] in REDES}
    filas_v = [dict(graph=c["graph"], semilla=c["semilla"], core_46=c["core"], core_19=ref[(c["graph"], c["semilla"])],
                    coincide=c["core"] == ref[(c["graph"], c["semilla"])])
               for c in filas_c if (c["graph"], c["semilla"]) in ref]
    escribir(os.path.join(a.salida, "particiones.csv"), filas_p)
    escribir(os.path.join(a.salida, "cadenas.csv"), filas_c)
    escribir(os.path.join(a.salida, "verificacion_19.csv"), filas_v)
    entradas = [os.path.join(a.entrada, f) for f in ("index.csv", "edges.csv", "nodes.csv")] + [a.ref19]
    escribir(os.path.join(a.salida, "entradas_md5.csv"),
             [dict(path=p.replace("\\", "/"), md5=md5(p)) for p in entradas] +
             [dict(path="analisis/python/46_SBM_F_secuencia_cadenas.py", md5=md5(os.path.abspath(__file__)))])
    cfg = dict(vars(a)); cfg.update(redes=REDES, commit_libreria=COMMIT_LIB, segundos_totales=round(time.time() - t0, 1))
    with open(os.path.join(a.salida, "config.json"), "w", encoding="utf-8") as f:
        json.dump(cfg, f, indent=2, ensure_ascii=False)
    with open(os.path.join(a.salida, "entorno.txt"), "w", encoding="utf-8") as f:
        f.write(f"python {sys.version}\nplataforma {platform.platform()}\nnumpy {np.__version__}\n"
                f"networkx {nx.__version__}\nscipy {scipy.__version__}\ncore_periphery_sbm {COMMIT_LIB}\n")
    ok = sum(v["coincide"] for v in filas_v)
    with open(os.path.join(a.salida, "LEEME.md"), "w", encoding="utf-8") as f:
        f.write("# 46: particiones por cadena de F secuencia ON/OFF (D-030)\n\n"
                f"Mismo procedimiento que 19; {a.cadenas} cadenas por red (semillas 1–{a.cadenas}).\n"
                f"Se guardan el consenso de cada cadena (últimas {a.ultimas} muestras) y {len(range(a.cada - 1, a.ultimas, a.cada))} muestras posteriores "
                f"individuales por cadena (una de cada {a.cada}).\n"
                f"Verificación: semillas 1–10 reproducen el tamaño de core de 19 en {ok}/{len(filas_v)} cadenas.\n"
                "Uso: entrada del script 47 (contrastes de Tabla 2 en F secuencia bajo incertidumbre de partición).\n")
    print("hecho", round(time.time() - t0), "s; verificación", ok, "/", len(filas_v))


if __name__ == "__main__":
    main()
