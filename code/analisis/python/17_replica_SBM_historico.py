# -*- coding: utf-8 -*-
# --- English ---
# Replication of the R2 Bayesian core-periphery SBM (Gallagher et al. 2021; vendored library) on the
# 18 R2 networks: 10 hub-and-spoke chains per network plus layered models for MDL comparison (E1, E2, E5).
# Inputs: ../parkinson/networks and R2 partitions. Outputs: analisis/output/17_replica_SBM_historico/ (overwritten).
"""
17_replica_SBM_historico.py

Propósito
    Replicar, sin el código de Virginia, el ajuste core-periphery (SBM bayesiano
    de Gallagher et al., 2021) de la R2 sobre las 18 redes históricas, y
    comprobar si las particiones guardadas en `parkinson/hub and spoke model/`
    se recuperan con la implementación oficial de la librería. Es una réplica de
    VALIDACIÓN: no sustituye ni sobrescribe ninguna partición histórica.

Puntos editoriales
    Prepara E1/E2/E5 (nuevas particiones sobre redes limpias/igualadas) y P5/P12.

Entradas (rutas por defecto relativas a revision_R3/)
    ../parkinson/networks/graph<Tarea>[_directed]<grupo>.graphml           (18 redes R2)
    ../parkinson/hub and spoke model/{directed,undirected}/csvs/*.graphml.csv (particiones R2)
    analisis/python/vendor/core_periphery_sbm/   (librería oficial, commit fijado)

Salidas  analisis/output/17_replica_SBM_historico/
    redes.csv                 descripción de cada red y variante de entrada
    particiones_cadenas.csv   etiqueta y P(core) por nodo, red, variante y semilla
    comparacion_cadenas.csv   cada cadena frente a la R2 (concordancia, ARI, MDL)
    resumen_hubspoke.csv      resumen por red/variante (R2 vs réplica vs entre cadenas)
    consenso.csv              partición de consenso por nodo frente a la R2
    mdl_modelos.csv           longitudes de descripción hub-and-spoke y capas (S1/S2)
    entorno.txt, entradas_md5.csv, config.json

Decisiones verificadas en el código de la librería (commit 90fb289, 15-11-2021)
    - El modelo es no dirigido y no ponderado: solo usa presencia de arista.
    - `infer()` elimina bucles; `mdl_*()` NO los elimina (la docstring lo
      presupone). Aquí se eliminan antes de calcular la MDL.
    - Con un DiGraph, las estadísticas iniciales cuentan todas las aristas, pero
      la actualización incremental usa solo sucesores (G[i]) y el número máximo
      de aristas es el no dirigido. Por eso, en redes dirigidas se prueban dos
      variantes: grafo dirigido tal cual y grafo convertido a no dirigido.
    - Parámetros por defecto del README: n_gibbs=100, n_mcmc=10*N,
      last_n_samples=50, MDL con n_samples=100000. La configuración de Virginia
      es desconocida (P5).
    - moves="neighbor" falla en este commit (firma de propose_neigh_move); solo
      se puede usar moves="random".
    - La partición R2 se lee del CSV (lo que usó la Tabla 2), no de los
      *Partition.graphml: en Acciones-secuencia-control el nodo n202
      ('partir<TAB>llegar') quedó allí con el valor centinela -2147483648.
    - Cada cadena fija np.random.seed(semilla); la MDL usa semillas derivadas.
      Resultados reproducibles con las mismas versiones de numpy/networkx.

Dependencias: Python >= 3.9, numpy, scipy, networkx (versiones en entorno.txt).

Uso (desde revision_R3/)
    python -I analisis/python/17_replica_SBM_historico.py [--cadenas 10] [--procesos 2]
    bash analisis/python/17_explorar_configuracion.sh   (sensibilidad del muestreador)
    Ejecución completa: ~1 h con 2 procesos (270 cadenas hub-and-spoke + 162 de capas).
"""
import argparse
import csv
import hashlib
import itertools
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

TAREAS = ["Acciones", "Animales", "Fpalabra"]
GRUPOS = ["control", "off", "on"]
TIPOS = {"": "coocurrencia", "_directed": "secuencia"}
COMMIT_LIB = "90fb2896172d0489490c3574b4761af7e7f778f9"


def args_cli():
    p = argparse.ArgumentParser()
    p.add_argument("--redes", default="../parkinson/networks")
    p.add_argument("--particiones", default="../parkinson/hub and spoke model")
    p.add_argument("--lib", default="analisis/python/vendor/core_periphery_sbm")
    p.add_argument("--salida", default="analisis/output/17_replica_SBM_historico")
    p.add_argument("--cadenas", type=int, default=10)
    p.add_argument("--cadenas-capas", type=int, default=3)
    p.add_argument("--capas", default="2,3,4")
    p.add_argument("--capas-variantes", default="no_dirigida",
                   help="variantes de entrada en las que se ajusta el modelo de capas")
    p.add_argument("--n-gibbs", type=int, default=100)
    p.add_argument("--mcmc-por-nodo", type=int, default=10)
    p.add_argument("--ultimas", type=int, default=50)
    p.add_argument("--n-mdl", type=int, default=100000)
    p.add_argument("--mcmc-init", default="random", choices=["random", "previous"])
    # moves="neighbor" falla en el commit fijado (firma de propose_neigh_move); no se ofrece
    p.add_argument("--moves", default="random", choices=["random"])
    p.add_argument("--procesos", type=int, default=2)
    p.add_argument("--solo", default="", help="filtrar redes por subcadena (pruebas)")
    return p.parse_args()


A = None  # argumentos globales para los procesos hijos


def init_hijo(a):
    global A, cp, mf
    A = a
    sys.path.insert(0, os.path.abspath(a.lib))
    from core_periphery_sbm import core_periphery as cp_  # noqa
    from core_periphery_sbm import model_fit as mf_  # noqa
    cp, mf = cp_, mf_


def md5(ruta):
    h = hashlib.md5()
    with open(ruta, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def catalogo(a):
    out = []
    for t in TAREAS:
        for suf, tipo in TIPOS.items():
            for g in GRUPOS:
                base = f"graph{t}{suf}{g}"
                if a.solo and a.solo not in base:
                    continue
                out.append(dict(red=base, tarea=t, tipo=tipo, grupo=g,
                                f_red=os.path.join(a.redes, base + ".graphml"),
                                f_csv=os.path.join(a.particiones, "directed" if suf else "undirected",
                                                   "csvs", base + ".graphml.csv")))
    return out


def variantes(G):
    """Variantes de entrada al SBM. Siempre sin bucles para la MDL."""
    if G.is_directed():
        return {"dirigida_tal_cual": G.copy(), "no_dirigida": nx.Graph(G)}
    return {"no_dirigida": G.copy()}


def sin_bucles(G):
    H = G.copy()
    H.remove_edges_from(list(nx.selfloop_edges(H)))
    return H


def etiquetas_hist(r):
    """Partición R2 desde el CSV histórico (node_id, hs_partition), que es el que
    usó la Tabla 2. Los nombres se toman del GraphML de la red (UTF-8). El CSV
    tiene codificación mixta y una fila entrecomillada entera (n202, celda
    'partir<TAB>llegar'), que se recupera aquí; solo se leen campos ASCII."""
    G = nx.read_graphml(r["f_red"])
    nombres = {n: d.get("name", "") for n, d in G.nodes(data=True)}
    with open(r["f_csv"], "rb") as f:
        texto = f.read().decode("latin-1")
    et = {}
    for k, fila in enumerate(csv.reader(texto.splitlines())):
        if k == 0:
            continue
        if len(fila) == 1:  # fila entrecomillada completa
            fila = fila[0].split(",")
        et[fila[1]] = int(float(fila[-1]))
    return et, nombres


def vector(G, d):
    return np.array([d[n] for n in sorted(G.nodes())], dtype=int)


def ari(x, y):
    x = np.asarray(x); y = np.asarray(y)
    cx = {v: i for i, v in enumerate(np.unique(x))}
    cy = {v: i for i, v in enumerate(np.unique(y))}
    M = np.zeros((len(cx), len(cy)))
    for a_, b_ in zip(x, y):
        M[cx[a_], cy[b_]] += 1
    c2 = lambda v: v * (v - 1) / 2
    s = c2(M).sum(); a_s = c2(M.sum(1)).sum(); b_s = c2(M.sum(0)).sum(); n = c2(len(x))
    esp = a_s * b_s / n
    mx = (a_s + b_s) / 2
    return float((s - esp) / (mx - esp)) if mx != esp else 1.0


def tarea_hubspoke(job):
    r, var, semilla = job
    G = nx.read_graphml(r["f_red"])
    Gv = variantes(G)[var]
    t0 = time.time()
    m = cp.HubSpokeCorePeriphery(n_gibbs=A.n_gibbs, n_mcmc=A.mcmc_por_nodo * len(Gv), seed=semilla,
                                 mcmc_init=A.mcmc_init, moves=A.moves)
    m.infer(Gv)
    lab = m.get_labels(last_n_samples=A.ultimas, prob=False, return_dict=True)
    prob = m.get_labels(last_n_samples=A.ultimas, prob=True, return_dict=True)
    seg = time.time() - t0
    H = sin_bucles(Gv)
    np.random.seed(10_000 + semilla)
    mdl = mf.mdl_hubspoke(H, vector(H, lab), n_samples=A.n_mdl)
    return dict(red=r["red"], variante=var, semilla=semilla, etiquetas=lab,
                p_core={n: float(p[0]) for n, p in prob.items()}, mdl=float(mdl),
                segundos=seg, aceptacion=float(np.mean(m.accepts)))


def tarea_capas(job):
    r, var, L, semilla = job
    G = nx.read_graphml(r["f_red"])
    Gv = variantes(G)[var]
    m = cp.LayeredCorePeriphery(n_layers=L, n_gibbs=A.n_gibbs, n_mcmc=A.mcmc_por_nodo * len(Gv), seed=semilla,
                                mcmc_init=A.mcmc_init, moves=A.moves)
    ok = m.infer(Gv)
    if ok is False:
        return dict(red=r["red"], variante=var, capas=L, semilla=semilla, mdl=np.nan, tamanos="")
    lab = m.get_labels(last_n_samples=A.ultimas, prob=False, return_dict=False)
    H = sin_bucles(Gv)
    np.random.seed(20_000 + 100 * L + semilla)
    mdl = mf.mdl_layered(H, lab, n_layers=L, n_samples=A.n_mdl)
    tam = np.bincount(lab, minlength=L)
    return dict(red=r["red"], variante=var, capas=L, semilla=semilla, mdl=float(mdl),
                tamanos="/".join(map(str, tam)))


def tarea_mdl_hist(job):
    r, var = job
    G = nx.read_graphml(r["f_red"])
    H = sin_bucles(variantes(G)[var])
    eh, _ = etiquetas_hist(r)
    np.random.seed(30_000)
    return dict(red=r["red"], variante=var, mdl=float(mf.mdl_hubspoke(H, vector(H, eh), n_samples=A.n_mdl)))


def escribir(ruta, filas, campos=None):
    if not filas:
        return
    campos = campos or list(filas[0].keys())
    with open(ruta, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        w.writerows(filas)


def main():
    a = args_cli()
    os.makedirs(a.salida, exist_ok=True)
    init_hijo(a)
    cat = catalogo(a)
    capas = [int(x) for x in a.capas.split(",") if x]

    # --- Descripción de entradas y comprobaciones de identidad -----------------
    redes, md5s = [], []
    for r in cat:
        G = nx.read_graphml(r["f_red"])
        eh, nh = etiquetas_hist(r)
        nombres = {n: d.get("name", "") for n, d in G.nodes(data=True)}
        assert set(eh) == set(G.nodes()), f"nodos distintos en {r['red']}"
        for var, Gv in variantes(G).items():
            redes.append(dict(red=r["red"], tarea=r["tarea"], tipo=r["tipo"], grupo=r["grupo"],
                              variante=var, clase=type(Gv).__name__, nodos=Gv.number_of_nodes(),
                              aristas=Gv.number_of_edges(), bucles=nx.number_of_selfloops(Gv),
                              aristas_sin_bucles=sin_bucles(Gv).number_of_edges(),
                              core_R2=sum(v == 0 for v in eh.values())))
        md5s.append(dict(fichero="parkinson/networks/" + os.path.basename(r["f_red"]), md5=md5(r["f_red"])))
        md5s.append(dict(fichero="parkinson/hub and spoke model/" + ("directed" if r["tipo"] == "secuencia" else "undirected")
                         + "/csvs/" + os.path.basename(r["f_csv"]), md5=md5(r["f_csv"])))
    escribir(os.path.join(a.salida, "redes.csv"), redes)
    escribir(os.path.join(a.salida, "entradas_md5.csv"), md5s)

    porred = {r["red"]: r for r in cat}
    jobs_hs = [(porred[x["red"]], x["variante"], s) for x in redes for s in range(1, a.cadenas + 1)]
    jobs_hist = [(porred[x["red"]], x["variante"]) for x in redes]
    jobs_cap = [(porred[x["red"]], x["variante"], L, s) for x in redes for L in capas
                for s in range(1, a.cadenas_capas + 1) if x["variante"] in a.capas_variantes.split(",")]
    t0 = time.time()
    with Pool(a.procesos, initializer=init_hijo, initargs=(a,)) as pool:
        res_hs = []
        for k, x in enumerate(pool.imap(tarea_hubspoke, jobs_hs, chunksize=1), 1):
            res_hs.append(x)
            print(f"  hs {k}/{len(jobs_hs)} {x['red']} {x['variante']} s{x['semilla']} {time.time()-t0:.0f}s", flush=True)
        print(f"hub-and-spoke: {len(res_hs)} cadenas, {time.time()-t0:.0f} s", flush=True)
        res_hist = pool.map(tarea_mdl_hist, jobs_hist, chunksize=1)
        res_cap = []
        for k, x in enumerate(pool.imap(tarea_capas, jobs_cap, chunksize=1), 1):
            res_cap.append(x)
            print(f"  capas {k}/{len(jobs_cap)} {x['red']} L{x['capas']} s{x['semilla']} {time.time()-t0:.0f}s", flush=True)
        print(f"capas: {len(res_cap)} cadenas, {time.time()-t0:.0f} s", flush=True)

    mdl_hist = {(x["red"], x["variante"]): x["mdl"] for x in res_hist}

    # --- Particiones por cadena, comparación con R2 -----------------------------
    filas_par, filas_cmp = [], []
    agrup = {}
    for x in res_hs:
        r = porred[x["red"]]
        eh, nh = etiquetas_hist(r)
        nodos = sorted(eh)
        v_rep = np.array([x["etiquetas"][n] for n in nodos])
        v_his = np.array([eh[n] for n in nodos])
        for n in nodos:
            filas_par.append(dict(red=x["red"], variante=x["variante"], semilla=x["semilla"], node_id=n,
                                  name=nh[n], hs_R2=eh[n], hs_replica=x["etiquetas"][n],
                                  p_core=round(x["p_core"][n], 4)))
        filas_cmp.append(dict(red=x["red"], variante=x["variante"], semilla=x["semilla"], nodos=len(nodos),
                              core_R2=int((v_his == 0).sum()), core_replica=int((v_rep == 0).sum()),
                              concordancia=round(float((v_rep == v_his).mean()), 4),
                              ARI=round(ari(v_rep, v_his), 4),
                              identica=bool((v_rep == v_his).all()),
                              mdl_replica=round(x["mdl"], 3),
                              mdl_R2=round(mdl_hist[(x["red"], x["variante"])], 3),
                              segundos=round(x["segundos"], 1), aceptacion=round(x["aceptacion"], 4)))
        agrup.setdefault((x["red"], x["variante"]), []).append((x, v_rep, v_his, nodos, nh))
    escribir(os.path.join(a.salida, "particiones_cadenas.csv"), filas_par)
    escribir(os.path.join(a.salida, "comparacion_cadenas.csv"), filas_cmp)

    # --- Resumen por red y consenso --------------------------------------------
    resumen, consenso = [], []
    for (red, var), lst in agrup.items():
        r = porred[red]
        v_his = lst[0][2]; nodos = lst[0][3]; nh = lst[0][4]
        V = np.vstack([z[1] for z in lst])
        conc = (V == v_his).mean(1)
        entre = [float((V[i] == V[j]).mean()) for i, j in itertools.combinations(range(len(V)), 2)]
        ari_entre = [ari(V[i], V[j]) for i, j in itertools.combinations(range(len(V)), 2)]
        pc = np.vstack([[z[0]["p_core"][n] for n in nodos] for z in lst]).mean(0)
        v_con = np.where(pc >= 0.5, 0, 1)
        frac_core = (V == 0).mean(0)
        inest = int(((frac_core > 0) & (frac_core < 1)).sum())
        mdls = [z[0]["mdl"] for z in lst]
        resumen.append(dict(red=red, tarea=r["tarea"], tipo=r["tipo"], grupo=r["grupo"], variante=var,
                            nodos=len(nodos), core_R2=int((v_his == 0).sum()),
                            core_replica_mediana=float(np.median((V == 0).sum(1))),
                            core_replica_min=int((V == 0).sum(1).min()), core_replica_max=int((V == 0).sum(1).max()),
                            conc_R2_media=round(float(conc.mean()), 4), conc_R2_min=round(float(conc.min()), 4),
                            conc_R2_max=round(float(conc.max()), 4),
                            ARI_R2_medio=round(float(np.mean([ari(v, v_his) for v in V])), 4),
                            cadenas_identicas_R2=int((V == v_his).all(1).sum()),
                            conc_entre_cadenas_media=round(float(np.mean(entre)), 4),
                            ARI_entre_cadenas_medio=round(float(np.mean(ari_entre)), 4),
                            nodos_inestables_entre_cadenas=inest,
                            conc_consenso_R2=round(float((v_con == v_his).mean()), 4),
                            discrepancias_consenso_R2=int((v_con != v_his).sum()),
                            mdl_R2=round(mdl_hist[(red, var)], 3),
                            mdl_replica_min=round(min(mdls), 3), mdl_replica_mediana=round(float(np.median(mdls)), 3)))
        for k, n in enumerate(nodos):
            consenso.append(dict(red=red, variante=var, node_id=n, name=nh[n], hs_R2=int(v_his[k]),
                                 p_core_media=round(float(pc[k]), 4), frac_cadenas_core=round(float(frac_core[k]), 4),
                                 hs_consenso=int(v_con[k]), discrepa=bool(v_con[k] != v_his[k])))
    resumen.sort(key=lambda z: (z["tipo"], z["tarea"], z["grupo"], z["variante"]))
    escribir(os.path.join(a.salida, "resumen_hubspoke.csv"), resumen)
    escribir(os.path.join(a.salida, "consenso.csv"), consenso)

    # --- Longitudes de descripción: hub-and-spoke frente a capas (S1/S2) --------
    filas_mdl = []
    for x in res_cap:
        filas_mdl.append(dict(red=x["red"], variante=x["variante"], modelo=f"capas_{x['capas']}",
                              semilla=x["semilla"], mdl=round(x["mdl"], 3), tamanos=x["tamanos"]))
    for x in res_hs:
        filas_mdl.append(dict(red=x["red"], variante=x["variante"], modelo="hubspoke",
                              semilla=x["semilla"], mdl=round(x["mdl"], 3),
                              tamanos=f"{sum(v == 0 for v in x['etiquetas'].values())}/"
                                      f"{sum(v == 1 for v in x['etiquetas'].values())}"))
    escribir(os.path.join(a.salida, "mdl_modelos.csv"), filas_mdl)

    cfg = {k: v for k, v in vars(a).items()}
    cfg.update(commit_libreria=COMMIT_LIB, segundos_totales=round(time.time() - t0, 1))
    with open(os.path.join(a.salida, "config.json"), "w", encoding="utf-8") as f:
        json.dump(cfg, f, indent=2, ensure_ascii=False)
    with open(os.path.join(a.salida, "entorno.txt"), "w", encoding="utf-8") as f:
        f.write(f"python {sys.version}\nplataforma {platform.platform()}\nnumpy {np.__version__}\n"
                f"networkx {nx.__version__}\nscipy {scipy.__version__}\ncore_periphery_sbm {COMMIT_LIB}\n")
    print("hecho", round(time.time() - t0), "s")


if __name__ == "__main__":
    main()
