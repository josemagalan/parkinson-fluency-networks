# -*- coding: utf-8 -*-
# --- English ---
# Hub-and-spoke SBM on the 210 R3 networks from script 18 (own, joint, size-matched): consensus of
# 10 chains (3 for size-matched), P(core), stability and MDL vs 'no core' and layered models (E1, E2, E5).
# Inputs: analisis/output/18_exportar_redes_SBM/. Outputs: analisis/output/19_SBM_redes_limpias/ (overwritten).
"""
19_SBM_redes_limpias.py

Propósito
    Ajustar el SBM core-periphery (Gallagher et al., 2021; D-021) sobre las redes R3
    exportadas por el script 18:
      - 18 redes propias completas y normalizadas (tarea × grupo × tipo);
      - 12 redes conjuntas (HC+ON+OFF y ON+OFF por tarea y tipo), que dan una partición
        de referencia común a los grupos comparados;
      - 180 redes igualadas HC/ON/OFF (extracciones condicionales 1–10 del script 07).
    Mismo procedimiento que la réplica validada (17): grafo no dirigido y binario, sin
    bucles; n_gibbs=100, n_mcmc=10·N, last_n_samples=50; varias cadenas con semilla fija;
    consenso por P(core) media y estabilidad entre cadenas.

Puntos editoriales: E1 (igualadas), E2/E5 (particiones para la nueva Tabla 2).

Entradas (desde revision_R3/)
    analisis/output/18_exportar_redes_SBM/{index,edges,nodes}.csv
    analisis/output/14_preparar_entrada_SBM/mapa_historico.csv   (comparación con la R2)
    analisis/python/vendor/core_periphery_sbm/

Salidas  analisis/output/19_SBM_redes_limpias/
    particiones_consenso.csv  por red (incluidas las igualadas, columna draw) y palabra: P(core)
                              media, fracción de cadenas en el core, etiqueta de consenso
                              (0 core / 1 periferia), estable
    cadenas.csv               por red y cadena: tamaño del core, MDL, aceptación
    resumen_redes.csv         por red: estabilidad, MDL hub-and-spoke, MDL sin core,
                              mejor modelo de capas, soporte de la estructura
    mdl_modelos.csv           MDL por modelo (hub-and-spoke y capas 2–4), redes completas
    comparacion_R2.csv        partición propia R3 frente a la histórica (por palabra)
    propia_vs_conjunta.csv    acuerdo entre la partición propia y la conjunta
    igualadas_resumen.csv     por extracción igualada: tamaño y proporción del core,
                              acuerdo con la red completa del mismo grupo
    entorno.txt, config.json, entradas_md5.csv, LEEME.md

Definiciones
    estable: el nodo cae en el mismo bloque en todas las cadenas.
    MDL sin core: longitud de descripción del hub-and-spoke con todos los nodos en la
        periferia (equivale a un grafo aleatorio sin estructura). «Soporte» = la mejor
        MDL hub-and-spoke es menor que la MDL sin core.

Uso
    python -I analisis/python/19_SBM_redes_limpias.py [--cadenas 10] [--cadenas-igualadas 3] [--procesos 2]
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
from collections import defaultdict
from multiprocessing import Pool

import numpy as np
import networkx as nx
import scipy

warnings.filterwarnings("ignore")
COMMIT_LIB = "90fb2896172d0489490c3574b4761af7e7f778f9"
A = None
G_CACHE = {}


def args_cli():
    p = argparse.ArgumentParser()
    p.add_argument("--entrada", default="analisis/output/18_exportar_redes_SBM")
    p.add_argument("--mapa", default="analisis/output/14_preparar_entrada_SBM/mapa_historico.csv")
    p.add_argument("--lib", default="analisis/python/vendor/core_periphery_sbm")
    p.add_argument("--salida", default="analisis/output/19_SBM_redes_limpias")
    p.add_argument("--cadenas", type=int, default=10)
    p.add_argument("--cadenas-igualadas", type=int, default=3)
    p.add_argument("--cadenas-capas", type=int, default=3)
    p.add_argument("--capas", default="2,3,4")
    p.add_argument("--n-gibbs", type=int, default=100)
    p.add_argument("--mcmc-por-nodo", type=int, default=10)
    p.add_argument("--ultimas", type=int, default=50)
    p.add_argument("--n-mdl", type=int, default=100000)
    p.add_argument("--n-mdl-igualadas", type=int, default=20000)
    p.add_argument("--procesos", type=int, default=2)
    p.add_argument("--solo", default="", help="subcadena de nombre de red (pruebas)")
    p.add_argument("--sin-igualadas-en-particiones", action="store_true",
                   help="no escribir en particiones_consenso.csv las redes igualadas (comportamiento anterior)")
    return p.parse_args()


def init_hijo(a, grafos):
    global A, cp, mf, G_CACHE
    A = a
    G_CACHE = grafos
    sys.path.insert(0, os.path.abspath(a.lib))
    from core_periphery_sbm import core_periphery as cp_  # noqa
    from core_periphery_sbm import model_fit as mf_  # noqa
    cp, mf = cp_, mf_


def leer_csv(ruta):
    with open(ruta, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def escribir(ruta, filas):
    if not filas:
        return
    with open(ruta, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0].keys()), lineterminator="\n")
        w.writeheader(); w.writerows(filas)


def md5(ruta):
    with open(ruta, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


def cargar(a):
    idx = leer_csv(os.path.join(a.entrada, "index.csv"))
    if a.solo:
        idx = [r for r in idx if a.solo in r["graph"]]
    nombres = {r["graph"] for r in idx}
    grafos = {g: nx.Graph() for g in nombres}
    for r in leer_csv(os.path.join(a.entrada, "nodes.csv")):
        if r["graph"] in grafos:
            grafos[r["graph"]].add_node(r["word_norm"])
    for r in leer_csv(os.path.join(a.entrada, "edges.csv")):
        if r["graph"] in grafos and r["from"] != r["to"]:
            grafos[r["graph"]].add_edge(r["from"], r["to"])  # dirigida → no dirigida (D-021)
    for r in idx:
        assert nx.number_of_isolates(grafos[r["graph"]]) == 0, r["graph"]
    return idx, grafos


def vec(G, d):
    return np.array([d[n] for n in sorted(G.nodes())], dtype=int)


def ari(x, y):
    x = np.asarray(x); y = np.asarray(y)
    cx = {v: i for i, v in enumerate(np.unique(x))}; cy = {v: i for i, v in enumerate(np.unique(y))}
    M = np.zeros((len(cx), len(cy)))
    for a_, b_ in zip(x, y):
        M[cx[a_], cy[b_]] += 1
    c2 = lambda v: v * (v - 1) / 2
    s = c2(M).sum(); a_s = c2(M.sum(1)).sum(); b_s = c2(M.sum(0)).sum(); n = c2(len(x))
    esp = a_s * b_s / n; mx = (a_s + b_s) / 2
    return float((s - esp) / (mx - esp)) if mx != esp else 1.0


def t_hubspoke(job):
    g, semilla, n_mdl = job
    G = G_CACHE[g]
    m = cp.HubSpokeCorePeriphery(n_gibbs=A.n_gibbs, n_mcmc=A.mcmc_por_nodo * len(G), seed=semilla)
    t0 = time.time(); m.infer(G)
    lab = m.get_labels(last_n_samples=A.ultimas, prob=False, return_dict=True)
    prob = m.get_labels(last_n_samples=A.ultimas, prob=True, return_dict=True)
    np.random.seed(10_000 + semilla)
    mdl = mf.mdl_hubspoke(G, vec(G, lab), n_samples=n_mdl)
    return dict(graph=g, semilla=semilla, etiquetas=lab, p_core={n: float(p[0]) for n, p in prob.items()},
                mdl=float(mdl), segundos=time.time() - t0, aceptacion=float(np.mean(m.accepts)))


def t_sin_core(job):
    g, n_mdl = job
    G = G_CACHE[g]
    np.random.seed(30_000)
    return g, float(mf.mdl_hubspoke(G, np.ones(len(G), dtype=int), n_samples=n_mdl))


def t_capas(job):
    g, L, semilla = job
    G = G_CACHE[g]
    m = cp.LayeredCorePeriphery(n_layers=L, n_gibbs=A.n_gibbs, n_mcmc=A.mcmc_por_nodo * len(G), seed=semilla)
    if m.infer(G) is False:
        return dict(graph=g, capas=L, semilla=semilla, mdl=float("nan"), tamanos="")
    lab = m.get_labels(last_n_samples=A.ultimas, prob=False, return_dict=False)
    np.random.seed(20_000 + 100 * L + semilla)
    mdl = mf.mdl_layered(G, lab, n_layers=L, n_samples=A.n_mdl)
    return dict(graph=g, capas=L, semilla=semilla, mdl=float(mdl), tamanos="/".join(map(str, np.bincount(lab, minlength=L))))


def mapa_imap(pool, f, jobs, etiqueta, t0):
    out = []
    for k, x in enumerate(pool.imap_unordered(f, jobs, chunksize=1), 1):
        out.append(x)
        if k % 25 == 0 or k == len(jobs):
            print(f"  {etiqueta} {k}/{len(jobs)} {time.time()-t0:.0f}s", flush=True)
    return out


A_GLOBAL = None


def main():
    global A_GLOBAL
    a = args_cli()
    A_GLOBAL = a
    os.makedirs(a.salida, exist_ok=True)
    idx, grafos = cargar(a)
    info = {r["graph"]: r for r in idx}
    completas = [r["graph"] for r in idx if r["family"] != "igualada"]
    igualadas = [r["graph"] for r in idx if r["family"] == "igualada"]
    capas = [int(x) for x in a.capas.split(",") if x]
    jobs_hs = [(g, s, a.n_mdl) for g in completas for s in range(1, a.cadenas + 1)] + \
              [(g, s, a.n_mdl_igualadas) for g in igualadas for s in range(1, a.cadenas_igualadas + 1)]
    # Primero las redes grandes para equilibrar los procesos
    jobs_hs.sort(key=lambda j: -grafos[j[0]].number_of_edges())
    jobs_sc = [(g, a.n_mdl if g in completas else a.n_mdl_igualadas) for g in completas + igualadas]
    jobs_cap = sorted([(g, L, s) for g in completas for L in capas for s in range(1, a.cadenas_capas + 1)],
                      key=lambda j: -grafos[j[0]].number_of_edges())
    t0 = time.time()
    with Pool(a.procesos, initializer=init_hijo, initargs=(a, grafos)) as pool:
        res_hs = mapa_imap(pool, t_hubspoke, jobs_hs, "hub-and-spoke", t0)
        res_sc = dict(mapa_imap(pool, t_sin_core, jobs_sc, "sin core", t0))
        res_cap = mapa_imap(pool, t_capas, jobs_cap, "capas", t0)

    por_red = defaultdict(list)
    for x in res_hs:
        por_red[x["graph"]].append(x)
    filas_cad, filas_part, filas_res = [], [], []
    consenso = {}
    for g in sorted(por_red):
        lst = sorted(por_red[g], key=lambda x: x["semilla"])
        G = grafos[g]; nodos = sorted(G.nodes())
        V = np.vstack([[x["etiquetas"][n] for n in nodos] for x in lst])
        pc = np.vstack([[x["p_core"][n] for n in nodos] for x in lst]).mean(0)
        frac = (V == 0).mean(0)
        v_con = np.where(pc >= 0.5, 0, 1)
        estable = (frac == 0) | (frac == 1)
        consenso[g] = dict(zip(nodos, v_con)); consenso[g + "#p"] = dict(zip(nodos, pc))
        for x, v in zip(lst, V):
            filas_cad.append(dict(graph=g, semilla=x["semilla"], nodos=len(nodos), core=int((v == 0).sum()),
                                  mdl=round(x["mdl"], 3), segundos=round(x["segundos"], 1), aceptacion=round(x["aceptacion"], 4)))
        # Desde el 22-10-2026 también se guardan las particiones de las redes igualadas (columna draw), que
        # necesita la comprobación E1 de niveles (script 22). Las redes completas no cambian.
        if info[g]["family"] != "igualada" or not A_GLOBAL.sin_igualadas_en_particiones:
            for k, n in enumerate(nodos):
                filas_part.append(dict(graph=g, family=info[g]["family"], task=info[g]["task"], group=info[g]["group"],
                                       network=info[g]["network"], draw=info[g]["draw"], word_norm=n, p_core=round(float(pc[k]), 4),
                                       frac_cadenas_core=round(float(frac[k]), 4), hs_consenso=int(v_con[k]),
                                       estable=bool(estable[k]), grado=G.degree(n)))
        pares = list(itertools.combinations(range(len(V)), 2))
        mdls = [x["mdl"] for x in lst]
        filas_res.append(dict(graph=g, family=info[g]["family"], task=info[g]["task"], group=info[g]["group"],
                              network=info[g]["network"], draw=info[g]["draw"], nodos=len(nodos), aristas=G.number_of_edges(),
                              cadenas=len(lst), core_consenso=int((v_con == 0).sum()),
                              core_min=int((V == 0).sum(1).min()), core_max=int((V == 0).sum(1).max()),
                              conc_entre_cadenas=round(float(np.mean([(V[i] == V[j]).mean() for i, j in pares])), 4) if pares else "",
                              ARI_entre_cadenas=round(float(np.mean([ari(V[i], V[j]) for i, j in pares])), 4) if pares else "",
                              nodos_inestables=int((~estable).sum()),
                              mdl_hubspoke_min=round(min(mdls), 3), mdl_sin_core=round(res_sc[g], 3),
                              soporte_core=bool(min(mdls) < res_sc[g])))

    # Modelos de capas (redes completas) y comparación de MDL
    filas_mdl = [dict(graph=x["graph"], modelo=f"capas_{x['capas']}", semilla=x["semilla"], mdl=round(x["mdl"], 3),
                      tamanos=x["tamanos"]) for x in res_cap]
    for x in res_hs:
        if info[x["graph"]]["family"] != "igualada":
            filas_mdl.append(dict(graph=x["graph"], modelo="hubspoke", semilla=x["semilla"], mdl=round(x["mdl"], 3),
                                  tamanos=f"{sum(v == 0 for v in x['etiquetas'].values())}/{sum(v == 1 for v in x['etiquetas'].values())}"))
    for g in completas:
        filas_mdl.append(dict(graph=g, modelo="sin_core", semilla="", mdl=round(res_sc[g], 3), tamanos=f"0/{len(grafos[g])}"))
    filas_mdl.sort(key=lambda r: (r["graph"], r["modelo"], str(r["semilla"])))
    mejor_capas = {}
    for x in res_cap:
        if not np.isnan(x["mdl"]):
            k = x["graph"]
            if k not in mejor_capas or x["mdl"] < mejor_capas[k][0]:
                mejor_capas[k] = (x["mdl"], x["capas"])
    for r in filas_res:
        if r["graph"] in mejor_capas:
            r["mdl_capas_min"] = round(mejor_capas[r["graph"]][0], 3)
            r["capas_mejor"] = mejor_capas[r["graph"]][1]
            r["modelo_ganador"] = "hubspoke" if r["mdl_hubspoke_min"] <= r["mdl_capas_min"] else f"capas_{r['capas_mejor']}"
            if r["mdl_sin_core"] < min(r["mdl_hubspoke_min"], r["mdl_capas_min"]):
                r["modelo_ganador"] = "sin_core"
        else:
            r["mdl_capas_min"] = r["capas_mejor"] = r["modelo_ganador"] = ""

    # Comparación con la partición histórica R2 (por palabra normalizada; colisiones excluidas)
    mapa = leer_csv(a.mapa)
    hist = defaultdict(set)
    for r in mapa:
        tipo = "directed" if r["network"] == "directed" else "undirected"
        hist[(r["task"], r["condition"], tipo, r["word_norm"])].add(int(r["historical_partition"]))
    filas_r2 = []
    for g in completas:
        if info[g]["family"] != "propia":
            continue
        t, c, tipo = info[g]["task"], info[g]["group"], info[g]["network"]
        n_ok = n_tot = n_col = 0
        for w, v in consenso[g].items():
            h = hist.get((t, c, tipo, w))
            if not h:
                continue
            if len(h) > 1:
                n_col += 1; continue
            n_tot += 1; n_ok += int(next(iter(h)) == v)
        filas_r2.append(dict(graph=g, palabras_comparadas=n_tot, coinciden=n_ok,
                             concordancia=round(n_ok / n_tot, 4) if n_tot else "", colisiones_excluidas=n_col))

    # Propia frente a conjunta
    filas_pc = []
    for g in completas:
        if info[g]["family"] != "propia":
            continue
        t, c, tipo = info[g]["task"], info[g]["group"], info[g]["network"]
        for fam in ("conjunta_3", "conjunta_PD"):
            if fam == "conjunta_PD" and c == "HC":
                continue
            gj = f"{t}_{fam}_{tipo}"
            if gj not in consenso:
                continue
            ws = sorted(consenso[g])
            x = np.array([consenso[g][w] for w in ws]); y = np.array([consenso[gj][w] for w in ws])
            filas_pc.append(dict(graph=g, referencia=gj, palabras=len(ws), concordancia=round(float((x == y).mean()), 4),
                                 ARI=round(ari(x, y), 4), core_propia=int((x == 0).sum()), core_en_conjunta=int((y == 0).sum())))

    # Igualadas: tamaño del core y acuerdo con la red completa del mismo grupo
    filas_ig = []
    for g in igualadas:
        t, c, tipo, dr = info[g]["task"], info[g]["group"], info[g]["network"], info[g]["draw"]
        gfull = f"{t}_{c}_{tipo}"
        ws = sorted(consenso[g]); x = np.array([consenso[g][w] for w in ws])
        r = next(rr for rr in filas_res if rr["graph"] == g)
        fila = dict(graph=g, task=t, group=c, network=tipo, draw=dr, nodos=len(ws), core=int((x == 0).sum()),
                    prop_core=round(float((x == 0).mean()), 4), conc_entre_cadenas=r["conc_entre_cadenas"],
                    soporte_core=r["soporte_core"])
        if gfull in consenso:
            y = np.array([consenso[gfull][w] for w in ws])
            fila.update(conc_con_completa=round(float((x == y).mean()), 4), ARI_con_completa=round(ari(x, y), 4))
        filas_ig.append(fila)

    filas_res.sort(key=lambda r: (r["family"] != "propia", r["family"], r["task"], r["network"], r["group"], str(r["draw"])))
    escribir(os.path.join(a.salida, "particiones_consenso.csv"), filas_part)
    escribir(os.path.join(a.salida, "cadenas.csv"), sorted(filas_cad, key=lambda r: (r["graph"], r["semilla"])))
    escribir(os.path.join(a.salida, "resumen_redes.csv"), filas_res)
    escribir(os.path.join(a.salida, "mdl_modelos.csv"), filas_mdl)
    escribir(os.path.join(a.salida, "comparacion_R2.csv"), filas_r2)
    escribir(os.path.join(a.salida, "propia_vs_conjunta.csv"), filas_pc)
    escribir(os.path.join(a.salida, "igualadas_resumen.csv"), sorted(filas_ig, key=lambda r: (r["task"], r["network"], r["group"], int(r["draw"]))))
    entradas = [os.path.join(a.entrada, f) for f in ("index.csv", "edges.csv", "nodes.csv")] + [a.mapa]
    escribir(os.path.join(a.salida, "entradas_md5.csv"),
             [dict(path=p.replace("\\", "/"), md5=md5(p)) for p in entradas] +
             [dict(path="analisis/python/19_SBM_redes_limpias.py", md5=md5(os.path.abspath(__file__)))])
    cfg = dict(vars(a)); cfg.update(commit_libreria=COMMIT_LIB, segundos_totales=round(time.time() - t0, 1))
    with open(os.path.join(a.salida, "config.json"), "w", encoding="utf-8") as f:
        json.dump(cfg, f, indent=2, ensure_ascii=False)
    with open(os.path.join(a.salida, "entorno.txt"), "w", encoding="utf-8") as f:
        f.write(f"python {sys.version}\nplataforma {platform.platform()}\nnumpy {np.__version__}\n"
                f"networkx {nx.__version__}\nscipy {scipy.__version__}\ncore_periphery_sbm {COMMIT_LIB}\n")
    print("hecho", round(time.time() - t0), "s")


if __name__ == "__main__":
    main()
