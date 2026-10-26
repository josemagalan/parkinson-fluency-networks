# -*- coding: utf-8 -*-
# --- English ---
# Figures 1-3: co-occurrence and sequence networks (HC, ON, OFF) for action, phonological and semantic
# fluency, drawn in the style of the previous Gephi figures and coloured with the consensus core-periphery
# partitions of script 19 (neutral grey where the structure has no support). Deterministic ForceAtlas2
# layout and label adjustment written here with numpy (no Gephi).
# Inputs: outputs 14, 18 and 19. Outputs: analisis/output/51_figuras_redes/.
"""
51_figuras_redes.py

Propósito
    Regenerar las Figuras 1–3 del manuscrito (redes de co-ocurrencia y de secuencia de los tres
    grupos en cada tarea) desde las redes limpias de la R3 y las particiones de consenso del
    script 19, con el aspecto de las figuras de Gephi de la R2:
      - disposición ForceAtlas2 (Jacomy et al., 2014), implementada aquí con numpy y semilla fija,
        con los valores por defecto de Gephi (escala 2, gravedad 1, atracción lineal, peso 1;
        peso = frecuencia de co-ocurrencia o de transición del output 14);
      - componentes desconectadas dispuestas por separado y colocadas junto a la principal, y colas
        largas (percentil 95 del radio) acortadas radialmente, para que no fijen la escala;
      - ajuste de etiquetas para que no se solapen (equivalente al «Label Adjust» de Gephi);
      - aristas curvas en sentido horario, con el color degradado del nodo origen al destino;
      - tamaño de nodo y de etiqueta según el grado (escala dentro de cada panel);
      - core en azul y periferia en magenta (colores muestreados de las figuras de la R2);
      - nodos en gris neutro en las redes sin estructura core-periphery con soporte en 19
        (soporte_core = False: F secuencia ON y OFF), en lugar de colorear una partición inestable.
    La escala (puntos por unidad de disposición) es común a los tres grupos de cada fila, de modo
    que el tamaño relativo de las redes HC, ON y OFF se puede comparar a simple vista.

Punto editorial: coherencia de las Figuras 1–3 con las particiones de la R3 (E2/E5); sustituye las
imágenes de Gephi de la R2, que usaban las particiones históricas (revisión pre-Oviedo, 26-10-2026).

Entradas (desde revision_R3/)
    analisis/output/14_preparar_entrada_SBM/edges.csv      aristas con peso (frecuencia) y bucles
    analisis/output/18_exportar_redes_SBM/edges.csv        aristas binarias de las redes del SBM
    analisis/output/18_exportar_redes_SBM/index.csv        nodos y aristas por red
    analisis/output/19_SBM_redes_limpias/particiones_consenso.csv
    analisis/output/19_SBM_redes_limpias/resumen_redes.csv

Salidas  analisis/output/51_figuras_redes/
    Figura1_acciones.{pdf,png}, Figura2_F.{pdf,png}, Figura3_animales.{pdf,png}
        PDF: texto y nodos vectoriales, aristas rasterizadas a 600 ppp; PNG a 600 ppp.
    *_300ppp.png         versión ligera (300 ppp) para insertar en el DOCX
    posiciones.csv       coordenadas finales (puntos), grado, partición y P(core) de cada nodo
    verificaciones.csv   comprobaciones de nodos, aristas y tamaño del core frente a 14, 18 y 19
    leyendas_borrador.txt, entorno.txt, entradas_md5.csv, LEEME.md

Uso
    python -I analisis/python/51_figuras_redes.py
    Requiere numpy y matplotlib (con Pillow, que matplotlib ya instala) (referencia: Python 3.13, numpy 2.5.3, matplotlib 3.11.2).
    Las fuentes son DejaVu Sans/Serif, incluidas en matplotlib, para que el resultado no dependa
    de las fuentes instaladas. Tarda unos 5 minutos.
"""
import csv
import hashlib
import os
import platform
import sys
from collections import defaultdict

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.font_manager import FontProperties  # noqa: E402
from matplotlib.textpath import TextPath  # noqa: E402

SEMILLA = 20261026
SALIDA = "analisis/output/51_figuras_redes"
E14 = "analisis/output/14_preparar_entrada_SBM/edges.csv"
E18 = "analisis/output/18_exportar_redes_SBM/edges.csv"
I18 = "analisis/output/18_exportar_redes_SBM/index.csv"
P19 = "analisis/output/19_SBM_redes_limpias/particiones_consenso.csv"
R19 = "analisis/output/19_SBM_redes_limpias/resumen_redes.csv"

FIGURAS = [("Figura1_acciones", "actions"), ("Figura2_F", "F"), ("Figura3_animales", "animals")]
GRUPOS = [("HC", "Control"), ("ON", "On"), ("OFF", "Off")]
FILAS = [("undirected", "Co-occurrence network"), ("directed", "Sequence network")]

# Colores muestreados de las figuras de la R2 (Gephi)
COLOR = {"core": "#62A8F8", "periphery": "#E0187A", "ni": "#A6A6A6"}
BORDE = {"core": "#3D7FD0", "periphery": "#9E0E52", "ni": "#787878"}

# Dibujo (unidades en puntos tipográficos, figura de 18 × 11,6 pulgadas)
FIG_W, FIG_H = 18.0, 11.6
R_MIN, R_MAX = 1.3, 5.5          # radio del nodo
FS_MIN, FS_MAX = 2.2, 9.5        # tamaño de la etiqueta
CURVA = 0.20                     # desplazamiento del punto de control / longitud de la arista
SEGMENTOS = 16                   # tramos por arista para el degradado
ALFA = {"undirected": 0.30, "directed": 0.60}
ANCHO = {"undirected": 0.28, "directed": 0.40}

FUENTE = FontProperties(family="DejaVu Sans")
FUENTE_TIT = "DejaVu Serif"


# ----------------------------------------------------------------------------- utilidades
def md5(ruta):
    h = hashlib.md5()
    with open(ruta, "rb") as f:
        for bloque in iter(lambda: f.read(1 << 20), b""):
            h.update(bloque)
    return h.hexdigest()


def leer_csv(ruta):
    with open(ruta, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


# ----------------------------------------------------------------------------- datos
def cargar():
    nodos = defaultdict(dict)                         # graph -> palabra -> fila 19
    for r in leer_csv(P19):
        if r["family"] == "propia":
            nodos[r["graph"]][r["word_norm"]] = r
    resumen = {r["graph"]: r for r in leer_csv(R19) if r["family"] == "propia"}
    indice = {r["graph"]: r for r in leer_csv(I18)}
    aristas14 = defaultdict(list)
    bucles14 = defaultdict(int)
    for r in leer_csv(E14):
        if r["is_loop"].upper() == "TRUE":
            bucles14[r["graph"]] += 1
            continue
        aristas14[r["graph"]].append((r["from"], r["to"], float(r["weight"])))
    aristas18 = defaultdict(set)
    for r in leer_csv(E18):
        if r["graph"] in nodos:
            aristas18[r["graph"]].add((r["from"], r["to"]))
    return nodos, resumen, indice, aristas14, bucles14, aristas18


# ----------------------------------------------------------------------------- ForceAtlas2
def forceatlas2(n, pares, pesos, iteraciones, rng, escala=2.0, gravedad=1.0,
                influencia_peso=1.0, tolerancia=1.0):
    """ForceAtlas2 exacto (sin Barnes-Hut) con la velocidad adaptativa de Gephi 0.9.

    Atracción lineal w^influencia·d, repulsión kr·(g_i+1)(g_j+1)/d, gravedad kg·(g_i+1).
    """
    i, j = pares[:, 0], pares[:, 1]
    w = pesos ** influencia_peso
    grado = np.bincount(np.concatenate([i, j]), minlength=n).astype(float)
    masa = grado + 1.0
    pos = rng.uniform(-1, 1, size=(n, 2)) * np.sqrt(n) * 10
    fuerza_prev = np.zeros((n, 2))
    velocidad, eficiencia = 1.0, 1.0
    mm = escala * np.outer(masa, masa)
    np.fill_diagonal(mm, 0.0)
    for _ in range(iteraciones):
        delta = pos[:, None, :] - pos[None, :, :]
        d2 = (delta ** 2).sum(-1)
        np.fill_diagonal(d2, 1.0)
        f = (mm / d2)[:, :, None] * delta                     # repulsión
        fuerza = f.sum(1)
        dist = np.sqrt((pos ** 2).sum(1))
        dist[dist == 0] = 1e-9
        fuerza -= (gravedad * masa / dist)[:, None] * pos     # gravedad
        dv = pos[i] - pos[j]                                  # atracción lineal
        fa = w[:, None] * dv
        np.add.at(fuerza, i, -fa)
        np.add.at(fuerza, j, fa)
        # velocidad adaptativa (Gephi 0.9)
        swing = masa * np.sqrt(((fuerza - fuerza_prev) ** 2).sum(1))
        tracc = masa * 0.5 * np.sqrt(((fuerza + fuerza_prev) ** 2).sum(1))
        tot_s, tot_t = swing.sum(), tracc.sum()
        est = 0.05 * np.sqrt(n)
        jt = tolerancia * max(np.sqrt(est), min(10.0, est * tot_t / n ** 2))
        if tot_t > 0 and tot_s / tot_t > 2.0:
            if eficiencia > 0.05:
                eficiencia *= 0.5
            jt = max(jt, tolerancia)
        objetivo = jt * eficiencia * tot_t / tot_s if tot_s > 0 else velocidad
        if tot_s > jt * tot_t:
            if eficiencia > 0.05:
                eficiencia *= 0.7
        elif velocidad < 1000:
            eficiencia *= 1.3
        velocidad = velocidad + min(objetivo - velocidad, 0.5 * velocidad)
        factor = velocidad / (1.0 + np.sqrt(velocidad * swing))
        pos = pos + fuerza * factor[:, None]
        fuerza_prev = fuerza
    return pos - pos.mean(0)


# ----------------------------------------------------------------------------- etiquetas
def rotulo(palabra):
    """Etiqueta en mayúsculas; el tabulador de la celda ambigua partir<TAB>llegar (P11) se muestra como «/»."""
    return palabra.replace("\t", "/").upper()


def caja_etiqueta(texto, fs):
    ext = TextPath((0, 0), texto, size=fs, prop=FUENTE).get_extents()
    return ext.width, fs * 0.95


def ajustar_etiquetas(xy, ancho, alto, iteraciones=4000, margen=0.8, paso=0.4):
    """Separa las cajas de etiqueta solapadas (en puntos), como el «Label Adjust» de Gephi.

    Cada par solapado se aleja a lo largo de la recta que une los centros, la distancia justa
    para separar las cajas; así la red se expande de forma isótropa y conserva su forma.
    Los pares con el mismo centro (palabras con vecindarios idénticos) usan un ángulo fijo.
    """
    xy = xy.copy()
    n = len(xy)
    idx = np.arange(n)
    ang = (idx[:, None] * 2.399963 + idx[None, :] * 0.618034) % (2 * np.pi)
    ang = np.where(idx[:, None] > idx[None, :], ang, ang.T + np.pi)   # antisimétrico
    for _ in range(iteraciones):
        dx = xy[:, None, 0] - xy[None, :, 0]
        dy = xy[:, None, 1] - xy[None, :, 1]
        ox = (ancho[:, None] + ancho[None, :]) / 2 + margen - np.abs(dx)
        oy = (alto[:, None] + alto[None, :]) / 2 + margen - np.abs(dy)
        sol = (ox > 0) & (oy > 0)
        np.fill_diagonal(sol, False)
        if not sol.any():
            break
        d = np.sqrt(dx ** 2 + dy ** 2)
        cero = d < 1e-9
        ux = np.where(cero, np.cos(ang), dx / np.where(cero, 1.0, d))
        uy = np.where(cero, np.sin(ang), dy / np.where(cero, 1.0, d))
        t = np.minimum(ox / np.maximum(np.abs(ux), 1e-6), oy / np.maximum(np.abs(uy), 1e-6))
        t = np.where(sol, t, 0.0)
        cuenta = np.maximum(sol.sum(1), 1)
        mov = np.stack([(0.5 * t * ux).sum(1), (0.5 * t * uy).sum(1)], 1) / np.sqrt(cuenta)[:, None]
        xy += paso * mov
    return xy


def solapes(xy, ancho, alto):
    """Número de pares de cajas de etiqueta que se solapan (sin margen)."""
    dx = np.abs(xy[:, None, 0] - xy[None, :, 0])
    dy = np.abs(xy[:, None, 1] - xy[None, :, 1])
    sol = ((ancho[:, None] + ancho[None, :]) / 2 - dx > 1e-6) & ((alto[:, None] + alto[None, :]) / 2 - dy > 1e-6)
    np.fill_diagonal(sol, False)
    return int(sol.sum() // 2)


def componentes(n, pares):
    """Etiqueta de componente conexa de cada nodo (unión-búsqueda)."""
    padre = list(range(n))

    def raiz(a):
        while padre[a] != a:
            padre[a] = padre[padre[a]]
            a = padre[a]
        return a
    for a, b in pares:
        ra, rb = raiz(int(a)), raiz(int(b))
        if ra != rb:
            padre[max(ra, rb)] = min(ra, rb)
    return np.array([raiz(a) for a in range(n)])


def disposicion(n, pares, pesos, rng, iteraciones=1500):
    """ForceAtlas2 por componente; las componentes pequeñas se colocan alrededor de la mayor.

    En ForceAtlas2 las componentes desconectadas se alejan hasta que la gravedad compensa la
    repulsión, lo que encoge la red principal al dibujarla; Gephi lo resuelve a mano y aquí se
    hace de forma reproducible: anillo exterior, de mayor a menor, junto al radio típico.
    """
    comp = componentes(n, pares)
    grupos = sorted((np.flatnonzero(comp == c) for c in np.unique(comp)), key=lambda v: (-len(v), v[0]))
    pos = np.zeros((n, 2))
    radios = []
    for nodos_c in grupos:
        if len(nodos_c) == 1:
            radios.append(0.0)
            continue
        mapa = {int(v): q for q, v in enumerate(nodos_c)}
        sel = np.isin(pares[:, 0], nodos_c)
        sub = np.array([[mapa[int(a)], mapa[int(b)]] for a, b in pares[sel]], dtype=int)
        p = forceatlas2(len(nodos_c), sub, pesos[sel], iteraciones, rng)
        pos[nodos_c] = p
        radios.append(float(np.sqrt((p ** 2).sum(1)).max()))
    principal = grupos[0]
    r0 = np.percentile(np.sqrt((pos[principal] ** 2).sum(1)), 90) if len(principal) > 1 else 1.0
    angulo = np.pi / 4
    hueco = max(r0 * 0.08, 1.0)
    for nodos_c, rc in zip(grupos[1:], radios[1:]):
        dist = r0 * 1.10 + rc + hueco
        pos[nodos_c] += dist * np.array([np.cos(angulo), np.sin(angulo)])
        angulo += 2 * (rc + hueco) / dist + 0.15
    pos = pos - pos[principal].mean(0)
    # Cadenas largas (frecuentes en las redes de secuencia) fijarían la escala de toda la fila y
    # encogerían el resto: se acortan radialmente los nodos más allá del percentil 95 del radio,
    # para que el más lejano quede a 1,3 veces ese radio. Solo cambia la longitud de esas colas.
    r = np.sqrt((pos ** 2).sum(1))
    r95 = np.percentile(r, 95)
    if r.max() > 1.3 * r95 > 0:
        k = 0.3 * r95 / (r.max() - r95)
        fuera = r > r95
        nuevo = r95 + (r[fuera] - r95) * k
        pos[fuera] *= (nuevo / r[fuera])[:, None]
    return pos, len(grupos)


# ----------------------------------------------------------------------------- un panel
def preparar_panel(g, nodos, resumen, aristas14, rng):
    palabras = sorted(nodos[g])
    k = {p: a for a, p in enumerate(palabras)}
    sin_soporte = resumen[g]["soporte_core"].strip().lower() != "true"
    clase = []
    for p in palabras:
        if sin_soporte:
            clase.append("ni")
        else:
            clase.append("core" if nodos[g][p]["hs_consenso"] == "0" else "periphery")
    dirigidas = [(k[a], k[b], w) for a, b, w in aristas14[g]]
    # disposición: aristas no dirigidas con peso sumado (como trata Gephi un par recíproco)
    acum = defaultdict(float)
    for a, b, w in dirigidas:
        acum[(min(a, b), max(a, b))] += w
    pares = np.array(list(acum.keys()), dtype=int)
    pesos = np.array(list(acum.values()), dtype=float)
    pos, n_comp = disposicion(len(palabras), pares, pesos, rng)
    grado = np.bincount(pares.ravel(), minlength=len(palabras))
    return dict(graph=g, palabras=palabras, clase=clase, dirigidas=dirigidas, pos=pos,
                grado=grado, sin_soporte=sin_soporte, componentes=n_comp)


def tamanos(panel):
    gr = panel["grado"].astype(float)
    t = (gr - gr.min()) / max(gr.max() - gr.min(), 1.0)
    radio = R_MIN + (R_MAX - R_MIN) * t
    fs = FS_MIN + (FS_MAX - FS_MIN) * t
    cajas = [caja_etiqueta(rotulo(p), f) for p, f in zip(panel["palabras"], fs)]
    ancho = np.array([c[0] for c in cajas])
    alto = np.array([c[1] for c in cajas])
    return radio, fs, ancho, alto


def dibujar_panel(ax, panel, xy, medio_w, medio_h, red):
    clase, dirigidas = panel["clase"], panel["dirigidas"]
    radio, fs = panel["radio"], panel["fs"]
    rgba = {c: np.array(matplotlib.colors.to_rgba(COLOR[c])) for c in COLOR}
    orden = {"periphery": 0, "ni": 0, "core": 1}
    dirigidas = sorted(dirigidas, key=lambda e: orden[clase[e[0]]] + orden[clase[e[1]]])
    t = np.linspace(0, 1, SEGMENTOS + 1)[:, None]
    segs, cols, anchos = [], [], []
    for a, b, w in dirigidas:
        p0, p1 = xy[a], xy[b]
        d = p1 - p0
        c = (p0 + p1) / 2 + CURVA * np.array([d[1], -d[0]])     # arco en sentido horario
        pts = (1 - t) ** 2 * p0 + 2 * (1 - t) * t * c + t ** 2 * p1
        f = (np.arange(SEGMENTOS) + 0.5) / SEGMENTOS
        col = rgba[clase[a]][None, :] * (1 - f[:, None]) + rgba[clase[b]][None, :] * f[:, None]
        col[:, 3] = ALFA[red]
        segs.extend(np.stack([pts[:-1], pts[1:]], axis=1))
        cols.extend(col)
        anchos.extend([ANCHO[red] * (1 + 0.5 * np.log(w))] * SEGMENTOS)
    lc = LineCollection(segs, colors=cols, linewidths=anchos, capstyle="butt", zorder=1)
    lc.set_rasterized(True)
    ax.add_collection(lc)
    ordn = np.argsort(radio)
    ax.scatter(xy[ordn, 0], xy[ordn, 1], s=(2 * radio[ordn]) ** 2,
               c=[COLOR[clase[q]] for q in ordn], edgecolors=[BORDE[clase[q]] for q in ordn],
               linewidths=0.35, zorder=2)
    for q in ordn:
        ax.text(xy[q, 0], xy[q, 1], rotulo(panel["palabras"][q]), fontsize=fs[q],
                family="DejaVu Sans", ha="center", va="center", color="#111122", zorder=3)
    ax.set_xlim(-medio_w, medio_w)
    ax.set_ylim(-medio_h, medio_h)
    ax.set_aspect("equal")
    ax.axis("off")


# ----------------------------------------------------------------------------- figura
def figura(nombre, tarea, datos, rng, filas_pos, verif):
    nodos, resumen, indice, aristas14, bucles14, aristas18 = datos
    izq, der, arr, aba = 0.035, 0.995, 0.995, 0.005
    sep_col, sep_fila = 0.005, 0.05                   # hueco central para los rótulos de grupo
    ancho_p = (der - izq - 2 * sep_col) / 3
    alto_p = (arr - aba - sep_fila) / 2
    medio_w = ancho_p * FIG_W * 72 / 2                # semiancho del panel en puntos
    medio_h = alto_p * FIG_H * 72 / 2
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    for f, (red, rotulo) in enumerate(FILAS):
        paneles = []
        for c, (grp, _) in enumerate(GRUPOS):
            g = f"{tarea}_{grp}_{red}"
            p = preparar_panel(g, nodos, resumen, aristas14, rng)
            p["radio"], p["fs"], p["ancho"], p["alto"] = tamanos(p)
            paneles.append(p)
            verif.extend(verificar(g, p, nodos, resumen, indice, aristas14, bucles14, aristas18))
        # escala común en la fila: la red más extensa ocupa el panel; luego ajuste de etiquetas
        extension = max(max(np.abs(p["pos"][:, 0]).max() / medio_w,
                            np.abs(p["pos"][:, 1]).max() / medio_h) for p in paneles)
        escala = 0.94 / extension
        for _ in range(8):
            finales = []
            for p in paneles:
                xy = ajustar_etiquetas(p["pos"] * escala, p["ancho"], p["alto"])
                finales.append(xy)
            necesaria = max(max((np.abs(xy[:, 0]) + p["ancho"] / 2).max() / medio_w,
                                (np.abs(xy[:, 1]) + p["alto"] / 2).max() / medio_h)
                            for xy, p in zip(finales, paneles))
            if necesaria <= 0.99:
                break
            escala *= 0.97 / necesaria
        for p, xy in zip(paneles, finales):
            verif.append(dict(graph=p["graph"], comprobacion="etiquetas_solapadas_tras_ajuste",
                              valor=solapes(xy, p["ancho"], p["alto"]), referencia=0,
                              ok=solapes(xy, p["ancho"], p["alto"]) == 0))
        for c, (p, xy) in enumerate(zip(paneles, finales)):
            x0 = izq + c * (ancho_p + sep_col)
            y0 = aba + (1 - f) * (alto_p + sep_fila)
            ax = fig.add_axes([x0, y0, ancho_p, alto_p])
            dibujar_panel(ax, p, xy, medio_w, medio_h, red)
            for q, palabra in enumerate(p["palabras"]):
                fila19 = nodos[p["graph"]][palabra]
                filas_pos.append(dict(
                    figura=nombre, graph=p["graph"], task=tarea, group=GRUPOS[c][0], network=red,
                    word_norm=palabra, x_pt=round(float(xy[q, 0]), 3), y_pt=round(float(xy[q, 1]), 3),
                    grado=int(p["grado"][q]), clase=p["clase"][q], p_core=fila19["p_core"],
                    estable=fila19["estable"], escala_pt_por_unidad=round(float(escala), 6)))
        fig.text(0.012, aba + (1 - f) * (alto_p + sep_fila) + alto_p / 2, rotulo, rotation=90,
                 ha="center", va="center", fontsize=22, family=FUENTE_TIT)
    for c, (_, etiqueta) in enumerate(GRUPOS):
        fig.text(izq + c * (ancho_p + sep_col) + ancho_p / 2, aba + alto_p + sep_fila / 2,
                 etiqueta, ha="center", va="center", fontsize=22, family=FUENTE_TIT)
    base = os.path.join(SALIDA, nombre)
    meta = {"Creator": "51_figuras_redes.py", "CreationDate": None}
    fig.savefig(base + ".pdf", dpi=600, metadata=meta)
    fig.savefig(base + ".png", dpi=600, metadata={"Software": "51_figuras_redes.py"})
    fig.savefig(base + "_300ppp.png", dpi=300, metadata={"Software": "51_figuras_redes.py"})
    plt.close(fig)
    # PNG en RGB sobre fondo blanco (sin canal alfa): ≈30 % menos de tamaño, mismo aspecto
    from PIL import Image  # Pillow, dependencia de matplotlib
    Image.MAX_IMAGE_PIXELS = None
    for ruta, ppp in ((base + ".png", 600), (base + "_300ppp.png", 300)):
        with Image.open(ruta) as im:
            rgb = Image.new("RGB", im.size, "white")
            rgb.paste(im, mask=im.getchannel("A"))
        rgb.save(ruta, optimize=True, dpi=(ppp, ppp))


def verificar(g, p, nodos, resumen, indice, aristas14, bucles14, aristas18):
    res = resumen[g]
    a14 = {(a, b) for a, b, _ in aristas14[g]}
    if g.endswith("_undirected"):
        norm = lambda s: {tuple(sorted(e)) for e in s}  # noqa: E731
        iguales = norm(a14) == norm(aristas18[g])
    else:
        iguales = a14 == aristas18[g]
    n_core = sum(1 for c in p["clase"] if c == "core")
    filas = [
        ("nodos_19_igual_index_18", len(p["palabras"]), indice[g]["nodes"],
         len(p["palabras"]) == int(indice[g]["nodes"])),
        ("aristas_14_sin_bucles_igual_18", len(a14), len(aristas18[g]), iguales),
        ("aristas_18_igual_index_18", len(aristas18[g]), indice[g]["edges"],
         len(aristas18[g]) == int(indice[g]["edges"])),
        ("nodos_aislados_dibujados", int((p["grado"] == 0).sum()), "", True),
        ("componentes_conexas", p["componentes"], "", True),
        ("bucles_14_no_dibujados", bucles14[g], "", True),
        ("core_dibujado_igual_core_consenso_19",
         n_core, res["core_consenso"] if not p["sin_soporte"] else "sin soporte (gris)",
         (n_core == int(res["core_consenso"])) if not p["sin_soporte"] else n_core == 0),
    ]
    return [dict(graph=g, comprobacion=a, valor=b, referencia=c, ok=d) for a, b, c, d in filas]


def leyendas():
    tareas = {"actions": "Action", "F": "Phonological", "animals": "Semantic"}
    nota = ("Networks built from the words produced in the {t} fluency task. The top row shows the "
            "co-occurrence networks (undirected) and the bottom row the sequence networks (directed links "
            "between consecutively produced words, drawn as clockwise arcs from the earlier to the later "
            "word); edge width increases with the frequency of the co-occurrence or transition. "
            "Columns correspond to the control group and to the PD group on and off medication. Each node is a word; node and label size "
            "increase with degree. Nodes are coloured by the consensus hub-and-spoke core–periphery "
            "partition (blue, core; pink, periphery), estimated with the stochastic block model on the "
            "binary undirected version of each network.{extra} Layouts were computed with ForceAtlas2 "
            "(Jacomy et al., 2014) and "
            "drawn at a common scale within each row.")
    extra_F = (" The sequence networks of the PD group (on and off) showed no supported core–periphery "
               "structure and are shown in grey.")
    out = []
    for n, (nombre, tarea) in enumerate(FIGURAS, 1):
        t = tareas[tarea]
        out.append(f"Figure {n}. {t} fluency co-occurrence and sequence (directed) networks for the three "
                   f"groups.\nNote. " + nota.format(t=t.lower(), extra=extra_F if tarea == "F" else "") + "\n")
    return "\n".join(out)


def main():
    if not os.path.isdir("analisis"):
        sys.exit("Ejecutar desde revision_R3/")
    os.makedirs(SALIDA, exist_ok=True)
    datos = cargar()
    rng = np.random.default_rng(SEMILLA)
    filas_pos, verif = [], []
    for nombre, tarea in FIGURAS:
        figura(nombre, tarea, datos, rng, filas_pos, verif)
        print("ok", nombre)
    with open(os.path.join(SALIDA, "posiciones.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(filas_pos[0]))
        w.writeheader()
        w.writerows(filas_pos)
    with open(os.path.join(SALIDA, "verificaciones.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(verif[0]))
        w.writeheader()
        w.writerows(verif)
    with open(os.path.join(SALIDA, "leyendas_borrador.txt"), "w", encoding="utf-8") as f:
        f.write(leyendas())
    with open(os.path.join(SALIDA, "entorno.txt"), "w", encoding="utf-8") as f:
        f.write(f"python {sys.version}\nplataforma {platform.platform()}\nnumpy {np.__version__}\n"
                f"matplotlib {matplotlib.__version__}\nsemilla {SEMILLA}\n")
    with open(os.path.join(SALIDA, "entradas_md5.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["path", "md5"])
        for ruta in (E14, E18, I18, P19, R19, "analisis/python/51_figuras_redes.py"):
            w.writerow([ruta, md5(ruta)])
    fallos = [v for v in verif if not v["ok"]]
    print(f"verificaciones: {len(verif)} ({len(fallos)} fallos)")
    for v in fallos:
        print("  FALLO", v)


if __name__ == "__main__":
    main()
