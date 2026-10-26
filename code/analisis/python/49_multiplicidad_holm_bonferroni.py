# -*- coding: utf-8 -*-
# --- English ---
# Multiplicity sensitivity: re-adjusts the stored raw p-values of every primary family with Holm and
# Bonferroni and compares decisions with the Benjamini-Hochberg adjustment used in the paper (E5).
# Inputs: outputs 04_v2, 06, 20_v2, 21_v2, 48. Outputs: analisis/output/49_multiplicidad_holm_bonferroni/.
"""
49_multiplicidad_holm_bonferroni.py

Propósito
    Comprobar si alguna conclusión depende del método de corrección por comparaciones múltiples.
    Se vuelven a ajustar los p crudos ya guardados de cada familia (las mismas familias, sin
    cambiar modelos ni p) con Holm y Bonferroni, y se comparan las decisiones a 0,05 con BH,
    que es el método del manuscrito. Primero se comprueba que BH reproduce el p_BH guardado.
    No sustituye a BH ni cambia la Tabla 2: es una comprobación interna preparada el 25-10-2026,
    después de conocer los resultados, para responder si un revisor pregunta por el método.
    Se guarda en el repositorio; el manuscrito solo justifica BH (D-033).

Punto editorial: E5 (cómo se controlaron las comparaciones múltiples).

Entradas (desde revision_R3/)
    analisis/output/04_produccion_sesion_v2/resultados.csv      (p, p_BH6, family)
    analisis/output/06_psicolinguistica_sesion/resultados.csv   (p, p_BH, family)
    analisis/output/20_tabla2_core_periferia_v2/resultados.csv  (p, p_BH, family)
    analisis/output/21_tabla2_niveles_R2_v2/resultados.csv      (p, p_BH, family)
    analisis/output/48_gradiente_por_grupo/gradientes.csv       (p, p_BH, family)
    Solo filas con p ajustado guardado (las n.i. y el intercepto quedan fuera, como en BH).

Salidas  analisis/output/49_multiplicidad_holm_bonferroni/
    resumen.csv        por fuente y familia: m, significativos nominales, BH, Holm, Bonferroni y
                       si BH reproduce el p ajustado guardado
    detalle.csv        cada test con p, p BH, p Holm, p Bonferroni
    cambios.csv        tests significativos con BH que dejan de serlo con Holm o Bonferroni
    entradas_md5.csv, LEEME.md

Uso
    python -I analisis/python/49_multiplicidad_holm_bonferroni.py
    Solo biblioteca estándar de Python (>= 3.8).
"""
import csv
import hashlib
import os
from collections import OrderedDict

ALFA = 0.05
FUENTES = [
    ("04_produccion", "analisis/output/04_produccion_sesion_v2/resultados.csv", "p_BH6"),
    ("06_medias", "analisis/output/06_psicolinguistica_sesion/resultados.csv", "p_BH"),
    ("20_brechas", "analisis/output/20_tabla2_core_periferia_v2/resultados.csv", "p_BH"),
    ("21_tabla2", "analisis/output/21_tabla2_niveles_R2_v2/resultados.csv", "p_BH"),
    ("48_gradiente_grupo", "analisis/output/48_gradiente_por_grupo/gradientes.csv", "p_BH"),
]
SALIDA = "analisis/output/49_multiplicidad_holm_bonferroni"
ID_COLS = ("task", "network", "variable", "term", "effect")


def bh(p):
    m = len(p)
    orden = sorted(range(m), key=lambda i: p[i])
    aj = [0.0] * m
    previo = 1.0
    for rango in range(m, 0, -1):
        i = orden[rango - 1]
        previo = min(previo, p[i] * m / rango)
        aj[i] = min(previo, 1.0)
    return aj


def holm(p):
    m = len(p)
    orden = sorted(range(m), key=lambda i: p[i])
    aj = [0.0] * m
    previo = 0.0
    for k, i in enumerate(orden):
        previo = max(previo, p[i] * (m - k))
        aj[i] = min(previo, 1.0)
    return aj


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def md5(ruta):
    h = hashlib.md5()
    with open(ruta, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def main():
    os.makedirs(SALIDA, exist_ok=True)
    resumen, detalle, cambios, entradas = [], [], [], []
    for fuente, ruta, col_bh in FUENTES:
        entradas.append({"path": ruta, "md5": md5(ruta)})
        familias = OrderedDict()
        with open(ruta, newline="", encoding="utf-8") as f:
            for fila in csv.DictReader(f):
                if num(fila.get(col_bh)) is None or num(fila.get("p")) is None:
                    continue
                familias.setdefault(fila["family"], []).append(fila)
        for fam, filas in familias.items():
            p = [num(r["p"]) for r in filas]
            guardado = [num(r[col_bh]) for r in filas]
            p_bh, p_holm = bh(p), holm(p)
            p_bonf = [min(x * len(p), 1.0) for x in p]
            reproduce = all(abs(a - b) <= 1e-12 for a, b in zip(p_bh, guardado))
            resumen.append({
                "fuente": fuente, "family": fam, "m": len(p),
                "nominal": sum(x < ALFA for x in p), "BH": sum(x < ALFA for x in p_bh),
                "Holm": sum(x < ALFA for x in p_holm), "Bonferroni": sum(x < ALFA for x in p_bonf),
                "BH_reproduce_guardado": reproduce})
            for r, a, b, c, d in zip(filas, p, p_bh, p_holm, p_bonf):
                ident = {k: r.get(k, "") for k in ID_COLS}
                fila = dict(fuente=fuente, family=fam, **ident, p=a, p_BH=b, p_Holm=c, p_Bonferroni=d)
                detalle.append(fila)
                if b < ALFA and (c >= ALFA or d >= ALFA):
                    cambios.append(fila)
    def escribir(nombre, filas):
        with open(os.path.join(SALIDA, nombre), "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(filas[0].keys()))
            w.writeheader()
            w.writerows(filas)
    escribir("resumen.csv", resumen)
    escribir("detalle.csv", detalle)
    escribir("cambios.csv", cambios)
    escribir("entradas_md5.csv", entradas)
    if not all(r["BH_reproduce_guardado"] for r in resumen):
        raise SystemExit("BH no reproduce algún p ajustado guardado: revisar antes de usar.")
    with open(os.path.join(SALIDA, "LEEME.md"), "w", encoding="utf-8") as f:
        f.write(
            "# 49 Multiplicidad: Holm y Bonferroni frente a BH\n\n"
            "Sensibilidad interna (D-033), preparada el 25-10-2026 después de conocer los resultados. "
            "Mismas familias y p crudos que los scripts 04, 06, 20, 21 y 48; solo cambia el ajuste. "
            "BH reproduce todos los p ajustados guardados.\n\n"
            "- resumen.csv: significativos a 0,05 por familia con cada método.\n"
            "- cambios.csv: tests significativos con BH que no lo son con Holm o Bonferroni.\n\n"
            "Lectura: ningún contraste directo paciente-control ni ON-OFF es significativo con ningún método; "
            "las señales aisladas (interacción OFF x Periphery en concreción de acciones secuencia, sesión en "
            "familiaridad F co-ocurrencia) y la producción ON-OFF se mantienen con Bonferroni. Solo cambia "
            "cuántas celdas del gradiente core-periferia son significativas; la frecuencia lo es en todos los "
            "grupos y redes con cualquier método. El manuscrito mantiene BH; esto no se publica como análisis "
            "adicional.\n")
    print("Familias:", len(resumen), "| tests:", len(detalle), "| cambian:", len(cambios))


if __name__ == "__main__":
    main()
