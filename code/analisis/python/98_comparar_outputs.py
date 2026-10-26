# -*- coding: utf-8 -*-
# --- English ---
# Reproducibility check: compares a run with the reference fingerprints (exact MD5, MD5 ignoring
# CRLF/LF, CSV rounded to 10 significant digits). Exit code 1 if any file differs or is missing.
"""
98_comparar_outputs.py

Propósito
    Comprobar si una ejecución nueva de la cadena (01–17) reproduce los outputs de
    referencia, sin necesitar los ficheros de referencia: basta con sus huellas.
    Para cada fichero se guardan tres huellas:
      md5      bytes exactos
      md5_lf   ignorando la diferencia de fin de línea CRLF/LF (Windows/Linux)
      md5_r10  solo CSV: números redondeados a 10 cifras significativas, para
               absorber el ruido de coma flotante entre plataformas (~1e-14)
    Con --csv-ref (carpeta con los CSV de referencia) se hace además la
    comparación celda a celda con tolerancia relativa --rtol.

Uso (desde revision_R3/)
    Comparar:   python -I analisis/python/98_comparar_outputs.py [--informe f.csv]
    Generar la referencia en la máquina de referencia:
                python -I analisis/python/98_comparar_outputs.py --generar \\
                    --salida ../referencia/outputs_md5_windows.csv [--carpetas REGEX]

Estados: idéntico · idéntico (salvo CRLF/LF) · igual (redondeo 1e-10) · igual
(tolerancia) · binario distinto (.rds/.png/.svg; esperable entre plataformas) ·
omitido · no ejecutado (no existe la carpeta del script) · distinto · falta.
Sale con 1 si hay «distinto» o «falta».

Solo biblioteca estándar. Las carpetas con sufijo _vN se normalizan (p. ej.,
07_redes_HC_bootstrap_v2 → 07_redes_HC_bootstrap); si hay varias versiones se usa
la de sufijo mayor.
"""
import argparse
import csv
import hashlib
import math
import os
import re
import sys

TEXTO = (".csv", ".graphml", ".txt", ".json", ".md")
BINARIOS = (".rds", ".png", ".svg")
EXCLUIR_GENERAR = ("sessionInfo.txt", "entorno.txt", "config.json")
# Ficheros que no genera la cadena o que dependen del entorno
OMITIR = {
    "09_auditoria_secuencia/figura4_R2_extraida.png": "extraída aparte del .docx R2; no la genera el script 09",
    "16_verificacion_integrada/huellas.csv": "huellas de manifiestos de la máquina que ejecuta",
    "16_verificacion_integrada/comprobaciones.csv": "el nº de huellas depende de los outputs presentes; "
                                                    "16 se detiene si alguna de sus 60 comprobaciones falla",
}
REF_DEFECTO = ["../referencia/outputs_md5_windows.csv", "../referencia/outputs_md5_linux.csv"]


def md5_bytes(b):
    return hashlib.md5(b).hexdigest()


def leer(ruta):
    with open(ruta, "rb") as f:
        return f.read()


def md5_lf(b):
    return md5_bytes(b.replace(b"\r\n", b"\n"))


def _redondear(x, sig=10):
    try:
        v = float(x)
    except ValueError:
        return x
    if math.isnan(v) or math.isinf(v):
        return x
    return format(v, f".{sig}g")


def md5_r10(b):
    texto = b.replace(b"\r\n", b"\n").decode("utf-8", errors="replace")
    filas = csv.reader(texto.splitlines())
    s = "\n".join("\x1f".join(_redondear(c) for c in fila) for fila in filas)
    return md5_bytes(s.encode("utf-8"))


def norm(rel):
    partes = rel.replace("\\", "/").split("/")
    partes[0] = re.sub(r"_v\d+$", "", partes[0])
    return "/".join(partes)


def indice(carpeta, patron=None):
    out = {}
    for raiz, _, fs in os.walk(carpeta):
        for f in fs:
            p = os.path.join(raiz, f)
            rel = os.path.relpath(p, carpeta).replace("\\", "/")
            if patron and not re.match(patron, rel.split("/")[0]):
                continue
            clave = norm(rel)
            v = re.match(r"^[^/]*_v(\d+)(/|$)", rel)
            n = int(v.group(1)) if v else 1
            if clave not in out or n > out[clave][1]:
                out[clave] = (p, n, rel.split("/")[0])
    return out


def num(x):
    try:
        return float(x)
    except ValueError:
        return None


def comparar_csv(a, b, rtol):
    with open(a, encoding="utf-8", newline="") as fa, open(b, encoding="utf-8", newline="") as fb:
        A = list(csv.reader(fa)); B = list(csv.reader(fb))
    if len(A) != len(B):
        return False, f"filas {len(A)} vs {len(B)}"
    peor = 0.0
    for i, (ra, rb) in enumerate(zip(A, B)):
        if len(ra) != len(rb):
            return False, f"columnas distintas en fila {i}"
        for x, y in zip(ra, rb):
            if x == y:
                continue
            nx_, ny_ = num(x), num(y)
            if nx_ is None or ny_ is None:
                return False, f"texto distinto en fila {i}: {x[:30]!r} vs {y[:30]!r}"
            if math.isnan(nx_) and math.isnan(ny_):
                continue
            peor = max(peor, abs(nx_ - ny_) / max(1.0, abs(nx_), abs(ny_)))
    return peor <= rtol, f"máx. diferencia relativa {peor:.3g}"


def generar(a):
    idx = indice(a.output, a.carpetas)
    filas = []
    for k in sorted(idx):
        p, _, carpeta = idx[k]
        nombre = os.path.basename(p)
        if nombre in EXCLUIR_GENERAR or nombre.startswith("LEEME") or nombre == "procedencia_figura.md":
            continue
        b = leer(p)
        filas.append(dict(fichero=k, md5=md5_bytes(b), md5_lf=md5_lf(b),
                          md5_r10=md5_r10(b) if k.endswith(".csv") else "", carpeta_origen=carpeta))
    with open(a.salida, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["fichero", "md5", "md5_lf", "md5_r10", "carpeta_origen"], lineterminator="\n")
        w.writeheader(); w.writerows(filas)
    print(f"{len(filas)} huellas escritas en {a.salida}")


def comparar(a):
    ref = {}
    for fr in a.referencia:
        if os.path.isfile(fr):
            with open(fr, encoding="utf-8") as f:
                for r in csv.DictReader(f):
                    ref[norm(r["fichero"])] = r
    if not ref:
        sys.exit("No se encontró ningún fichero de referencia")
    nuevo = {k: v[0] for k, v in indice(a.output).items()}
    refcsv = {k: v[0] for k, v in indice(a.csv_ref).items()} if a.csv_ref else {}
    carpetas_nuevas = {k.split("/")[0] for k in nuevo}
    filas = []
    for k in sorted(ref):
        r = ref[k]
        if k.endswith("entradas_md5.csv") or k in OMITIR:
            estado, det = "omitido", OMITIR.get(k, "contiene rutas/huellas del entorno")
        elif k.split("/")[0] not in carpetas_nuevas:
            estado, det = "no ejecutado", "no existe la carpeta de este script"
        elif k not in nuevo:
            estado, det = "falta", ""
        else:
            b = leer(nuevo[k])
            if md5_bytes(b) == r["md5"]:
                estado, det = "idéntico", ""
            elif k.endswith(TEXTO) and r.get("md5_lf") and md5_lf(b) == r["md5_lf"]:
                estado, det = "idéntico", "salvo fin de línea CRLF/LF"
            elif k.endswith(".csv") and r.get("md5_r10") and md5_r10(b) == r["md5_r10"]:
                estado, det = "igual (redondeo 1e-10)", "diferencias solo de coma flotante"
            elif k.endswith(".csv") and k in refcsv:
                ok, det = comparar_csv(refcsv[k], nuevo[k], a.rtol)
                estado = "igual (tolerancia)" if ok else "distinto"
            elif k.endswith(BINARIOS):
                estado, det = "binario distinto", "esperable entre versiones de R/sistemas (ver ENTORNO.md)"
            else:
                estado, det = "distinto", "huellas distintas"
        filas.append(dict(fichero=k, estado=estado, detalle=det))
    resumen = {}
    for r in filas:
        resumen[r["estado"]] = resumen.get(r["estado"], 0) + 1
        print(f"{r['estado']:24} {r['fichero']}  {r['detalle']}")
    print("RESUMEN:", resumen)
    if a.informe:
        with open(a.informe, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["fichero", "estado", "detalle"], lineterminator="\n")
            w.writeheader(); w.writerows(filas)
    sys.exit(1 if any(r["estado"] in ("distinto", "falta") for r in filas) else 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--generar", action="store_true")
    ap.add_argument("--salida", default="../referencia/outputs_md5_windows.csv")
    ap.add_argument("--carpetas", default=r"^(0[1-9]|1[0-6])_", help="regex de carpetas al generar")
    ap.add_argument("--referencia", nargs="+", default=REF_DEFECTO)
    ap.add_argument("--output", default="analisis/output")
    ap.add_argument("--csv-ref", default="")
    ap.add_argument("--rtol", type=float, default=1e-8)
    ap.add_argument("--informe", default="")
    a = ap.parse_args()
    generar(a) if a.generar else comparar(a)


if __name__ == "__main__":
    main()
