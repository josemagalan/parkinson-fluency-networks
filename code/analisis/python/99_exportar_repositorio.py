# -*- coding: utf-8 -*-
# --- English ---
# Maintainer tool: copies the analysis code into repositorio_reproducible/ and regenerates
# MANIFEST.csv and datos_entrada_md5.csv. Never copies data or outputs.
"""
99_exportar_repositorio.py

Propósito
    Sincronizar el código de los análisis nuevos (R3) con la carpeta
    `repositorio_reproducible/`, que es el borrador de un repositorio público y
    reproducible. La fuente de verdad sigue siendo `analisis/`: este script solo
    COPIA código (nunca datos ni outputs) y regenera dos inventarios.

Entradas (desde revision_R3/)
    analisis/R/*.R, analisis/python/*.py, analisis/python/*.sh,
    analisis/python/vendor/**                (código; se excluye analisis/R/library)
    analisis/output/*/entradas_md5.csv       (para listar los datos de entrada usados)
    ../parkinson/networks/*.graphml           (entradas del script 17)

Salidas
    repositorio_reproducible/code/analisis/...          copia exacta del código
    (en el repositorio público la carpeta de trabajo revision_R3/ se llama code/)
    repositorio_reproducible/MANIFEST.csv               ruta, bytes y md5 del código copiado
    repositorio_reproducible/datos_entrada_md5.csv      datos que hay que aportar, con md5
    (README.md, DATOS.md, ENTORNO.md, scripts de ejecución, etc. se editan a mano
     en repositorio_reproducible/ y este script no los toca)

Uso
    python -I analisis/python/99_exportar_repositorio.py [--comprobar]
    --comprobar: no copia nada; sale con código 1 si el repositorio está desfasado.

Solo usa la biblioteca estándar de Python (>= 3.8).
"""
import argparse
import csv
import glob
import hashlib
import os
import shutil
import sys

DEST = "repositorio_reproducible"
SUB = "code"  # nombre público de revision_R3/ (24-10-2026)
CODIGO = ["analisis/R/*.R", "analisis/python/*.py", "analisis/python/*.sh",
          "analisis/python/vendor/**/*"]
EXCLUIR = ("__pycache__", ".pyc")
# Código que no se publica: lee correspondencia privada (no forma parte de la cadena de análisis)
NO_PUBLICAR = ("analisis/python/00d_extraer_correo_elena.py",)
# Ficheros de entrada que no se distribuyen (correo privado; solo se registra su huella)
PRIVADOS = ("../RE_ Major revision requested_ TLCD-2025-0015_elena.msg",)


def md5(ruta):
    h = hashlib.md5()
    with open(ruta, "rb") as f:
        for bloque in iter(lambda: f.read(1 << 20), b""):
            h.update(bloque)
    return h.hexdigest()


def ficheros_codigo():
    out = []
    for patron in CODIGO:
        for f in glob.glob(patron, recursive=True):
            f = f.replace("\\", "/")
            if os.path.isfile(f) and not any(x in f for x in EXCLUIR) and f not in NO_PUBLICAR:
                out.append(f)
    return sorted(set(out))


def entradas():
    rutas = set()
    for f in glob.glob("analisis/output/*/entradas_md5.csv"):
        with open(f, encoding="utf-8") as h:
            for fila in csv.DictReader(h):
                p = fila.get("path") or fila.get("fichero") or ""
                p = p.replace("\\", "/")
                if p.startswith("parkinson/"):
                    p = "../" + p
                if p.startswith("../") or p.startswith("historial/"):
                    rutas.add(p)
    rutas.update(f.replace("\\", "/") for f in glob.glob("../parkinson/networks/graph*.graphml")
                 if not f.endswith("Partition.graphml"))
    filas = []
    for p in sorted(rutas):
        existe = os.path.isfile(p)
        # El nombre del correo privado no se publica: solo se indica que existe y su huella
        filas.append(dict(ruta_desde_code="(correo privado; no se distribuye)" if p in PRIVADOS else p,
                          distribuible="no" if p in PRIVADOS else "por decidir (datos clínicos/derechos)",
                          bytes=os.path.getsize(p) if existe else "",
                          md5=md5(p) if existe else "NO ENCONTRADO"))
    return filas


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--comprobar", action="store_true")
    a = ap.parse_args()
    if not os.path.isdir("analisis") or not os.path.isfile("AGENTS.md"):
        sys.exit("Ejecutar desde revision_R3/")
    cod = ficheros_codigo()
    desfase = []
    manifest = []
    for f in cod:
        destino = os.path.join(DEST, SUB, f)
        h = md5(f)
        manifest.append(dict(ruta=f"{SUB}/{f}", bytes=os.path.getsize(f), md5=h))
        if not os.path.isfile(destino) or md5(destino) != h:
            desfase.append(f)
            if not a.comprobar:
                os.makedirs(os.path.dirname(destino), exist_ok=True)
                shutil.copy2(f, destino)
    # Código en el repositorio que ya no existe en analisis/
    sobrantes = []
    base = os.path.join(DEST, SUB)
    for f in glob.glob(os.path.join(base, "analisis", "**", "*"), recursive=True):
        rel = os.path.relpath(f, base).replace("\\", "/")
        if os.path.isfile(f) and rel not in cod:
            sobrantes.append(rel)
    if a.comprobar:
        print(f"desfasados: {len(desfase)}; sobrantes: {len(sobrantes)}")
        for f in desfase + sobrantes:
            print("  ", f)
        sys.exit(1 if (desfase or sobrantes) else 0)
    with open(os.path.join(DEST, "MANIFEST.csv"), "w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=["ruta", "bytes", "md5"], lineterminator="\n")
        w.writeheader(); w.writerows(manifest)
    filas = entradas()
    with open(os.path.join(DEST, "datos_entrada_md5.csv"), "w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=list(filas[0].keys()), lineterminator="\n")
        w.writeheader(); w.writerows(filas)
    print(f"copiados/actualizados: {len(desfase)} de {len(cod)} ficheros de código")
    print(f"entradas inventariadas: {len(filas)}")
    if sobrantes:
        print("AVISO: ficheros en el repositorio que ya no están en analisis/ (no se borran):")
        for f in sobrantes:
            print("  ", f)


if __name__ == "__main__":
    main()
