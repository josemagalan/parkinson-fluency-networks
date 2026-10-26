# --- English ---
# Optional: one-off download of the pinned TabPFN v2 classifier weights and licence, with SHA-256
# check. Writes analisis/python/.models_tabpfn_v2/ (not distributed). Uses no study data.
"""E3/D-026: descargar únicamente software/pesos oficiales para inferencia local.

Entradas: repositorio público Prior-Labs/TabPFN-v2-clf de Hugging Face.
Salidas: analisis/python/.models_tabpfn_v2 (pesos originales, licencia y huellas).
No abre datos del estudio ni utiliza credenciales. Ejecutar desde revision_R3/.
Python estándar; descarga explícita separada de41 (que trabaja sin conexión).
"""
import hashlib
import json
import urllib.request
from pathlib import Path

REVISION = "f851f2a3c941544733b712d8c0f96dfae9b28862"
SHA256 = "f65a35685aeef42e31b796d9bfa34e68d6fc780bc98e7bff7763802964cf435f"


def fetch(url):
    return urllib.request.urlopen(urllib.request.Request(
        url, headers={"User-Agent": "Parkinson-R3-local-TabPFN-v2"}), timeout=120)


def main():
    out = Path("analisis/python/.models_tabpfn_v2")
    out.mkdir(parents=True, exist_ok=True)
    metadata = out / "weights_manifest.json"
    if metadata.exists():
        info = json.loads(metadata.read_text(encoding="utf-8"))
        path = out / info["filename"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == info["sha256"] == SHA256
        assert info["revision"] == REVISION
        print("Pesos existentes verificados:", info["sha256"], flush=True)
        return
    repo = "Prior-Labs/TabPFN-v2-clf"
    with fetch(f"https://huggingface.co/api/models/{repo}/revision/{REVISION}?blobs=true") as r:
        model = json.load(r)
    commit = model["sha"]
    assert commit == REVISION
    filename = "tabpfn-v2-classifier.ckpt"
    licence_url = "https://raw.githubusercontent.com/PriorLabs/TabPFN/49394b053a6759cfe68e90c21a2d51c31b396768/LICENSE"
    for name in [filename, "README.md", "LICENSE"]:
        target = out / name
        if target.exists():
            # Reanudar una descarga interrumpida solo tras verificar contra el hub.
            entry = next(x for x in model["siblings"] if x["rfilename"] == name)
            if "lfs" in entry:
                assert hashlib.sha256(target.read_bytes()).hexdigest() == entry["lfs"]["sha256"]
            else:
                with fetch(f"https://huggingface.co/{repo}/resolve/{commit}/{name}") as src:
                    assert target.read_bytes() == src.read()
            print("Reutilizado tras verificación:", name, flush=True)
            continue
        url = licence_url if name == "LICENSE" else f"https://huggingface.co/{repo}/resolve/{commit}/{name}?download=true"
        partial = out / (name + ".partial")
        with fetch(url) as src, partial.open("wb") as dst:
            while block := src.read(1 << 20):
                dst.write(block)
        partial.rename(target)
        print("Descargado", name, target.stat().st_size, "bytes", flush=True)
    path = out / filename
    assert hashlib.sha256(path.read_bytes()).hexdigest() == SHA256
    info = {"repo": repo, "revision": commit, "filename": filename,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "bytes": path.stat().st_size,
            "url": f"https://huggingface.co/{repo}/resolve/{commit}/{filename}",
            "licence_file": "LICENSE", "licence_url": licence_url,
            "inference": "local; sin datos del estudio"}
    metadata.write_text(json.dumps(info, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print("Checkpoint original v2 identificado:", info["sha256"], flush=True)


if __name__ == "__main__":
    main()
