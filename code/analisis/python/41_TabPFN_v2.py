# --- English ---
# Optional TabPFN v2 add-on (E3, D-026): fixed configuration, local CUDA GPU, same participants, folds
# and bootstrap plans as scripts 15/40; no tuning, no new p-values.
# Inputs: 40 outputs, weights from 42. Outputs: analisis/output/41_TabPFN_v2[_vN]/.
"""E3/D-026: sensibilidad TabPFN v2, sin ajuste, con validación por persona.

Entradas: output/15_ML_agrupado y planes CSV exportados por40; pesos de42.
Salidas: output/41_TabPFN_v2[_vN], predicciones LOPO/CV, métricas, IC y diferencias.
Python3.12: tabpfn2.2.1, torch2.8.0/cu128, numpy, pandas, scipy, sklearn (lock).
Ejecución desde revision_R3/: python -I analisis/python/41_TabPFN_v2.py
--phase lopo/repeated/report/all; --resume <output> reanuda sin repetir folds.
No conexión a servicios, selección de variables por resultado, tuning ni nuevos p.
"""
from __future__ import annotations

import os

# Antes de importar bibliotecas: inferencia local, telemetría desactivada.
for key, value in {"TABPFN_DISABLE_TELEMETRY": "1", "DO_NOT_TRACK": "1",
                   "HF_HUB_DISABLE_TELEMETRY": "1", "HF_HUB_OFFLINE": "1",
                   "CUBLAS_WORKSPACE_CONFIG": ":4096:8", "OMP_NUM_THREADS": "4",
                   "MKL_NUM_THREADS": "4"}.items():
    os.environ[key] = value

import argparse
import dataclasses
import hashlib
import importlib.metadata
import json
import platform
import random
import sys
import time
import warnings
from pathlib import Path

RUNTIME = Path("analisis/python/.runtime_tabpfn_v2").resolve()
RUNTIME.mkdir(parents=True, exist_ok=True)
for key in ("TMP", "TEMP", "TMPDIR", "CUDA_CACHE_PATH", "XDG_CACHE_HOME"):
    os.environ[key] = str(RUNTIME)

import numpy as np
import pandas as pd
import torch
from tabpfn import TabPFNClassifier
from tabpfn.preprocessing import default_classifier_preprocessor_configs

BASE = Path("analisis/output/15_ML_agrupado")
PREP = Path("analisis/output/40_preparar_TabPFN")
MODELS = Path("analisis/python/.models_tabpfn_v2")
ALG = "tabpfn_v2"
ENDPOINTS = ("three_class", "ON_OFF")
PANELS = ("counts", "distributions", "distributions_AoA")
CLASSES = ("HC", "OFF", "ON")


def digest(path, name="md5"):
    h = hashlib.new(name)
    with Path(path).open("rb") as f:
        while block := f.read(1 << 20):
            h.update(block)
    return h.hexdigest()


def serialise(value):
    if dataclasses.is_dataclass(value):
        return dataclasses.asdict(value)
    if isinstance(value, (Path, torch.dtype, torch.device)):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)


def save_json(path, value):
    path = Path(path)
    if path.exists():
        raise FileExistsError(path)
    temporary = path.with_suffix(path.suffix + ".partial")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False,
                                    default=serialise) + "\n", encoding="utf-8")
    temporary.rename(path)


def read_csv(path):
    return pd.read_csv(path, float_precision="round_trip")


def save_csv(path, frame):
    if Path(path).exists():
        raise FileExistsError(path)
    frame.to_csv(path, index=False, encoding="utf-8", lineterminator="\n")


def prepare(X, train, test):
    """R15: constantes, mediana, media y DE muestral aprendidas solo en train."""
    a, b = X[train].copy(), X[test].copy()
    keep = np.array([len(np.unique(v[np.isfinite(v)])) > 1 for v in a.T])
    a, b = a[:, keep], b[:, keep]
    med = np.array([np.median(v[np.isfinite(v)]) for v in a.T])
    a = np.where(np.isfinite(a), a, med)
    b = np.where(np.isfinite(b), b, med)
    mu, sd = a.mean(axis=0), a.std(axis=0, ddof=1)
    sd[~np.isfinite(sd) | (sd == 0)] = 1
    return (a - mu) / sd, (b - mu) / sd, keep


def validate_preprocessing(features, sets):
    checks = []
    for endpoint in ENDPOINTS:
        dat = features if endpoint == "three_class" else features.query("condition != 'HC'")
        dat = dat.reset_index(drop=True)
        first_PD = dat.loc[dat.participant.str.startswith("PD_"), "participant"].iloc[0]
        test = np.flatnonzero(dat.participant.to_numpy() == first_PD)
        train = np.setdiff1d(np.arange(len(dat)), test)
        for panel in PANELS:
            ref = read_csv(PREP / f"reference_{endpoint}_{panel}.csv")
            a, b, keep = prepare(dat[sets[panel]].to_numpy(float), train, test)
            expected_names = list(np.array(sets[panel])[keep])
            assert list(ref.columns[2:]) == expected_names
            assert np.array_equal(ref.row, np.r_[train, test] + 1)
            gap = float(np.max(np.abs(ref.iloc[:, 2:].to_numpy() - np.vstack([a, b]))))
            assert gap < 1e-11, (endpoint, panel, gap)
            checks.append(dict(endpoint=endpoint, panel=panel, n_features=int(keep.sum()),
                               max_abs_R_Python=gap, passed=True))
    return pd.DataFrame(checks)


def scores(truth, predicted, weights=None):
    truth, predicted = np.asarray(truth), np.asarray(predicted)
    w = np.ones(len(truth)) if weights is None else np.asarray(weights)
    recalls = {g: np.sum(w[truth == g] * (predicted[truth == g] == g)) /
               np.sum(w[truth == g]) for g in sorted(np.unique(truth))}
    return {"accuracy": float(np.sum(w * (truth == predicted)) / np.sum(w)),
            "balanced_accuracy": float(np.mean(list(recalls.values()))),
            **{f"recall_{g}": float(v) for g, v in recalls.items()}}


def make_model(seed, device):
    return TabPFNClassifier(
        n_estimators=8, categorical_features_indices=None, softmax_temperature=0.9,
        balance_probabilities=False, average_before_softmax=False,
        model_path=MODELS / "tabpfn-v2-classifier.ckpt", device=device,
        ignore_pretraining_limits=False, inference_precision=torch.float32,
        fit_mode="fit_preprocessors", memory_saving_mode="auto",
        random_state=seed, n_jobs=1, inference_config=None)


def predict_cv(dat, cols, folds, endpoint, panel, key, seed, out, device):
    X = dat[cols].to_numpy(float)
    y = dat.condition.to_numpy(str)
    rows = []
    start = time.perf_counter()
    for fold in sorted(np.unique(folds)):
        file = out / "folds" / f"{key}_fold{int(fold):02d}.json"
        test, train = np.flatnonzero(folds == fold), np.flatnonzero(folds != fold)
        assert not set(dat.participant.iloc[test]) & set(dat.participant.iloc[train])
        assert set(y[train]) == set(y)
        if file.exists():
            value = json.loads(file.read_text(encoding="utf-8"))
            assert value["test_rows"] == (test + 1).tolist()
            assert value["seed"] == seed + int(fold)
            rows.extend(value["predictions"])
            continue
        xy, zz, keep = prepare(X, train, test)
        s = seed + int(fold)
        random.seed(s); np.random.seed(s); torch.manual_seed(s)
        torch.cuda.manual_seed_all(s)
        model = make_model(s, device)
        tick = time.perf_counter()
        model.fit(xy, y[train])
        prob = model.predict_proba(zz)
        torch.cuda.synchronize()
        elapsed = time.perf_counter() - tick
        assert np.isfinite(prob).all() and np.allclose(prob.sum(axis=1), 1, atol=1e-6)
        pred = model.classes_[prob.argmax(axis=1)]
        records = []
        for j, row in enumerate(test):
            p = {f"p_{c}": 0.0 for c in CLASSES}
            p.update({f"p_{c}": float(prob[j, k]) for k, c in enumerate(model.classes_)})
            records.append(dict(row=int(row + 1), participant=dat.participant.iloc[row],
                                truth=y[row], predicted=pred[j], fold=int(fold),
                                n_train=len(train), n_features=int(keep.sum()), **p,
                                endpoint=endpoint, panel=panel, algorithm=ALG))
        save_json(file, {"seed": s, "test_rows": (test + 1).tolist(),
                        "seconds": elapsed, "features": list(np.array(cols)[keep]),
                        "classes": list(model.classes_), "predictions": records})
        internal = out / "internal_config.json"
        if not internal.exists():
            save_json(internal, {"interface_config": model.interface_config_,
                                "architecture_config": model.config_,
                                "default_preprocessors": default_classifier_preprocessor_configs(),
                                "first_fit_categorical_indices": model.inferred_categorical_indices_,
                                "inference_dtype": model.forced_inference_dtype_,
                                "devices": model.devices_})
        rows.extend(records)
        del model
        if int(fold) % 10 == 0 or int(fold) == max(folds):
            print(f"{key}: fold {fold}/{max(folds)}, último {elapsed:.1f}s, "
                  f"acumulado {time.perf_counter()-start:.1f}s", flush=True)
    result = pd.DataFrame(rows).sort_values("row").reset_index(drop=True)
    assert result.row.tolist() == list(range(1, len(dat) + 1))
    assert result.participant.tolist() == dat.participant.tolist()
    assert result.truth.tolist() == dat.condition.tolist()
    return result


def bootstrap_summary(predictions, boot_plans, people):
    weights = np.array([[np.sum(ids == person) for person in people]
                        for _, ids in boot_plans.groupby("replicate").participant])
    output, draws = [], {}
    for key, p in predictions.groupby(["endpoint", "panel", "algorithm"], sort=True):
        p = p.sort_values("row")
        assert not p.row.duplicated().any()
        lookup = {person: j for j, person in enumerate(people)}
        wi = weights[:, [lookup[person] for person in p.participant]]
        point = scores(p.truth, p.predicted)
        vals = {metric: [] for metric in point}
        for w in wi:
            value = scores(p.truth, p.predicted, w)
            for metric in vals:
                vals[metric].append(value[metric])
        draws[key] = {metric: np.array(value) for metric, value in vals.items()}
        for metric, estimate in point.items():
            lo, hi = np.quantile(vals[metric], [0.025, 0.975])
            output.append(dict(endpoint=key[0], panel=key[1], algorithm=key[2],
                               metric=metric, estimate=estimate, lo=lo, hi=hi))
    return pd.DataFrame(output), draws


def report(out, features):
    lopo_path, repeated_path = out / "predictions.csv", out / "repeated.csv"
    assert lopo_path.exists() and repeated_path.exists(), "Completar LOPO y CV antes del informe"
    new = read_csv(lopo_path)
    old = read_csv(BASE / "predictions.csv")
    combined = pd.concat([old, new], ignore_index=True)
    assert len(new) == 324 and len(combined) == 972
    summary, draws = bootstrap_summary(combined, read_csv(PREP / "bootstrap_plans.csv"),
                                      features.participant.unique())
    ref = read_csv(BASE / "summary.csv")
    check = summary[summary.algorithm != ALG].merge(
        ref, on=["endpoint", "panel", "algorithm", "metric"], suffixes=("_new", "_old"))
    assert len(check) == len(ref)
    gap = max(np.max(np.abs(check[f"{c}_new"] - check[f"{c}_old"])) for c in ["estimate", "lo", "hi"])
    assert gap < 1e-12, f"Bootstrap R/Python no idéntico: {gap}"
    save_csv(out / "summary_all_algorithms.csv", summary)
    delta = []
    # Comparaciones fijas: TabPFN frente a cada algoritmo15 en cada panel;
    # distribuciones (±AoA) frente a conteos de TabPFN. Sin selección de ganador.
    comparisons = []
    for e in ENDPOINTS:
        for panel in PANELS:
            for alg in ["ridge_multinomial", "random_forest"]:
                comparisons.append(((e, panel, ALG), (e, panel, alg)))
        for panel in ["distributions", "distributions_AoA"]:
            comparisons.append(((e, panel, ALG), (e, "counts", ALG)))
    for a, b in comparisons:
        def point(key):
            return summary.loc[(summary.endpoint==key[0]) & (summary.panel==key[1]) &
                               (summary.algorithm==key[2]) & (summary.metric=="balanced_accuracy"), "estimate"].iloc[0]
        point_a, point_b = point(a), point(b)
        lo, hi = np.quantile(draws[a]["balanced_accuracy"] - draws[b]["balanced_accuracy"], [0.025, 0.975])
        delta.append(dict(endpoint=a[0], panel=a[1], algorithm=a[2], reference_panel=b[1],
                          reference_algorithm=b[2], difference=point_a-point_b, lo=lo, hi=hi))
    delta = pd.DataFrame(delta)
    save_csv(out / "paired_differences.csv", delta)
    confusion = combined.groupby(["endpoint", "panel", "algorithm", "truth", "predicted"]).size().reset_index(name="n")
    save_csv(out / "confusion_all_algorithms.csv", confusion)
    repeated = pd.concat([read_csv(BASE / "repeated.csv"), read_csv(repeated_path)], ignore_index=True)
    save_csv(out / "repeated_all_algorithms.csv", repeated)
    ranges = repeated.groupby(["endpoint", "panel", "algorithm"]).balanced_accuracy.agg(
        n="size", mean="mean", median="median", minimum="min", maximum="max").reset_index()
    save_csv(out / "repeated_summary.csv", ranges)
    checks = [dict(check="OOF records TabPFN",value=len(new),passed=len(new)==324),
              dict(check="Bootstrap métricas15 R/Python: máximo error",value=gap,passed=gap<1e-12),
              dict(check="20 CV × dos objetivos × dos paneles",value=len(read_csv(repeated_path)),passed=len(read_csv(repeated_path))==80),
              dict(check="Folds individuales calculados",value=len(list((out/"folds").glob("*.json"))),passed=len(list((out/"folds").glob("*.json")))==598)]
    save_csv(out / "verification.csv", pd.DataFrame(checks))
    lines = ["# TabPFN v2 — sensibilidad E3 (D-026)", "",
             "Análisis autorizado por Josema; protocolo fijado antes de resultados TabPFN, después de conocer15.",
             "LOPO por persona:66 evaluaciones/45 personas; ON/OFF42/21. Mismos datos, predictores y folds15.",
             "Checkpoint originalv2, biblioteca2.2.1,8 miembros, configuración fija; sin tuning ni ajuste de pesos en el estudio.",
             "Intervalos95%:2000 bootstrap estratificados PD/HC de personas sobre predicciones OOF fijas;",
             "condicionales a los modelos ajustados, no incorporan toda incertidumbre de entrenamiento.",
             "20 CV agrupadas de5 folds: variación de entrenamiento/partición; mínimo-máximo no es IC.",
             "AoA provisional (P10). No hay p nuevo, nuevo nulo ni ampliación de BH2 de15.",
             "Clasificación no estima efecto clínico ni causal de medicación. No hay validación externa.", "",
             "## Balanced accuracy LOPO", "",
             "| Objetivo | Panel | Algoritmo | BA | IC95% condicional |",
             "|---|---|---|---:|---| "]
    for r in summary[summary.metric == "balanced_accuracy"].itertuples():
        lines.append(f"| {r.endpoint} | {r.panel} | {r.algorithm} | {r.estimate:.4f} | {r.lo:.4f}–{r.hi:.4f} |")
    lines += ["", "## Diferencias emparejadas de BA", "",
              "Todas las comparaciones fijadas, sin escoger el mejor resultado; IC descriptivos sin ajuste simultáneo.", "",
              "| Objetivo | TabPFN panel | Referencia | Diferencia | IC95% condicional |",
              "|---|---|---|---:|---|"]
    for r in delta.itertuples():
        lines.append(f"| {r.endpoint} | {r.panel} | {r.reference_algorithm}/{r.reference_panel} | {r.difference:.4f} | {r.lo:.4f}–{r.hi:.4f} |")
    lines += ["", "## Trazabilidad", "",
              "config.json, internal_config.json, requirements-lock.txt, entradas_md5.csv y pesos_sha256.csv identifican la ejecución.",
              "folds/*.json guarda cada ajuste/test, semilla, columnas retenidas y tiempo; reanudación solo con idéntica configuración/entradas.",
              "preprocessing_verification.csv verifica seis referencias de R; verification.csv reproduce también métricas/IC de15.",
              "Datos y planes de15 preservados; pesos/entorno no se exportan al repositorio público. Inferencia local sin telemetría."]
    (out / "RESULTADOS.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
    print(summary.query("metric == 'balanced_accuracy' and algorithm == @ALG").to_string(index=False), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=["all", "lopo", "repeated", "report"], default="all")
    parser.add_argument("--resume", type=Path)
    args = parser.parse_args()
    if not BASE.exists() or not PREP.exists():
        sys.exit("Ejecutar desde revision_R3/ después de15 y40.")
    assert importlib.metadata.version("tabpfn") == "2.2.1"
    assert torch.cuda.is_available(), "Protocolo de esta ejecución: GPU CUDA local"
    torch.set_num_threads(4)
    torch.set_float32_matmul_precision("highest")
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
    torch.use_deterministic_algorithms(True)
    device = "cuda"
    weights = json.loads((MODELS / "weights_manifest.json").read_text(encoding="utf-8"))
    assert digest(MODELS / weights["filename"], "sha256") == weights["sha256"]
    features = read_csv(BASE / "features.csv")
    panel_table = read_csv(PREP / "panels.csv")
    sets = {panel: panel_table[panel_table.panel==panel].sort_values("position").feature.tolist() for panel in PANELS}
    assert [len(sets[p]) for p in PANELS] == [3,107,131]
    assert len(features)==66 and features.participant.nunique()==45
    assert features.condition.value_counts().to_dict()=={"HC":24,"OFF":21,"ON":21}
    files = [BASE/f for f in ["features.csv", "predictions.csv", "summary.csv", "repeated.csv", "planes.rds"]]
    files += sorted(PREP.glob("*.csv"))
    files += [Path(__file__), Path("analisis/python/42_descargar_pesos_TabPFN.py"),
              Path("analisis/R/40_preparar_TabPFN.R"),
              Path("notas/propuesta_TabPFN_v2_ML_20261022.md")]
    inputs = [{"path":str(p).replace("\\","/"),"md5":digest(p)} for p in files]
    packages = {d.metadata["Name"]: d.version for d in importlib.metadata.distributions()}
    config = {"python":sys.version,"platform":platform.platform(),"packages":packages,
              "weights":weights,"device":torch.cuda.get_device_name(),"cuda":torch.version.cuda,
              "inputs":inputs,"n_estimators":8,"softmax_temperature":0.9,
              "balance_probabilities":False,"average_before_softmax":False,
              "fit_mode":"fit_preprocessors","inference_precision":"float32",
              "memory_saving_mode":"auto","categorical_features_indices":None,"n_jobs":1,
              "inference_config":None,"torch_threads":4,"lopo_seed":"7321+fold",
              "repeated_seed":"5000+10*repeat+fold","deterministic_algorithms":True,
              "bootstrap":"2000 same people plans15; conditional OOF",
              "permutation":"not extended; no new p", "telemetry":"disabled; offline"}
    if args.resume:
        out = args.resume
        assert json.loads((out/"config.json").read_text(encoding="utf-8")) == config, "Configuración/entradas cambiadas"
    else:
        if args.phase == "report":
            sys.exit("--phase report requiere --resume")
        out = Path("analisis/output/41_TabPFN_v2")
        version = 2
        while out.exists():
            out = Path(f"analisis/output/41_TabPFN_v2_v{version}"); version += 1
        (out/"folds").mkdir(parents=True)
        save_json(out/"config.json",config)
        save_csv(out/"entradas_md5.csv",pd.DataFrame(inputs))
        save_csv(out/"pesos_sha256.csv",pd.DataFrame([weights]))
        save_csv(out/"preprocessing_verification.csv",validate_preprocessing(features,sets))
        (out/"requirements-lock.txt").write_text("\n".join(f"{k}=={v}" for k,v in sorted(packages.items()))+"\n",encoding="utf-8")
    print("Salida:",out, "GPU:",torch.cuda.get_device_name(),flush=True)
    folds_table = read_csv(PREP/"folds.csv")
    if args.phase in ["all","lopo"] and not (out/"predictions.csv").exists():
        results = []
        for e in ENDPOINTS:
            dat = (features if e=="three_class" else features.query("condition != 'HC'")).reset_index(drop=True)
            ft = folds_table[folds_table.endpoint==e].sort_values("row")
            assert ft.participant.tolist()==dat.participant.tolist() and ft.truth.tolist()==dat.condition.tolist()
            for panel in PANELS:
                results.append(predict_cv(dat,sets[panel],ft.fold.to_numpy(),e,panel,
                                          f"LOPO_{e}_{panel}",7321,out,device))
        save_csv(out/"predictions.csv",pd.concat(results,ignore_index=True))
    if args.phase in ["all","repeated"] and not (out/"repeated.csv").exists():
        plans = read_csv(PREP/"repeated_plans.csv")
        metrics, predictions = [], []
        for r in range(1,21):
            plan = plans[plans.repeat_id==r].set_index("participant").fold
            for e in ENDPOINTS:
                dat = (features if e=="three_class" else features.query("condition != 'HC'")).reset_index(drop=True)
                fold = plan.loc[dat.participant].to_numpy()
                for panel in ["counts","distributions"]:
                    p = predict_cv(dat,sets[panel],fold,e,panel,f"CV{r:02d}_{e}_{panel}",
                                   5000+r*10,out,device)
                    p["repeat_id"] = r
                    predictions.append(p)
                    metrics.append(dict(repeat_id=r,endpoint=e,panel=panel,algorithm=ALG,
                                        balanced_accuracy=scores(p.truth,p.predicted)["balanced_accuracy"]))
            print(f"CV repetida {r}/20 completada",flush=True)
        save_csv(out/"repeated_predictions.csv",pd.concat(predictions,ignore_index=True))
        save_csv(out/"repeated.csv",pd.DataFrame(metrics))
    if args.phase in ["all","report"] and not (out/"RESULTADOS.md").exists():
        report(out,features)
    print("Fase completada:",args.phase,"Salida:",out,flush=True)


if __name__ == "__main__":
    main()
