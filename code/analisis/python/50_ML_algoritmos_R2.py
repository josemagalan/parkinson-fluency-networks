"""E3: reevalúa con validación por persona los clasificadores de la versión anterior (XGBoost,
rotation forest y super learner), con los mismos datos, paneles, folds y bootstrap que 15/40/41.

Entradas: output/15_ML_agrupado (features, predictions, summary, repeated), output/40_preparar_TabPFN
(panels, folds, bootstrap_plans, repeated_plans), output/41_TabPFN_v2 (predictions, repeated).
Salidas: output/50_ML_algoritmos_R2[_vN]: predicciones LOPO y de 20 CV agrupadas, métricas con
IC bootstrap de las seis familias de modelos, diferencias emparejadas frente a la regresión logística.
Python >= 3.10: numpy, pandas, scipy, scikit-learn, xgboost (versiones en entorno.txt).
Ejecución desde revision_R3/: python -I analisis/python/50_ML_algoritmos_R2.py
Configuraciones fijas antes de ver resultados; sin tuning ni selección por rendimiento.
--- English ---
Purpose (E3): evaluate the classifiers of the previous version (XGBoost, rotation forest and super learner)
with participant-grouped validation, using exactly the data, predictor panels, folds and bootstrap plans of
scripts 15/40/41. Fixed configurations, no tuning. Outputs: analisis/output/50_ML_algoritmos_R2[_vN]/.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import os
import platform
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from xgboost import XGBClassifier

BASE = Path("analisis/output/15_ML_agrupado")
PREP = Path("analisis/output/40_preparar_TabPFN")
TABPFN = Path("analisis/output/41_TabPFN_v2")
ENDPOINTS = ("three_class", "ON_OFF")
PANELS = ("counts", "distributions", "distributions_AoA")
CLASSES = ("HC", "OFF", "ON")
NEW = ("xgboost", "rotation_forest", "super_learner")
SEED = 7321  # misma semilla base por fold que 15 y 41


def digest(path):
    h = hashlib.md5()
    with Path(path).open("rb") as f:
        while block := f.read(1 << 20):
            h.update(block)
    return h.hexdigest()


def read_csv(path):
    return pd.read_csv(path, float_precision="round_trip")


def prepare(X, train, test):
    """Igual que 15/41: constantes, mediana, media y DE muestral aprendidas solo en train."""
    a, b = X[train].copy(), X[test].copy()
    keep = np.array([len(np.unique(v[np.isfinite(v)])) > 1 for v in a.T])
    a, b = a[:, keep], b[:, keep]
    med = np.array([np.median(v[np.isfinite(v)]) for v in a.T])
    a = np.where(np.isfinite(a), a, med)
    b = np.where(np.isfinite(b), b, med)
    mu, sd = a.mean(axis=0), a.std(axis=0, ddof=1)
    sd[~np.isfinite(sd) | (sd == 0)] = 1
    return (a - mu) / sd, (b - mu) / sd, keep


def scores(truth, predicted, weights=None):
    truth, predicted = np.asarray(truth), np.asarray(predicted)
    w = np.ones(len(truth)) if weights is None else np.asarray(weights)
    recalls = {g: np.sum(w[truth == g] * (predicted[truth == g] == g)) / np.sum(w[truth == g])
               for g in sorted(np.unique(truth))}
    return {"accuracy": float(np.sum(w * (truth == predicted)) / np.sum(w)),
            "balanced_accuracy": float(np.mean(list(recalls.values()))),
            **{f"recall_{g}": float(v) for g, v in recalls.items()}}


# ------------------------------------------------------------------ models (fixed configurations)
def fit_xgboost(x, y, seed):
    k = len(np.unique(y))
    m = XGBClassifier(n_estimators=100, learning_rate=0.1, max_depth=3, subsample=0.8,
                      colsample_bytree=0.8, min_child_weight=1, reg_lambda=1.0, tree_method="exact",
                      objective="multi:softprob" if k > 2 else "binary:logistic",
                      n_jobs=1, random_state=seed, verbosity=0)
    m.fit(x, y)
    return m


class RotationForest:
    """Rotation forest (Rodríguez, Kuncheva & Alonso, 2006), con los valores por defecto de Weka:
    10 árboles, grupos de 3 variables, PCA sobre una submuestra bootstrap del 75 % de las
    observaciones de un subconjunto aleatorio no vacío de clases; árbol de decisión (entropía,
    hoja mínima 2) sobre los datos rotados. Probabilidades: media de los árboles."""

    def __init__(self, n_trees=10, group_size=3, seed=0):
        self.n_trees, self.group_size, self.rng = n_trees, group_size, np.random.default_rng(seed)

    def fit(self, x, y):
        self.classes_ = np.unique(y)
        p = x.shape[1]
        self.members = []
        for _ in range(self.n_trees):
            perm = self.rng.permutation(p)
            groups = [perm[i:i + self.group_size] for i in range(0, p, self.group_size)]
            rot = np.zeros((p, p))
            for g in groups:
                cls = self.classes_[self.rng.random(len(self.classes_)) < 0.5]
                if len(cls) == 0:
                    cls = self.classes_[[self.rng.integers(len(self.classes_))]]
                idx = np.flatnonzero(np.isin(y, cls))
                size = max(2, int(round(0.75 * len(idx))))
                sample = idx[self.rng.integers(len(idx), size=size)]
                xs = x[np.ix_(sample, g)]
                if np.all(xs.std(axis=0) == 0):
                    comp = np.eye(len(g))
                else:
                    pca = PCA(n_components=len(g), svd_solver="full").fit(xs)
                    comp = pca.components_.T
                    if comp.shape[1] < len(g):
                        comp = np.eye(len(g))
                rot[np.ix_(g, g)] = comp
            tree = DecisionTreeClassifier(criterion="entropy", min_samples_leaf=2,
                                          random_state=int(self.rng.integers(2**31 - 1)))
            tree.fit(x @ rot, y)
            self.members.append((rot, tree))
        return self

    def predict_proba(self, z):
        out = np.zeros((z.shape[0], len(self.classes_)))
        for rot, tree in self.members:
            pr = tree.predict_proba(z @ rot)
            out[:, np.searchsorted(self.classes_, tree.classes_)] += pr
        return out / len(self.members)


class XGBWrap:
    """XGBoost con etiquetas de texto (codificadas internamente)."""

    def __init__(self, seed):
        self.seed = seed

    def fit(self, x, y):
        self.classes_ = np.unique(y)
        self.m = fit_xgboost(x, np.searchsorted(self.classes_, y), self.seed)
        return self

    def predict_proba(self, z):
        return self.m.predict_proba(z)


def base_learners(seed):
    """Biblioteca del super learner: regresión logística L2, random forest, XGBoost y media de clases."""
    return {
        "logistic_L2": lambda x, y: LogisticRegression(C=1.0, max_iter=5000).fit(x, y),
        "random_forest": lambda x, y: RandomForestClassifier(n_estimators=200, max_features="sqrt",
                                                             min_samples_leaf=3, random_state=seed,
                                                             n_jobs=1).fit(x, y),
        "xgboost": lambda x, y: XGBWrap(seed).fit(x, y),
    }


def proba(model, z, classes):
    pr = model.predict_proba(z)
    out = np.zeros((z.shape[0], len(classes)))
    out[:, np.searchsorted(classes, model.classes_)] = pr
    return out


def inner_folds(groups, cohort, k, rng):
    """Folds internos por persona, estratificados PD/HC como en 15."""
    people = pd.unique(groups)
    fold = {}
    for c in (True, False):
        ppl = [p for p in people if cohort[p] == c]
        lab = rng.permutation(np.resize(np.arange(k), len(ppl)))
        fold.update(dict(zip(ppl, lab)))
    return np.array([fold[g] for g in groups])


class SuperLearner:
    """Super learner (van der Laan, Polley & Hubbard, 2007): combinación convexa de las probabilidades
    de los modelos base y de la proporción de clases, con pesos que minimizan la log-verosimilitud
    negativa en una validación cruzada interna de 5 folds agrupada por persona dentro del train."""

    def __init__(self, seed=0, k=5):
        self.seed, self.k = seed, k

    def fit(self, x, y, groups):
        self.classes_ = np.unique(y)
        rng = np.random.default_rng(self.seed)
        cohort = {g: g.startswith("PD_") for g in groups}
        f = inner_folds(groups, cohort, self.k, rng)
        libs = base_learners(self.seed)
        names = list(libs) + ["class_mean"]
        Z = np.zeros((len(y), len(self.classes_), len(names)))
        for j in range(self.k):
            te, tr = f == j, f != j
            if len(np.unique(y[tr])) < len(self.classes_) or not te.any():
                continue
            for m, (nm, fit) in enumerate(libs.items()):
                Z[te, :, m] = proba(fit(x[tr], y[tr]), x[te], self.classes_)
            prior = np.array([np.mean(y[tr] == c) for c in self.classes_])
            Z[te, :, -1] = prior
        Y = (y[:, None] == self.classes_[None, :]).astype(float)

        def loss(theta):
            w = np.exp(theta - theta.max()); w /= w.sum()
            p = np.clip(np.tensordot(Z, w, axes=([2], [0])), 1e-12, 1)
            return -np.mean(np.sum(Y * np.log(p), axis=1))
        res = minimize(loss, np.zeros(len(names)), method="Nelder-Mead",
                       options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 20000})
        w = np.exp(res.x - res.x.max()); self.weights = w / w.sum()
        self.names = names
        self.models = {nm: fit(x, y) for nm, fit in libs.items()}
        self.prior = np.array([np.mean(y == c) for c in self.classes_])
        return self

    def predict_proba(self, z):
        ps = [proba(self.models[nm], z, self.classes_) for nm in self.names[:-1]]
        ps.append(np.tile(self.prior, (z.shape[0], 1)))
        return np.tensordot(np.stack(ps, axis=2), self.weights, axes=([2], [0]))


def fit_predict(algorithm, x, y, z, seed, groups):
    if algorithm == "xgboost":
        classes = np.unique(y)
        enc = np.searchsorted(classes, y)
        m = fit_xgboost(x, enc, seed)
        pr = m.predict_proba(z)
        return classes, pr, None
    if algorithm == "rotation_forest":
        m = RotationForest(seed=seed).fit(x, y)
        return m.classes_, m.predict_proba(z), None
    if algorithm == "super_learner":
        m = SuperLearner(seed=seed).fit(x, y, groups)
        return m.classes_, m.predict_proba(z), dict(zip(m.names, m.weights))
    raise ValueError(algorithm)


def predict_cv(dat, cols, folds, endpoint, panel, algorithm, seed):
    X = dat[cols].to_numpy(float)
    y = dat.condition.to_numpy(str)
    groups = dat.participant.to_numpy(str)
    rows, weights = [], []
    for fold in sorted(np.unique(folds)):
        test, train = np.flatnonzero(folds == fold), np.flatnonzero(folds != fold)
        assert not set(groups[test]) & set(groups[train])
        assert set(y[train]) == set(y)
        a, b, keep = prepare(X, train, test)
        s = seed + int(fold)
        classes, prob, w = fit_predict(algorithm, a, y[train], b, s, groups[train])
        assert np.isfinite(prob).all() and np.allclose(prob.sum(axis=1), 1, atol=1e-6)
        pred = classes[prob.argmax(axis=1)]
        if w is not None:
            weights.append(dict(endpoint=endpoint, panel=panel, fold=int(fold), **w))
        for j, row in enumerate(test):
            p = {f"p_{c}": 0.0 for c in CLASSES}
            p.update({f"p_{c}": float(prob[j, k]) for k, c in enumerate(classes)})
            rows.append(dict(row=int(row + 1), participant=groups[row], truth=y[row], predicted=pred[j],
                             fold=int(fold), n_train=len(train), n_features=int(keep.sum()), **p,
                             endpoint=endpoint, panel=panel, algorithm=algorithm))
    out = pd.DataFrame(rows).sort_values("row").reset_index(drop=True)
    assert out.participant.tolist() == dat.participant.tolist()
    return out, weights


def bootstrap_summary(predictions, boot_plans, people):
    weights = np.array([[np.sum(ids == person) for person in people]
                        for _, ids in boot_plans.groupby("replicate").participant])
    lookup = {person: j for j, person in enumerate(people)}
    output, draws = [], {}
    for key, p in predictions.groupby(["endpoint", "panel", "algorithm"], sort=True):
        p = p.sort_values("row")
        wi = weights[:, [lookup[person] for person in p.participant]]
        point = scores(p.truth, p.predicted)
        vals = {metric: [] for metric in point}
        for w in wi:
            value = scores(p.truth, p.predicted, w)
            for metric in vals:
                vals[metric].append(value[metric])
        draws[key] = {metric: np.array(v) for metric, v in vals.items()}
        for metric, estimate in point.items():
            lo, hi = np.quantile(vals[metric], [0.025, 0.975])
            output.append(dict(endpoint=key[0], panel=key[1], algorithm=key[2], metric=metric,
                               estimate=estimate, lo=lo, hi=hi))
    return pd.DataFrame(output), draws


def main():
    if not BASE.exists() or not PREP.exists() or not TABPFN.exists():
        sys.exit("Ejecutar desde revision_R3/ después de 15, 40 y 41.")
    features = read_csv(BASE / "features.csv")
    pt = read_csv(PREP / "panels.csv")
    sets = {p: pt[pt.panel == p].sort_values("position").feature.tolist() for p in PANELS}
    assert [len(sets[p]) for p in PANELS] == [3, 107, 131]
    assert len(features) == 66 and features.participant.nunique() == 45
    folds_table = read_csv(PREP / "folds.csv")
    out = Path("analisis/output/50_ML_algoritmos_R2")
    v = 2
    while out.exists():
        out = Path(f"analisis/output/50_ML_algoritmos_R2_v{v}"); v += 1
    out.mkdir(parents=True)
    inputs = [BASE / f for f in ("features.csv", "predictions.csv", "summary.csv", "repeated.csv")]
    inputs += [PREP / f for f in ("panels.csv", "folds.csv", "bootstrap_plans.csv", "repeated_plans.csv")]
    inputs += [TABPFN / f for f in ("predictions.csv", "repeated.csv")] + [Path(__file__)]
    pd.DataFrame([dict(path=str(p).replace("\\", "/"), md5=digest(p)) for p in inputs]).to_csv(
        out / "entradas_md5.csv", index=False)
    pk = ["numpy", "pandas", "scipy", "scikit-learn", "xgboost"]
    (out / "entorno.txt").write_text(f"python {sys.version}\n{platform.platform()}\n" + "\n".join(
        f"{p} {importlib.metadata.version(p)}" for p in pk) + "\n", encoding="utf-8")

    preds, sl_weights = [], []
    for e in ENDPOINTS:
        dat = (features if e == "three_class" else features.query("condition != 'HC'")).reset_index(drop=True)
        ft = folds_table[folds_table.endpoint == e].sort_values("row")
        assert ft.participant.tolist() == dat.participant.tolist()
        for panel in PANELS:
            for alg in NEW:
                p, w = predict_cv(dat, sets[panel], ft.fold.to_numpy(), e, panel, alg, SEED)
                preds.append(p); sl_weights += w
                print(f"LOPO {e} {panel} {alg}: BA {scores(p.truth, p.predicted)['balanced_accuracy']:.3f}", flush=True)
    new = pd.concat(preds, ignore_index=True)
    new.to_csv(out / "predictions.csv", index=False)
    pd.DataFrame(sl_weights).to_csv(out / "super_learner_weights.csv", index=False)

    plans = read_csv(PREP / "repeated_plans.csv")
    rep = []
    for r in range(1, 21):
        plan = plans[plans.repeat_id == r].set_index("participant").fold
        for e in ENDPOINTS:
            dat = (features if e == "three_class" else features.query("condition != 'HC'")).reset_index(drop=True)
            fold = plan.loc[dat.participant].to_numpy()
            for panel in ("counts", "distributions"):
                for alg in NEW:
                    p, _ = predict_cv(dat, sets[panel], fold, e, panel, alg, 5000 + r * 10)
                    rep.append(dict(repeat_id=r, endpoint=e, panel=panel, algorithm=alg,
                                    balanced_accuracy=scores(p.truth, p.predicted)["balanced_accuracy"]))
        print(f"CV repetida {r}/20", flush=True)
    rep = pd.DataFrame(rep)
    rep.to_csv(out / "repeated.csv", index=False)

    combined = pd.concat([read_csv(BASE / "predictions.csv"), read_csv(TABPFN / "predictions.csv"), new],
                         ignore_index=True)
    summary, draws = bootstrap_summary(combined, read_csv(PREP / "bootstrap_plans.csv"),
                                      features.participant.unique())
    ref = read_csv(BASE / "summary.csv")
    chk = summary.merge(ref, on=["endpoint", "panel", "algorithm", "metric"], suffixes=("_new", "_old"))
    gap = max(np.max(np.abs(chk[f"{c}_new"] - chk[f"{c}_old"])) for c in ("estimate", "lo", "hi"))
    assert len(chk) == len(ref) and gap < 1e-12, gap
    summary.to_csv(out / "summary_all_algorithms.csv", index=False)
    delta = []
    for e in ENDPOINTS:
        for panel in PANELS:
            for alg in ("random_forest", "tabpfn_v2") + NEW:
                a, b = (e, panel, alg), (e, panel, "ridge_multinomial")
                pa = summary.query("endpoint==@e and panel==@panel and algorithm==@alg and metric=='balanced_accuracy'").estimate.iloc[0]
                pb = summary.query("endpoint==@e and panel==@panel and algorithm=='ridge_multinomial' and metric=='balanced_accuracy'").estimate.iloc[0]
                lo, hi = np.quantile(draws[a]["balanced_accuracy"] - draws[b]["balanced_accuracy"], [.025, .975])
                delta.append(dict(endpoint=e, panel=panel, algorithm=alg, reference="ridge_multinomial",
                                  difference=pa - pb, lo=lo, hi=hi))
    pd.DataFrame(delta).to_csv(out / "paired_differences.csv", index=False)
    allrep = pd.concat([read_csv(BASE / "repeated.csv"), read_csv(TABPFN / "repeated.csv"), rep], ignore_index=True)
    allrep.groupby(["endpoint", "panel", "algorithm"]).balanced_accuracy.agg(
        n="size", mean="mean", median="median", minimum="min", maximum="max").reset_index().to_csv(
        out / "repeated_summary.csv", index=False)
    combined.groupby(["endpoint", "panel", "algorithm", "truth", "predicted"]).size().reset_index(
        name="n").to_csv(out / "confusion_all_algorithms.csv", index=False)
    pd.DataFrame([dict(check="Bootstrap de 15 reproducido (máx. error)", value=gap, passed=gap < 1e-12),
                  dict(check="Predicciones LOPO nuevas", value=len(new), passed=len(new) == 3 * 3 * (66 + 42)),
                  dict(check="CV repetidas nuevas", value=len(rep), passed=len(rep) == 20 * 2 * 2 * 3)]).to_csv(
        out / "verification.csv", index=False)
    lines = ["# 50 Clasificadores de la versión anterior con validación por persona (E3, D-034)", "",
             "Clasificadores de la versión anterior del análisis evaluados con el protocolo de 15.",
             "Mismos datos, paneles (3/107/131), folds LOPO, 2000 planes bootstrap y 20 CV agrupadas que 15/40/41.",
             "XGBoost: 100 árboles, eta 0,1, profundidad 3, submuestreo 0,8 de filas y columnas, sin tuning.",
             "Rotation forest (Rodríguez et al., 2006): implementación propia con los valores por defecto de Weka (10 árboles,",
             "  grupos de 3 variables, PCA sobre bootstrap del 75 % de un subconjunto aleatorio de clases; árbol con entropía).",
             "Super learner (van der Laan et al., 2007): pesos convexos que minimizan la log-verosimilitud en CV interna de",
             "  5 folds por persona dentro de cada train; biblioteca: logística L2, random forest, XGBoost y proporción de clases.",
             "La regresión logística es la multinomial con penalización L2 de 15: sin penalización no es estimable con 107 predictores.",
             "summary_all_algorithms.csv reúne los seis modelos (15, 41 y 50) con los mismos planes bootstrap; 15 se reproduce exactamente.",
             "Intervalos condicionados a las predicciones OOF; no hay nuevos p ni permutaciones para estos modelos.",
             "AoA provisional. Sin validación externa; no estima efecto clínico ni causal.", ""]
    ba = summary[summary.metric == "balanced_accuracy"]
    for r in ba.itertuples():
        lines.append(f"- {r.endpoint} | {r.panel} | {r.algorithm}: {r.estimate:.4f} [{r.lo:.4f}, {r.hi:.4f}]")
    (out / "LEEME.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(ba.to_string(index=False))
    print("Salida:", out)


if __name__ == "__main__":
    main()
