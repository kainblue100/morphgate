"""Leave-one-patient fit that refuses genes the image cannot see."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

from morphgate.features import section_features
from morphgate.her2st import load_section, section_names

GENES = 250
GATE = 0.10
RIDGE = 10.0


def _log_cp10k(counts: np.ndarray) -> np.ndarray:
    library = np.maximum(counts.sum(axis=1, keepdims=True), 1.0)
    return np.log1p(counts / library * 1e4)


def _pearson(truth: np.ndarray, pred: np.ndarray) -> np.ndarray:
    scores = np.full(truth.shape[1], np.nan, dtype=np.float64)
    for index in range(truth.shape[1]):
        y = truth[:, index].astype(np.float64)
        p = pred[:, index].astype(np.float64)
        y = y - y.mean()
        p = p - p.mean()
        denom = np.sqrt((y * y).sum() * (p * p).sum())
        if denom < 1e-8:
            continue
        scores[index] = float((y * p).sum() / denom)
    return scores


def _prepare(root: Path, name: str, cache: Path) -> tuple[np.ndarray, np.ndarray, list[str]]:
    path = cache / f"{name}.npz"
    if path.is_file():
        stored = np.load(path, allow_pickle=False)
        return stored["features"], stored["expression"], stored["genes"].astype(str).tolist()
    section = load_section(root, name)
    genes = [str(gene) for gene in section.expression.columns]
    counts = section.expression.to_numpy(dtype=np.float64)
    features = section_features(section.image, section.pixel_x, section.pixel_y).astype(np.float32)
    expression = _log_cp10k(counts).astype(np.float32)
    cache.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, features=features, expression=expression, genes=np.array(genes))
    print("SECTION", name, expression.shape[0], expression.shape[1], flush=True)
    return features, expression, genes


def _stack(prepared: list[tuple[np.ndarray, np.ndarray, list[str]]], genes: list[str]) -> tuple[np.ndarray, np.ndarray]:
    features = []
    expression = []
    for feat, expr, names in prepared:
        index = {gene: position for position, gene in enumerate(names)}
        columns = [index.get(gene, -1) for gene in genes]
        block = np.zeros((len(expr), len(genes)), dtype=np.float32)
        for column, source in enumerate(columns):
            if source >= 0:
                block[:, column] = expr[:, source]
        features.append(feat)
        expression.append(block)
    return np.vstack(features), np.vstack(expression)


def _variable(expression: np.ndarray, names: list[str], count: int) -> list[str]:
    order = np.argsort(expression.var(axis=0))[::-1][:count]
    return [names[index] for index in order]


def _scale(train: np.ndarray, test: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    scaler = StandardScaler()
    return scaler.fit_transform(train), scaler.transform(test)


def _gate_scores(prepared: dict, train_names: list[str], genes: list[str]) -> tuple[np.ndarray, str]:
    """Score genes on one training patient that the fit has not seen."""
    gate_patient = sorted({name[0] for name in train_names})[0]
    fit_names = [name for name in train_names if name[0] != gate_patient]
    gate_names = [name for name in train_names if name[0] == gate_patient]
    fit_x, fit_y = _stack([prepared[name] for name in fit_names], genes)
    gate_x, gate_y = _stack([prepared[name] for name in gate_names], genes)
    fit_x, gate_x = _scale(fit_x, gate_x)
    model = Ridge(alpha=RIDGE)
    model.fit(fit_x, fit_y)
    return _pearson(gate_y, model.predict(gate_x)), gate_patient


def fit_her2st(root: Path, out: Path) -> dict:
    root = Path(root)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    names = section_names(root)
    cache = out / "cache"
    prepared = {name: _prepare(root, name, cache) for name in names}
    patients = sorted({name[0] for name in names})
    folds = []
    for patient in patients:
        test_names = [name for name in names if name[0] == patient]
        train_names = [name for name in names if name[0] != patient]
        print("FOLD", patient, "train", len(train_names), "test", len(test_names), flush=True)
        gate_patient = sorted({name[0] for name in train_names})[0]
        screen_names = [name for name in train_names if name[0] != gate_patient]
        all_genes = sorted(set.intersection(*(set(prepared[name][2]) for name in screen_names)))
        _, screen_y = _stack([prepared[name] for name in screen_names], all_genes)
        chosen = _variable(screen_y, all_genes, GENES)
        link, gate_patient = _gate_scores(prepared, train_names, chosen)
        positive = [gene for gene, score in zip(chosen, link) if np.isfinite(score) and score >= GATE]
        print("GATE", patient, "via", gate_patient, len(positive), "of", GENES, flush=True)
        train_x_all, train_y = _stack([prepared[name] for name in train_names], chosen)
        mean = train_y.mean(axis=0)
        test_x, test_y = _stack([prepared[name] for name in test_names], chosen)
        train_x, test_x = _scale(train_x_all, test_x)
        model = None
        if positive:
            columns = [chosen.index(gene) for gene in positive]
            model = Ridge(alpha=RIDGE)
            model.fit(train_x, train_y[:, columns])
        pred = np.repeat(mean[None, :], len(test_y), axis=0)
        if model is not None:
            pred[:, [chosen.index(gene) for gene in positive]] = model.predict(test_x)
        scores = _pearson(test_y, pred)
        positive_index = [chosen.index(gene) for gene in positive]
        negative_index = [index for index in range(GENES) if index not in set(positive_index)]
        def average(indexes: list[int]) -> float | None:
            if not indexes:
                return None
            values = scores[indexes].astype(np.float64)
            # A constant fill has no correlation. Score that gene as 0 so it stays in the mean.
            values = np.where(np.isfinite(values), values, 0.0)
            return round(float(values.mean()), 4)

        folds.append(
            {
                "patient": patient,
                "sections": test_names,
                "spots": int(len(test_y)),
                "morphology_positive": positive,
                "mean_r_positive": average(positive_index),
                "mean_r_withheld": average(negative_index),
                "mean_r_all_250": average(list(range(GENES))),
            }
        )
    def across(key: str) -> float | None:
        values = [fold[key] for fold in folds if fold[key] is not None]
        if not values:
            return None
        return round(float(np.mean(values)), 4)

    summary = {
        "dataset": "HER2ST",
        "accession": "10.5281/zenodo.4751624",
        "split": "leave one patient out",
        "genes_screened": GENES,
        "gate": "alphabetically first training patient held out, Pearson >= 0.10; a constant fill scores 0",
        "features": "spot stain summary plus the mean of its six nearest spots",
        "mean_r_positive": across("mean_r_positive"),
        "mean_r_withheld": across("mean_r_withheld"),
        "mean_r_all_250": across("mean_r_all_250"),
        "mean_positive_genes": round(float(np.mean([len(fold["morphology_positive"]) for fold in folds])), 1),
        "folds": folds,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("FIT_DONE", flush=True)
    return summary
