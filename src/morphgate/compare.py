"""Wang et al. Nat Commun 2025 protocol for MorphGate and a spot-level GHIST model.

The published GHIST network is not rerun. GHIST-spot is the spot-resolution
form of its cell-program mixture: expression is a weighted sum of programs
learned on the training patients, and the weights come from the H&E features.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from sklearn.linear_model import Ridge
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

from morphgate.features import section_features
from morphgate.her2st import _count_path, _spot_path

Image.MAX_IMAGE_PIXELS = None

HER2_GENES = 785
CSCC_GENES = 997
HVG_COUNT = 30
SVG_COUNT = 20
PROGRAMS = 6
GATE = 0.10
RIDGE = 10.0
MORPH = "#1f4e79"
GHIST = "#c4704b"


def _log_cp10k(counts: np.ndarray) -> np.ndarray:
    library = np.maximum(counts.sum(axis=1, keepdims=True), 1.0)
    return np.log1p(counts / library * 1e4)


def _pearson(truth: np.ndarray, pred: np.ndarray) -> np.ndarray:
    left = truth - truth.mean(axis=0)
    right = pred - pred.mean(axis=0)
    denom = np.sqrt((left * left).sum(axis=0) * (right * right).sum(axis=0))
    scores = np.zeros(truth.shape[1], dtype=np.float64)
    ok = denom >= 1e-8
    scores[ok] = (left[:, ok] * right[:, ok]).sum(axis=0) / denom[ok]
    return scores


def _ssim(truth: np.ndarray, pred: np.ndarray) -> np.ndarray:
    low = truth.min(axis=0)
    high = truth.max(axis=0)
    span = high - low
    ok = span >= 1e-8
    scaled_truth = np.zeros_like(truth, dtype=np.float64)
    scaled_pred = np.zeros_like(pred, dtype=np.float64)
    scaled_truth[:, ok] = (truth[:, ok] - low[ok]) / span[ok]
    scaled_pred[:, ok] = (pred[:, ok] - low[ok]) / span[ok]
    mu_truth = scaled_truth.mean(axis=0)
    mu_pred = scaled_pred.mean(axis=0)
    var_truth = scaled_truth.var(axis=0)
    var_pred = scaled_pred.var(axis=0)
    cov = ((scaled_truth - mu_truth) * (scaled_pred - mu_pred)).mean(axis=0)
    c1 = 0.01 ** 2
    c2 = 0.03 ** 2
    scores = ((2 * mu_truth * mu_pred + c1) * (2 * cov + c2)) / (
        (mu_truth ** 2 + mu_pred ** 2 + c1) * (var_truth + var_pred + c2)
    )
    scores = np.where(ok, scores, 0.0)
    return scores.astype(np.float64)


def _mutual(truth: np.ndarray, pred: np.ndarray, bins: int = 20) -> np.ndarray:
    scores = np.zeros(truth.shape[1], dtype=np.float64)
    for index in range(truth.shape[1]):
        observed = truth[:, index]
        predicted = pred[:, index]
        if np.std(observed) < 1e-8 or np.std(predicted) < 1e-8:
            continue
        joint, _, _ = np.histogram2d(observed, predicted, bins=bins)
        total = joint.sum()
        if total <= 0:
            continue
        probability = joint / total
        px = probability.sum(axis=1, keepdims=True)
        py = probability.sum(axis=0, keepdims=True)
        mask = probability > 0
        scores[index] = (probability[mask] * np.log(probability[mask] / (px @ py)[mask])).sum()
    return scores


def _auc(truth: np.ndarray, pred: np.ndarray) -> np.ndarray:
    scores = np.full(truth.shape[1], np.nan, dtype=np.float64)
    for index in range(truth.shape[1]):
        binary = truth[:, index] > 0
        if binary.all() or not binary.any():
            continue
        scores[index] = float(roc_auc_score(binary, pred[:, index]))
    return scores


def _nrmse(truth: np.ndarray, pred: np.ndarray) -> np.ndarray:
    span = truth.max(axis=0) - truth.min(axis=0)
    error = np.sqrt(((pred - truth) ** 2).mean(axis=0))
    scores = np.full(truth.shape[1], np.nan, dtype=np.float64)
    ok = span >= 1e-8
    scores[ok] = error[ok] / span[ok]
    return scores


def _section_metrics(truth: np.ndarray, pred: np.ndarray) -> dict[str, np.ndarray]:
    return {
        "pcc": _pearson(truth, pred),
        "ssim": _ssim(truth, pred),
        "mi": _mutual(truth, pred),
        "auc": _auc(truth, pred),
        "nrmse": _nrmse(truth, pred),
    }


def _her2_coords(root: Path, name: str) -> tuple[np.ndarray, np.ndarray]:
    identifiers = pd.read_csv(_count_path(root, name), sep="\t", usecols=[0], dtype=str).iloc[:, 0].astype(str)
    spots = pd.read_csv(_spot_path(root, name), sep="\t")
    if "selected" in spots.columns:
        spots = spots[spots["selected"].astype(int) == 1]
    array_x = np.rint(spots["x"].to_numpy()).astype(int)
    array_y = np.rint(spots["y"].to_numpy()).astype(int)
    spots = spots.copy()
    spots["spot_id"] = [f"{x}x{y}" for x, y in zip(array_x, array_y)]
    spots = spots.drop_duplicates("spot_id").set_index("spot_id")
    shared = [spot for spot in identifiers if spot in spots.index]
    aligned = spots.loc[shared]
    return aligned["x"].to_numpy(dtype=float), aligned["y"].to_numpy(dtype=float)


def _load_her2(root: Path, cache: Path) -> dict[str, dict]:
    sections = {}
    for path in sorted(cache.glob("*.npz")):
        name = path.stem
        stored = np.load(path, allow_pickle=False)
        array_x, array_y = _her2_coords(root, name)
        features = stored["features"]
        if len(array_x) != len(features):
            raise RuntimeError(f"{name}: {len(array_x)} coordinates and {len(features)} cached spots")
        sections[name] = {
            "patient": name[0],
            "features": features,
            "expression": stored["expression"],
            "genes": stored["genes"].astype(str).tolist(),
            "x": np.rint(array_x).astype(int),
            "y": np.rint(array_y).astype(int),
        }
        print("HER2", name, features.shape[0], flush=True)
    return sections


def _cscc_samples(root: Path) -> list[tuple[str, str, Path, Path, Path]]:
    samples = []
    for counts in sorted(root.rglob("*_stdata.tsv")):
        token = counts.name.split("_", 1)[1].replace("_stdata.tsv", "")
        patient = token.split("_")[0]
        image = next(root.rglob(f"*{token}.jpg"))
        spots = next(root.rglob(f"*selection*{token}.tsv"))
        samples.append((token, patient, counts, spots, image))
    return samples


def _load_cscc(root: Path, cache: Path) -> dict[str, dict]:
    cache.mkdir(parents=True, exist_ok=True)
    sections = {}
    for token, patient, counts_path, spots_path, image_path in _cscc_samples(root):
        stored_path = cache / f"{token}.npz"
        if stored_path.is_file():
            stored = np.load(stored_path, allow_pickle=False)
            sections[token] = {
                "patient": patient,
                "features": stored["features"],
                "expression": stored["expression"],
                "genes": stored["genes"].astype(str).tolist(),
                "x": stored["x"],
                "y": stored["y"],
            }
            print("CSCC cache", token, stored["features"].shape[0], flush=True)
            continue
        counts = pd.read_csv(counts_path, sep="\t", index_col=0)
        spots = pd.read_csv(spots_path, sep="\t")
        spots = spots.copy()
        spots["spot_id"] = [f"{int(x)}x{int(y)}" for x, y in zip(spots["x"], spots["y"])]
        spots = spots.drop_duplicates("spot_id").set_index("spot_id")
        shared = counts.index.astype(str).intersection(spots.index.astype(str))
        counts = counts.loc[shared]
        spots = spots.loc[shared]
        image = np.asarray(Image.open(image_path).convert("RGB"))
        features = section_features(image, spots["pixel_x"].to_numpy(float), spots["pixel_y"].to_numpy(float))
        expression = _log_cp10k(counts.to_numpy(dtype=np.float64)).astype(np.float32)
        genes = np.array([str(gene) for gene in counts.columns])
        array_x = np.rint(spots["x"].to_numpy(dtype=float)).astype(np.int32)
        array_y = np.rint(spots["y"].to_numpy(dtype=float)).astype(np.int32)
        np.savez_compressed(
            stored_path,
            features=features.astype(np.float32),
            expression=expression,
            genes=genes,
            x=array_x,
            y=array_y,
        )
        sections[token] = {
            "patient": patient,
            "features": features.astype(np.float32),
            "expression": expression,
            "genes": genes.tolist(),
            "x": array_x,
            "y": array_y,
        }
        print("CSCC", token, expression.shape[0], expression.shape[1], flush=True)
    return sections


def _panel(sections: dict[str, dict], count: int) -> list[str]:
    shared = set.intersection(*(set(section["genes"]) for section in sections.values()))
    names = sorted(shared)
    total = 0
    first = np.zeros(len(names), dtype=np.float64)
    second = np.zeros(len(names), dtype=np.float64)
    lookup = {gene: index for index, gene in enumerate(names)}
    for section in sections.values():
        index = [lookup[gene] for gene in section["genes"] if gene in lookup]
        columns = [position for position, gene in enumerate(section["genes"]) if gene in lookup]
        values = section["expression"][:, columns].astype(np.float64)
        first[index] += values.sum(axis=0)
        second[index] += (values * values).sum(axis=0)
        total += len(values)
    variance = second / total - (first / total) ** 2
    order = np.argsort(variance)[::-1][:count]
    return [names[index] for index in order]


def _take(section: dict, genes: list[str]) -> np.ndarray:
    index = {gene: position for position, gene in enumerate(section["genes"])}
    columns = [index[gene] for gene in genes]
    return section["expression"][:, columns].astype(np.float64)


def _stack(sections: dict[str, dict], names: list[str], genes: list[str]) -> tuple[np.ndarray, np.ndarray]:
    features = [sections[name]["features"] for name in names]
    expression = [_take(sections[name], genes) for name in names]
    return np.vstack(features), np.vstack(expression)


def _kmeans(expression: np.ndarray, k: int, seed: int = 42) -> np.ndarray:
    rng = np.random.default_rng(seed)
    centers = expression[rng.choice(len(expression), size=k, replace=False)].copy()
    for _ in range(25):
        distance = ((expression[:, None, :] - centers[None, :, :]) ** 2).sum(axis=2)
        labels = distance.argmin(axis=1)
        for program in range(k):
            members = expression[labels == program]
            if len(members):
                centers[program] = members.mean(axis=0)
    return centers


def _weights(expression: np.ndarray, centers: np.ndarray) -> np.ndarray:
    distance = ((expression[:, None, :] - centers[None, :, :]) ** 2).sum(axis=2)
    nearest = np.median(distance.min(axis=1))
    scale = nearest if nearest > 1e-8 else 1.0
    shifted = -distance / scale
    shifted -= shifted.max(axis=1, keepdims=True)
    raw = np.exp(shifted)
    return raw / raw.sum(axis=1, keepdims=True)


def _predict_morphgate(train_x, train_y, test_x, patients: list[str], train_names: list[str], sections, genes):
    gate_patient = sorted(set(patients))[0]
    fit_names = [name for name in train_names if sections[name]["patient"] != gate_patient]
    gate_names = [name for name in train_names if sections[name]["patient"] == gate_patient]
    fit_x, fit_y = _stack(sections, fit_names, genes)
    gate_x, gate_y = _stack(sections, gate_names, genes)
    scaler = StandardScaler().fit(fit_x)
    model = Ridge(alpha=RIDGE).fit(scaler.transform(fit_x), fit_y)
    link = _pearson(gate_y, model.predict(scaler.transform(gate_x)))
    positive = [gene for gene, score in zip(genes, link) if score >= GATE]
    full = StandardScaler().fit(train_x)
    mean = train_y.mean(axis=0)
    pred_test_shape = (len(test_x), len(genes))
    prediction = np.repeat(mean[None, :], pred_test_shape[0], axis=0)
    if positive:
        columns = [genes.index(gene) for gene in positive]
        ridge = Ridge(alpha=RIDGE).fit(full.transform(train_x), train_y[:, columns])
        prediction[:, columns] = ridge.predict(full.transform(test_x))
    return prediction, positive


def _predict_ghist(train_x, train_y, test_x) -> np.ndarray:
    centers = _kmeans(train_y, PROGRAMS)
    target = _weights(train_y, centers)
    scaler = StandardScaler().fit(train_x)
    model = Ridge(alpha=RIDGE).fit(scaler.transform(train_x), target)
    raw = model.predict(scaler.transform(test_x))
    raw = np.maximum(raw, 0.0)
    total = raw.sum(axis=1, keepdims=True)
    total[total[:, 0] == 0] = 1.0
    empty = raw.sum(axis=1) == 0
    raw[empty] = 1.0
    total = raw.sum(axis=1, keepdims=True)
    return (raw / total) @ centers


def _predict_open(train_x, train_y, test_x) -> np.ndarray:
    scaler = StandardScaler().fit(train_x)
    model = Ridge(alpha=RIDGE).fit(scaler.transform(train_x), train_y)
    return model.predict(scaler.transform(test_x))


def _folds(sections: dict[str, dict]) -> list[list[str]]:
    patients = sorted({section["patient"] for section in sections.values()}, key=lambda patient: (len(patient), patient))
    if len(patients) == 8:
        groups = [patients[index : index + 2] for index in range(0, 8, 2)]
    else:
        groups = [[patient] for patient in patients]
    folds = []
    for group in groups:
        folds.append([name for name, section in sections.items() if section["patient"] in group])
    return folds


def _moran_rank(sections: dict[str, dict], genes: list[str]) -> list[str]:
    ranks = []
    for section in sections.values():
        expression = _take(section, genes)
        points = np.column_stack([section["x"], section["y"]]).astype(np.float64)
        distance = np.sqrt(((points[:, None, :] - points[None, :, :]) ** 2).sum(axis=2))
        np.fill_diagonal(distance, np.inf)
        neighbors = min(6, len(points) - 1)
        chosen = np.argpartition(distance, neighbors, axis=1)[:, :neighbors]
        weights = np.zeros((len(points), len(points)), dtype=np.float64)
        rows = np.repeat(np.arange(len(points)), neighbors)
        weights[rows, chosen.ravel()] = 1.0
        weights = np.maximum(weights, weights.T)
        centered = expression - expression.mean(axis=0)
        numerator = (centered * (weights @ centered)).sum(axis=0)
        denominator = (centered * centered).sum(axis=0)
        score = np.zeros(len(genes))
        ok = denominator >= 1e-8
        score[ok] = numerator[ok] / denominator[ok]
        order = np.argsort(score)[::-1]
        rank = np.empty(len(genes), dtype=int)
        rank[order] = np.arange(len(genes))
        ranks.append(rank)
    stacked = np.vstack(ranks)
    in_top = (stacked < SVG_COUNT).sum(axis=0)
    mean_rank = stacked.mean(axis=0)
    chosen = np.argsort(np.where(in_top > 0, mean_rank, mean_rank + len(genes)))[:SVG_COUNT]
    return [genes[index] for index in chosen]


def _evaluate(sections: dict[str, dict], genes: list[str], label: str) -> dict:
    folds = _folds(sections)
    records = {method: [] for method in ("morphgate", "ghist_spot", "morphgate_open")}
    elapsed = {method: 0.0 for method in records}
    spatial = {}
    for fold_index, test_names in enumerate(folds):
        train_names = [name for index, names in enumerate(folds) if index not in {fold_index, (fold_index + 1) % len(folds)} for name in names]
        print(label, "FOLD", fold_index, "train", len(train_names), "test", test_names, flush=True)
        train_patients = [sections[name]["patient"] for name in train_names]
        train_x, train_y = _stack(sections, train_names, genes)
        started = time.perf_counter()
        # Predictions are produced per test section so coordinates stay aligned.
        test_matrix = np.vstack([sections[name]["features"] for name in test_names])
        morph_rows, positive = _predict_morphgate(train_x, train_y, test_matrix, train_patients, train_names, sections, genes)
        elapsed["morphgate"] += time.perf_counter() - started
        started = time.perf_counter()
        ghist_rows = _predict_ghist(train_x, train_y, test_matrix)
        elapsed["ghist_spot"] += time.perf_counter() - started
        started = time.perf_counter()
        open_rows = _predict_open(train_x, train_y, test_matrix)
        elapsed["morphgate_open"] += time.perf_counter() - started
        cursor = 0
        blocks = {"morphgate": morph_rows, "ghist_spot": ghist_rows, "morphgate_open": open_rows}
        for name in test_names:
            stop = cursor + len(sections[name]["features"])
            truth = _take(sections[name], genes)
            for method, block in blocks.items():
                kept = positive if method == "morphgate" else genes
                metrics = _section_metrics(truth, block[cursor:stop])
                metrics.update({"section": name, "patient": sections[name]["patient"], "fold": fold_index, "positive": kept})
                records[method].append(metrics)
            if name == "B1":
                gene_index = genes.index("FASN") if "FASN" in genes else 0
                spatial = {
                    "gene": genes[gene_index],
                    "x": sections[name]["x"],
                    "y": sections[name]["y"],
                    "truth": truth[:, gene_index],
                    "morphgate": morph_rows[cursor:stop, gene_index],
                    "ghist_spot": ghist_rows[cursor:stop, gene_index],
                }
            cursor = stop
    return {"records": records, "elapsed": elapsed, "spatial": spatial, "folds": folds}


def _gene_table(records: list[dict], genes: list[str]) -> dict[str, np.ndarray]:
    """Average section metrics within a patient, then across patients and folds."""
    metrics = {}
    for key in ("pcc", "ssim", "mi", "auc", "nrmse"):
        per_gene = []
        folds = sorted({record["fold"] for record in records})
        for fold in folds:
            fold_records = [record for record in records if record["fold"] == fold]
            patients = sorted({record["patient"] for record in fold_records})
            patient_means = []
            for patient in patients:
                chosen = [record[key] for record in fold_records if record["patient"] == patient]
                patient_means.append(np.nanmean(np.vstack(chosen), axis=0))
            per_gene.append(np.nanmean(np.vstack(patient_means), axis=0))
        metrics[key] = np.nanmean(np.vstack(per_gene), axis=0)
    metrics["genes"] = np.array(genes)
    return metrics


def _mean(values: np.ndarray) -> float | None:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return None
    return round(float(finite.mean()), 4)


def _median(values: np.ndarray) -> float | None:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return None
    return round(float(np.median(finite)), 4)


def _subset_mean(table: dict, names: list[str], key: str) -> float | None:
    index = {gene: position for position, gene in enumerate(table["genes"])}
    columns = [index[gene] for gene in names if gene in index]
    if not columns:
        return None
    return _mean(table[key][columns])


def _summarize(dataset: str, genes: list[str], result: dict, hvgs: list[str], svgs: list[str]) -> dict:
    summary = {"dataset": dataset, "genes": len(genes), "hvg": hvgs, "svg": svgs, "methods": {}}
    for method, records in result["records"].items():
        table = _gene_table(records, genes)
        positive = set.intersection(*(set(record["positive"]) for record in records)) if method == "morphgate" else set(genes)
        summary["methods"][method] = {
            "seconds": round(result["elapsed"][method], 2),
            "pcc_mean": _mean(table["pcc"]),
            "pcc_median": _median(table["pcc"]),
            "ssim_mean": _mean(table["ssim"]),
            "mi_mean": _mean(table["mi"]),
            "auc_mean": _mean(table["auc"]),
            "nrmse_mean": _mean(table["nrmse"]),
            "hvg_pcc": _subset_mean(table, hvgs, "pcc"),
            "hvg_ssim": _subset_mean(table, hvgs, "ssim"),
            "svg_pcc": _subset_mean(table, svgs, "pcc"),
            "svg_ssim": _subset_mean(table, svgs, "ssim"),
            "gated_gene_count": len(positive) if method == "morphgate" else len(genes),
            "top_genes": [genes[index] for index in np.argsort(table["pcc"])[::-1][:8]],
            "top_pcc": [round(float(table["pcc"][index]), 4) for index in np.argsort(table["pcc"])[::-1][:8]],
        }
        np.savez_compressed(
            result["out"] / f"{dataset}_{method}.npz",
            genes=np.array(genes),
            **{key: table[key] for key in ("pcc", "ssim", "mi", "auc", "nrmse")},
        )
    summary["published_ghist_her2"] = {
        "source": "Fu et al. Nature Methods 2025, Figure 3, their network, not this run",
        "pcc": 0.16,
        "ssim": 0.10,
        "hvg_pcc": 0.20,
        "svg_pcc": 0.27,
        "hvg_ssim": 0.17,
        "svg_ssim": 0.26,
    }
    return summary


OPEN = "#7aa0c4"


def _one_violin(ax, position: float, values: np.ndarray, color: str, width: float = 0.28) -> None:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return
    body = ax.violinplot(finite, positions=[position], widths=width, showmedians=False, showextrema=False)
    for part in body["bodies"]:
        part.set_facecolor(color)
        part.set_alpha(0.8)
        part.set_edgecolor("white")
    quartile = np.quantile(finite, [0.25, 0.5, 0.75])
    ax.scatter([position], [quartile[1]], color="white", s=12, zorder=3, edgecolor=color)
    ax.vlines(position, quartile[0], quartile[2], color="white", linewidth=1.6, zorder=2)


def _panel_style(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _columns(genes: list[str], names: list[str]) -> list[int]:
    index = {gene: position for position, gene in enumerate(genes)}
    return [index[gene] for gene in names if gene in index]


def _draw(dataset: str, genes: list[str], result: dict, hvgs: list[str], svgs: list[str], out: Path) -> None:
    tables = {method: _gene_table(records, genes) for method, records in result["records"].items()}
    _draw_tables(dataset, genes, tables, hvgs, svgs, result.get("spatial", {}), out)


def _draw_tables(dataset: str, genes: list[str], tables: dict, hvgs: list[str], svgs: list[str], spatial: dict, out: Path) -> None:
    methods = [("MorphGate", "morphgate", MORPH), ("No gate", "morphgate_open", OPEN), ("GHIST-spot", "ghist_spot", GHIST)]
    figure, axes = plt.subplots(1, 4, figsize=(11.4, 3.6))
    for axis, key, title in zip(axes, ("pcc", "ssim", "mi", "auc"), ("PCC", "SSIM", "Mutual information", "AUC")):
        for offset, (_, method, color) in zip((-0.28, 0.0, 0.28), methods):
            _one_violin(axis, 1 + offset, tables[method][key], color, width=0.22)
        axis.set_xlim(0.4, 1.6)
        axis.set_xticks([1])
        axis.set_xticklabels([title])
        _panel_style(axis)
    handles = [plt.Line2D([0], [0], color=color, linewidth=6) for _, _, color in methods]
    figure.legend(handles, [label for label, _, _ in methods], frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    figure.tight_layout()
    figure.savefig(out / f"{dataset}_metrics.png", dpi=160, bbox_inches="tight")
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(8.8, 3.8))
    sets = [("All", genes), ("HVG", hvgs), ("SVG", svgs)]
    for axis, key, title in zip(axes, ("pcc", "ssim"), ("PCC", "SSIM")):
        for index, (_, names) in enumerate(sets):
            columns = _columns(genes, names)
            for offset, (_, method, color) in zip((-0.24, 0.0, 0.24), methods):
                _one_violin(axis, index + 1 + offset, tables[method][key][columns], color, width=0.22)
        axis.set_xticks([1, 2, 3])
        axis.set_xticklabels([name for name, _ in sets])
        axis.set_title(title)
        _panel_style(axis)
    figure.legend(handles, [label for label, _, _ in methods], frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    figure.tight_layout()
    figure.savefig(out / f"{dataset}_sets.png", dpi=160, bbox_inches="tight")
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(6.6, 3.5))
    order = np.argsort(tables["morphgate_open"]["pcc"])[::-1][:5]
    labels = [genes[index] for index in order]
    y = np.arange(len(labels))
    axis.barh(y + 0.24, [tables["morphgate"]["pcc"][index] for index in order], height=0.22, color=MORPH, label="MorphGate")
    axis.barh(y, [tables["morphgate_open"]["pcc"][index] for index in order], height=0.22, color=OPEN, label="No gate")
    axis.barh(y - 0.24, [tables["ghist_spot"]["pcc"][index] for index in order], height=0.22, color=GHIST, label="GHIST-spot")
    axis.set_yticks(y)
    axis.set_yticklabels(labels)
    axis.invert_yaxis()
    axis.set_xlabel("PCC")
    axis.set_title(f"{dataset}: five genes with the highest ungated PCC")
    axis.legend(frameon=False, loc="lower right")
    _panel_style(axis)
    figure.tight_layout()
    figure.savefig(out / f"{dataset}_top.png", dpi=160)
    plt.close(figure)

    if dataset == "HER2ST" and spatial:
        figure, axes = plt.subplots(1, 3, figsize=(9.2, 3.3))
        low = np.quantile(spatial["truth"], 0.02)
        high = np.quantile(spatial["truth"], 0.98)
        for axis, key, title in zip(axes, ("truth", "morphgate", "ghist_spot"), ("Measured", "MorphGate", "GHIST-spot")):
            points = axis.scatter(spatial["x"], -spatial["y"], c=spatial[key], s=16, cmap="viridis", vmin=low, vmax=high)
            axis.set_title(title)
            axis.set_aspect("equal")
            axis.set_xticks([])
            axis.set_yticks([])
            figure.colorbar(points, ax=axis, fraction=0.046)
        figure.suptitle(f"B1 {spatial['gene']}, test fold", fontsize=11)
        figure.tight_layout()
        figure.savefig(out / "her2_fasn.png", dpi=160)
        plt.close(figure)


def run(her2_root: Path, her2_cache: Path, cscc_root: Path, out: Path) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    her2 = _load_her2(her2_root, her2_cache)
    her2_genes = _panel(her2, HER2_GENES)
    her2_hvg = her2_genes[:HVG_COUNT]
    her2_svg = _moran_rank(her2, her2_genes)
    her2_result = _evaluate(her2, her2_genes, "HER2ST")
    her2_result["out"] = out
    cscc = _load_cscc(cscc_root, out / "cscc_cache")
    cscc_genes = _panel(cscc, CSCC_GENES)
    cscc_hvg = cscc_genes[:HVG_COUNT]
    cscc_svg = _moran_rank(cscc, cscc_genes)
    cscc_result = _evaluate(cscc, cscc_genes, "cSCC")
    cscc_result["out"] = out
    summary = {
        "protocol": "Wang et al. Nat Commun 2025, 4-fold CV, adjacent sections from one patient stay in one fold; two folds train, one validates, one tests. Metrics are computed per section and gene, averaged over a patient's sections, then over patients, then over folds.",
        "ghist_spot": "Six expression programs from k-means on the training patients. H&E features predict the program weights. This is not the published GHIST network.",
        "published_ghist": {
            "source": "Fu et al. Nature Methods 2025 Figure 3, HER2ST, their trained network",
            "pcc": 0.16,
            "ssim": 0.10,
            "hvg_pcc": 0.20,
            "hvg_ssim": 0.17,
            "svg_pcc": 0.27,
            "svg_ssim": 0.26,
        },
        "HER2ST": _summarize("HER2ST", her2_genes, her2_result, her2_hvg, her2_svg),
        "cSCC": _summarize("cSCC", cscc_genes, cscc_result, cscc_hvg, cscc_svg),
    }
    _draw("HER2ST", her2_genes, her2_result, her2_hvg, her2_svg, out)
    _draw("cSCC", cscc_genes, cscc_result, cscc_hvg, cscc_svg, out)
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("COMPARE_DONE", flush=True)
    return summary


def redraw(out: Path) -> None:
    summary = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    for dataset in ("HER2ST", "cSCC"):
        block = summary[dataset]
        tables = {}
        genes: list[str] = []
        for method in ("morphgate", "morphgate_open", "ghist_spot"):
            stored = np.load(out / f"{dataset}_{method}.npz")
            genes = stored["genes"].astype(str).tolist()
            tables[method] = {key: stored[key] for key in ("pcc", "ssim", "mi", "auc", "nrmse")}
        _draw_tables(dataset, genes, tables, block["hvg"], block["svg"], {}, out)
    print("REDRAW_DONE", flush=True)


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--her2-root", type=Path)
    parser.add_argument("--her2-cache", type=Path)
    parser.add_argument("--cscc-root", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--redraw", action="store_true")
    args = parser.parse_args()
    if args.redraw:
        redraw(args.out)
        return
    run(args.her2_root, args.her2_cache, args.cscc_root, args.out)


if __name__ == "__main__":
    main()
