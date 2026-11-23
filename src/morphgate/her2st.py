"""Read the HER2ST sections used by most histology-to-expression models."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

Image.MAX_IMAGE_PIXELS = None


@dataclass
class Section:
    name: str
    patient: str
    image: np.ndarray
    expression: pd.DataFrame
    pixel_x: np.ndarray
    pixel_y: np.ndarray


def _read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0)


def _find(root: Path, pattern: str) -> list[Path]:
    return sorted(root.rglob(pattern))


def gene_columns(root: Path, name: str) -> list[str]:
    table = pd.read_csv(_count_path(root, name), sep="\t", nrows=0, index_col=0)
    return [str(column) for column in table.columns]


def section_names(root: Path) -> list[str]:
    counts = _find(root, "A1.tsv*") + _find(root, "*.tsv") + _find(root, "*.tsv.gz")
    names = []
    for path in counts:
        stem = path.name.split(".")[0]
        if len(stem) == 2 and stem[0].isalpha() and stem[1].isdigit():
            names.append(stem)
    return sorted(set(names))


def _count_path(root: Path, name: str) -> Path:
    found = _find(root, f"{name}.tsv") + _find(root, f"{name}.tsv.gz")
    counts = [path for path in found if "selection" not in path.name and "label" not in path.name]
    if not counts:
        raise FileNotFoundError(name)
    return counts[0]


def _spot_path(root: Path, name: str) -> Path:
    found = _find(root, f"{name}_selection.tsv") + _find(root, f"{name}_selection.tsv.gz")
    if not found:
        raise FileNotFoundError(name + " spots")
    return found[0]


def _image_path(root: Path, name: str) -> Path:
    found = [path for path in _find(root, f"{name}.jpg") if path.is_file()]
    found += [path for path in _find(root, f"{name}.jpeg") if path.is_file()]
    plain = [path for path in found if path.parent.name.upper() == "HE"]
    if not plain and not found:
        raise FileNotFoundError(name + " image")
    return (plain or found)[0]


def load_section(root: Path, name: str) -> Section:
    root = Path(root)
    counts = _read_table(_count_path(root, name))
    spots = pd.read_csv(_spot_path(root, name), sep="\t")
    if "selected" in spots.columns:
        spots = spots[spots["selected"].astype(int) == 1]
    array_x = np.rint(spots["x"].to_numpy()).astype(int)
    array_y = np.rint(spots["y"].to_numpy()).astype(int)
    spots = spots.copy()
    spots["spot_id"] = [f"{x}x{y}" for x, y in zip(array_x, array_y)]
    spots = spots.drop_duplicates("spot_id").set_index("spot_id")
    shared = counts.index.intersection(spots.index)
    if len(shared) == 0:
        raise RuntimeError(f"{name}: spot ids {list(spots.index[:3])} did not match count ids {list(counts.index[:3])}")
    counts = counts.loc[shared]
    spots = spots.loc[shared]
    image = np.asarray(Image.open(_image_path(root, name)).convert("RGB"))
    return Section(
        name=name,
        patient=name[0],
        image=image,
        expression=counts.apply(pd.to_numeric, errors="coerce").fillna(0.0),
        pixel_x=spots["pixel_x"].to_numpy(dtype=float),
        pixel_y=spots["pixel_y"].to_numpy(dtype=float),
    )
