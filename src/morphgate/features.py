"""Stain features for one spot and for its immediate neighbors."""

from __future__ import annotations

import numpy as np

# Ruifrok and Johnston, optical-density stain matrix, rows are H, E, DAB.
_STAINS = np.array(
    [
        [0.650, 0.704, 0.286],
        [0.072, 0.990, 0.105],
        [0.268, 0.570, 0.776],
    ],
    dtype=np.float64,
)
_MIX = np.linalg.pinv(_STAINS.T)
HALF = 112


def _crop(image: np.ndarray, x: float, y: float) -> np.ndarray:
    height, width = image.shape[:2]
    cx = int(np.clip(round(x), 0, width - 1))
    cy = int(np.clip(round(y), 0, height - 1))
    patch = image[max(0, cy - HALF) : cy + HALF, max(0, cx - HALF) : cx + HALF]
    if patch.size == 0:
        patch = image[cy : cy + 1, cx : cx + 1]
    return patch


def _describe(patch: np.ndarray) -> np.ndarray:
    pixels = patch.reshape(-1, 3).astype(np.float64)
    rgb_mean = pixels.mean(axis=0) / 255.0
    rgb_std = pixels.std(axis=0) / 255.0
    optical = -np.log10((pixels + 1.0) / 256.0)
    stain = optical @ _MIX
    hematoxylin = stain[:, 0]
    eosin = stain[:, 1]
    gray = pixels.mean(axis=1)
    # A one-pixel shift approximates how much the tissue edge moves inside the patch.
    shifted = gray[1:] if gray.size > 1 else gray
    gradient = np.abs(gray[: shifted.size] - shifted).mean() / 255.0
    return np.concatenate(
        [
            rgb_mean,
            rgb_std,
            [hematoxylin.mean(), hematoxylin.std(), eosin.mean(), eosin.std(), gradient],
        ]
    )


def section_features(image: np.ndarray, pixel_x: np.ndarray, pixel_y: np.ndarray, neighbors: int = 6) -> np.ndarray:
    local = np.vstack([_describe(_crop(image, x, y)) for x, y in zip(pixel_x, pixel_y)])
    if len(local) == 1:
        return np.hstack([local, local])
    points = np.column_stack([pixel_x, pixel_y])
    distance = np.sqrt(((points[:, None, :] - points[None, :, :]) ** 2).sum(axis=2))
    np.fill_diagonal(distance, np.inf)
    chosen = np.argpartition(distance, min(neighbors, len(local) - 1), axis=1)[:, :neighbors]
    context = local[chosen].mean(axis=1)
    return np.hstack([local, context])
