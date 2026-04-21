from __future__ import annotations

import math
from typing import List

import numpy as np
import pandas as pd


def make_quantile_grid(start: float, stop: float, step: float) -> List[float]:
    vals = []
    x = start
    while x <= stop + 1e-9:
        vals.append(round(float(x), 2))
        x += step
    return vals


def doy_noleap(dt: pd.Series) -> pd.Series:
    doy = dt.dt.dayofyear.astype(int)
    leap = dt.dt.is_leap_year
    after_feb28 = (dt.dt.month > 2) | ((dt.dt.month == 2) & (dt.dt.day == 29))
    return doy - ((leap & after_feb28).astype(int))


def circular_day_distance(days: np.ndarray, center: int, max_day: int = 365) -> np.ndarray:
    raw = np.abs(days - center)
    return np.minimum(raw, max_day - raw)


def moving_block_bootstrap(arr: np.ndarray, block_length: int, rng: np.random.Generator) -> np.ndarray:
    arr = np.asarray(arr, dtype=float)
    n = len(arr)
    if n == 0:
        return arr.copy()
    block_length = max(1, min(int(block_length), n))
    n_blocks = int(math.ceil(n / block_length))
    starts = rng.integers(0, n, size=n_blocks)
    offsets = np.arange(block_length)
    indices = (starts[:, None] + offsets[None, :]) % n
    sampled = arr[indices].reshape(-1)
    return sampled[:n]


def iid_bootstrap(arr: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    arr = np.asarray(arr, dtype=float)
    n = len(arr)
    if n == 0:
        return arr.copy()
    return arr[rng.integers(0, n, size=n)]
