from __future__ import annotations

import math
from typing import List

import numpy as np
import pandas as pd
import statsmodels.api as sm


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


def select_block_length(
    n_obs: int,
    block_length: int | str | None = "auto",
    *,
    rule: str = "cube_root",
    min_block_length: int = 2,
    max_block_length: int | None = None,
) -> int:
    n_obs = int(n_obs)
    if n_obs <= 0:
        return int(min_block_length)

    if isinstance(block_length, str) and block_length.lower() == "auto":
        rule_name = str(rule).lower()
        if rule_name == "sqrt":
            chosen = int(math.ceil(math.sqrt(n_obs)))
        else:
            chosen = int(math.ceil(n_obs ** (1.0 / 3.0)))
    elif block_length is None:
        chosen = int(math.ceil(n_obs ** (1.0 / 3.0)))
    else:
        chosen = int(block_length)

    lower = max(2, int(min_block_length))
    upper = n_obs if max_block_length is None else min(int(max_block_length), n_obs)
    upper = max(lower, upper)
    return int(np.clip(chosen, lower, upper))


def moving_block_bootstrap(arr: np.ndarray, block_length: int, rng: np.random.Generator) -> np.ndarray:
    arr = np.asarray(arr, dtype=float)
    n = len(arr)
    if n == 0:
        return arr.copy()
    if n < 3:
        return arr.copy()
    block_length = int(np.clip(block_length, 2, n))
    starts = np.arange(0, n - block_length + 1)
    n_blocks = int(math.ceil(n / block_length))
    pieces = []
    for _ in range(n_blocks):
        s = int(rng.choice(starts))
        pieces.append(arr[s : s + block_length])
    return np.concatenate(pieces)[:n]


def moving_block_bootstrap_indices(n: int, block_length: int, rng: np.random.Generator) -> np.ndarray:
    n = int(n)
    if n <= 0:
        return np.array([], dtype=int)
    if n < 3:
        return np.arange(n, dtype=int)
    block_length = int(np.clip(block_length, 2, n))
    starts = np.arange(0, n - block_length + 1)
    n_blocks = int(math.ceil(n / block_length))
    pieces = []
    for _ in range(n_blocks):
        s = int(rng.choice(starts))
        pieces.append(np.arange(s, s + block_length, dtype=int))
    return np.concatenate(pieces)[:n]


def iid_bootstrap(arr: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    arr = np.asarray(arr, dtype=float)
    n = len(arr)
    if n == 0:
        return arr.copy()
    return arr[rng.integers(0, n, size=n)]


def maximum_entropy_bootstrap(arr: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    x = np.asarray(arr, dtype=float)
    n = len(x)
    if n < 3:
        return x.copy()

    order = np.argsort(x)
    xs = x[order]
    mids = 0.5 * (xs[:-1] + xs[1:])

    left_width = mids[0] - xs[0]
    right_width = xs[-1] - mids[-1]
    z = np.empty(n + 1, dtype=float)
    z[0] = xs[0] - left_width
    z[1:-1] = mids
    z[-1] = xs[-1] + right_width

    u = np.sort(rng.uniform(size=n))
    scaled = np.clip(u * n, 0, n - 1e-12)
    idx = np.floor(scaled).astype(int)
    frac = scaled - idx
    ys_sorted = z[idx] + frac * (z[idx + 1] - z[idx])

    y = np.empty(n, dtype=float)
    y[order] = ys_sorted
    return y - y.mean() + x.mean()


def residual_bootstrap(x: np.ndarray, y: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(y) < 3:
        return y.copy()
    X = sm.add_constant(x)
    ols = sm.OLS(y, X).fit()
    y_hat = np.asarray(ols.predict(X), dtype=float)
    resid = np.asarray(ols.resid, dtype=float)
    resid = resid - np.mean(resid)
    sampled_resid = rng.choice(resid, size=len(resid), replace=True)
    return y_hat + sampled_resid
