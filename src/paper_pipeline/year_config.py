from __future__ import annotations

from typing import Iterable

import pandas as pd


def _normalize_year_pair(years: Iterable[int] | None) -> tuple[int, int] | None:
    if years is None:
        return None
    y0, y1 = [int(x) for x in years]
    return (y0, y1) if y0 <= y1 else (y1, y0)


def get_analysis_years(cfg: dict) -> tuple[int, int] | None:
    return _normalize_year_pair(cfg.get("data", {}).get("analysis_years"))


def get_data_year_range(df: pd.DataFrame, cfg: dict) -> tuple[int, int] | None:
    if df.empty:
        return None
    year_col = cfg["data"]["year_col"]
    years = pd.to_numeric(df[year_col], errors="coerce").dropna()
    if years.empty:
        return None
    return int(years.min()), int(years.max())


def get_effective_year_range(cfg: dict, df: pd.DataFrame | None = None) -> tuple[int, int] | None:
    configured = get_analysis_years(cfg)
    if configured is not None:
        return configured
    if df is None:
        return None
    return get_data_year_range(df, cfg)


def format_year_range_label(year_range: tuple[int, int] | None, sep: str = "-") -> str:
    if year_range is None:
        return "unknown"
    return f"{year_range[0]}{sep}{year_range[1]}"


def filter_to_analysis_years(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    year_range = get_analysis_years(cfg)
    if year_range is None or df.empty:
        return df.copy()
    year_col = cfg["data"]["year_col"]
    y0, y1 = year_range
    return df.loc[pd.to_numeric(df[year_col], errors="coerce").between(y0, y1)].copy()


def resolve_reference_years(cfg: dict, df: pd.DataFrame | None = None) -> tuple[tuple[int, int] | None, str]:
    index_cfg = cfg.get("index_construction", {})
    if "reference_years" not in index_cfg:
        year_range = get_effective_year_range(cfg, df)
        if year_range is None:
            return None, "all_available_years"
        return year_range, f"{year_range[0]}_{year_range[1]}"
    ref_years = index_cfg.get("reference_years")
    if ref_years is None:
        return None, "all_available_years"
    if isinstance(ref_years, str):
        key = ref_years.strip().lower()
        if key in {"all_available", "all_available_years"}:
            return None, "all_available_years"
        if key in {"analysis_years", "analysis_period", "analysis_range"}:
            year_range = get_effective_year_range(cfg, df)
            if year_range is None:
                return None, "all_available_years"
            return year_range, f"{year_range[0]}_{year_range[1]}"
        raise ValueError(f"Unsupported reference_years setting: {ref_years}")
    year_range = _normalize_year_pair(ref_years)
    assert year_range is not None
    return year_range, f"{year_range[0]}_{year_range[1]}"


def build_split_periods(year_range: tuple[int, int] | None) -> list[tuple[str, int, int]]:
    if year_range is None:
        return []
    start_year, end_year = year_range
    n_years = end_year - start_year + 1
    if n_years <= 1:
        label = format_year_range_label(year_range)
        return [(label, start_year, end_year)]
    first_end = start_year + (n_years // 2) - 1
    second_start = first_end + 1
    first = (format_year_range_label((start_year, first_end)), start_year, first_end)
    second = (format_year_range_label((second_start, end_year)), second_start, end_year)
    return [first, second]
