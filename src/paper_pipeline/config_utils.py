from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from .math_utils import make_quantile_grid


class ConfigError(ValueError):
    """Raised when required analysis configuration is missing or invalid."""


def require_config_value(cfg: dict, path: str) -> Any:
    current: Any = cfg
    for key in path.split("."):
        if not isinstance(current, dict) or key not in current:
            raise ConfigError(f"Missing required config value: {path}")
        current = current[key]
    if current is None:
        raise ConfigError(f"Config value cannot be null: {path}")
    return current


def _normalize_quantile_list(values: Any, path: str) -> list[float]:
    if not isinstance(values, list) or not values:
        raise ConfigError(f"Config value must be a non-empty list: {path}")
    quantiles = sorted({round(float(x), 2) for x in values})
    for tau in quantiles:
        if not (0.0 < tau < 1.0):
            raise ConfigError(f"Invalid quantile {tau} in {path}; quantiles must be between 0 and 1.")
    return quantiles


def get_time_scale_years(cfg: dict) -> float:
    scale = float(require_config_value(cfg, "quantile_regression.time_scale_years"))
    if scale <= 0:
        raise ConfigError("quantile_regression.time_scale_years must be positive.")
    return scale


def get_time_unit_label(cfg: dict) -> str:
    return str(require_config_value(cfg, "quantile_regression.time_unit_label"))


def get_plot_dpi(cfg: dict) -> int:
    dpi = int(require_config_value(cfg, "plots.dpi"))
    if dpi <= 0:
        raise ConfigError("plots.dpi must be positive.")
    return dpi


def get_default_figure_dpi(cfg: dict) -> int:
    dpi = int(require_config_value(cfg, "plots.default_figure_dpi"))
    if dpi <= 0:
        raise ConfigError("plots.default_figure_dpi must be positive.")
    return dpi


def get_homogeneity_alpha(cfg: dict) -> float:
    alpha = float(require_config_value(cfg, "quality_control.homogeneity_alpha"))
    if not (0.0 < alpha < 1.0):
        raise ConfigError("quality_control.homogeneity_alpha must be between 0 and 1.")
    return alpha


def get_homogeneity_permutations(cfg: dict) -> int:
    n_perm = int(require_config_value(cfg, "quality_control.homogeneity_permutations"))
    if n_perm <= 0:
        raise ConfigError("quality_control.homogeneity_permutations must be positive.")
    return n_perm


def get_full_quantiles(cfg: dict) -> list[float]:
    full_cfg = require_config_value(cfg, "quantile_regression.full_quantiles")
    if not isinstance(full_cfg, dict):
        raise ConfigError("quantile_regression.full_quantiles must be a mapping.")
    for key in ("start", "stop", "step"):
        if key not in full_cfg:
            raise ConfigError(f"Missing required config value: quantile_regression.full_quantiles.{key}")
    return make_quantile_grid(float(full_cfg["start"]), float(full_cfg["stop"]), float(full_cfg["step"]))


def get_focus_quantiles(cfg: dict) -> list[float]:
    return _normalize_quantile_list(require_config_value(cfg, "quantile_regression.focus_quantiles"), "quantile_regression.focus_quantiles")


def get_sensitivity_quantiles(cfg: dict) -> list[float]:
    return _normalize_quantile_list(
        require_config_value(cfg, "quantile_regression.sensitivity_check_quantiles"),
        "quantile_regression.sensitivity_check_quantiles",
    )


def get_plot_quantiles(cfg: dict, key: str) -> list[float]:
    return _normalize_quantile_list(require_config_value(cfg, f"plots.quantile_selections.{key}"), f"plots.quantile_selections.{key}")


def format_tau_suffix(tau: float) -> str:
    return f"{float(tau):0.2f}"


def slope_col(tau: float) -> str:
    return f"slope_{format_tau_suffix(tau)}"


def ci_low_col(tau: float) -> str:
    return f"ci_low_{format_tau_suffix(tau)}"


def ci_high_col(tau: float) -> str:
    return f"ci_high_{format_tau_suffix(tau)}"


def boot_mean_col(metric: str) -> str:
    return f"boot_mean_{metric}"


def boot_sd_col(metric: str) -> str:
    return f"boot_sd_{metric}"


def boot_ci_low_col(metric: str) -> str:
    return f"boot_ci_low_{metric}"


def boot_ci_high_col(metric: str) -> str:
    return f"boot_ci_high_{metric}"


def tau_label(tau: float) -> str:
    return f"q{format_tau_suffix(tau)}"


def metric_label(metric: str) -> str:
    if metric.startswith("slope_"):
        return tau_label(float(metric.replace("slope_", "")))
    return metric


def get_delta_definitions(cfg: dict) -> dict[str, list[str]]:
    delta_defs = require_config_value(cfg, "feature_engineering.delta_definitions")
    if not isinstance(delta_defs, dict) or not delta_defs:
        raise ConfigError("feature_engineering.delta_definitions must be a non-empty mapping.")
    out: dict[str, list[str]] = {}
    for name, operands in delta_defs.items():
        if not isinstance(operands, list) or len(operands) != 2:
            raise ConfigError(f"Delta definition '{name}' must contain exactly two operand columns.")
        out[str(name)] = [str(operands[0]), str(operands[1])]
    return out


def get_primary_delta(cfg: dict) -> str:
    primary = str(require_config_value(cfg, "feature_engineering.primary_delta"))
    delta_defs = get_delta_definitions(cfg)
    if primary not in delta_defs:
        raise ConfigError(f"feature_engineering.primary_delta='{primary}' is not defined in delta_definitions.")
    return primary


def compute_defined_deltas(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    out = df.copy()
    for name, (left_col, right_col) in get_delta_definitions(cfg).items():
        left = pd.to_numeric(out.get(left_col), errors="coerce")
        right = pd.to_numeric(out.get(right_col), errors="coerce")
        out[name] = left - right
    return out


def validate_analysis_config(cfg: dict) -> None:
    focus_quantiles = get_focus_quantiles(cfg)
    sensitivity_quantiles = get_sensitivity_quantiles(cfg)
    full_quantiles = get_full_quantiles(cfg)
    get_time_scale_years(cfg)
    get_time_unit_label(cfg)
    get_plot_dpi(cfg)
    get_default_figure_dpi(cfg)
    get_homogeneity_alpha(cfg)
    get_homogeneity_permutations(cfg)
    get_primary_delta(cfg)
    for key in (
        "station_timeseries",
        "station_comparison",
        "paper2_maps",
        "bootstrap_distributions",
        "paper1_dendrograms",
        "split_period_bars",
    ):
        get_plot_quantiles(cfg, key)

    missing_focus_refs: list[str] = []
    focus_cols = {slope_col(tau) for tau in focus_quantiles}
    sensitivity_cols = {format_tau_suffix(tau) for tau in sensitivity_quantiles}
    for name, operands in get_delta_definitions(cfg).items():
        for operand in operands:
            if operand not in focus_cols:
                missing_focus_refs.append(f"{name}:{operand}")
    if missing_focus_refs:
        joined = ", ".join(missing_focus_refs)
        raise ConfigError(f"All delta operands must reference configured focus quantile slopes. Invalid operands: {joined}")

    missing_sensitivity = [tau for tau in sensitivity_quantiles if tau not in focus_quantiles]
    if missing_sensitivity:
        joined = ", ".join(format_tau_suffix(tau) for tau in missing_sensitivity)
        raise ConfigError(f"Sensitivity quantiles must be included in focus_quantiles. Missing: {joined}")

    if not set(get_plot_quantiles(cfg, "bootstrap_distributions")).issubset(set(focus_quantiles)):
        raise ConfigError("plots.quantile_selections.bootstrap_distributions must be a subset of quantile_regression.focus_quantiles.")
    if not set(get_plot_quantiles(cfg, "paper1_dendrograms")).issubset(set(focus_quantiles)):
        raise ConfigError("plots.quantile_selections.paper1_dendrograms must be a subset of quantile_regression.focus_quantiles.")
    if not set(get_plot_quantiles(cfg, "split_period_bars")).issubset(set(focus_quantiles)):
        raise ConfigError("plots.quantile_selections.split_period_bars must be a subset of quantile_regression.focus_quantiles.")

    if len(full_quantiles) == 0:
        raise ConfigError("quantile_regression.full_quantiles produced an empty grid.")
