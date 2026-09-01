"""Modular pipeline for quantile-regression climate-extreme analysis."""


def run_pipeline(*args, **kwargs):
    """Import the full pipeline only when execution is requested."""
    from .pipeline import run_pipeline as _run_pipeline

    return _run_pipeline(*args, **kwargs)


__all__ = ["run_pipeline"]
