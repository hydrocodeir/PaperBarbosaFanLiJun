# PaperBarbosaFanLiJun

Modular Python pipeline for Q1-grade analysis of climate extreme indices using:

- Quantile Regression (full quantile grid + focus quantiles 0.05, 0.50, 0.95)
- Bootstrap uncertainty analysis
- Feature engineering from quantile slopes (Delta1/Delta2/Delta3)
- Station clustering (hierarchical or k-means)
- Publication-style figures and summary tables

## Run

```bash
python run_analysis.py --config config.yaml
```

## Structure

- `run_analysis.py`: CLI entrypoint
- `src/paper_pipeline/indices.py`: extreme-index construction
- `src/paper_pipeline/quantile.py`: quantile/OLS fitting + bootstrap
- `src/paper_pipeline/clustering.py`: feature table + clustering
- `src/paper_pipeline/plotting.py`: figure generation
- `src/paper_pipeline/reporting.py`: report generation
- `src/paper_pipeline/pipeline.py`: end-to-end orchestrator
