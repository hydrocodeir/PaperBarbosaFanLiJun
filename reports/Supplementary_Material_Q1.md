# Supplementary material

## S1. Scope and reproducibility

This supplementary material documents the station-level estimates, diagnostic checks, sensitivity analyses, climate-regime comparisons, clustering diagnostics, and compound-event results supporting the main manuscript. All numerical tables are generated from the analysis outputs; no supplementary value is intended to replace the primary results in the manuscript.

The analysis uses 124 stations and 34 annual observations (1991–2024). Station counts are network counts, not area-weighted estimates of affected land. Empirical return periods are within-record rarity classes and should not be interpreted as stationary design return periods.

## S2. Supporting tables

The following machine-readable tables are available under `../outputs/tables/` and `../outputs/compound_dry_hot/tables/`:

- `annual_extreme_indices.csv`: annual warm- and cool-index counts.
- `qr_focus_slopes_and_bootstrap_summary.csv`: station-level quantile slopes, bootstrap intervals, and significance fields.
- `qr_all_quantiles_long.csv`: the full 0.10–0.90 quantile grid.
- `station_significance_fdr.csv`: Benjamini–Hochberg results at the median quantile.
- `spatial_autocorrelation_moran.csv`: Moran's *I* diagnostics for focal quantiles.
- `koppen_geiger_station_assignments.csv` and `koppen_geiger_regime_summary.csv`: climate-regime assignments and summaries.
- `homogeneity_flag_exclusion_sensitivity.csv`: exclusion sensitivity for detrended homogeneity flags.
- `bootstrap_depth_sensitivity_summary.csv` and `bootstrap_method_sensitivity_summary.csv`: uncertainty-method checks.
- `temperature_internal_consistency_screening.csv`, `temperature_internal_consistency_quantile_sensitivity.csv`, and `temperature_internal_consistency_compound_sensitivity.csv`: internal temperature-consistency sensitivity.
- `compound_dry_hot_station_year.csv`, `compound_dry_hot_trend_summary.csv`, `compound_dry_hot_distribution_shift.csv`, and `compound_dry_hot_serial_dependence_sensitivity.csv`: compound-event diagnostics and dependence-aware inference.

Additional tables in the output directory provide exploratory clustering, driver-association, signal-emergence, and fixed-baseline summaries. These outputs are descriptive unless a corresponding inferential procedure is explicitly reported in the main text.

The exploratory climate-fingerprint composite is retained in the output directory for auditability but is not used as headline evidence. Its components are dependent and its equal weighting has no external calibration; readers should therefore interpret the underlying trend, uncertainty, field-significance, and robustness components separately.

## S3. Supporting figures

The publication-oriented figures used in the manuscript are linked from `../outputs/figures/`:

1. `ijoc_study_area_regional_context_new.png` — station network and regional context.
2. `ijoc_regional_quantile_panels.png` — network-mean quantile slope profiles.
3. `ijoc_main_delta1_maps.png` — station-level tail-asymmetry maps.
4. `advanced_climate_regimes/climate_regime_quantile_profiles.png` — climate-regime profiles.
5. `advanced_climate_change_signal/fixed_baseline_period_change_summary.png` — fixed-threshold early–late changes.
6. `ijoc_robustness_synthesis.png` — robustness diagnostics.
7. `ijoc_split_period_comparison.png` — split-period quantile comparison.
8. `compound_dry_hot/figures/compound_dry_hot_extent_timeseries.png` — compound-event extent.
9. `compound_dry_hot/figures/compound_dry_hot_driver_shift.png` — marginal component dominance.
10. `compound_dry_hot/figures/compound_dry_hot_station_frequency_maps.png` — station-level compound-event frequency.
11. `compound_dry_hot/figures/compound_dry_hot_moran_connectedness.png` — spatial connectedness diagnostics.

The data-completeness and homogeneity figure (`ijoc_data_quality_homogeneity.png`) is retained as a diagnostic supplementary figure.

Diagnostic figures include data coverage, homogeneity, bootstrap depth and method sensitivity, alternative clustering, spatial inference, station comparisons, and compound-event connectedness. They should be used to audit the reported results, not as additional confirmatory tests.

## S4. Interpretation limits

The supplementary outputs do not overcome the 34-year record length, the lack of metadata-supported daily homogenization, or uneven station density. The Köppen–Geiger comparison is an association with present-day climate classes, and compound-event connectedness is a property of the station network rather than a gridded estimate of physical area. Exploratory cluster assignments are retained for transparency but are not treated as fixed natural regions.
