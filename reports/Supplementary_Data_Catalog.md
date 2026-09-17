# Supplementary data catalog

This catalog distinguishes numerical evidence from display exports. Files remain in their canonical locations so scripts can reuse them without duplicate copies. The supplement embeds eight selected tables and eleven figures; the CSV files below provide complete station-level and sensitivity evidence. Raw observations remain in `data/` and are not redistributed by this catalog.

Verification means reproducibility from the stated inputs, not proof that observations are error-free or homogenized. See [the audit](Output_Audit_2026.md) for bootstrap scope and limitations.

## Current publication evidence

| Dataset | Rows × columns | Role and interpretation |
| --- | --- | --- |
| [bootstrap_network_annual.csv](../outputs/publication_v2/tables/bootstrap_network_annual.csv) | 4,999 × 4 | Reproducibility: synchronized bootstrap draws for primary/group intervals. |
| [bootstrap_network_warm_season.csv](../outputs/publication_v2/tables/bootstrap_network_warm_season.csv) | 4,999 × 4 | Reproducibility: synchronized bootstrap draws for primary/group intervals. |
| [climate_regime_compound_bootstrap.csv](../outputs/publication_v2/tables/climate_regime_compound_bootstrap.csv) | 59,988 × 7 | Reproducibility: synchronized bootstrap draws for primary/group intervals. |
| [climate_regime_compound_partition.csv](../outputs/publication_v2/tables/climate_regime_compound_partition.csv) | 48 × 10 | Main/S8 evidence: climate-regime summaries or four-component partitions. |
| [compound_fixed_threshold_extent.csv](../outputs/publication_v2/tables/compound_fixed_threshold_extent.csv) | 68 × 6 | Main/S1–S2 evidence: fixed-threshold event definitions, station partitions or sensitivity scenarios. |
| [compound_partition_all_scenario_stations.csv](../outputs/publication_v2/tables/compound_partition_all_scenario_stations.csv) | 2,017 × 24 | Main/S1–S2 evidence: fixed-threshold event definitions, station partitions or sensitivity scenarios. |
| [compound_partition_primary.csv](../outputs/publication_v2/tables/compound_partition_primary.csv) | 8 × 20 | Main/S1–S2 evidence: fixed-threshold event definitions, station partitions or sensitivity scenarios. |
| [compound_partition_sensitivity.csv](../outputs/publication_v2/tables/compound_partition_sensitivity.csv) | 80 × 20 | Main/S1–S2 evidence: fixed-threshold event definitions, station partitions or sensitivity scenarios. |
| [compound_partition_stations.csv](../outputs/publication_v2/tables/compound_partition_stations.csv) | 207 × 24 | Main/S1–S2 evidence: fixed-threshold event definitions, station partitions or sensitivity scenarios. |
| [network_quantile_profiles_recomputed.csv](../outputs/publication_v2/tables/network_quantile_profiles_recomputed.csv) | 324 × 4 | Main evidence: independently recomputed network thermal profiles/trends. |
| [screened_annual_seasonal_aggregates.csv](../outputs/publication_v2/tables/screened_annual_seasonal_aggregates.csv) | 8,432 × 10 | Supplementary data: screened annual/seasonal aggregates and calendar coverage (S1). |
| [table01_thermal_trends.csv](../outputs/publication_v2/tables/table01_thermal_trends.csv) | 4 × 7 | Main evidence: independently recomputed network thermal profiles/trends. |
| [table02_climate_regime_thermal.csv](../outputs/publication_v2/tables/table02_climate_regime_thermal.csv) | 6 × 8 | Main/S8 evidence: climate-regime summaries or four-component partitions. |

## Historical thermal and diagnostic evidence

| Dataset | Rows × columns | Role and interpretation |
| --- | --- | --- |
| [alternative_clustering_assignments.csv](../outputs/tables/alternative_clustering_assignments.csv) | 496 × 8 | Exploratory regionalization (S6–S8): features, assignments, stability or representative selection. |
| [alternative_clustering_sensitivity_summary.csv](../outputs/tables/alternative_clustering_sensitivity_summary.csv) | 16 × 8 | Exploratory regionalization (S6–S8): features, assignments, stability or representative selection. |
| [annual_extreme_indices.csv](../outputs/tables/annual_extreme_indices.csv) | 4,216 × 12 | Reproducibility: full-record thermal annual station indices rebuilt from daily data. |
| [bootstrap_depth_sensitivity_station_comparison.csv](../outputs/tables/bootstrap_depth_sensitivity_station_comparison.csv) | 496 × 56 | Sensitivity (S2): saved alternative station estimates/aggregation checked; alternative random ensembles not regenerated. |
| [bootstrap_depth_sensitivity_summary.csv](../outputs/tables/bootstrap_depth_sensitivity_summary.csv) | 32 × 9 | Sensitivity (S2): saved alternative station estimates/aggregation checked; alternative random ensembles not regenerated. |
| [bootstrap_distributions_long.csv](../outputs/tables/bootstrap_distributions_long.csv) | 99,200 × 10 | Reproducibility: all 99,200 saved station bootstrap draws; all summary columns recalculated. |
| [bootstrap_method_sensitivity_station_level.csv](../outputs/tables/bootstrap_method_sensitivity_station_level.csv) | 496 × 17 | Sensitivity (S2): saved alternative station estimates/aggregation checked; alternative random ensembles not regenerated. |
| [bootstrap_method_sensitivity_summary.csv](../outputs/tables/bootstrap_method_sensitivity_summary.csv) | 24 × 9 | Sensitivity (S2): saved alternative station estimates/aggregation checked; alternative random ensembles not regenerated. |
| [climate_regime_difference_tests.csv](../outputs/tables/climate_regime_difference_tests.csv) | 32 × 13 | Climate classification and group diagnostics; descriptive partitions with multiplicity/size limitations. |
| [climate_regime_emergence_summary.csv](../outputs/tables/climate_regime_emergence_summary.csv) | 96 × 10 | Precision diagnostics (S4); legacy filename does not denote an emergence date or calibrated detection. |
| [climate_regime_fixed_baseline_summary.csv](../outputs/tables/climate_regime_fixed_baseline_summary.csv) | 24 × 8 | Climate classification and group diagnostics; descriptive partitions with multiplicity/size limitations. |
| [climate_regime_quantile_summary.csv](../outputs/tables/climate_regime_quantile_summary.csv) | 24 × 20 | Climate classification and group diagnostics; descriptive partitions with multiplicity/size limitations. |
| [climate_regime_warming_response_summary.csv](../outputs/tables/climate_regime_warming_response_summary.csv) | 24 × 15 | Exploratory association (S9); shared observations/trends, not external-forcing attribution. Station fits were reproduced with the historical iteration limit; no independent station-level convergence certification. |
| [climate_signal_emergence_station_level.csv](../outputs/tables/climate_signal_emergence_station_level.csv) | 1,984 × 10 | Precision diagnostics (S4); legacy filename does not denote an emergence date or calibrated detection. |
| [climate_signal_emergence_summary.csv](../outputs/tables/climate_signal_emergence_summary.csv) | 16 × 8 | Precision diagnostics (S4); legacy filename does not denote an emergence date or calibrated detection. |
| [cluster_assignments.csv](../outputs/tables/cluster_assignments.csv) | 496 × 4 | Exploratory regionalization (S6–S8): features, assignments, stability or representative selection. |
| [cluster_assignments_reduced_features.csv](../outputs/tables/cluster_assignments_reduced_features.csv) | 496 × 4 | Exploratory regionalization (S6–S8): features, assignments, stability or representative selection. |
| [cluster_robustness_summary.csv](../outputs/tables/cluster_robustness_summary.csv) | 4 × 4 | Exploratory regionalization (S6–S8): features, assignments, stability or representative selection. |
| [clustering_feature_screening_summary.csv](../outputs/tables/clustering_feature_screening_summary.csv) | 52 × 7 | Exploratory regionalization (S6–S8): features, assignments, stability or representative selection. |
| [clustering_feature_table.csv](../outputs/tables/clustering_feature_table.csv) | 496 × 59 | Exploratory regionalization (S6–S8): features, assignments, stability or representative selection. |
| [data_homogeneity_tests_station_summary.csv](../outputs/tables/data_homogeneity_tests_station_summary.csv) | 124 × 25 | Quality/sensitivity (S1–S3); diagnostics do not constitute station homogenization. |
| [data_quality_homogeneity_overview.csv](../outputs/tables/data_quality_homogeneity_overview.csv) | 9 × 2 | Quality/sensitivity (S1–S3); diagnostics do not constitute station homogenization. |
| [data_quality_station_summary.csv](../outputs/tables/data_quality_station_summary.csv) | 124 × 10 | Quality/sensitivity (S1–S3); diagnostics do not constitute station homogenization. |
| [driver_analysis_summary.csv](../outputs/tables/driver_analysis_summary.csv) | 48 × 9 | Geographical association (S11); legacy driver terminology does not imply causality. |
| [fixed_baseline_annual_extreme_indices.csv](../outputs/tables/fixed_baseline_annual_extreme_indices.csv) | 4,216 × 12 | Main fixed-baseline sensitivity: annual indices and station/network trends or period contrasts. |
| [fixed_baseline_period_change_station_level.csv](../outputs/tables/fixed_baseline_period_change_station_level.csv) | 496 × 10 | Main fixed-baseline sensitivity: annual indices and station/network trends or period contrasts. |
| [fixed_baseline_period_change_summary.csv](../outputs/tables/fixed_baseline_period_change_summary.csv) | 4 × 6 | Main fixed-baseline sensitivity: annual indices and station/network trends or period contrasts. |
| [fixed_baseline_qr_summary.csv](../outputs/tables/fixed_baseline_qr_summary.csv) | 496 × 19 | Main fixed-baseline sensitivity: annual indices and station/network trends or period contrasts. |
| [homogeneity_flag_exclusion_sensitivity.csv](../outputs/tables/homogeneity_flag_exclusion_sensitivity.csv) | 4 × 22 | Quality/sensitivity (S1–S3); diagnostics do not constitute station homogenization. |
| [interpolation_method_sensitivity_summary.csv](../outputs/tables/interpolation_method_sensitivity_summary.csv) | 24 × 7 | Reproducibility only: interpolation comparison; no interpolated surface is used as inferential evidence. |
| [koppen_geiger_regime_summary.csv](../outputs/tables/koppen_geiger_regime_summary.csv) | 6 × 10 | Climate classification and group diagnostics; descriptive partitions with multiplicity/size limitations. |
| [koppen_geiger_station_assignments.csv](../outputs/tables/koppen_geiger_station_assignments.csv) | 124 × 14 | Climate classification and group diagnostics; descriptive partitions with multiplicity/size limitations. |
| [publication_summary_table.csv](../outputs/tables/publication_summary_table.csv) | 496 × 43 | Thermal evidence (main/S3–S4/S8): station quantiles, slopes, intervals and feature summaries. |
| [qr_all_quantiles_long.csv](../outputs/tables/qr_all_quantiles_long.csv) | 40,176 × 6 | Thermal evidence (main/S3–S4/S8): station quantiles, slopes, intervals and feature summaries. |
| [qr_focus_slopes_and_bootstrap_summary.csv](../outputs/tables/qr_focus_slopes_and_bootstrap_summary.csv) | 496 × 57 | Thermal evidence (main/S3–S4/S8): station quantiles, slopes, intervals and feature summaries. |
| [regional_cluster_composites.csv](../outputs/tables/regional_cluster_composites.csv) | 64 × 9 | Exploratory regionalization (S6–S8): features, assignments, stability or representative selection. |
| [regional_cluster_spatial_validation.csv](../outputs/tables/regional_cluster_spatial_validation.csv) | 4 × 8 | Exploratory spatial diagnostics (S5/S7); nominal permutation probabilities. |
| [regional_temperature_anomaly.csv](../outputs/tables/regional_temperature_anomaly.csv) | 34 × 4 | Exploratory association (S9); shared observations/trends, not external-forcing attribution. Station fits were reproduced with the historical iteration limit; no independent station-level convergence certification. |
| [regional_temperature_anomaly_summary.csv](../outputs/tables/regional_temperature_anomaly_summary.csv) | 1 × 7 | Exploratory association (S9); shared observations/trends, not external-forcing attribution. Station fits were reproduced with the historical iteration limit; no independent station-level convergence certification. |
| [representative_station_selection.csv](../outputs/tables/representative_station_selection.csv) | 16 × 11 | Exploratory regionalization (S6–S8): features, assignments, stability or representative selection. |
| [spatial_autocorrelation_moran.csv](../outputs/tables/spatial_autocorrelation_moran.csv) | 12 × 5 | Exploratory spatial diagnostics (S5/S7); nominal permutation probabilities. |
| [station_significance_fdr.csv](../outputs/tables/station_significance_fdr.csv) | 1,488 × 12 | Spatial diagnostics (S5): local analytic FDR; missing tail tests remain NA, not zero. |
| [temperature_internal_consistency_compound_sensitivity.csv](../outputs/tables/temperature_internal_consistency_compound_sensitivity.csv) | 6 × 18 | Quality/sensitivity (S1–S3); diagnostics do not constitute station homogenization. |
| [temperature_internal_consistency_quantile_sensitivity.csv](../outputs/tables/temperature_internal_consistency_quantile_sensitivity.csv) | 4 × 19 | Quality/sensitivity (S1–S3); diagnostics do not constitute station homogenization. |
| [temperature_internal_consistency_screening.csv](../outputs/tables/temperature_internal_consistency_screening.csv) | 2 × 6 | Quality/sensitivity (S1–S3); diagnostics do not constitute station homogenization. |
| [warming_link_network_quantile_response.csv](../outputs/tables/warming_link_network_quantile_response.csv) | 16 × 7 | Exploratory association (S9); shared observations/trends, not external-forcing attribution. Station fits were reproduced with the historical iteration limit; no independent station-level convergence certification. |
| [warming_link_station_quantile_response.csv](../outputs/tables/warming_link_station_quantile_response.csv) | 496 × 17 | Exploratory association (S9); shared observations/trends, not external-forcing attribution. Station fits were reproduced with the historical iteration limit; no independent station-level convergence certification. |

## Historical event-definition sensitivity

| Dataset | Rows × columns | Role and interpretation |
| --- | --- | --- |
| [compound_dry_hot_distribution_shift_tests.csv](../outputs/compound_dry_hot/tables/compound_dry_hot_distribution_shift_tests.csv) | 6 × 13 | Historical definition sensitivity (S10); varying network and empirical joint rarity, not fixed AND events or causal drivers. |
| [compound_dry_hot_driver_summary.csv](../outputs/compound_dry_hot/tables/compound_dry_hot_driver_summary.csv) | 12 × 10 | Historical definition sensitivity (S10); varying network and empirical joint rarity, not fixed AND events or causal drivers. |
| [compound_dry_hot_moran_trend_summary.csv](../outputs/compound_dry_hot/tables/compound_dry_hot_moran_trend_summary.csv) | 6 × 9 | Historical definition sensitivity (S10); varying network and empirical joint rarity, not fixed AND events or causal drivers. |
| [compound_dry_hot_moran_yearly.csv](../outputs/compound_dry_hot/tables/compound_dry_hot_moran_yearly.csv) | 204 × 7 | Historical definition sensitivity (S10); varying network and empirical joint rarity, not fixed AND events or causal drivers. |
| [compound_dry_hot_serial_dependence_sensitivity.csv](../outputs/compound_dry_hot/tables/compound_dry_hot_serial_dependence_sensitivity.csv) | 6 × 16 | Historical definition sensitivity (S10); varying network and empirical joint rarity, not fixed AND events or causal drivers. |
| [compound_dry_hot_station_frequency.csv](../outputs/compound_dry_hot/tables/compound_dry_hot_station_frequency.csv) | 738 × 20 | Historical definition sensitivity (S10); varying network and empirical joint rarity, not fixed AND events or causal drivers. |
| [compound_dry_hot_station_year.csv](../outputs/compound_dry_hot/tables/compound_dry_hot_station_year.csv) | 8,309 × 26 | Historical definition sensitivity (S10); varying network and empirical joint rarity, not fixed AND events or causal drivers. |
| [compound_dry_hot_trend_summary.csv](../outputs/compound_dry_hot/tables/compound_dry_hot_trend_summary.csv) | 6 × 17 | Historical definition sensitivity (S10); varying network and empirical joint rarity, not fixed AND events or causal drivers. |
| [compound_dry_hot_yearly_extent.csv](../outputs/compound_dry_hot/tables/compound_dry_hot_yearly_extent.csv) | 204 × 7 | Historical definition sensitivity (S10); varying network and empirical joint rarity, not fixed AND events or causal drivers. |

## Provenance and validation

- [Numerical audit](../outputs/audit_cleanup/numerical_checks.json), [dependency checks](../check_output_dependencies.py), and [independent network warming solver check](../outputs/audit_cleanup/network_warming_solver_check.csv).
- [Raw-index and compound validation](../outputs/publication_v2/validation.json).
- [Source hashes for supplementary figures](../outputs/publication_v2/supplementary_figure_sources.json).
- [Original run metadata](../outputs/run_metadata.json) and [publication configuration](../publication_config.yaml).
- [Reference audit](Reference_Audit_2026.md) and [curated reference metadata](verified_references.json).
- [File-by-file curation decisions](../outputs/audit_cleanup/curation_manifest.csv) and [cleanup summary](../outputs/audit_cleanup/cleanup_summary.json).

Archived historical reports may contain obsolete paths and claims. The current manuscript, supplement and this catalog define the active publication package. Running the full legacy pipeline can recreate superseded exports; rerun curation review before treating those exports as publication evidence.
