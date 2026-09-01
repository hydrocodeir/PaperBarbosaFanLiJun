# Supplementary material

## Quantile-dependent changes in thermal extremes and compound dry–hot conditions across Iran, 1991–2024

This file indexes the derived evidence supporting `Manuscript_Q1_Revised.md`. The linked CSV files are machine-readable and retain greater numerical precision than the manuscript. Before journal submission, the items should be renumbered and packaged according to the selected journal's supplementary-material policy.

## S1. Data quality and homogeneity

- [Table S1. Data-quality and homogeneity overview](../outputs/tables/data_quality_homogeneity_overview.csv)
- [Table S2. Station-level homogeneity diagnostics](../outputs/tables/data_homogeneity_tests_station_summary.csv)
- [Table S3. Sensitivity to excluding detrended homogeneity flags](../outputs/tables/homogeneity_flag_exclusion_sensitivity.csv)
- [Table S4. Internal temperature-consistency screening](../outputs/tables/temperature_internal_consistency_screening.csv)
- [Table S5. Quantile sensitivity after internal-consistency masking](../outputs/tables/temperature_internal_consistency_quantile_sensitivity.csv)
- [Table S6. Compound-event sensitivity after internal-consistency masking](../outputs/tables/temperature_internal_consistency_compound_sensitivity.csv)

The internal-consistency screen is non-destructive: it masks contradictory daily temperature combinations only in the sensitivity calculation and does not alter the source archive or primary outputs.

## S2. Quantile-regression and bootstrap results

- [Table S7. Focal station and regional quantile slopes with bootstrap summaries](../outputs/tables/qr_focus_slopes_and_bootstrap_summary.csv)
- [Table S8. Bootstrap-depth comparison](../outputs/tables/bootstrap_depth_sensitivity_summary.csv)
- [Table S9. Station-level bootstrap-depth comparison](../outputs/tables/bootstrap_depth_sensitivity_station_comparison.csv)
- [Table S10. Moving-block versus maximum-entropy bootstrap summary](../outputs/tables/bootstrap_method_sensitivity_summary.csv)
- [Table S11. Station-level bootstrap-method comparison](../outputs/tables/bootstrap_method_sensitivity_station_level.csv)

Primary moving-block inference uses 200 replicates. The 400-replicate rerun diagnoses Monte Carlo stability, while the maximum-entropy comparison assesses sensitivity to resampling method. These comparisons support the direction of the principal findings but show greater method dependence in some upper-quantile interval estimates.

## S3. Climate regimes and exploratory regionalization

- [Table S12. Köppen–Geiger quantile summaries](../outputs/tables/climate_regime_quantile_summary.csv)
- [Table S13. Climate-regime fixed-baseline summaries](../outputs/tables/climate_regime_fixed_baseline_summary.csv)
- [Table S14. Climate-regime response to regional warming](../outputs/tables/climate_regime_warming_response_summary.csv)
- [Table S15. Permutation tests for climate-regime contrasts](../outputs/tables/climate_regime_difference_tests.csv)
- [Table S16. Cluster robustness summary](../outputs/tables/cluster_robustness_summary.csv)
- [Table S17. Alternative clustering sensitivity](../outputs/tables/alternative_clustering_sensitivity_summary.csv)
- [Table S18. Spatial validation of cluster assignments](../outputs/tables/regional_cluster_spatial_validation.csv)

Climate-regime contrasts were corrected across 32 comparisons. Cluster assignments are retained only as exploratory descriptions because stability and spatial compactness differ among indices and analytical choices.

## S4. Compound dry–hot analysis

- [Table S19. Annual affected-station extent](../outputs/compound_dry_hot/tables/compound_dry_hot_yearly_extent.csv)
- [Table S20. Raw trend summaries](../outputs/compound_dry_hot/tables/compound_dry_hot_trend_summary.csv)
- [Table S21. Dependence-aware trend sensitivity](../outputs/compound_dry_hot/tables/compound_dry_hot_serial_dependence_sensitivity.csv)
- [Table S22. Early–late distribution-shift tests](../outputs/compound_dry_hot/tables/compound_dry_hot_distribution_shift_tests.csv)
- [Table S23. Event-component classification](../outputs/compound_dry_hot/tables/compound_dry_hot_driver_summary.csv)
- [Table S24. Moran's I trend summaries](../outputs/compound_dry_hot/tables/compound_dry_hot_moran_trend_summary.csv)
- [Table S25. Station-level event frequencies](../outputs/compound_dry_hot/tables/compound_dry_hot_station_frequency.csv)

The 5-, 10-, and 20-year labels are empirical rarity classes estimated within the 34-year record. They should not be interpreted as stable design return periods. The dependence-aware sensitivity uses 4-year residual moving blocks and 4,999 replicates.

## S5. Additional figures

- [Figure S1. Homogeneity diagnostics](../outputs/figures/ijoc_data_quality_homogeneity.png)
- [Figure S2. Sensitivity to homogeneity-flag exclusion](../outputs/figures/ijoc_homogeneity_sensitivity.png)
- [Figure S3. Bootstrap-depth sensitivity](../outputs/figures/ijoc_bootstrap_depth_sensitivity.png)
- [Figure S4. Bootstrap-method sensitivity](../outputs/figures/advanced_method_sensitivity/bootstrap_method_sensitivity.png)
- [Figure S5. Alternative clustering sensitivity](../outputs/figures/ijoc_alternative_clustering_sensitivity.png)
- [Figure S6. Köppen–Geiger station regimes](../outputs/figures/advanced_climate_regimes/koppen_geiger_station_regimes.png)
- [Figure S7. Compound-event frequency changes](../outputs/compound_dry_hot/figures/compound_dry_hot_frequency_change_maps.png)
- [Figure S8. Compound-event connectedness](../outputs/compound_dry_hot/figures/compound_dry_hot_moran_connectedness.png)

## S6. Scope of the supporting evidence

The supplementary outputs support observational claims about the analyzed station network and period. They do not add evidence for anthropogenic attribution, causal physical mechanisms, health impacts, or stable century-scale return periods. Exploratory composite “fingerprint scores” remain available among historical pipeline outputs but are intentionally excluded from this evidentiary package because the score weights are not externally calibrated and several components are dependent.
