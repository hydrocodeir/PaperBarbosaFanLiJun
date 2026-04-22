# IJoC Upgrade Log

This note records how the 18 high-priority upgrades for raising the manuscript toward a stronger `Q1` / `International Journal of Climatology` standard were implemented in the current project outputs.

## Status overview

The main upgraded manuscript is:

- [IJoC_submission_refined.md](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/IJoC_submission_refined.md)

The supporting polished figures are:

- [ijoc_study_area.png](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/figures/ijoc_study_area.png)
- [ijoc_regional_quantile_panels.png](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/figures/ijoc_regional_quantile_panels.png)
- [main_figure_representative_stations.png](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/figures/ijoc_station_comparisons/main_figure_representative_stations.png)
- [ijoc_main_delta1_maps.png](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/figures/ijoc_main_delta1_maps.png)
- [ijoc_robustness_synthesis.png](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/figures/ijoc_robustness_synthesis.png)
- [ijoc_split_period_comparison.png](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/figures/ijoc_split_period_comparison.png)

The new manuscript tables are:

- [ijoc_table_main_summary.csv](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/tables/ijoc_table_main_summary.csv)
- [ijoc_table_day_night_asymmetry.csv](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/tables/ijoc_table_day_night_asymmetry.csv)
- [ijoc_table_representative_stations.csv](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/tables/ijoc_table_representative_stations.csv)
- [ijoc_split_period_regional_summary.csv](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/tables/ijoc_split_period_regional_summary.csv)
- [ijoc_elevation_class_summary.csv](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/tables/ijoc_elevation_class_summary.csv)

## Step-by-step implementation

1. Central claim sharpened.
Implemented in the Abstract, Introduction, Section `3.1`, and Conclusions of [IJoC_submission_refined.md](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/IJoC_submission_refined.md). The manuscript is now explicitly framed around the idea that Iran's thermal extremes reflect a `quantile-dependent reorganization`, not a simple mean shift.

2. Discussion strengthened from description to explanation.
Implemented mainly in Sections `3.3`, `3.5`, `3.6`, and `3.8`, where the text now interprets day-night differences, upper-tail amplification, upper-tail collapse, elevation effects, and wider dryland / Mediterranean parallels.

3. International comparison strengthened.
Implemented in the Introduction and Section `3.8` through explicit connection to Mediterranean, Middle Eastern, and global extremes literature, rather than keeping the argument purely national.

4. Limitations written proactively.
Implemented in Section `3.9`, including finite station network, tail uncertainty, interpolation-as-visualization, and descriptive rather than formal breakpoint interpretation.

5. Novelty statement made explicit.
Implemented in the Introduction through a five-part contribution statement: full-distribution trends, tail-contrast metrics, uncertainty integration, cluster-based regionalization, and robustness testing.

6. Day-night asymmetry added as a dedicated synthesis result.
Implemented in Section `3.3`, supported by [ijoc_table_day_night_asymmetry.csv](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/tables/ijoc_table_day_night_asymmetry.csv).

7. Upper-tail amplification versus upper-tail collapse made explicit.
Implemented across Sections `3.2` to `3.5`, especially the contrast between warm-day / warm-night amplification and cool-day / cool-night contraction.

8. Statistical significance translated into climatological relevance.
Implemented in Section `3.2` and Table 1 by translating slopes into approximate cumulative change over the 1961-2024 period.

9. Driver analysis strengthened.
Implemented in Section `3.6` through stronger interpretation of geographic controls and the addition of an elevation-class synthesis table:
[ijoc_elevation_class_summary.csv](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/tables/ijoc_elevation_class_summary.csv)

10. Temporal nonlinearity tested.
Implemented using split-period comparison for `1961-1990` versus `1991-2024`, reported in Section `3.7`, [ijoc_split_period_comparison.png](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/figures/ijoc_split_period_comparison.png), and [ijoc_split_period_regional_summary.csv](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/tables/ijoc_split_period_regional_summary.csv).

11. Main-text figure strategy made selective and memorable.
Implemented by restructuring the narrative around six main figures instead of many scattered figures. This is reflected directly in [IJoC_submission_refined.md](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/IJoC_submission_refined.md).

12. Study-area / station-network figure added.
Implemented as [ijoc_study_area.png](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/figures/ijoc_study_area.png) and integrated as Figure 1.

13. Main multi-panel map figure added.
Implemented as [ijoc_main_delta1_maps.png](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/figures/ijoc_main_delta1_maps.png), used to summarize the spatial expression of `Delta1` across all four indices.

14. Robustness figure consolidated.
Implemented as [ijoc_robustness_synthesis.png](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/figures/ijoc_robustness_synthesis.png), combining reference-period, bootstrap-method, interpolation, and clustering diagnostics.

15. Figure clutter reduced.
Implemented by moving many detailed outputs to the supplementary strategy section and replacing multiple separate station-comparison figures with one polished multi-panel figure:
[main_figure_representative_stations.png](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/figures/ijoc_station_comparisons/main_figure_representative_stations.png)

16. Main summary table created.
Implemented as [ijoc_table_main_summary.csv](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/tables/ijoc_table_main_summary.csv) and embedded as Table 1 in the refined manuscript.

17. Representative-station table created.
Implemented as [ijoc_table_representative_stations.csv](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/tables/ijoc_table_representative_stations.csv) and embedded as Table 3 in the refined manuscript.

18. Main-text vs supplementary-material separation clarified.
Implemented in Section `3.9` and the `Supplementary material strategy` section, which now explicitly identifies which figures and tables should move out of the main text.

## Code changes supporting the upgrades

The plotting and pipeline changes supporting these outputs were implemented in:

- [plotting.py](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/src/paper_pipeline/plotting.py)
- [pipeline.py](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/src/paper_pipeline/pipeline.py)

Key additions include:

- study-area plotting
- regional quantile multi-panel plotting
- representative-station multi-panel plotting
- main `Delta1` map synthesis
- robustness synthesis plotting
- split-period comparison plotting
- updated `region_quantile_slopes_*.png` display range to `q = 0.05-0.95`

## Recommended next manuscript file

For the strongest journal-facing draft, use:

- [IJoC_submission_refined.md](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/IJoC_submission_refined.md)

For the earlier fuller draft with broader background material and the existing reference block, keep:

- [IJoC_submission_ready.md](/d:/Pooya/w/GitHub/HydroCodeIR/PaperBarbosaFanLiJun/outputs/IJoC_submission_ready.md)
