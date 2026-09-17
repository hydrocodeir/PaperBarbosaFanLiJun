# Output audit and curation — 17 September 2026

## Inventory and decision

All 705 original files (721.7 MiB) were inventoried, hashed and inspected for file readability. CSV schemas, missing values, exact row duplicates and numerical infinities were recorded. There were 704 readable files and one PNG with an invalid IDAT checksum. Every original file has a decision and reason in the [curation manifest](../outputs/audit_cleanup/curation_manifest.csv).

583 obsolete, redundant or excluded files (667.2 MiB) were selected for removal from the active tree. An additional old output index was backed up before replacement. The archive is [output_cleanup_20260917.zip](../archives/output_cleanup_20260917.zip); every archived source file was verified against its original SHA-256. The archive's SHA-256 is `b100424322fda6f4837e81b03cbc395e1a9b2c38c1c412aba01b1c0b9457a4ed`. Deletion completed: **True**. The archive preserves recovery, including the originally corrupt image; it does not repair that image.

## Numerical verification

- All 40,176 station quantile estimates and 496 focal station–index summaries were rerun from the annual indices and matched the stored estimates.
- All 32 stored bootstrap-summary fields across 496 station–index groups were recalculated from 99,200 saved draws; maximum discrepancy was floating-point rounding. Saved 200-replicate draws were not generated again.
- Quality and detrended homogeneity diagnostics, exclusions, baseline/alternative clustering, representative selection, spatial tests, fixed-baseline indices, warming associations, climate-raster assignments and group summaries, historical joint-rarity analysis and screening sensitivities were reproduced in an isolated work directory. See the 64 [numerical records](../outputs/audit_cleanup/numerical_checks.json) for exact table coverage.
- The compound extension's eight primary estimates and pointwise/family intervals were independently recalculated from station values and saved synchronized bootstrap draws. All 48 climate-group intervals, eight group-weighted network closures and 80 sensitivity component means/sample sizes were checked.
- Both thermal annual index definitions were independently rebuilt from daily observations by [validate_publication.py](../validate_publication.py); its nine checks passed. This validation's historical bootstrap/clustering scope statement pertains to that script alone; the separate audit additionally reran clustering.
- The 12 network warming-response slopes displayed in S9 were also solved independently by linear programming; the largest absolute difference was below 0.000005 days per degree Celsius. Station-level warming fits were reproduced with the historical iteration limit; some iterative fits emitted convergence warnings and their convergence is not independently certified.

No unresolved numerical mismatch was found in these checks. This conclusion establishes consistency with recorded inputs and implementations, not observational truth or the validity of every historical statistical interpretation. Alternative 400-replicate and maximum-entropy bootstrap draws were unavailable: only their saved station results and aggregate summaries were checked, not new ensembles. The sensitivity-scenario intervals were checked for ordering; this audit did not independently rerun all ten 4,999-draw extension scenarios.

## Corrections and scientific selection

1. Tail analytic probabilities/intervals are absent in the stored station fits. The curated FDR display uses **NA**, replacing the misleading zero produced by summing missing values. Median retained counts are 115, 104, 97 and 86. Bootstrap tail intervals remain separately available.
2. “Signal emergence” maps were relabeled as descriptive bootstrap precision. The ratio is not an emergence date or calibrated detection probability.
3. Historical composite fingerprint outputs were removed from active evidence: overlapping components and independent shifts of related indices do not provide an external-forcing attribution test. Five numerical score tables and their graphics remain recoverable in the archive.
4. Interpolated surfaces and hundreds of repeated station exports were replaced by station-point panels and selected, reproducible representative profiles. All underlying retained station quantiles and bootstrap draws remain available.
5. Historical joint-rarity outputs are explicitly separated from fixed marginal AND events. Their varying station network, available-row denominator and inverse-probability definition are disclosed; “return period” and causal “driver” interpretations are not carried into the supplement.
6. Duplicate display tables were removed; canonical partition, sensitivity and regime tables supply the manuscript and supplement directly. Climate classification is retained in main Figures 5/9, main tables and supplementary Table S8.
7. One broken PNG (`40790_robat_e_poshtebadam_figure4.png`) was excluded with the superseded station exports.

All nonstation legacy figure families were visually inspected using contact sheets, and the complete curated supplementary set was visually reviewed. All original station images underwent integrity checks and their shared numerical sources were checked; this was not a separate full-resolution visual inspection of every station image. The final set has **10 main figures, 11 supplementary figures, 5 main tables and 8 supplementary tables**, plus the linked [machine-readable data catalog](Supplementary_Data_Catalog.md). Each figure is available as PDF, SVG, PNG and TIFF.

## Reproducibility and recovery

Use `python audit_output_data.py recompute`, then `python check_output_dependencies.py` and `python check_warming_solver.py` for the numerical audit. Use `python build_supplementary.py` and `python build_publication_docs.py` to rebuild the curated supplement and atlases. The initial inventory is immutable; do not rerun the inventory phase over this cleaned release. Historical score comparisons after deletion require restoring their archived inputs into a separate review copy.

The recovery ZIP stores project-relative paths. Extract selected files into a separate folder for inspection before choosing to restore them; a blind extraction over the current output tree would reintroduce superseded material. The raw `data/`, reference `assets/` and original manuscript were not modified by cleanup. The archive remains outside `outputs/` and therefore reduces clutter, not total project disk usage.

Data units, source-data and boundary redistribution permissions, station metadata/homogenization and submission-specific requirements remain author responsibilities as described in the manuscript. These limitations are not resolved by a successful numerical audit.
