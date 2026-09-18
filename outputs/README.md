# Curated publication outputs

Start with the [revised manuscript](../reports/Manuscript_Q1_2026.md) and [supplementary material](../reports/Supplementary_Q1_2026.md).

- [Main figure atlas](publication_v2/Figure_Atlas.pdf): Figures 1–10.
- [Supplementary figure atlas](publication_v2/Supplementary_Figure_Atlas.pdf): Figures S1–S15.
- [Complete supplementary data catalog](../reports/Supplementary_Data_Catalog.md): canonical CSV files, roles and limitations.
- [Output audit](../reports/Output_Audit_2026.md): numerical checks, corrections and archive location.
- [Curation manifest](audit_cleanup/curation_manifest.csv): a decision and reason for every original output.

The active publication figures are defined by the two manifests in `publication_v2/`; each has PDF, SVG, PNG and TIFF exports. Numerical sources remain in `tables/`, `compound_dry_hot/tables/` and `publication_v2/tables/`. Historical station graphics, superseded map exports and unsupported composite scores were archived before deletion. Original observations and manuscript files were not changed by cleanup.

Build curated figures with `python build_supplementary.py`, then documents and atlases with `python build_publication_docs.py`, from the project root. A complete legacy `run_analysis.py` run can regenerate superseded outputs and requires renewed curation.

Thermal network trends now use a common 108-station set, with direct uncertainty for slope asymmetry and day–night contrasts. See the [revision report](../reports/Thermal_Network_Revision_2026_FA.md). Reproduce with `python run_thermal_network.py` and verify with `python validate_thermal_network.py`.

Index-construction and structural-zero comparisons are documented in the [second revision report](../reports/Index_Zero_Revision_2026_FA.md), with [independent checks](publication_v2/index_zero_validation.json). The active supplement contains 15 figures and 14 tables.
