# Curated publication outputs

Start with the [revised manuscript](../reports/Manuscript_Q1_2026.md) and [supplementary material](../reports/Supplementary_Q1_2026.md).

- [Main figure atlas](publication_v2/Figure_Atlas.pdf): Figures 1–10.
- [Supplementary figure atlas](publication_v2/Supplementary_Figure_Atlas.pdf): Figures S1–S11.
- [Complete supplementary data catalog](../reports/Supplementary_Data_Catalog.md): canonical CSV files, roles and limitations.
- [Output audit](../reports/Output_Audit_2026.md): numerical checks, corrections and archive location.
- [Curation manifest](audit_cleanup/curation_manifest.csv): a decision and reason for every original output.

The active publication figures are defined by the two manifests in `publication_v2/`; each has PDF, SVG, PNG and TIFF exports. Numerical sources remain in `tables/`, `compound_dry_hot/tables/` and `publication_v2/tables/`. Historical station graphics, superseded map exports and unsupported composite scores were archived before deletion. Original observations and manuscript files were not changed by cleanup.

Build curated figures with `python build_supplementary.py`, then documents and atlases with `python build_publication_docs.py`, from the project root. A complete legacy `run_analysis.py` run can regenerate superseded outputs and requires renewed curation.
