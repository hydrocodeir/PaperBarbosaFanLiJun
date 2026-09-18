# Reference and comparison audit

Initial audit: 17 September 2026; targeted novelty/comparison update: 18 September 2026. This is a targeted bibliographic and claim audit, not a systematic literature review. The manuscript now contains 27 references, including Zhao and Xiong (2026) for literature positioning and Fitzenberger (1998) for block-bootstrap regression inference. No claim of being the first Iranian quantile-regression study is retained.

## Bibliographic verification

The [metadata audit](../outputs/publication_v2/references/reference_metadata_audit.csv) records the DOI, registered title, authors, publication dates, journal, volume and pages for every reference. [Raw Crossref responses](../outputs/publication_v2/references/crossref_records.json) are retained. Twenty-four DOIs returned metadata in the initial audit; the latest update verifies Zhao and Xiong as the twenty-fifth. Fan's original DOI returned HTTP 404 from Crossref; its title, author, 2014 issue, pages and original DOI were instead verified against the supplied `assets/Paper 2.pdf` and the [publisher-hosted paper](https://www.tandfonline.com/doi/pdf/10.3878/j.issn.1674-2834.13.0102). A failed registry query is not recorded as successful validation.

Corrections made:

- Dunn et al. (2020): the inherited 34-name list omitted authors and included incorrect names. It was replaced with the 56 registered authors in source order. [HadEX3 source](https://doi.org/10.1029/2019JD032263).
- Alexander et al. (2006): corrected the author name from “Kumar KR” to “Rupa Kumar K.” [Source](https://doi.org/10.1029/2005JD006290).
- Künsch (1989): restored the surname diacritic consistently.
- Retained Soltani's 2016 volume year despite online publication in 2015, Fan's original 2014 issue year despite later hosting metadata, and Koenker's 2005 book year rather than its later electronic publication date.
- Distinguished Zhang et al. (2005a), the regional trends paper, from Zhang et al. (2005b), the percentile-index construction paper.

The curated bibliography is stored in [verified_references.json](verified_references.json); rebuilding the documents reads this file so new references and corrections are preserved. Metadata matching alone does not verify a scientific claim. The comparisons below were checked separately against supplied full papers or primary publisher text.

## New references and their roles

| Addition | Role in the manuscript | Evidence accessed |
| --- | --- | --- |
| [Ghasemi (2026)](https://doi.org/10.1029/2025EA004860) | Direct Iranian precedent for PCA and quantile regression; distinguishes monthly temperature intensity from annual event counts | Publisher full text, especially abstract, data section and Figure 3 discussion |
| [Jamali et al. (2026)](https://doi.org/10.1016/j.jaridenv.2026.105606) | Comparison with changing Iranian climate types; distinguishes dynamic transitions from fixed grouping | Publisher abstract and accessible article text; no unreported numerical results inferred |
| [Beck et al. (2023)](https://doi.org/10.1038/s41597-023-02549-6) | Documents the availability of updated climate maps and a future sensitivity option | Publisher/author-repository metadata and abstract; this raster was not used in our analysis |
| [Zhang et al. (2005b)](https://doi.org/10.1175/JCLI3366.1) | Explains potential baseline-related inhomogeneity in percentile indices | Publisher abstract and registered bibliographic record |

## Traceability of comparative statements

| Comparison | Source location and check | Interpretation constraint |
| --- | --- | --- |
| Barbosa et al. (2011) | Supplied `assets/Paper1.pdf`; [publisher](https://doi.org/10.5194/nhess-11-3227-2011) | Daily-temperature quantiles are a methodological precedent, not annual-count effect sizes. |
| Fan (2014) | Supplied `assets/Paper2.pdf`, results for warm nights / Figure 1c; publisher PDF | Three warm-night slopes were transcribed in their original quantile order; different baselines, homogenization and periods preclude a regional sensitivity inference. |
| Soltani et al. (2016) | [Publisher abstract and Section 3.1.1](https://link.springer.com/article/10.1007/s00704-015-1609-5) | The count trends cited are for 1995–2010, not the entire 1975–2010 archive. Temperature-unit figures elsewhere in the source were not used for count comparisons. |
| Ghasemi (2026) | Publisher abstract and Figure 3 discussion | Reported cold/hot temperature-tail rates concern different quantiles and units from our annual-count analysis. |
| Alizadeh et al. (2020) | Supplied `assets/Paper3.pdf`; [publisher](https://doi.org/10.1126/sciadv.aaz4571), data, methods and event-composition results | Record period, climate-division support and annual/seasonal windows differ; joint rarity is not a marginal AND threshold. |
| Bevacqua et al. (2022) | [Publisher](https://doi.org/10.1038/s41558-022-01309-5), abstract and study framing | Projected compound occurrence cannot be pooled with an observed two-period change. |
| Schmutz et al. (2026) | [Publisher](https://nhess.copernicus.org/articles/26/881/2026/index.html), Sections 3–5 | Full copula components and emergence timing differ from our exact binary probability partition. |
| Jamali et al. (2026) | [Publisher abstract](https://www.sciencedirect.com/science/article/abs/pii/S0140196326000583) | Published transition percentages refer to record-length subsets, not all 279 stations or our fixed classes. |

The original method references remain for trend estimation, bootstrap design, multiplicity, index construction and homogeneity diagnostics. They were checked bibliographically; this audit does not claim an independent replication of all their methods or results. There is no evidence synthesis based on secondary summaries or invented citations.

## Reviewer implications

The discussion now compares agreement, disagreement and non-comparability, rather than simply listing citations. Climate classes remain in the main paper because they expose heterogeneous point estimates and the uneven impact of zero precipitation thresholds. No time-trend difference between classes is promoted to an FDR-supported result. A 2023 classification citation does not imply that the older, actually used 2018 raster was replaced. Zhang's baseline-bias correction is also distinguished from resampling to estimate uncertainty; the former was not implemented by the latter.

## Targeted novelty update — 18 September 2026

Table 5 now compares nine studies through their published finding, present evidence, added value and comparability constraint. The introduction and Section 4.5 identify an empirical and interpretive contribution; they do not claim that quantile regression or probability decomposition is new. The conclusion reflects the same scope. The previous manuscript is preserved in [the dated revision snapshot](archive/Manuscript_Q1_2026_before_novelty_20260918.md).

- **Fan (2014):** the supplied `assets/Paper 2.pdf`, Section 4.1 and Figure 1b, reports cool-day q10 and q90 slopes of +0.19 and +0.07 days per decade, with little evidence of tail trends. The present −5.60 and −15.37 values were checked against `table01_thermal_trends.csv`. This adds a concrete difference in the reported tail pattern. It is not a harmonized between-country significance test. The source's 1961–1990 threshold reference and homogenized records are now explicit in the comparison.
- **Ghasemi (2026):** the [publisher text](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2025EA004860), abstract and Section 3.3, supports the distinction between monthly temperature intensity and annual event counts. No count/intensity effect-size ratio is inferred.
- **Schmutz et al. (2026):** the [publisher full text](https://nhess.copernicus.org/articles/26/881/2026/) supports the clarification that dependence can affect emergence timing. Our unresolved indicator-covariance difference cannot establish unchanged copula dependence.
- **Zhao and Xiong (2026):** added from the [publisher abstract and bibliographic record](https://rmets.onlinelibrary.wiley.com/doi/abs/10.1002/joc.70324). Their decomposition of projected probability ratios into distribution-parameter contributions is directly relevant to positioning our observed frequency partition. Only that methodological scope is used. A numerical inconsistency visible in one abstract interval was not carried into our comparison; no numerical effect-size comparison is made. Full-text methods were not independently replicated.
- **Jamali et al. (2026):** the [publisher abstract](https://www.sciencedirect.com/science/article/abs/pii/S0140196326000583) again supports the distinction between changing classifications and our fixed classes.
- **Soltani (2016) and Bevacqua (2022):** direct page retrieval was inconsistent during this update. Existing comparisons were retained from the earlier recorded source check, without extracting additional quantitative claims. Their findings are not described as falsified by the present observations.

All new present-study numerical entries in Table 5 are taken from Tables 1, 3 and 4 and the canonical CSV sources. The comparisons remain descriptive where periods, definitions or spatial support differ. A [reviewer report](Q1_Reviewer_Report_2026_FA.md) identifies additional analyses that could strengthen the contribution; those recommendations are not represented as completed analyses.

## Thermal network uncertainty update — 18 September 2026

[Fitzenberger (1998)](https://doi.org/10.1016/S0304-4076(97)00058-4) was added as a methodological reference for moving-block bootstrap inference in least-squares and quantile regression. DOI metadata and the publisher abstract support this scope; full-text assumptions were not independently replicated. The reference is not used to claim exact validity for the present 34-year nonstationary record. That limitation is explicit in Section 2.3. Table 5 now uses the common 108-station estimates and distinguishes unresolved thermal asymmetry from supported trend directions. See the [implementation report](Thermal_Network_Revision_2026_FA.md).

## Index-construction method check — 18 September 2026

Zhang et al. (2005b) was already cited and now supports an implemented in-base replacement correction. The publisher landing page could not be retrieved during this update; the algorithm and type-8 convention were cross-checked against the maintainers’ [climdex.pcic implementation](https://raw.githubusercontent.com/pacificclimate/climdex.pcic/master/src/zhang_running_quantile.cc) and [ETCCDI definitions](https://etccdi.pacificclimate.org/indices_def.shtml). The method excludes the target year, duplicates each other baseline year, evaluates target-year counts separately, and averages the resulting indices. Independent tests use literal duplicated samples and NumPy quantiles. The revised manuscript explicitly distinguishes the study’s 17-year baseline/annual coverage rule from full ETCCDI compliance. No extra paper was needed and the reference count remains 27. The [implementation report](Index_Zero_Revision_2026_FA.md) records scope and results.
