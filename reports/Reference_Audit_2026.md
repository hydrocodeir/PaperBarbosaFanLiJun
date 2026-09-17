# Reference and comparison audit

Checked: 17 September 2026. This was a targeted bibliographic and claim audit, not a systematic literature review. The manuscript contains 25 references, including four additions in this expansion. No claim of being the first Iranian quantile-regression study is retained.

## Bibliographic verification

The [metadata audit](../outputs/publication_v2/references/reference_metadata_audit.csv) records the DOI, registered title, authors, publication dates, journal, volume and pages for every reference. [Raw Crossref responses](../outputs/publication_v2/references/crossref_records.json) are retained. Twenty-four DOIs returned metadata. Fan's original DOI returned HTTP 404 from Crossref; its title, author, 2014 issue, pages and original DOI were instead verified against the supplied `assets/Paper2.pdf` and the [publisher-hosted paper](https://www.tandfonline.com/doi/pdf/10.3878/j.issn.1674-2834.13.0102). A failed registry query is not recorded as successful validation.

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
