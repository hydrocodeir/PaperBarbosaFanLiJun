# Diagnostic review and revision decisions

Date: 17 September 2026. Source manuscript: `reports/Manuscript.md`; preserved without modification. The originally referenced `DONT READ/Manuscript.md` is absent from the working tree; the author confirmed the reports copy.

## Article-specific adaptation

- **Type:** observational hydroclimatology; 124-station Iranian archive, 1991–2024.
- **Target:** international climate/environmental journal; no journal or current quartile has been selected or verified.
- **Core contribution:** connect distributional changes in annual thermal-event counts with fixed-baseline dry–hot concurrence, while distinguishing changing marginal event frequencies from excess joint occurrence.
- **Strongest existing evidence:** warm-day quantile contrast; widespread direction of warm/cool station trends; fixed-baseline thermal changes; small effects of temperature consistency screening; dry–hot network changes under the historical joint-rarity definition.
- **Strengths:** daily source archive, reproducible modules, bootstrap and sensitivity outputs, three relevant supplied reference articles.
- **Vulnerabilities:** short baseline and analysis period; unhomogenized daily observations; sparse tails; uncertain data provenance; ambiguous compound-event definition; dependence and multiplicity; many competing figures and repeated discussion/conclusions.

## Decisions before rewriting

1. Do not sell quantile regression, clustering, or empirical copulas as new methods: Barbosa et al. (2011), Fan (2014), and Alizadeh et al. (2020) already supply those precedents.
2. Replace joint-rarity classes as the headline compound result with an explicit AND event: precipitation below an early-period quantile AND temperature above its early-period quantile. Low joint survival probability alone does not ensure both margins cross an extreme threshold.
3. Add a symmetric, exact probability partition, calculated separately at each station. The residual is a change in indicator covariance; it must not be called isolated copula change, physical feedback, or causal attribution.
4. Rebuild compound aggregates from screened daily observations. Count coverage against calendar days, including days absent from the CSV. Use balanced stations to prevent changing network membership from driving the period comparison.
5. Resample year blocks synchronously across stations and variables. Include baseline-threshold re-estimation, coverage, block-length, quantile, zero-rainfall, tie-rule, and common-network sensitivities.
6. Following the author’s expansion request, organize the article into ten main figures and one supplementary figure. Use point maps, common centered color scales, vector exports, explicit units, and captions separating station spread from uncertainty.
7. Retain historical analyses in the supplement. Distinguish their original run provenance from the new raw-data extension and validation.

## Literature positioning

The supplied PDFs were inspected directly. Paper 1 is Barbosa et al. (2011); Paper 2 is Fan (2014); Paper 3 is Alizadeh et al. (2020). Their scientific questions and estimands inform the framing; their prose and figure artwork are not copied.

Recent primary literature also rules out a blanket novelty claim for marginal/dependence decomposition: [Bevacqua et al. (2022)](https://www.nature.com/articles/s41558-022-01309-5) and [Schmutz et al. (2026)](https://nhess.copernicus.org/articles/26/881/2026/index.html). The contribution is a transparent station-based application and a sharper test of the interpretation, not invention of an established family of methods. This was a targeted positioning search, not a systematic proof that no Iranian study has used the approach.

## Exact limitation wording

“The excess-joint term measures a change in covariance between threshold indicators. Because indicator covariance can change when marginal frequencies change even under an unchanged copula, this term does not isolate a change in the dependence structure or identify a physical feedback. The 17-year baseline, discrete event counts, and zero summer precipitation further limit the precision and generality of the partition.”

## Remaining submission requirements

The author confirmed the Iran Meteorological Organization (IRIMO) as data provider. Confirm measurement units and access/license, and boundary-data provenance/license; supply author declarations and an archived code/data identifier; choose the target journal and verify its current scope and quartile in the intended ranking system. The complete historical station-bootstrap/clustering pipeline has not been rerun by this revision. New compound analyses and both full-network thermal index reconstructions are covered by the separate validation report.

## Final reviewer assessment

The new analysis and publication graphics have been produced, and the English manuscript has been rewritten around the resulting evidence. All 20 scenario–definition combinations were run with 4,999 replicates. Nine scientific validation checks passed, including complete reconstruction of both thermal index sets. Ten main and one supplementary figures were rendered in four formats; the principal plots and maps were visually inspected. Figure and supplement links and numerical tables were checked separately.

The principal new finding is increasing explicit AND-event frequency, with a larger summer hot-frequency point component and no resolved network excess-joint change. This is a defensible observational contribution, with no claim of a new statistical method or physical attribution. A reviewer can still question baseline length, bootstrap coverage, unhomogenized daily data, nonrandom missing precipitation, and thermal analytic probabilities. These limitations are explicit in the revised Methods and Discussion. Administrative declarations, journal-specific formatting and a permanent archive remain necessary before submission.


## Expansion following author review

The climate classification has been restored to the main methods, results and discussion, with two new figures and two numerical tables. Compound uncertainty was summarized within six fixed climate groups using the original synchronized primary bootstrap fields; matching network draws and independently reproduced thermal group means were verified. No between-group significance test was inferred from interval overlap.

The conclusion now develops the findings, regime context, robustness and limits in five paragraphs. Table 5 compares eight relevant studies. Four references were added after checking primary sources, including two direct Iranian studies published in 2026. Their existence narrows the novelty claim and strengthens the comparison. The inherited HadEX3 author list was corrected. See `Reference_Audit_2026.md` for access limits and bibliographic evidence. The curated reference library prevents later document builds from reverting these changes.
