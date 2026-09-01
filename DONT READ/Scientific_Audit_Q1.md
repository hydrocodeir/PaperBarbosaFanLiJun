# Scientific audit and revision record

## Scope

This audit covers `reports/Manuscript.md`, the analysis code and configuration, the available tabular and graphical outputs, and the rewritten manuscript `reports/Manuscript_Q1_Revised.md`. The review followed the project-level Q1 protocol, the local AI-writing rules, and the reference-based rewriting instructions. Alizadeh et al. (2020), supplied as `assets/Paper 3.pdf`, was used as a benchmark for results-led organization, empirical compound-event framing, and concise scientific prose; no text was copied from that article.

## Article-specific assessment

- **Manuscript type:** observational hydroclimatology and environmental statistics.
- **Target tier:** international Q1, journal not yet selected.
- **Core contribution:** station-scale evidence that recent changes in Iranian thermal extremes differ across the annual count distribution and coincide with a widening empirical warm-season dry–hot footprint.
- **Strongest evidence:** (1) warming-consistent station signs and median FDR counts; (2) marked upper-versus-lower quantile differences for warm and cool days; (3) persistence under fixed early-period thresholds and multiple sensitivity checks; (4) positive dependence-aware trends in the number of stations affected by empirical dry–hot conditions; and (5) a late-period shift toward heat as the rarer marginal component of joint events.
- **Principal vulnerabilities:** short 34-year record; no metadata-supported daily homogenization; incomplete data provenance; many station/quantile comparisons; uneven station geography; rarity estimates defined within the same short record; and no direct physical-driver or attribution data.

## Major findings of the audit

### Scientific strengths

1. The 124-station network provides broad coverage of Iran's principal topographic and climatic settings.
2. Quantile regression answers a question that a mean-only trend cannot: whether changes differ between low- and high-count years.
3. Multiple-testing control, block resampling, fixed-threshold analysis, homogeneity exclusion, and method comparisons provide several independent robustness diagnostics.
4. The empirical compound analysis complements the thermal indices without claiming that precipitation deficit is equivalent to every form of drought.
5. The available outputs are sufficiently rich to support a focused article without relying on the exploratory composite score.

### Major weaknesses in the original manuscript

1. The title and several passages used promotional “fingerprint” framing that exceeded what the observational design could establish.
2. Results and Discussion were combined, long, and repetitive; primary evidence was difficult to distinguish from secondary diagnostics.
3. The compound affected-station trends were initially reported without an explicit dependence-aware trend sensitivity.
4. The methods described a “5-day moving window,” whereas the implementation uses all observations within ±5 calendar days, an 11-day centered window.
5. The equal-weight composite fingerprint score combined dependent quantities without external calibration and was not suitable as headline evidence.
6. Internal daily inconsistencies among minimum, mean, and maximum temperature were not quantified in the manuscript.
7. Several mechanism-oriented statements went beyond the variables analyzed; no circulation, soil moisture, radiation, humidity, land-use, or attribution data are present.
8. Data provenance, licensing, author contributions, funding, conflicts, and permanent code/data availability were absent or incomplete.

## New analyses and outputs added

### Dependence-aware compound-trend sensitivity

A residual moving-block analysis was added for the annual affected-station series, using 4-year blocks and 4,999 replicates. For warm-season events with empirical joint return periods of at least 10 years, the trend remained positive (Kendall τ = 0.583, block-null p = 0.0002; slope = 19.24 stations decade⁻¹, 95% block-bootstrap interval 13.94–24.53). The corresponding 20-year-class result was also positive (τ = 0.539, p = 0.0002; slope = 8.72, interval 5.66–11.68 stations decade⁻¹). Annual-definition results were weaker but remained positive at the 10- and 20-year thresholds.

Output: `outputs/compound_dry_hot/tables/compound_dry_hot_serial_dependence_sensitivity.csv`.

### Internal temperature-consistency sensitivity

The archive contains 45 daily records at 32 stations where minimum temperature exceeds maximum temperature (0.0029% of daily rows). After those overlapping cases were excluded, 2,165 records at 61 stations had mean temperature outside the minimum–maximum interval (0.1406%). A non-destructive sensitivity rerun masked these entries and rebuilt the thermal indices and compound summaries. Across the focal thermal results, the maximum absolute slope change was 0.026 days decade⁻¹ and the maximum absolute change in Δ1 was 0.024 days decade⁻¹. Compound affected-station means changed by no more than 0.53 station, and slopes by no more than 0.38 station decade⁻¹.

Outputs:

- `outputs/tables/temperature_internal_consistency_screening.csv`
- `outputs/tables/temperature_internal_consistency_quantile_sensitivity.csv`
- `outputs/tables/temperature_internal_consistency_compound_sensitivity.csv`

### Reproducibility changes

- The pipeline now registers and produces the added sensitivity tables.
- The compound configuration records the 4,999-replicate dependence check.
- A project `requirements.txt` records the Python dependency set.
- Package initialization and one utility import were made lazy so that light-weight checks do not unnecessarily load the full analysis stack.

## Rewriting decisions

1. The title now describes the variables, method-relevant contribution, region, and period without unsupported novelty language.
2. The abstract is quantitative and states the principal limits.
3. Results and Discussion are separated.
4. The main Results retain seven figures and three compact tables; exploratory clustering and detailed sensitivity outputs are assigned to supplementary material.
5. Composite fingerprint scores are explicitly excluded from primary evidence.
6. Causal and attribution claims are replaced by observational, mechanism-aware wording.
7. The 5-, 10-, and 20-year compound thresholds are described as within-record empirical rarity classes, not stable long-return-period hazard estimates.
8. The supplied reference article influenced the argument architecture: direct research questions, results-led subsections, explicit empirical probability definition, and a restrained comparison between event components.

## Reviewer-style verdict

### Recommendation at the current stage: major revision before submission

The revised scientific narrative is substantially stronger and the added tests address two important reviewer concerns. The manuscript should not yet be submitted because the data source and access conditions are unknown in the project text, the declaration fields are incomplete, and no target-journal format has been selected. These are submission blockers rather than evidence that the central numerical findings are invalid.

### Required before submission

1. Replace every bracketed declaration placeholder with verified information.
2. Name the station-data provider, product/version or acquisition route, original units, access date, license, and redistribution limits.
3. Archive the exact code/configuration release and cite a permanent DOI or commit.
4. Select the target journal and conform structure, word count, reference style, figure resolution, and supplementary-file conventions.
5. Run the complete pipeline in the declared environment and preserve its log. The targeted additions were executed against the existing derived data, but the entire multi-hour pipeline was not rerun during this revision because `statsmodels` was unavailable in the current environment and network installation was unavailable.
6. Confirm author names, affiliations, CRediT roles, funding, conflicts, and acknowledgements.

### Strongly recommended

1. If station histories or a suitable reference network can be obtained, perform formal relative homogenization of the daily series; otherwise retain the present limitation language.
2. Deposit station-level outputs and a data dictionary even if raw observations cannot be redistributed.
3. Prepare a supplementary inventory mapping every reported result to a file and code stage.
4. Avoid extending the discussion to health impacts, physical mechanisms, or anthropogenic attribution unless corresponding observations or model experiments are added.

## Reproducibility status

- Source-data description: **incomplete**.
- Preprocessing and analytical code: **present**.
- Derived tables and figures: **present**.
- Added sensitivity outputs: **present**.
- Dependency declaration: **added**.
- Full clean-environment rerun: **still required**.
- Permanent repository/archive identifier: **still required**.
- Manuscript-to-output consistency: **checked for the principal reported results and all linked figure paths**.

## Final claim calibration

- **Strongly supported:** widespread warming-consistent changes in percentile-based station indices over 1991–2024; positive median station trends after FDR control; positive warm-season affected-station trends under the block sensitivity.
- **Moderately supported:** stronger distributional change in daytime than nighttime indices; meaningful variation among climate settings; an increasing role of heat rarity among empirical joint events.
- **Exploratory:** cluster-defined regions, split-period acceleration/deceleration, and specific physical mechanisms.
- **Not supported by this design:** anthropogenic attribution, causal effects of Köppen–Geiger class or elevation, long-return-period hazard estimation, and impact or policy effectiveness claims.
