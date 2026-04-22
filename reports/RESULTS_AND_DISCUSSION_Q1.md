# Q1-Oriented Results, Discussion, and Manuscript Positioning

## 1. Q1-readiness assessment

### Overall assessment

This study already has several elements that are compatible with a strong Q1 submission: a long observational period (1961-2024), national-scale station coverage, distribution-aware trend analysis through quantile regression, uncertainty quantification through bootstrap procedures, and a second layer of spatial/regional interpretation through clustering and driver analysis. In its current form, the study is methodologically stronger than a routine ETCCDI-style mean-trend paper.

The present manuscript appears closest to a **solid upper-mid-tier manuscript with real Q1 potential**, provided that the framing, synthesis, and physical interpretation are tightened. The main opportunity is not merely adding more figures or statistics, but turning the existing outputs into a sharper scientific argument about **how and where the distribution of thermal extremes is being reshaped across Iran**.

### Current strengths

1. The analysis goes beyond mean trends and explicitly resolves quantile-dependent change.
2. The warm and cool indices are examined jointly, which helps reveal asymmetry across the full distribution.
3. Station-level heterogeneity is preserved rather than masked by regional averaging.
4. Robustness checks are already available for reference period, bootstrap method, interpolation method, and clustering configuration.
5. The dataset spans 64 years and 30 stations, which is adequate for a national observational synthesis.

### Current weaknesses that still separate the manuscript from a top Q1 paper

1. The central scientific claim needs to be stated much more sharply.
2. The discussion must connect statistical patterns to plausible climate mechanisms more explicitly.
3. The strongest results are currently dispersed across many outputs rather than synthesized into a small number of memorable take-home findings.
4. The spatial interpretation should be presented carefully because the Moran's I results suggest weak field-scale spatial autocorrelation.
5. The international relevance must be strengthened by linking the findings to broader arid, semi-arid, Mediterranean, and Middle Eastern warming literature.

### Practical verdict

If the manuscript is rewritten around a clear claim such as:

> Observed warming in Iran is not a uniform shift in mean temperature extremes, but a quantile-dependent reorganization of warm and cool event distributions, with especially strong amplification in the upper tail of warm days and warm nights at selected stations and a widespread contraction of cool extremes.

then the paper can credibly target a Q1 journal in climate dynamics, regional climate, or climate impacts/extremes.

---

## 2. Recommended full-paper structure

## Suggested title direction

**Quantile-dependent changes in warm and cool temperature extremes across Iran during 1961-2024: station-level trends, uncertainty, and regionalization**

Alternative:

**Beyond mean warming: quantile regression reveals asymmetric changes in warm and cool temperature extremes across Iran**

## Recommended manuscript architecture

### 1. Introduction

1. Explain why mean-based trend analysis is insufficient for extremes.
2. Position Iran as a climatically diverse, water-stressed, and societally exposed region.
3. Review warm/cool day-night asymmetry from global and regional literature.
4. Identify the gap: limited full-distribution analyses with uncertainty and regionalization.
5. End with explicit contributions:
   - full-distribution quantile trend analysis for four thermal-extreme indices
   - bootstrap-supported uncertainty characterization
   - spatial and regional clustering of quantile signatures
   - sensitivity testing across methodological choices

### 2. Data and Methods

1. Stations, temporal coverage, and data quality.
2. Index construction.
3. Quantile regression and OLS trend estimation.
4. Bootstrap uncertainty framework.
5. Spatial significance and multiple-testing control.
6. Clustering and driver-analysis framework.
7. Sensitivity analyses.

### 3. Results and Discussion

Use the chapter drafted below.

### 4. Conclusions

1. State the 3-4 strongest findings only.
2. Emphasize why distribution-aware diagnostics matter for climate-risk interpretation.
3. Clarify the national and international implications.
4. End with limitations and future work.

---

## 3. Draft abstract

Thermal extremes are often summarized through mean changes, although societal and environmental impacts are frequently governed by distributional shifts that occur unevenly across quantiles. Here we investigate quantile-dependent trends in four annual temperature-extreme indices across Iran for 1961-2024 using observations from 30 meteorological stations. We analyze warm days, cool days, warm nights, and cool nights using quantile regression across the full conditional distribution, complement the estimates with bootstrap-based uncertainty diagnostics, and examine spatial structure through clustering, field-significance assessment, and geographic driver analysis. The results show that warming-related indices increase strongly but non-uniformly across stations and quantiles. Warm days exhibit clear upper-tail amplification at many stations, with especially large positive Delta1 values in Esfahan, Ahvaz, and Zanjan, indicating that the most rapidly increasing part of the distribution is often its warmest tail. Warm nights also increase substantially, but with stronger station-to-station heterogeneity and several cases of extreme upper-tail intensification, most notably at Shiraz. In contrast, cool-day and cool-night indices generally decline, with the strongest negative slopes concentrated in upper quantiles, implying a rapid collapse of years historically characterized by many cool events. Field significance is strongest for intermediate quantiles, whereas Moran's I indicates limited national-scale spatial autocorrelation, highlighting the importance of local and regional controls rather than a single coherent national pattern. Sensitivity analyses show that the principal conclusions are robust to alternative reference periods, bootstrap methods, and interpolation choices, although the most extreme tails remain inherently uncertain. These findings demonstrate that thermal-extreme change across Iran is better understood as a quantile-dependent reorganization of the temperature-extremes distribution than as a uniform shift in the mean, with direct implications for climate-risk assessment in arid and topographically complex regions.

## Suggested contribution statement

1. This study provides one of the few national-scale, station-based assessments of thermal extremes in Iran that explicitly resolves quantile-dependent trends rather than relying on mean changes alone.
2. It integrates quantile regression, bootstrap uncertainty estimation, clustering, and methodological sensitivity analysis in a single observational workflow.
3. It shows that warm and cool extremes are evolving asymmetrically across both quantiles and stations, revealing tail-specific changes that are obscured in conventional trend analyses.

---

## 4. Results and Discussion

## 4.1. Data completeness and analytical scope

The observational basis of this analysis is strong enough for publication-grade inference. The final dataset includes **30 stations**, **64 years** of annual records, and a common analysis window spanning **1961-2024**. No station-index series fell below either the minimum requirement for quantile regression or the more conservative publication threshold of 20 years, so the reported spatial contrasts are not driven by obvious differences in record length.

![Data coverage](figures/data_coverage_by_station.png)

This uniform temporal depth is important because the central aim of the study is not simply to detect monotonic change, but to compare how the **shape of the trend distribution** varies across indices, stations, and quantiles. In that sense, the data structure supports a genuinely distribution-aware analysis rather than a collection of short, weakly comparable station records.

## 4.2. Regional-scale quantile behavior of thermal extremes

At the regional scale, averaging annual indices across stations reveals a highly asymmetric climate signal. The two warming-related indices show persistent positive trends, while the two cooling-related indices show persistent negative trends. However, the magnitude of change is quantile-dependent, and the dependence differs sharply by index.

### Table 1. Regional mean and selected quantile trend estimates

| Index | OLS slope | OLS 95% CI | QR q0.01 | QR q0.50 | QR q0.95 | QR q0.99 |
|---|---:|---|---:|---:|---:|---:|
| Warm Days | +9.53 | [ +7.54, +11.52 ] | +15.38 | +8.15 | +12.05 | +11.40 |
| Warm Nights | +13.58 | [ +12.05, +15.10 ] | +18.88 | +13.65 | +12.55 | +12.09 |
| Cool Days | -2.55 | [ -4.01, -1.10 ] | -1.12 | -2.01 | -8.05 | -16.58 |
| Cool Nights | -2.98 | [ -4.03, -1.93 ] | -2.04 | -2.47 | -6.78 | -10.49 |

The regional warm-day signal is clearly positive throughout the distribution, but it is not uniform. The quantile-coefficient profile shows that upper quantiles tend to exceed the median trend, indicating that years already characterized by relatively high counts of warm days have intensified faster than the center of the distribution. This is a strong sign of **upper-tail amplification**, which is climatologically more consequential than a simple mean increase because it implies accelerated growth in already warm years.

![Regional quantile coefficients for warm days](figures/region_quantile_slopes_warm_days.png)

Warm nights also exhibit strong positive change across almost the full distribution, and the regional OLS slope is the largest among the four indices. However, unlike warm days, the warm-night signal is somewhat more uniform across central and upper quantiles, while very low quantiles also show strong positive slopes. This suggests that nighttime warming is not confined to the most extreme subset of years, but instead reflects a broader shift in the baseline distribution of warm-night occurrence.

![Regional quantile coefficients for warm nights](figures/region_quantile_slopes_warm_nights.png)

The two cool indices tell a different story. Both cool days and cool nights decline regionally, but the strongest negative changes emerge in the upper quantiles. Because high quantiles of a cool-extremes index correspond to years with unusually many cool events, the steeply negative slopes at q0.95-q0.99 indicate that those historically cool years are disappearing fastest. In other words, the contraction of cool extremes is not merely a reduction in average frequency; it is a disproportionate erosion of the upper end of the cool-event distribution.

![Regional quantile coefficients for cool days](figures/region_quantile_slopes_cool_days.png)

![Regional quantile coefficients for cool nights](figures/region_quantile_slopes_cool_nights.png)

Taken together, the regional analysis supports a first central conclusion: **Iran's thermal-extreme climate has not shifted as a rigid block**. Instead, warm and cool indices have been restructured differently across quantiles, with warm extremes generally increasing and cool extremes contracting, but with the strongest responses often concentrated away from the center of the distribution.

## 4.3. Station-scale heterogeneity and tail amplification

Regional means are informative, but they conceal substantial station-scale heterogeneity. The station heatmaps show that the four indices do not share a single national pattern; rather, they consist of geographically mixed but statistically coherent groups with very different slope signatures.

![Warm-days station signatures](figures/station_focus_heatmap_warm_days.png)

![Warm-nights station signatures](figures/station_focus_heatmap_warm_nights.png)

![Cool-days station signatures](figures/station_focus_heatmap_cool_days.png)

![Cool-nights station signatures](figures/station_focus_heatmap_cool_nights.png)

The most interpretable heterogeneity appears in **Delta1 = slope(q0.95) - slope(q0.05)**, which measures whether the upper tail is intensifying faster than the lower tail.

### Table 2. Stations with the strongest upper-tail amplification or contraction

| Index | Strongest positive Delta1 stations | Strongest negative Delta1 stations |
|---|---|---|
| Warm Days | Esfahan (+16.58), Ahvaz (+15.63), Zanjan (+14.11) | Bushehr Airport (-6.33), Arak (-2.95), Khoy (-2.89) |
| Warm Nights | Shiraz (+25.90), Bam (+13.00), Tehran (+11.91) | Ramsar (-5.38), Orumiyeh (-4.89), Khorramabad (-4.00) |
| Cool Days | Birjand (-2.14), Kerman (-2.99), Zahedan (-3.30) | Orumiyeh (-16.95), Torbat-e Heydariyeh (-12.10), Qazvin (-11.00) |
| Cool Nights | Shahrekord (+11.96), Esfahan (+2.86), Tabriz (-2.26) | Khorramabad (-21.01), Rasht (-16.94), Qazvin (-13.76) |

For **warm days**, the positive Delta1 values at Esfahan, Ahvaz, and Zanjan indicate that the hottest portion of the warm-day distribution is intensifying faster than its lower tail. This is a particularly strong signal because these stations are not merely warming in mean terms; they are becoming more unevenly weighted toward the upper end of the warm-day distribution.

For **warm nights**, the station-level heterogeneity is even stronger. Shiraz stands out as an exceptional upper-tail outlier, with Delta1 approaching +26 days per decade, far above the rest of the network. Bam and Tehran also show marked upper-tail intensification. These results imply that the most extreme warm-night years are escalating at a faster pace than modestly warm-night years in selected locations, which is highly relevant for heat-stress exposure because nocturnal heat accumulation undermines physiological recovery.

The **cool-day** index is strikingly different: Delta1 is negative almost everywhere of practical importance, showing that upper quantiles decline faster than lower quantiles. This means that years once characterized by many cool days are collapsing more rapidly than already low-cool-day years. The same broad pattern is visible for **cool nights**, although with a few local reversals such as Shahrekord and Esfahan.

This station-scale evidence strengthens the second central conclusion of the manuscript: **the direction of change is only half the story; the internal distribution of change matters just as much**. Warm extremes tend to expand through upper-tail amplification at many stations, whereas cool extremes generally contract through upper-tail collapse.

## 4.4. Spatial distribution and field significance

The mapped results show that spatial organization exists, but it is not strong enough to justify treating Iran as a smoothly varying, nationally coherent field in all cases. This is an important methodological and scientific point.

Illustrative Delta1 maps are shown below.

![Warm-days Delta1 map](figures/map_warm_days_Delta1.png)

![Warm-nights Delta1 map](figures/map_warm_nights_Delta1.png)

![Cool-days Delta1 map](figures/map_cool_days_Delta1.png)

![Cool-nights Delta1 map](figures/map_cool_nights_Delta1.png)

The spatial significance analysis indicates that intermediate quantiles carry the strongest and most spatially persistent signal. After false-discovery-rate control:

### Table 3. Number of FDR-retained significant station results by index and quantile

| Quantile | Warm Days | Warm Nights | Cool Days | Cool Nights |
|---|---:|---:|---:|---:|
| 0.10 | 27 | 29 | 0 | 12 |
| 0.50 | 27 | 28 | 17 | 23 |
| 0.90 | 17 | 25 | 26 | 25 |

![Raw vs FDR-retained significance counts](figures/advanced_spatial_inference/raw_vs_fdr_counts.png)

These counts reveal three meaningful features. First, **warm indices are highly robust at lower and median quantiles**, especially warm nights. Second, **cool indices become most significant toward the upper quantiles**, consistent with the regional coefficient profiles showing that the loss of cool events is especially concentrated in years formerly rich in cool extremes. Third, the outermost tails at q0.05 and q0.95 remain analytically less stable, which is expected given the finite sample and the inherent difficulty of tail estimation.

At the same time, Moran's I values remain non-significant for all 20 station-slope fields tested, with none reaching p < 0.05. This suggests that although many stations share the same sign of change, the detailed spatial arrangement is not well described by a smooth, strongly autocorrelated national field.

![Moran's I summary](figures/advanced_spatial_inference/moran_i_heatmap.png)

This result should be treated as a strength, not a weakness. It implies that the study is correctly identifying a climate system shaped by **regional contrasts, local controls, and topographic complexity**, rather than forcing an overly smooth interpretation on a sparse observational network.

## 4.5. Regionalization and clustering of quantile signatures

The clustering analysis reinforces the idea that thermal-extreme change in Iran is structured by multiple regimes rather than one single national trajectory.

![Warm-days cluster composites](figures/advanced_regionalization/cluster_composite_heatmap_warm_days.png)

![Warm-nights cluster composites](figures/advanced_regionalization/cluster_composite_heatmap_warm_nights.png)

![Cool-days cluster composites](figures/advanced_regionalization/cluster_composite_heatmap_cool_days.png)

![Cool-nights cluster composites](figures/advanced_regionalization/cluster_composite_heatmap_cool_nights.png)

For **warm days**, four clusters are clearly meaningful. Cluster 1, although small, has the highest median Delta1 (+13.98), indicating very strong upper-tail amplification. Cluster 4 also shows strong amplification (median Delta1 = +8.22), while Cluster 2 is much weaker (median Delta1 = +2.26) and Cluster 3 is slightly negative (median Delta1 = -2.60), indicating stations where warm-day change is more uniform or even mildly inverted across tails.

For **warm nights**, the cluster structure is even more revealing. One singleton cluster centered on Shiraz has an extremely high Delta1 (+25.90), while another multi-station cluster shows consistently strong positive Delta1 (median around +10.85). By contrast, the coastal/high-baseline warm-night group exhibits high warming across all quantiles but only modest tail contrast, meaning that some stations are warming strongly without necessarily showing severe within-distribution amplification.

For **cool days** and **cool nights**, the network is dominated by one large cluster in each case, but the outlier clusters are scientifically valuable. They indicate that most stations share a broad common pattern of declining cool extremes, yet a few locations diverge strongly enough to form their own regimes, especially in the upper tail.

The clustering robustness metrics are sufficiently reassuring for interpretation. The adjusted Rand index is 0.632 for warm days, 0.934 for warm nights, 0.685 for cool days, and 0.821 for cool nights when comparing the configured and reduced-feature runs. This indicates that the regionalization is not an artifact of one arbitrary feature choice, even though exact label identities shift, as expected in unsupervised classification.

## 4.6. Geographic drivers and physical interpretation

The driver analysis is intentionally simple, but it still provides clues that help anchor the quantile patterns in physical geography.

### Table 4. Most interpretable driver relationships

| Index-metric | Predictor | Standardized beta | p-value | Spearman rho | Interpretation |
|---|---|---:|---:|---:|---|
| Warm Nights, q0.05 | Elevation | -0.507 | 0.006 | -0.521 | Lower-elevation stations tend to warm more strongly in the lower tail of warm nights |
| Warm Nights, q0.50 | Elevation | -0.479 | 0.010 | -0.550 | Median warm-night intensification is stronger at lower elevations |
| Warm Nights, q0.95 | Elevation | -0.318 | 0.094 | -0.432 | Upper-tail warm-night growth also weakens with elevation, though less strongly |
| Cool Days, q0.50 | Longitude | +0.471 | 0.021 | +0.355 | East-west geographic structure influences the median decline of cool days |
| Warm Days, Delta1 | Latitude | -0.281 | 0.185 | -0.378 | Tail amplification tends to strengthen southward, although the regression is not formally significant |

The clearest signal is the **negative relationship between warm-night slopes and elevation**. This implies that lower-elevation stations experience stronger increases in warm-night frequency, especially in lower and median quantiles. Such a pattern is physically plausible because low-elevation and often drier or more urbanized environments may retain and re-radiate heat more effectively, suppressing nighttime cooling.

The median cool-day slope also shows a significant positive association with longitude. Because cool-day slopes are generally negative, a positive beta suggests that the strongest declines tend to occur toward the western side of the network or, conversely, that eastern stations show somewhat weaker declines in the median cool-day regime. This east-west organization deserves further discussion in relation to continentality, regional circulation, and land-surface dryness in the full manuscript.

More broadly, the spatially mixed but statistically coherent patterns suggest that no single mechanism can explain all four indices. A defensible interpretation is that:

1. **Warm nights** are especially sensitive to lower-elevation and possibly urban or moisture-limited conditions.
2. **Warm days** show strong tail amplification at selected interior stations, consistent with intensification of the hottest years rather than uniform warming alone.
3. **Cool days and cool nights** decline widely, but their strongest contractions are concentrated in upper quantiles, indicating rapid disappearance of historically cool years.

## 4.7. Robustness and sensitivity of the main findings

One of the strongest aspects of this study is that the central conclusions do not rely on a single fragile analytical setting.

![Reference-period sensitivity](figures/advanced_method_sensitivity/reference_period_sensitivity.png)

![Bootstrap-method sensitivity](figures/advanced_method_sensitivity/bootstrap_method_sensitivity.png)

![Interpolation sensitivity at tau=0.05](figures/advanced_method_sensitivity/interpolation_sensitivity_tau_0.05.png)

![Interpolation sensitivity at tau=0.95](figures/advanced_method_sensitivity/interpolation_sensitivity_tau_0.95.png)

The sensitivity checks show that the broad directional conclusions remain stable across changes in reference period, bootstrap method, and interpolation strategy. The largest reported discrepancies occur for a limited subset of metrics, such as warm-day Delta1 under reference-period changes and cool-index bootstrap summaries at low quantiles. Importantly, these sensitivity signals modify the **magnitude** of some estimates more than their **qualitative interpretation**.

This is exactly the kind of robustness profile that strengthens a Q1 manuscript. It allows the discussion to distinguish between:

1. findings that are structurally robust and can be emphasized confidently, and
2. findings near the outermost tails that should be interpreted cautiously.

Accordingly, the manuscript should foreground the following as robust results:

1. widespread increase in warm-day and warm-night indices,
2. widespread decrease in cool-day and cool-night indices,
3. strong upper-tail amplification at selected warm-index stations,
4. strong upper-quantile contraction in cool indices,
5. meaningful but spatially non-smooth regional heterogeneity.

## 4.8. Synthesis and high-level interpretation

The combined evidence supports three major scientific conclusions.

First, **warming in Iran is distributionally asymmetric**. Warm-day and warm-night indices increase broadly, but their rate of increase depends on quantile and station. This means that the climate system is not merely shifting upward in a mean sense; it is being redistributed internally, often in ways that intensify the warmest years faster than the rest.

Second, **cool extremes are collapsing through upper-quantile contraction**. The largest negative slopes for cool days and cool nights occur in high quantiles, implying that years that were once unusually rich in cool events are disappearing fastest. This is a stronger and more policy-relevant statement than saying only that cool events are decreasing on average.

Third, **local and regional controls matter strongly**. The non-significant Moran's I results, the strong station-to-station variation in Delta1, the presence of outlier clusters, and the elevation signal in warm nights all suggest that thermal-extreme change across Iran is shaped by a combination of broad warming and geographically specific modifiers.

In publication terms, the core added value of this manuscript is therefore not just that it documents warming, which is already widely known, but that it shows **how the distribution of thermal extremes is being reorganized**, where that reorganization is strongest, and how robust those findings are to methodological choices.

---

## 5. Tables that should be cited in the manuscript even if not printed in full

The following tables are too detailed to embed fully in the narrative chapter, but they should be cited explicitly in the manuscript text or supplementary material:

1. `outputs/tables/publication_summary_table.csv`
   This is the main station-level table for slopes, confidence intervals, Delta metrics, bootstrap summaries, and cluster assignments.

2. `outputs/tables/regional_cluster_composites.csv`
   This contains the cluster-level composite statistics used to interpret regionalization.

3. `outputs/tables/station_significance_fdr.csv`
   This is the full field-significance table by station, index, and quantile.

4. `outputs/tables/driver_analysis_summary.csv`
   This should support the driver-analysis subsection and any supplementary discussion.

5. `outputs/tables/reference_period_sensitivity_summary.csv`
   This provides the compact summary for reference-period robustness.

6. `outputs/tables/bootstrap_method_sensitivity_summary.csv`
   This provides the compact summary for bootstrap-method robustness.

7. `outputs/tables/interpolation_method_sensitivity_summary.csv`
   This supports the discussion of spatial-interpolation sensitivity.

---

## 6. Suggested concluding sentence for the Results and Discussion chapter

Overall, the results show that thermal-extreme change across Iran is best understood not as a uniform rise or fall in annual indices, but as a quantile-dependent restructuring of warm and cool event distributions, with strong local heterogeneity, robust broad-scale directional signals, and important implications for interpreting climate risk in an arid and topographically diverse setting.
