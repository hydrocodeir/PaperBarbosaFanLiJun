# Quantile-dependent changes in warm and cool temperature extremes across Iran during 1961-2024

**Running title:** Quantile-dependent thermal extremes in Iran

## Abstract

Thermal extremes are often summarized using mean trends, although many of their climatic and societal impacts depend on changes across the full distribution. This study investigates quantile-dependent trends in four annual temperature-extreme indices across Iran during 1961-2024 using observations from 30 meteorological stations. Warm days, warm nights, cool days, and cool nights were analysed using quantile regression, with ordinary least squares estimates used as a mean-trend benchmark. The analysis was complemented by bootstrap uncertainty estimation, field-significance testing with false-discovery-rate control, cluster-based regionalization, and sensitivity analyses for reference period, bootstrap method, and interpolation strategy. The results show strong but non-uniform warming-related changes. Warm days and warm nights increased across much of the network, but the magnitude of the increase depended on quantile and station, with pronounced upper-tail amplification at several locations. At the regional scale, warm-day and warm-night mean trends were +9.53 and +13.58 days year^-1 per decade, respectively. In contrast, cool-day and cool-night indices declined regionally, with mean trends of -2.55 and -2.98 days year^-1 per decade, and with the strongest negative slopes concentrated in upper quantiles, indicating rapid loss of years historically characterized by many cool events. Station-level contrasts were substantial: the strongest warm-day upper-tail amplification occurred in Esfahan, Ahvaz, and Zanjan, whereas warm-night amplification was most pronounced in Shiraz. Field significance was strongest for intermediate quantiles, but Moran's I statistics showed little evidence of strong national-scale spatial autocorrelation, implying an important role for local and regional controls. Sensitivity analyses confirmed that the principal findings are robust despite larger uncertainty in the most extreme tails. These findings show that thermal-extreme change across Iran is better interpreted as a quantile-dependent restructuring of warm and cool event distributions than as a uniform shift in the mean, with implications for climatological understanding in arid and topographically complex regions.

**Keywords:** climate extremes; quantile regression; warm nights; cool days; regionalization; Iran; thermal extremes

## 1. Introduction

Changes in climate extremes are among the most consequential manifestations of contemporary climate change because their impacts are transmitted through agriculture, public health, water resources, ecosystems, and infrastructure `(Easterling et al., 2000; IPCC, 2021)`. While mean temperature trends provide a useful first-order description of climatic change, they are often insufficient for understanding extremes, since the practical significance of climate change depends not only on a shift in the centre of the distribution but also on changes in its tails `(Katz and Brown, 1992; Katz et al., 2002; Easterling et al., 2000)`. A system may exhibit moderate mean warming while simultaneously undergoing much stronger intensification in the warmest conditions or a pronounced contraction of historically cool conditions. For this reason, approaches that resolve the full conditional distribution of climate indicators are increasingly important in climatology `(Koenker and Bassett, 1978; Koenker, 2005)`.

This issue is especially relevant for temperature-extreme indices. Warm and cool day-night indices have been widely used to characterize climate variability and change `(Frich et al., 2002; Alexander et al., 2006; Donat et al., 2013; Zhang et al., 2011)`, yet their interpretation often relies on linear trends in the mean or, at most, the median. Such summaries can mask quantile-dependent behaviour, including upper-tail amplification, lower-tail stabilization, or asymmetric decline across the distribution `(Moberg and Jones, 2005; Kostopoulou and Jones, 2005)`. In physically diverse regions, these distributional effects may vary sharply across space because of topographic contrasts, continentality, land-atmosphere coupling, urban influence, and differences in moisture availability.

Iran provides an especially instructive setting for such an analysis. The country spans strong gradients in elevation, aridity, latitude, and marine influence, and includes climatic regimes ranging from humid Caspian lowlands to hot southern and interior lowlands, as well as high-elevation continental environments. These contrasts make Iran a useful natural laboratory for examining whether changes in warm and cool extremes are spatially coherent or whether they are strongly conditioned by local and regional controls. Previous studies in the broader climatological literature have documented warming and changes in thermal extremes across dryland, Mediterranean, and semi-arid regions `(Zhang et al., 2005; Kostopoulou and Jones, 2005; Soltani et al., 2016; Zittis et al., 2016)`, yet relatively few station-based studies for Iran have resolved these changes across the full distribution while simultaneously addressing uncertainty, field significance, and regionalization `(Soltani et al., 2016)`.

This distinction is important beyond the study area. The *International Journal of Climatology* emphasizes that regional studies should contribute clearly to wider climatological understanding, even when their empirical focus is geographically specific. In that spirit, the relevance of the present study lies not simply in documenting change over Iran, but in examining a broader question: whether observed thermal-extreme change is better represented as a uniform mean shift or as a quantile-dependent reorganization of the distribution. This issue is likely to be important for other arid and topographically complex regions where climatic gradients are sharp and local modifiers interact strongly with broader warming signals `(Kostopoulou and Jones, 2005; Zittis et al., 2016)`.

The present study addresses this issue using 30 meteorological stations with annual data spanning 1961-2024. Four annual indices are analysed: warm days, warm nights, cool days, and cool nights. Quantile regression is used to estimate trend behaviour from low to high quantiles, and these estimates are interpreted alongside ordinary least squares trends, bootstrap uncertainty measures, field-significance diagnostics, spatial analysis, and clustering. The specific objectives are:

1. to quantify how thermal-extreme trends vary across quantiles for each index;
2. to identify stations and regions exhibiting strong upper-tail amplification or contraction;
3. to assess whether the detected patterns are field-significant and spatially coherent; and
4. to evaluate the robustness of the main conclusions to methodological choices.

By doing so, the study aims to contribute not only a national assessment of thermal extremes in Iran, but also a more general demonstration of how quantile-based approaches can refine climatological interpretation in observational extreme-index studies.

## 2. Data and methods

### 2.1. Data and study design

The analysis is based on daily station observations from 30 meteorological stations distributed across Iran. After preprocessing, the common analysis period spans 1961-2024, yielding 64 annual values for each station-index series. The workflow uses station-level observations of daily minimum temperature (`tmin`) and daily maximum temperature (`tmax`) to construct annual thermal-extreme indices. Annual aggregation is performed after daily threshold exceedances are identified.

The analysis is station-based throughout. National and regional summaries are derived from station estimates rather than from gridded products. This design is appropriate for an observational quantile-trend study because it preserves local climatic heterogeneity and avoids imposing interpolation assumptions at the primary inference stage.

### 2.2. Construction of annual thermal-extreme indices

Four annual indices were derived from daily temperature observations: warm days, warm nights, cool days, and cool nights. Daily thresholds were calculated separately for each station using a reference period of 1961-1990, with a 5-day moving calendar window around each day of year. The 90th percentile of daily maximum temperature was used to define warm days, the 10th percentile of daily maximum temperature to define cool days, the 90th percentile of daily minimum temperature to define warm nights, and the 10th percentile of daily minimum temperature to define cool nights. Leap-day observations were removed before threshold estimation to ensure a consistent 365-day annual cycle.

For each station and day of year, thresholds were computed from reference-period samples. Where the moving-window sample size was insufficient, the available station-specific reference sample was used as fallback. Daily binary exceedance indicators were then summed within year to obtain annual counts of the four indices. In effect, each annual series represents the number of days per year meeting the corresponding warm or cool condition, following the logic commonly used in percentile-based temperature-extreme diagnostics `(Frich et al., 2002; Alexander et al., 2006; Zhang et al., 2011)`.

### 2.3. Quantile regression and mean-trend estimation

The principal analytical tool is quantile regression `(Koenker and Bassett, 1978; Koenker, 2005)`. For each station-index series, linear trends were estimated across the conditional quantiles of the annual counts. Time was expressed in decades using:

`x = (year - min(year)) / 10`

so that slopes are interpretable as changes in days year^-1 per decade. Quantile-specific slopes were estimated on the full quantile grid used by the workflow, and a set of focal quantiles was retained for summary and comparison.

Ordinary least squares (OLS) regression was used as a benchmark for mean-trend estimation. For both OLS and quantile regression, the fitted slope represents the linear rate of change in the annual index as a function of time. For quantile regression, where available, confidence intervals for the slope coefficients were extracted from the fitted model. When exact estimation was required for very short series, or when the regression procedure encountered instability, fallback logic in the workflow was used to obtain stable slope estimates from the same station series.

At the regional scale, annual indices were first averaged across stations by year and then analysed using the same quantile-regression framework. These regional estimates should be interpreted as distributional summaries of the station-network mean series rather than as spatially interpolated national fields.

### 2.4. Summary metrics of distributional asymmetry

To summarize whether the upper tail of the trend distribution was intensifying or contracting faster than the lower tail, the following contrast metric was computed:

`Delta1 = slope(q0.95) - slope(q0.05)`

Positive `Delta1` values indicate stronger upper-tail increases, whereas negative values indicate either upper-tail weakening or, for declining indices, stronger upper-tail contraction. Two auxiliary metrics were also available in the workflow:

`Delta2 = slope(q0.95) - slope(q0.50)`

`Delta3 = slope(q0.50) - slope(q0.05)`

Among these, `Delta1` is used most prominently because it captures the full contrast between upper- and lower-quantile trend behaviour and therefore offers a compact summary of distributional asymmetry.

### 2.5. Bootstrap uncertainty analysis

Bootstrap resampling was used to characterize uncertainty in the station-level quantile slopes. The configured workflow uses 200 bootstrap replicates per station series and applies a maximum-entropy bootstrap (`meboot`) as the default method, while also enabling sensitivity comparison with other bootstrap schemes. For each replicate, slopes were re-estimated at the focal quantiles and summarized using bootstrap means, standard deviations, and percentile-based confidence intervals.

This step serves two purposes. First, it provides an uncertainty estimate that does not rely exclusively on asymptotic analytical inference. Second, it allows direct evaluation of the stability of derived contrast metrics such as `Delta1`. In the present manuscript, bootstrap summaries are used primarily to support robustness and interpretation rather than to replace the main quantile-regression estimates `(Efron and Tibshirani, 1993; Vinod, 2006)`.

### 2.6. Field significance and spatial dependence

To assess whether statistically significant station results exceeded what might be expected under multiple testing, field significance was examined using false-discovery-rate (FDR) control `(Benjamini and Hochberg, 1995)`. Station-level analytic p-values for selected quantile slopes were adjusted across the network, and the number of FDR-retained significant results was summarized by quantile and index.

Because a network of stations may exhibit spatial dependence, Moran's I was also used to assess whether station-slope fields showed significant spatial autocorrelation `(Moran, 1950)`. The workflow evaluates Moran's I using permutation-based inference for selected quantiles and indices. This step is important because it helps distinguish widespread same-sign change from genuinely smooth spatial organization across the station field.

### 2.7. Spatial visualization and regionalization

The study distinguishes between **station-based inference** and **map-based visualization**. Station maps of slope metrics and cluster labels are generated directly from observed station coordinates. For selected quantiles, interpolated fields are also produced using radial-basis or gridded interpolation methods and masked to the national boundary where available. These interpolated surfaces are used primarily for visualization of broad spatial structure and are not treated as the main inferential basis, particularly given the finite station network.

Regionalization was performed using hierarchical clustering applied to station-level feature vectors derived from quantile slopes, `Delta` metrics, and uncertainty summaries. The configured workflow uses average linkage with Euclidean distance after feature standardization. Cluster robustness was evaluated through a reduced-feature rerun and comparison of solutions using adjusted Rand index and same-label fractions `(Hubert and Arabie, 1985)`.

### 2.8. Driver analysis and sensitivity assessment

To explore whether simple geographic controls were associated with thermal-extreme behaviour, station metrics were related to latitude, longitude, and elevation. Both standardized linear-regression coefficients and Spearman rank correlations were used. These analyses were intentionally exploratory and are interpreted as first-order geographic clues rather than as full mechanistic attribution.

Sensitivity testing was performed for three main methodological choices:

1. reference period used for threshold construction;
2. bootstrap method; and
3. interpolation approach used for visualization.

These checks were designed to determine whether the principal climatological conclusions were stable to reasonable analytical alternatives. In the present study, sensitivity analysis is used not merely as a technical appendix, but as part of the substantive interpretation of result robustness.

## 3. Results and discussion

### 3.1. Data coverage and analytical scope

The final dataset includes 30 stations and 64 years of annual records, covering the common period 1961-2024. No station-index series fell below either the minimum requirement for quantile regression or the more conservative publication threshold of 20 years. This matters because the study is fundamentally comparative: the aim is not only to detect temporal change, but also to compare the shape of change across indices, stations, and quantiles. The broadly uniform record length therefore provides an adequate basis for comparing distributional trend behaviour across the national network.

![Figure 1. Data coverage by station.](figures/data_coverage_by_station.png)

*Figure 1. Number of annual observations available for each station included in the analysis. The common temporal depth of the station series supports comparison of quantile-dependent trend behaviour across the national network.*

### 3.2. Regional-scale quantile behaviour

Regional averaging across stations reveals a clear asymmetry between warming-related and cooling-related indices. Warm days and warm nights both increased strongly over the study period, whereas cool days and cool nights generally declined. However, this broad contrast conceals marked quantile dependence.

**Table 1.** Regional mean and selected quantile trend estimates for the four annual thermal-extreme indices.

| Index | OLS slope | OLS 95% CI | QR q0.01 | QR q0.50 | QR q0.95 | QR q0.99 |
|---|---:|---|---:|---:|---:|---:|
| Warm Days | +9.53 | [+7.54, +11.52] | +15.38 | +8.15 | +12.05 | +11.40 |
| Warm Nights | +13.58 | [+12.05, +15.10] | +18.88 | +13.65 | +12.55 | +12.09 |
| Cool Days | -2.55 | [-4.01, -1.10] | -1.12 | -2.01 | -8.05 | -16.58 |
| Cool Nights | -2.98 | [-4.03, -1.93] | -2.04 | -2.47 | -6.78 | -10.49 |

For warm days, the regional signal is positive throughout the distribution, but the upper quantiles generally exceed the median slope. This indicates that years already characterized by relatively high numbers of warm days have intensified faster than the centre of the distribution. Such behaviour is consistent with upper-tail amplification and suggests that warming is strengthening the most anomalously warm years more rapidly than typical years.

![Figure 2. Regional quantile coefficients for warm days.](figures/region_quantile_slopes_warm_days.png)

*Figure 2. Quantile-regression coefficients for the regional warm-day series from q = 0.01 to q = 0.99. Black dots show quantile-specific slopes, grey shading shows available 95% confidence intervals, and red solid and dashed lines denote the OLS mean trend and its 95% confidence limits, respectively.*

Warm nights also increased strongly, and the regional OLS trend (+13.58 days year^-1 per decade) is the largest among the four indices. The quantile profile, however, is less strongly skewed toward the upper tail than for warm days. Instead, positive slopes occur across almost the full distribution, including low quantiles, indicating a broad upward shift in the frequency of warm nights. From a climatological perspective, this is important because widespread nighttime warming implies not only increased thermal extremes but also reduced nocturnal relief during hot periods.

![Figure 3. Regional quantile coefficients for warm nights.](figures/region_quantile_slopes_warm_nights.png)

*Figure 3. As in Figure 2, but for warm nights.*

The cool indices display the opposite sign and a different pattern of asymmetry. For cool days and cool nights, the strongest negative slopes are concentrated in the upper quantiles. Because high quantiles of these indices correspond to years with unusually many cool events, the steeply negative q0.95-q0.99 slopes imply that the historically coolest years are disappearing faster than the central tendency alone would suggest. This indicates a contraction of the cool-event distribution that is sharper in its upper tail than in its lower tail.

![Figure 4. Regional quantile coefficients for cool days.](figures/region_quantile_slopes_cool_days.png)

*Figure 4. As in Figure 2, but for cool days.*

![Figure 5. Regional quantile coefficients for cool nights.](figures/region_quantile_slopes_cool_nights.png)

*Figure 5. As in Figure 2, but for cool nights.*

Taken together, the regional results indicate that thermal-extreme change across Iran is not well described as a rigid translation of the mean. Instead, the warm and cool indices are being reshaped differently across quantiles, with warm extremes generally expanding and cool extremes contracting, but with the strongest responses often concentrated away from the centre of the distribution.

### 3.3. Station-scale heterogeneity and tail amplification

Regional means conceal substantial local heterogeneity. The station heatmaps show that the four indices do not share a single national signature; rather, they comprise multiple trend regimes with contrasting low-, median-, and high-quantile behaviour.

![Figure 6. Station heatmap for warm days.](figures/station_focus_heatmap_warm_days.png)

*Figure 6. Station-level heatmap for warm days showing q0.05 slope, q0.50 slope, q0.95 slope, and `Delta1` for each station.*

![Figure 7. Station heatmap for warm nights.](figures/station_focus_heatmap_warm_nights.png)

*Figure 7. As in Figure 6, but for warm nights.*

![Figure 8. Station heatmap for cool days.](figures/station_focus_heatmap_cool_days.png)

*Figure 8. As in Figure 6, but for cool days.*

![Figure 9. Station heatmap for cool nights.](figures/station_focus_heatmap_cool_nights.png)

*Figure 9. As in Figure 6, but for cool nights.*

An especially useful summary metric is `Delta1 = slope(q0.95) - slope(q0.05)`, which measures whether the upper tail is changing faster than the lower tail.

**Table 2.** Stations with the strongest upper-tail amplification or contraction for each index.

| Index | Strongest positive Delta1 stations | Strongest negative Delta1 stations |
|---|---|---|
| Warm Days | Esfahan (+16.58), Ahvaz (+15.63), Zanjan (+14.11) | Bushehr Airport (-6.33), Arak (-2.95), Khoy (-2.89) |
| Warm Nights | Shiraz (+25.90), Bam (+13.00), Tehran (+11.91) | Ramsar (-5.38), Orumiyeh (-4.89), Khorramabad (-4.00) |
| Cool Days | Birjand (-2.14), Kerman (-2.99), Zahedan (-3.30) | Orumiyeh (-16.95), Torbat-e Heydariyeh (-12.10), Qazvin (-11.00) |
| Cool Nights | Shahrekord (+11.96), Esfahan (+2.86), Tabriz (-2.26) | Khorramabad (-21.01), Rasht (-16.94), Qazvin (-13.76) |

For warm days, upper-tail amplification is especially pronounced in Esfahan, Ahvaz, and Zanjan. These stations are not simply warming in mean terms; the highest-quantile years are intensifying more rapidly than the low-quantile years. This indicates that the changing distribution of warm days is becoming increasingly weighted toward its upper end. In contrast, stations such as Bushehr Airport show negative Delta1 values, implying a more uniform trend or even relative moderation in the upper tail.

Warm nights show even stronger station-to-station contrasts. Shiraz is a particularly notable outlier, with a Delta1 of +25.90 days year^-1 per decade, far exceeding the rest of the network. Bam and Tehran also display marked upper-tail amplification. These findings are climatologically important because intensification in the upper tail of warm-night occurrence implies a disproportionate growth in years with especially high nocturnal heat burden.

The cool indices exhibit a different structure. Cool-day Delta1 is negative at nearly all stations of practical importance, indicating that the upper quantiles are declining faster than the lower quantiles. In other words, years once characterized by many cool days are collapsing faster than already less-cool years. Cool nights show a broadly similar tendency, although there are local reversals, most notably at Shahrekord and, to a lesser extent, Esfahan. These reversals underline the importance of local context and confirm that the decline of cool events is not spatially uniform.

Overall, the station-scale analysis shows that the direction of change alone is insufficient. What matters equally is how change is distributed within the index distribution. Warm extremes commonly expand through upper-tail amplification, whereas cool extremes tend to contract through upper-tail collapse.

### 3.3.1. Comparison of representative stations across cluster structure and focal quantiles

To make the cluster logic more transparent, a small set of representative stations was selected for direct comparison within each index. The selection was data-driven and designed to contrast both **cluster membership** and the relative behaviour of the three focal quantiles (`q0.05`, `q0.50`, and `q0.95`). For warm days and warm nights, one representative station was selected from each cluster so that the comparison spans the full range of observed quantile-structure regimes. For cool days and cool nights, where one dominant cluster coexists with a few highly distinctive secondary clusters, the selected stations were chosen to highlight both the dominant pattern and the most informative departures from it.

![Figure 9. Multi-panel comparison of representative stations across cluster-defined quantile structures.](figures/ijoc_station_comparisons/main_figure_representative_stations.png)

*Figure 9. Multi-panel comparison of representative stations selected to contrast cluster membership and focal-quantile behaviour for (a) warm days, (b) warm nights, (c) cool days, and (d) cool nights. Colours denote the selected representative stations, and markers identify the focal quantiles q = 0.05, 0.50, and 0.95. The small in-panel summaries report the corresponding focal-quantile slopes for each station. Together, the panels show that the clustering structure corresponds to visibly different quantile-slope geometries rather than only differences in average trend magnitude.*

For **warm days**, Esfahan was selected as the representative of Cluster 1 because it exhibits the strongest upper-tail amplification, Bushehr Airport represents Cluster 2 because it shows a rare negative `Delta1`, Khoy represents Cluster 3 because it combines a strong lower-tail slope with weaker upper-tail growth, and Zanjan represents Cluster 4 because it shows one of the clearest upward transitions from `q0.05` to `q0.95`. The resulting comparison shows that the warm-day regime is not merely a question of stronger or weaker warming; rather, the curvature of the quantile-slope profile differs fundamentally among clusters.

For **warm nights**, Ramsar was selected from Cluster 1 because it combines high warm-night trends with negative `Delta1`, Tabriz represents Cluster 2 as a moderate positive-amplification case, Tehran (Mehrabad Airport) represents Cluster 3 because it shows clear upper-tail intensification, and Shiraz represents Cluster 4 because it is the most extreme upper-tail warm-night outlier in the dataset. The comparison demonstrates that warm-night change spans at least three distinct behaviours: uniformly high warming with muted upper-tail gain, moderate warming with gradual amplification, and extremely strong upper-tail intensification.

For **cool days**, Sanandaj was selected to represent the dominant Cluster 1 because it shows the characteristic monotonic decline toward high quantiles, Zabol represents Cluster 2 because it combines positive `q0.05` behaviour with strongly negative `q0.95`, Orumiyeh represents Cluster 3 because it exhibits the most extreme upper-tail contraction, and Torbat-e Heydariyeh represents Cluster 4 because it is unusual in showing a positive median slope despite strongly negative upper-tail behaviour. This set demonstrates that the decline of cool days is widespread but not structurally uniform: some stations lose cool days most rapidly only in the upper tail, whereas others show a more monotonic decline through the distribution.

For **cool nights**, Khorramabad was selected from Cluster 1 because it exhibits the strongest negative `Delta1`, Qazvin represents the dominant Cluster 2 because it shows the characteristic transition from weak low-quantile change to strong upper-tail decline, Gorgan represents Cluster 3 because it retains positive slopes through much of the distribution, and Shahrekord represents Cluster 4 because it shows positive slopes at all three focal quantiles and stands out as the most positive cool-night case in the network. This comparison is especially instructive because it reveals that cool-night change is not solely a monotonic decline everywhere; a small subset of stations departs strongly from the national majority and therefore justifies separate cluster treatment.

These station comparisons strengthen the interpretation of the clustering results. They show that the clusters are not abstract statistical classes, but correspond to visibly distinct quantile-slope geometries. This is particularly useful in a climatological context because it links the station classification directly to interpretable changes in the lower, central, and upper parts of the thermal-extreme distribution.

### 3.4. Spatial expression and field significance

The mapped results indicate that spatial organization exists, but that it is not strong enough to justify a simple interpretation based on smooth national-scale gradients alone.

![Figure 10. Spatial distribution of warm-day Delta1.](figures/map_warm_days_Delta1.png)

*Figure 10. Station map of `Delta1` for warm days. Positive values indicate stronger upper-tail amplification relative to lower-tail change.*

![Figure 11. Spatial distribution of warm-night Delta1.](figures/map_warm_nights_Delta1.png)

*Figure 11. As in Figure 10, but for warm nights.*

![Figure 12. Spatial distribution of cool-day Delta1.](figures/map_cool_days_Delta1.png)

*Figure 12. As in Figure 10, but for cool days.*

![Figure 13. Spatial distribution of cool-night Delta1.](figures/map_cool_nights_Delta1.png)

*Figure 13. As in Figure 10, but for cool nights.*

After false-discovery-rate control, the number of field-significant station results varies markedly by quantile and index.

**Table 3.** Number of FDR-retained significant station results by index and quantile.

| Quantile | Warm Days | Warm Nights | Cool Days | Cool Nights |
|---|---:|---:|---:|---:|
| 0.10 | 27 | 29 | 0 | 12 |
| 0.50 | 27 | 28 | 17 | 23 |
| 0.90 | 17 | 25 | 26 | 25 |

![Figure 14. Raw and FDR-retained significance counts.](figures/advanced_spatial_inference/raw_vs_fdr_counts.png)

*Figure 14. Number of raw significant and FDR-retained station results across indices and selected quantiles. Intermediate quantiles show the strongest field-significant signal.*

These counts suggest three main features. First, warm indices are highly robust at lower and median quantiles, especially warm nights. Second, cool indices become most strongly significant toward the upper quantiles, consistent with the quantile profiles showing rapid loss of years with historically many cool events. Third, the most extreme tails remain less stable statistically, which is expected because tail estimation is inherently more uncertain in finite samples.

At the same time, Moran's I values are non-significant for all 20 station-slope fields tested. None reached p < 0.05, indicating weak evidence for strong national-scale spatial autocorrelation.

![Figure 15. Moran's I summary for station-slope fields.](figures/advanced_spatial_inference/moran_i_heatmap.png)

*Figure 15. Moran's I values and associated permutation-based significance for station-slope fields across indices and quantiles. The absence of significant field-scale autocorrelation indicates that same-sign trends do not necessarily translate into smooth national spatial structure.*

This is a useful result rather than a limitation. It implies that, although many stations share the same sign of change, the detailed arrangement of slopes is not well represented by a simple, smoothly varying field. In a climatologically diverse country such as Iran, this is plausible and suggests an important role for local physiography, continentality, and land-surface conditions.

### 3.5. Regionalization of quantile signatures

The clustering results reinforce the view that thermal-extreme change across Iran is structured by multiple regimes rather than a single national trajectory.

![Figure 16. Cluster composite heatmap for warm days.](figures/advanced_regionalization/cluster_composite_heatmap_warm_days.png)

*Figure 16. Cluster composite heatmap for warm days, summarizing cluster-level means or medians of the principal quantile metrics.*

![Figure 17. Cluster composite heatmap for warm nights.](figures/advanced_regionalization/cluster_composite_heatmap_warm_nights.png)

*Figure 17. As in Figure 16, but for warm nights.*

![Figure 18. Cluster composite heatmap for cool days.](figures/advanced_regionalization/cluster_composite_heatmap_cool_days.png)

*Figure 18. As in Figure 16, but for cool days.*

![Figure 19. Cluster composite heatmap for cool nights.](figures/advanced_regionalization/cluster_composite_heatmap_cool_nights.png)

*Figure 19. As in Figure 16, but for cool nights.*

For warm days, Cluster 1 has the highest median Delta1 (+13.98), indicating exceptionally strong upper-tail amplification, while Cluster 4 also shows pronounced amplification (median +8.22). Cluster 2 is much weaker (median +2.26), and Cluster 3 is slightly negative (median -2.60), implying that some stations experience more uniform change or even relatively weaker growth in the upper tail. This diversity is consistent with the mixed station-level results noted above.

Warm nights exhibit an even more differentiated clustering structure. One singleton cluster, centred on Shiraz, displays an extremely large Delta1 (+25.90), while another multi-station cluster shows strong and consistent amplification (median approximately +10.85). By contrast, the high-baseline warm-night cluster exhibits substantial warming across all quantiles but more modest within-distribution contrast. This distinction is useful because it separates stations experiencing strong warming overall from stations experiencing especially strong tail reorganization.

For cool days and cool nights, most stations belong to one dominant cluster in each case, indicating a broadly coherent signal of decline. However, the smaller outlier clusters are scientifically informative because they identify locations where the upper-tail contraction of cool events is exceptionally strong or where local behaviour diverges from the national majority. Such departures are valuable for interpretation because they suggest that local climatic controls are modifying the background warming signal.

The regionalization appears sufficiently robust for interpretation. The adjusted Rand index between the configured and reduced-feature clustering runs is 0.632 for warm days, 0.934 for warm nights, 0.685 for cool days, and 0.821 for cool nights. Although exact label identities shift, as expected in unsupervised analysis, the broad cluster structure remains stable enough to support a climatological interpretation.

### 3.6. Geographic drivers and physical interpretation

The geographic driver analysis is deliberately simple, but it offers several useful clues. The clearest and most physically interpretable relationship is the negative association between warm-night slopes and elevation.

**Table 4.** Selected driver relationships supporting geographic interpretation.

| Index-metric | Predictor | Standardized beta | p-value | Spearman rho |
|---|---:|---:|---:|---:|
| Warm Nights, q0.05 | Elevation | -0.507 | 0.006 | -0.521 |
| Warm Nights, q0.50 | Elevation | -0.479 | 0.010 | -0.550 |
| Warm Nights, q0.95 | Elevation | -0.318 | 0.094 | -0.432 |
| Cool Days, q0.50 | Longitude | +0.471 | 0.021 | +0.355 |
| Warm Days, Delta1 | Latitude | -0.281 | 0.185 | -0.378 |

Lower-elevation stations tend to exhibit stronger increases in warm-night frequency, particularly at low and median quantiles. This pattern is physically plausible. Lower-elevation, drier, and often more urbanized environments may suppress nocturnal cooling through greater heat storage and re-radiation, thereby strengthening nighttime warming. The persistence of the signal from q0.05 to q0.50 and its weaker continuation at q0.95 suggests that elevation influences not only the most extreme warm nights, but the broader warm-night regime.

The cool-day results suggest an east-west component as well. The positive standardized beta for longitude in the q0.50 cool-day slope indicates that the strongest median declines tend to occur toward the western part of the network, or conversely that eastern stations show somewhat weaker contraction in median cool-day behaviour. This may reflect differences in continentality, regional circulation influences, or land-surface dryness, although these mechanisms would require targeted follow-up analysis.

The broader implication is that no single mechanism explains all four indices. A more defensible interpretation is that warm nights are especially sensitive to lower-elevation settings, warm days show strong upper-tail amplification at selected interior and southern stations, and cool indices are undergoing broad contraction with local modulation by geography and regional climate context.

### 3.7. Robustness of the principal findings

One of the strongest aspects of the study is that the main conclusions remain stable across multiple sensitivity checks.

![Figure 20. Reference-period sensitivity.](figures/advanced_method_sensitivity/reference_period_sensitivity.png)

*Figure 20. Sensitivity of selected station metrics to alternative reference-period definitions used in threshold construction.*

![Figure 21. Bootstrap-method sensitivity.](figures/advanced_method_sensitivity/bootstrap_method_sensitivity.png)

*Figure 21. Sensitivity of bootstrap-based summary metrics to alternative resampling methods.*

![Figure 22. Interpolation sensitivity at tau = 0.05.](figures/advanced_method_sensitivity/interpolation_sensitivity_tau_0.05.png)

*Figure 22. Agreement among interpolation methods for the station field visualizations at q = 0.05.*

![Figure 23. Interpolation sensitivity at tau = 0.95.](figures/advanced_method_sensitivity/interpolation_sensitivity_tau_0.95.png)

*Figure 23. Agreement among interpolation methods for the station field visualizations at q = 0.95.*

The largest reported differences arise for a limited subset of metrics, including warm-day Delta1 under alternative reference periods and some low-quantile bootstrap summaries for the cool indices. However, these differences affect the magnitude of certain estimates more than their qualitative interpretation. The direction of change remains stable: warm indices increase, cool indices decrease, warm-tail amplification emerges at several stations, and cool-event loss is especially pronounced in upper quantiles.

This distinction is important in a climatological context. Sensitivity analyses are often presented defensively, but here they add substantive value because they help separate highly robust features of the climate signal from results that are more uncertain near the extreme tails. That is exactly the balance required for a strong regional climatology paper with wider relevance.

### 3.8. Synthesis in broader climatological context

The results support three principal conclusions with significance beyond the study area.

First, thermal-extreme change across Iran is distributionally asymmetric. Warm-day and warm-night indices increase broadly, but their rates of increase depend on quantile and station. This means that climate change is not simply shifting the distribution upward; it is altering the relative behaviour of central and tail conditions.

Second, cool extremes are disappearing through a process of upper-quantile contraction. The strongest negative slopes for cool days and cool nights occur in high quantiles, implying that years that were once unusually rich in cool events are declining fastest. This is a more informative climatological statement than a simple decline in average cool-event frequency.

Third, local and regional controls remain important. The weak Moran's I results, marked station-to-station differences in Delta1, and interpretable elevation signal in warm nights all suggest that the thermal-extreme response of Iran reflects a combination of broad-scale warming and geographically specific modifiers. This is likely to be relevant to other arid and topographically complex regions where strong climatic contrasts occur over relatively short distances `(Kostopoulou and Jones, 2005; Zittis et al., 2016)`.

In this sense, the broader contribution of the study lies not merely in documenting warming in Iran, but in showing that the response of temperature extremes is quantile-dependent, regionally differentiated, and robust to a range of methodological decisions. That contribution is of direct relevance to the international climatological literature, particularly for studies seeking to move beyond mean-trend analyses in observational records.

## 4. Conclusions

This study examined quantile-dependent trends in four annual temperature-extreme indices across Iran during 1961-2024 using observations from 30 stations. The results show that the dominant climate signal is one of increasing warm extremes and declining cool extremes, but that the strength and structure of change vary substantially across quantiles and stations.

Three main findings emerge. First, warm days and warm nights increased regionally, with warm nights showing the strongest mean trend. Second, the evolution of these indices is strongly distribution-dependent: warm indices frequently exhibit upper-tail amplification, whereas cool indices show the strongest declines in upper quantiles, indicating rapid loss of years historically characterized by many cool events. Third, the spatial expression of change is heterogeneous. Although field significance is strong for several quantile-index combinations, Moran's I reveals weak national-scale spatial autocorrelation, and cluster analysis identifies multiple regimes rather than a single coherent national pattern.

The study also shows that these conclusions are robust. Sensitivity analyses indicate that the main climatological interpretation is stable across alternative reference periods, bootstrap methods, and interpolation strategies, even though uncertainty increases near the most extreme tails. This gives confidence that the principal findings reflect meaningful climatic structure rather than a narrow modelling choice.

From a wider climatological perspective, the results suggest that thermal-extreme change in arid and topographically complex regions should not be interpreted solely through mean trends. Quantile-aware analysis reveals shifts in the structure of extremes that are invisible in conventional summaries. For Iran, this means that climate change is not only increasing the frequency of warm conditions and reducing cool conditions, but also reorganizing how those conditions are distributed through time and across space.

Several limitations should nevertheless be acknowledged. The station network is finite, the strongest tails remain the most uncertain part of the distribution, and the driver analysis is intentionally simple. In addition, interpolated surfaces should be interpreted primarily as visualization tools rather than as definitive spatial reconstructions. Future work could build on this framework by linking station-scale quantile behaviour to circulation patterns, land-surface processes, urbanization, humidity, and compound heat metrics. Even with these limitations, the present results provide a strong observational basis for understanding how the distribution of thermal extremes is being reshaped across Iran and for situating those changes within the broader international climatological literature.

## Acknowledgement placeholder

Acknowledgements, funding statements, and data-availability statements can be inserted here according to journal requirements.

## References

Alexander LV, Zhang X, Peterson TC, Caesar J, Gleason B, Klein Tank AMG, Haylock M, Collins D, Trewin B, Rahimzadeh F, Tagipour A, Kumar KR, Revadekar J, Griffiths G, Vincent L, Stephenson DB, Burn J, Aguilar E, Brunet M, Taylor M, New M, Zhai P, Rusticucci M, Vazquez-Aguirre JL (2006) Global observed changes in daily climate extremes of temperature and precipitation. *Journal of Geophysical Research: Atmospheres* 111:D05109. https://doi.org/10.1029/2005JD006290

Benjamini Y, Hochberg Y (1995) Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B* 57(1):289-300.

Donat MG, Alexander LV, Yang H, Durre I, Vose R, Caesar J, Gleason BE, Klein Tank AMG, Manton MJ, Molanus C, Peterson TC, Renom M, Brunet M, Caesar J, Zhang X, Kitching S, Haylock M (2013) Updated analyses of temperature and precipitation extreme indices since the beginning of the twentieth century: the HadEX2 dataset. *Journal of Geophysical Research: Atmospheres* 118(5):2098-2118. https://doi.org/10.1002/jgrd.50150

Easterling DR, Meehl GA, Parmesan C, Changnon SA, Karl TR, Mearns LO (2000) Climate extremes: observations, modeling, and impacts. *Science* 289(5487):2068-2074. https://doi.org/10.1126/science.289.5487.2068

Efron B, Tibshirani RJ (1993) *An Introduction to the Bootstrap*. Chapman & Hall, New York.

Frich P, Alexander LV, Della-Marta P, Gleason B, Haylock M, Klein Tank AMG, Peterson T (2002) Observed coherent changes in climatic extremes during the second half of the twentieth century. *Climate Research* 19:193-212. https://doi.org/10.3354/cr019193

Hubert L, Arabie P (1985) Comparing partitions. *Journal of Classification* 2:193-218. https://doi.org/10.1007/BF01908075

IPCC (2021) Weather and climate extreme events in a changing climate. In: Masson-Delmotte V, Zhai P, Pirani A, Connors SL, Péan C, Berger S, Caud N, Chen Y, Goldfarb L, Gomis MI, Huang M, Leitzell K, Lonnoy E, Matthews JBR, Maycock TK, Waterfield T, Yelekçi O, Yu R, Zhou B (eds) *Climate Change 2021: The Physical Science Basis*. Cambridge University Press, Cambridge, Chapter 11. https://www.ipcc.ch/report/ar6/wg1/chapter/chapter-11/

Katz RW, Brown BG (1992) Extreme events in a changing climate: variability is more important than averages. *Climatic Change* 21(3):289-302. https://doi.org/10.1007/BF00139728

Katz RW, Parlange MB, Naveau P (2002) Statistics of extremes in hydrology. *Advances in Water Resources* 25(8-12):1287-1304. https://doi.org/10.1016/S0309-1708(02)00056-8

Koenker R (2005) *Quantile Regression*. Cambridge University Press, Cambridge.

Koenker R, Bassett G Jr (1978) Regression quantiles. *Econometrica* 46(1):33-50. https://doi.org/10.2307/1913643

Kostopoulou E, Jones PD (2005) Assessment of climate extremes in the eastern Mediterranean. *Meteorology and Atmospheric Physics* 89:69-85. https://doi.org/10.1007/s00703-005-0122-2

Moberg A, Jones PD (2005) Trends in indices for extremes in daily temperature and precipitation in central and western Europe, 1901-99. *International Journal of Climatology* 25(9):1149-1171. https://doi.org/10.1002/joc.1163

Moran PAP (1950) Notes on continuous stochastic phenomena. *Biometrika* 37(1-2):17-23. https://doi.org/10.1093/biomet/37.1-2.17

Soltani M, Laux P, Kunstmann H, Stan K, Sohrabi MM, Molanejad M, Sabziparvar AA, Ranjbar SaadatAbadi A, Ranjbar F, Rousta I, Zawar-Reza P, Khoshakhlagh F, Soltanzadeh I, Babu CA, Azizi GH, Martin MV (2016) Assessment of climate variations in temperature and precipitation extreme events over Iran. *Theoretical and Applied Climatology* 126:775-795. https://doi.org/10.1007/s00704-015-1609-5

Vinod HD (2006) Maximum entropy ensembles for time series inference in economics. *Journal of Asian Economics* 17(6):955-978. https://doi.org/10.1016/j.asieco.2006.09.003

Zhang X, Aguilar E, Sensoy S, Melkonyan H, Tagiyeva U, Ahmed N, Kutaladze N, Rahimzadeh F, Taghipour A, Hantosh TH, Albert P, Semawi M, Karam Ali M, Al-Shabibi MHS, Al-Oulan Z, Zatari T, Al Dean Khelet I, Hamoud S, Sagir R, Demircan M, Eken M, Adiguzel M, Alexander LV, Peterson TC, Wallis T (2005) Trends in Middle East climate extreme indices from 1950 to 2003. *Journal of Geophysical Research: Atmospheres* 110:D22104. https://doi.org/10.1029/2005JD006181

Zhang X, Hegerl G, Zwiers FW, Kenyon J (2005) Avoiding inhomogeneity in percentile-based indices of temperature extremes. *Journal of Climate* 18(11):1641-1651. https://doi.org/10.1175/JCLI3366.1

Zhang X, Alexander L, Hegerl GC, Jones P, Klein Tank AMG, Peterson TC, Trewin B, Zwiers FW (2011) Indices for monitoring changes in extremes based on daily temperature and precipitation data. *WIREs Climate Change* 2(6):851-870. https://doi.org/10.1002/wcc.147

Zittis G, Hadjinicolaou P, Lelieveld J (2016) Strongly increasing heat extremes in the Middle East and North Africa (MENA) in the 21st century. *Climatic Change* 137:245-260. https://doi.org/10.1007/s10584-016-1665-6

## Supplementary material pointers

The following generated tables are suitable for citation in the main text or placement in supplementary material:

1. `outputs/tables/publication_summary_table.csv`
2. `outputs/tables/regional_cluster_composites.csv`
3. `outputs/tables/station_significance_fdr.csv`
4. `outputs/tables/driver_analysis_summary.csv`
5. `outputs/tables/reference_period_sensitivity_summary.csv`
6. `outputs/tables/bootstrap_method_sensitivity_summary.csv`
7. `outputs/tables/interpolation_method_sensitivity_summary.csv`
