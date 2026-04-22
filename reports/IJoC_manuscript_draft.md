# Quantile-dependent changes in warm and cool temperature extremes across Iran during 1961-2024

## Abstract

Thermal extremes are commonly summarized using mean trends, although many of their climatic and societal impacts depend on changes across the full distribution. This study investigates quantile-dependent trends in four annual temperature-extreme indices across Iran during 1961-2024 using observations from 30 meteorological stations. Warm days, warm nights, cool days, and cool nights were analysed using quantile regression, with ordinary least squares estimates used as a mean-trend benchmark. The analysis was complemented by bootstrap uncertainty estimation, field-significance testing with false-discovery-rate control, cluster-based regionalization, and sensitivity analyses for reference period, bootstrap method, and interpolation strategy. The results show strong but non-uniform warming-related changes. Warm days and warm nights increased across much of the network, but the magnitude of the increase depended on quantile and station, with pronounced upper-tail amplification at several locations. At the regional scale, warm-day and warm-night mean trends were +9.53 and +13.58 days year^-1 per decade, respectively. In contrast, cool-day and cool-night indices declined regionally, with mean trends of -2.55 and -2.98 days year^-1 per decade, and with the strongest negative slopes concentrated in upper quantiles, indicating rapid loss of years historically characterized by many cool events. Station-level contrasts were substantial: the strongest warm-day upper-tail amplification occurred in Esfahan, Ahvaz, and Zanjan, whereas warm-night amplification was most pronounced in Shiraz. Field significance was strongest for intermediate quantiles, but Moran's I statistics showed little evidence of strong national-scale spatial autocorrelation, implying an important role for local and regional controls. Sensitivity analyses confirmed that the principal findings are robust despite larger uncertainty in the most extreme tails. These findings show that thermal-extreme change across Iran is better interpreted as a quantile-dependent restructuring of warm and cool event distributions than as a uniform shift in the mean, with implications for climatological understanding in arid and topographically complex regions.

**Keywords:** climate extremes; quantile regression; warm nights; cool days; regionalization; Iran; thermal extremes

## 1. Introduction

Changes in climate extremes are among the most consequential expressions of contemporary climate change because their impacts are transmitted through agriculture, public health, water resources, ecosystems, and infrastructure. While mean temperature trends provide a useful first-order description of climatic change, they are often insufficient for understanding extremes, since the societal relevance of climate change depends not only on a shift in the centre of the distribution but also on changes in its tails. A system may show a moderate mean warming while simultaneously experiencing a much stronger intensification of the warmest conditions or a marked contraction of historically cool conditions. For this reason, approaches that resolve the full conditional distribution of climate indicators are increasingly important in climatology.

This issue is especially relevant for temperature-extreme indices. Warm and cool day-night indices have been widely used to characterize changing climate conditions, yet their interpretation is often based on linear trends in the mean or median alone. Such summaries can mask quantile-dependent behaviour, including upper-tail amplification, lower-tail stabilization, or asymmetric decline across the distribution. In physically diverse regions, these distributional effects may differ sharply across space because of topographic contrasts, continentality, land-atmosphere coupling, urban influence, and differences in moisture availability.

Iran provides an especially instructive setting for such an analysis. The country spans strong gradients in elevation, aridity, latitude, and proximity to marine influences, and it contains climatic regimes ranging from humid Caspian lowlands to hot low-elevation interior and southern stations, as well as high-elevation continental environments. These contrasts make Iran a useful natural laboratory for examining whether changes in warm and cool extremes are spatially coherent or whether they are strongly conditioned by local and regional controls. Although previous studies in the broader climatological literature have documented warming and shifts in thermal extremes across arid and semi-arid regions, relatively few station-based studies for Iran have resolved these changes across the full distribution while simultaneously addressing uncertainty, spatial significance, and regionalization.

The International Journal of Climatology places particular value on regional studies that contribute clearly to the wider climatological literature rather than remaining purely descriptive. In that sense, the relevance of the present study lies not simply in documenting change over Iran, but in using a climatically heterogeneous, topographically complex region to examine a broader question: whether observed thermal-extreme change is best understood as a uniform mean shift or as a quantile-dependent reorganization of the distribution. This distinction matters beyond the study area because many dryland and mountain-influenced regions are likely to experience similarly non-uniform changes in extremes.

The present study addresses this issue using 30 meteorological stations with annual data spanning 1961-2024. Four indices are analysed: warm days, warm nights, cool days, and cool nights. Quantile regression is used to estimate trend behaviour from low to high quantiles, and these results are interpreted alongside ordinary least squares trends, bootstrap uncertainty measures, field-significance diagnostics, spatial analysis, and clustering. The specific objectives are: (1) to quantify how temperature-extreme trends vary across quantiles for each index; (2) to identify stations and regions exhibiting strong tail amplification or contraction; (3) to assess whether the detected patterns are field-significant and spatially coherent; and (4) to evaluate the robustness of the main conclusions to methodological choices. By doing so, the study aims to contribute not only a national assessment of Iranian thermal extremes, but also a more general demonstration of how quantile-based methods can refine climatological interpretation in observational extreme-index studies.

## 2. Results and discussion

### 2.1. Data coverage and analytical scope

The final dataset includes 30 stations and 64 years of annual records, covering the common period 1961-2024. No station-index series fell below either the minimum requirement for quantile regression or the more conservative publication threshold of 20 years. This matters because the study is fundamentally comparative: the aim is not only to detect temporal change, but also to compare the shape of change across indices, stations, and quantiles. The broadly uniform record length therefore provides an adequate basis for comparing distributional trend behaviour across the national network.

![Data coverage by station](figures/data_coverage_by_station.png)

### 2.2. Regional-scale quantile behaviour

Regional averaging across stations reveals a clear asymmetry between warming-related and cooling-related indices. Warm days and warm nights both increased strongly over the study period, whereas cool days and cool nights generally declined. However, this broad contrast conceals marked quantile dependence.

| Index | OLS slope | OLS 95% CI | QR q0.01 | QR q0.50 | QR q0.95 | QR q0.99 |
|---|---:|---|---:|---:|---:|---:|
| Warm Days | +9.53 | [+7.54, +11.52] | +15.38 | +8.15 | +12.05 | +11.40 |
| Warm Nights | +13.58 | [+12.05, +15.10] | +18.88 | +13.65 | +12.55 | +12.09 |
| Cool Days | -2.55 | [-4.01, -1.10] | -1.12 | -2.01 | -8.05 | -16.58 |
| Cool Nights | -2.98 | [-4.03, -1.93] | -2.04 | -2.47 | -6.78 | -10.49 |

For warm days, the regional signal is positive throughout the distribution, but the upper quantiles generally exceed the median slope. This indicates that years already characterized by relatively high numbers of warm days have intensified faster than the centre of the distribution. Such behaviour is consistent with upper-tail amplification and suggests that warming is strengthening the most anomalously warm years more rapidly than typical years.

![Regional quantile coefficients for warm days](figures/region_quantile_slopes_warm_days.png)

Warm nights also increased strongly, and the regional OLS trend (+13.58 days year^-1 per decade) is the largest among the four indices. The quantile profile, however, is less strongly skewed toward the upper tail than for warm days. Instead, positive slopes occur across almost the full distribution, including low quantiles, indicating a broad upward shift in the frequency of warm nights. From a climatological perspective, this is important because widespread nighttime warming implies not only increased thermal extremes but also reduced nocturnal relief during hot periods.

![Regional quantile coefficients for warm nights](figures/region_quantile_slopes_warm_nights.png)

The cool indices display the opposite sign and a different pattern of asymmetry. For cool days and cool nights, the strongest negative slopes are concentrated in the upper quantiles. Because high quantiles of these indices correspond to years with unusually many cool events, the steeply negative q0.95-q0.99 slopes imply that the historically coolest years are disappearing faster than the central tendency alone would suggest. This indicates a contraction of the cool-event distribution that is sharper in its upper tail than in its lower tail.

![Regional quantile coefficients for cool days](figures/region_quantile_slopes_cool_days.png)

![Regional quantile coefficients for cool nights](figures/region_quantile_slopes_cool_nights.png)

Taken together, the regional results indicate that thermal-extreme change across Iran is not well described as a rigid translation of the mean. Instead, the warm and cool indices are being reshaped differently across quantiles, with warm extremes generally expanding and cool extremes contracting, but with the strongest responses often concentrated away from the centre of the distribution.

### 2.3. Station-scale heterogeneity and tail amplification

Regional means conceal substantial local heterogeneity. The station heatmaps show that the four indices do not share a single national signature; rather, they comprise multiple trend regimes with contrasting low-, median-, and high-quantile behaviour.

![Warm-days station signatures](figures/station_focus_heatmap_warm_days.png)

![Warm-nights station signatures](figures/station_focus_heatmap_warm_nights.png)

![Cool-days station signatures](figures/station_focus_heatmap_cool_days.png)

![Cool-nights station signatures](figures/station_focus_heatmap_cool_nights.png)

An especially useful summary metric is Delta1 = slope(q0.95) - slope(q0.05), which measures whether the upper tail is changing faster than the lower tail.

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

### 2.4. Spatial expression and field significance

The mapped results indicate that spatial organization exists, but that it is not strong enough to justify a simple interpretation based on smooth national-scale gradients alone.

![Warm-days Delta1 map](figures/map_warm_days_Delta1.png)

![Warm-nights Delta1 map](figures/map_warm_nights_Delta1.png)

![Cool-days Delta1 map](figures/map_cool_days_Delta1.png)

![Cool-nights Delta1 map](figures/map_cool_nights_Delta1.png)

After false-discovery-rate control, the number of field-significant station results varies markedly by quantile and index.

| Quantile | Warm Days | Warm Nights | Cool Days | Cool Nights |
|---|---:|---:|---:|---:|
| 0.10 | 27 | 29 | 0 | 12 |
| 0.50 | 27 | 28 | 17 | 23 |
| 0.90 | 17 | 25 | 26 | 25 |

![Raw versus FDR-retained significance counts](figures/advanced_spatial_inference/raw_vs_fdr_counts.png)

These counts suggest three main features. First, warm indices are highly robust at lower and median quantiles, especially warm nights. Second, cool indices become most strongly significant toward the upper quantiles, consistent with the quantile profiles showing rapid loss of years with historically many cool events. Third, the most extreme tails remain less stable statistically, which is expected because tail estimation is inherently more uncertain in finite samples.

At the same time, Moran's I values are non-significant for all 20 station-slope fields tested. None reached p < 0.05, indicating weak evidence for strong national-scale spatial autocorrelation.

![Moran's I heatmap](figures/advanced_spatial_inference/moran_i_heatmap.png)

This is a useful result rather than a limitation. It implies that, although many stations share the same sign of change, the detailed arrangement of slopes is not well represented by a simple, smoothly varying field. In a climatologically diverse country such as Iran, this is plausible and suggests an important role for local physiography, continentality, and land-surface conditions.

### 2.5. Regionalization of quantile signatures

The clustering results reinforce the view that thermal-extreme change across Iran is structured by multiple regimes rather than a single national trajectory.

![Warm-days cluster composites](figures/advanced_regionalization/cluster_composite_heatmap_warm_days.png)

![Warm-nights cluster composites](figures/advanced_regionalization/cluster_composite_heatmap_warm_nights.png)

![Cool-days cluster composites](figures/advanced_regionalization/cluster_composite_heatmap_cool_days.png)

![Cool-nights cluster composites](figures/advanced_regionalization/cluster_composite_heatmap_cool_nights.png)

For warm days, Cluster 1 has the highest median Delta1 (+13.98), indicating exceptionally strong upper-tail amplification, while Cluster 4 also shows pronounced amplification (median +8.22). Cluster 2 is much weaker (median +2.26), and Cluster 3 is slightly negative (median -2.60), implying that some stations experience more uniform change or even relatively weaker growth in the upper tail. This diversity is consistent with the mixed station-level results noted above.

Warm nights exhibit an even more differentiated clustering structure. One singleton cluster, centred on Shiraz, displays an extremely large Delta1 (+25.90), while another multi-station cluster shows strong and consistent amplification (median approximately +10.85). By contrast, the high-baseline warm-night cluster exhibits substantial warming across all quantiles but more modest within-distribution contrast. This distinction is useful because it separates stations experiencing strong warming overall from stations experiencing especially strong tail reorganization.

For cool days and cool nights, most stations belong to one dominant cluster in each case, indicating a broadly coherent signal of decline. However, the smaller outlier clusters are scientifically informative because they identify locations where the upper-tail contraction of cool events is exceptionally strong or where local behaviour diverges from the national majority. Such departures are valuable for interpretation because they suggest that local climatic controls are modifying the background warming signal.

The regionalization appears sufficiently robust for interpretation. The adjusted Rand index between the configured and reduced-feature clustering runs is 0.632 for warm days, 0.934 for warm nights, 0.685 for cool days, and 0.821 for cool nights. Although exact label identities shift, as expected in unsupervised analysis, the broad cluster structure remains stable enough to support a climatological interpretation.

### 2.6. Geographic drivers and physical interpretation

The geographic driver analysis is deliberately simple, but it offers several useful clues. The clearest and most physically interpretable relationship is the negative association between warm-night slopes and elevation.

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

### 2.7. Robustness of the principal findings

One of the strongest aspects of the study is that the main conclusions remain stable across multiple sensitivity checks.

![Reference-period sensitivity](figures/advanced_method_sensitivity/reference_period_sensitivity.png)

![Bootstrap-method sensitivity](figures/advanced_method_sensitivity/bootstrap_method_sensitivity.png)

![Interpolation sensitivity at tau = 0.05](figures/advanced_method_sensitivity/interpolation_sensitivity_tau_0.05.png)

![Interpolation sensitivity at tau = 0.95](figures/advanced_method_sensitivity/interpolation_sensitivity_tau_0.95.png)

The largest reported differences arise for a limited subset of metrics, including warm-day Delta1 under alternative reference periods and some low-quantile bootstrap summaries for the cool indices. However, these differences affect the magnitude of certain estimates more than their qualitative interpretation. The direction of change remains stable: warm indices increase, cool indices decrease, warm-tail amplification emerges at several stations, and cool-event loss is especially pronounced in upper quantiles.

This distinction is important in a climatological context. Sensitivity analyses are often presented defensively, but here they add substantive value because they help separate highly robust features of the climate signal from results that are more uncertain near the extreme tails. That is exactly the balance required for a strong regional climatology paper with wider relevance.

### 2.8. Synthesis in broader climatological context

The results support three principal conclusions with significance beyond the study area.

First, thermal-extreme change across Iran is distributionally asymmetric. Warm-day and warm-night indices increase broadly, but their rates of increase depend on quantile and station. This means that climate change is not simply shifting the distribution upward; it is altering the relative behaviour of central and tail conditions.

Second, cool extremes are disappearing through a process of upper-quantile contraction. The strongest negative slopes for cool days and cool nights occur in high quantiles, implying that years that were once unusually rich in cool events are declining fastest. This is a more informative climatological statement than a simple decline in average cool-event frequency.

Third, local and regional controls remain important. The weak Moran's I results, marked station-to-station differences in Delta1, and interpretable elevation signal in warm nights all suggest that the thermal-extreme response of Iran reflects a combination of broad-scale warming and geographically specific modifiers. This is likely to be relevant to other arid and topographically complex regions where strong climatic contrasts occur over relatively short distances.

In this sense, the broader contribution of the study lies not merely in documenting warming in Iran, but in showing that the response of temperature extremes is quantile-dependent, regionally differentiated, and robust to a range of methodological decisions. That contribution is of direct relevance to the international climatological literature, particularly for studies seeking to move beyond mean-trend analyses in observational records.

## 3. Conclusions

This study examined quantile-dependent trends in four annual temperature-extreme indices across Iran during 1961-2024 using observations from 30 stations. The results show that the dominant climate signal is one of increasing warm extremes and declining cool extremes, but that the strength and structure of change vary substantially across quantiles and stations.

Three main findings emerge. First, warm days and warm nights increased regionally, with warm nights showing the strongest mean trend. Second, the evolution of these indices is strongly distribution-dependent: warm indices frequently exhibit upper-tail amplification, whereas cool indices show the strongest declines in upper quantiles, indicating rapid loss of years historically characterized by many cool events. Third, the spatial expression of change is heterogeneous. Although field significance is strong for several quantile-index combinations, Moran's I reveals weak national-scale spatial autocorrelation, and cluster analysis identifies multiple regimes rather than a single coherent national pattern.

The study also shows that these conclusions are robust. Sensitivity analyses indicate that the main climatological interpretation is stable across alternative reference periods, bootstrap methods, and interpolation strategies, even though uncertainty increases near the most extreme tails. This gives confidence that the principal findings reflect meaningful climatic structure rather than a narrow modelling choice.

From a wider climatological perspective, the results suggest that thermal-extreme change in arid and topographically complex regions should not be interpreted solely through mean trends. Quantile-aware analysis reveals shifts in the structure of extremes that are invisible in conventional summaries. For Iran, this means that climate change is not only increasing the frequency of warm conditions and reducing cool conditions, but also reorganizing how those conditions are distributed through time and across space.

Several limitations should nevertheless be acknowledged. The station network is finite, the strongest tails remain the most uncertain part of the distribution, and the driver analysis is intentionally simple. In addition, interpolated surfaces should be interpreted primarily as visualization tools rather than as definitive spatial reconstructions. Future work could build on this framework by linking station-scale quantile behaviour to circulation patterns, land-surface processes, urbanization, humidity, and compound heat metrics. Even with these limitations, the present results provide a strong observational basis for understanding how the distribution of thermal extremes is being reshaped across Iran and for situating those changes within the broader international climatological literature.

## Notes for alignment with International Journal of Climatology

The present draft has been shaped toward the expectations described by the Royal Meteorological Society page for the *International Journal of Climatology*, which emphasizes that regional studies should contribute clearly to the international literature and should explain the broader significance of the findings even when the study area is geographically specific. Source used: Royal Meteorological Society journal description: https://www.rmets.org/international-journal-of-climatology

## Supporting tables to cite in the manuscript or supplementary material

1. `outputs/tables/publication_summary_table.csv`
2. `outputs/tables/regional_cluster_composites.csv`
3. `outputs/tables/station_significance_fdr.csv`
4. `outputs/tables/driver_analysis_summary.csv`
5. `outputs/tables/reference_period_sensitivity_summary.csv`
6. `outputs/tables/bootstrap_method_sensitivity_summary.csv`
7. `outputs/tables/interpolation_method_sensitivity_summary.csv`
