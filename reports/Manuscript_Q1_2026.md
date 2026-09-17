# Distributional thermal change and the components of increasing dry–hot concurrence across Iran, 1991–2024

**Running title:** Thermal extremes and dry–hot concurrence in Iran

## Abstract

Changes in thermal extremes and compound dry–hot conditions need not reflect the same statistical processes. We analyzed daily observations from 124 Iranian stations during 1991–2024 to connect quantile-dependent trends in annual temperature-extreme counts with changes in explicitly defined dry–hot concurrence. Quantile regression described warm-day, warm-night, cool-day, and cool-night distributions. A complementary analysis classified years using fixed 1991–2007 precipitation and temperature quartiles and partitioned joint-frequency changes into dry-frequency, hot-frequency, and excess-joint components. Warm-day trends in the annual network-mean series increased from 9.45 to 22.50 days decade⁻¹ between the 0.10 and 0.90 quantiles; the corresponding warm-night trends were 8.74 and 12.21 days decade⁻¹. Fixed thermal thresholds yielded 36.78 additional warm days and 28.41 additional warm nights per year in 2008–2024. Among stations meeting coverage requirements in every year, annual dry–hot frequency increased from 9.16% to 29.81% (104 stations), while June–September frequency increased from 5.54% to 18.39% (103 stations). The summer increase of 12.85 percentage points partitioned into 2.89 from dry-frequency change, 9.91 from hot-frequency change, and 0.06 from excess joint occurrence. Synchronized year-block resampling supported positive joint-frequency changes but did not resolve an excess-joint change. These results connect distributionally uneven thermal change with increasing compound concurrence without implying strengthened physical coupling. Inference remains conditional on the observed network, short baseline, precipitation-tie treatment, and absence of metadata-supported daily homogenization.

**Keywords:** quantile regression; temperature extremes; compound events; fixed thresholds; block bootstrap; Iran

## 1. Introduction

Climate change can alter both the center and shape of a temperature distribution, so a mean trend alone does not describe the behavior of extremes (Katz and Brown, 1992). Changes in the annual frequency of warm and cool days are particularly useful observational indicators because they translate daily temperature departures into interpretable counts (Alexander et al., 2006; Zhang et al., 2011; Dunn et al., 2020). Nevertheless, their annual distributions may change unevenly: years containing many warm days can exhibit a different trend from years containing relatively few. Distinguishing these responses requires an explicit distributional estimand rather than an assumption that the whole distribution follows its mean.

Quantile regression provides that distinction by estimating conditional trends across the response distribution (Koenker and Bassett, 1978; Koenker, 2005). Barbosa et al. (2011) combined quantile trends, bootstrap uncertainty, and clustering for European temperatures, while Fan (2014) examined quantile-dependent annual temperature-extreme counts in China. These studies establish the methodological precedent for the present thermal analysis. Our question concerns the regional expression of those distributional changes and their relationship to a separately defined compound hydroclimatic signal.

Concurrent low precipitation and high temperature create a second interpretive challenge. Their joint frequency depends on the frequencies of the individual conditions and on their concurrence. Greater joint occurrence can therefore accompany more frequent heat without requiring stronger dependence between temperature and precipitation (Zscheischler and Seneviratne, 2017). Alizadeh et al. (2020) used a century-scale United States archive to investigate joint rarity, spatial extent, and changing event composition. Bevacqua et al. (2022) highlighted the importance of precipitation changes and sampling uncertainty for compound hot–dry occurrence, and Schmutz et al. (2026) examined marginal and dependence contributions to emerging European events. Separating these statistical components is therefore an established scientific question, rather than a new method introduced here.

Event definition matters when applying that question to a short archive. An observation may have a small empirical joint survival probability because one variable is unusual, even when the other is not beyond a marginal extreme threshold. A joint-rarity classification and an explicit low-precipitation/high-temperature AND classification consequently describe different sets of years. In arid environments, the distinction is complicated further by zero precipitation and tied observations. Transparent definitions are necessary before changes in joint frequency can be interpreted as changes in concurrent dryness and heat.

Iran provides a heterogeneous setting for evaluating these distinctions. Its station network spans Caspian lowlands, the Alborz and Zagros mountains, interior plateaus, deserts, and southern coasts. Previous work documented regional changes in temperature and precipitation extremes (Zhang et al., 2005a; Soltani et al., 2016). More recently, Ghasemi (2026) combined PCA and quantile regression for monthly and seasonal Iranian temperatures, while Jamali et al. (2026) examined changes in station climate classifications. Quantile analysis and climatic regionalization therefore already have direct Iranian precedents. The remaining question addressed here is how changes in annual thermal-event count distributions align with explicit dry–hot concurrence and its marginal-frequency components within a common climate-regime framework.

We address three objectives: quantify distributional thermal change and its climate-regime contrasts across the 124-station archive; determine whether explicit fixed-baseline dry–hot concurrence increased on stable observing networks; and partition the period difference into changes associated with marginal event frequencies and excess joint occurrence, both nationally and within climate regimes. We assess sensitivity to daily screening, thresholds, data coverage, baseline uncertainty, block length, and precipitation ties. The contribution is a reproducible observational synthesis with clearly separated estimands. Neither the partition nor agreement among the diagnostics constitutes physical attribution.

## 2. Data and methods

### 2.1. Observations, coverage, and quality diagnostics

The archive contains daily minimum, maximum, and mean temperature and precipitation from 124 Iranian stations during 1991–2024 (Fig. 1), supplied by the Iran Meteorological Organization (IRIMO). Coordinates and elevation are supplied in the station metadata. The original measurement units, acquisition procedure, and redistribution conditions require author confirmation: **[UNITS, ACCESS AND LICENSE REQUIRED]**.

![Station network](../outputs/publication_v2/figures/fig01_station_network.png)

*Figure 1. Locations of the 124 stations in the 1991–2024 thermal archive. Color denotes station elevation. Coordinates are displayed in WGS84 longitude and latitude with latitude-adjusted aspect. Boundaries provide geographical context; symbols represent individual observations rather than a spatially continuous sample. Compound analyses use the balanced subsets described in Section 2.4.*

There were 1,539,956 station-day records and no duplicate station dates. Median completeness was 98.96% for minimum temperature, 99.42% for maximum temperature, 99.35% for mean temperature, and 98.94% for precipitation. Internal checks found 45 records at 32 stations with minimum temperature exceeding maximum temperature. After excluding those overlapping conflicts, a further 2,165 records at 61 stations had mean temperature outside the minimum–maximum interval. Historical thermal estimates use the supplied observations, with a separate consistency-screened sensitivity. The new compound analysis masks all three temperatures on minimum–maximum conflict days and masks only mean temperature for the second category. Negative precipitation, if present, is treated as missing; source files are never changed.

Pettitt, standard normal homogeneity, and Buishand-type tests applied to annual mean temperature flagged 119 raw and 39 detrended station series under at least one diagnostic (Pettitt, 1979; Buishand, 1982; Alexandersson, 1986). These screens do not constitute daily homogenization. We retain historical exclusion sensitivity for the 39 detrended flags, because station histories and reference-network information are insufficient to support automatic adjustments.

### 2.2. Annual thermal indices and baseline sensitivity

Day-of-year 10th- and 90th-percentile thresholds were estimated from the 1991–2024 record using an 11-day circular window, comprising each calendar day and five days on either side. Leap days were omitted. When a local reference sample contained fewer than 15 observations, the station-wide reference sample was used. Warm days and nights are strict exceedances of the upper maximum- and minimum-temperature thresholds; cool days and nights are strict departures below the corresponding lower thresholds. Annual values are counts of observed qualifying days, set to missing below 80% valid coverage. They are not annualized to compensate for missing days, so residual coverage differences remain a limitation.

A second set of daily thresholds was estimated using 1991–2007 alone and applied unchanged to both periods. We compared station mean annual counts in 2008–2024 with those in 1991–2007. This tests sensitivity to inclusion of recent observations in the reference distribution; it does not provide an independent climatological baseline. Sampling differences inside and outside a percentile-estimation period can introduce artificial discontinuities, for which Zhang et al. (2005b) proposed an in-base correction. That correction is not implemented here, and re-estimating thresholds during uncertainty resampling is not equivalent to correcting the original index series.

### 2.3. Quantile trends and spatial summaries

For annual count \(Y\) and time \(t\), measured in decades, we fitted

$$Q_Y(\tau\mid t)=\beta_0(\tau)+\beta_1(\tau)t.$$

The quantile coefficient minimizes the sum of check losses \(\rho_\tau(u)=u[\tau-\mathbf{1}(u<0)]\). Models were fitted at quantiles 0.10–0.90 in steps of 0.01, with 0.10, 0.50, and 0.90 as focal summaries. Ordinary least squares (OLS) provides a mean-trend comparison. We distinguish regression of the annual network-mean series from averaging individual-station regression coefficients; these operations generally yield different values. The asymmetry metric is

$$\Delta_1=\hat\beta_1(0.90)-\hat\beta_1(0.10).$$

It describes differences in slopes across annual count quantiles, not a change in daily temperature intensity. For declining cool indices, a negative value denotes a more negative slope at the upper count quantile.

Historical station uncertainty used 200 moving-block bootstrap replicates, with \(\lceil n^{1/3}\rceil\) years per block, bounded between two and eight; the 34-year series therefore used four-year blocks. A 400-replicate calculation and a maximum-entropy alternative provide sensitivity checks (Künsch, 1989; Vinod, 2006). Tail intervals are pointwise percentile intervals. At the median, analytic station probabilities were adjusted separately within each index using Benjamini–Hochberg FDR at 0.05 (Benjamini and Hochberg, 1995). These retained counts are local tests after multiplicity adjustment; they are not a separate test of field significance, and FDR does not repair misspecification of the underlying probabilities.

Station maps display the observed estimates without interpolation. Historical Moran diagnostics used five-nearest-neighbor weights and 499 label permutations; their nominal probabilities are exploratory across the 12 fields. Hierarchical clustering was examined through feature and algorithm sensitivities. Climate-regime stratification, described in Section 2.7, uses a separately established classification rather than clusters fitted to the same trends being compared.

### 2.4. Explicit dry–hot events on a balanced network

The compound analysis uses annual precipitation totals with annual mean temperature, and June–September precipitation totals with the seasonal mean of daily maximum temperature. Coverage is computed against the number of calendar days, including missing daily rows, separately for precipitation and temperature. A year requires at least 80% coverage for both variables. To keep the comparison independent of changes in station membership, a station must meet these conditions in all 34 years; this yields 104 annual and 103 warm-season stations. Precipitation totals are not scaled for missing days. A 90% coverage sensitivity tests the influence of this choice, without assuming that missing precipitation is zero.

For station \(s\), define dry and hot indicators

$$D_{sy}=\mathbf{1}(P_{sy}<q^{P}_{s,0.25}),\qquad H_{sy}=\mathbf{1}(T_{sy}>q^{T}_{s,0.75}),\qquad J_{sy}=D_{sy}H_{sy},$$

where both thresholds use the 17 baseline years, 1991–2007, with linearly interpolated sample quantiles. Quartiles provide a less sparse primary classification than more extreme cutoffs in this short baseline; these events are relative dry–hot conditions, not rare design hazards. Fixed thresholds are applied to both periods. Period frequencies are first calculated at each station and then averaged with equal station weights. Annual network extent is the percentage of those same stations with \(J=1\).

Strict inequalities prevent equality to a tied threshold from automatically qualifying as an extreme. In particular, a zero lower precipitation threshold permits no strictly drier year. We retain these stations in the primary network summary but report their number and examine both inclusive inequalities and exclusion of zero-threshold stations. Additional sensitivities use 20th/80th and 30th/70th percentiles, stricter coverage, and the intersection of annual and warm-season station sets.

### 2.5. Exact partition of the joint-frequency difference

Let \(d_k\), \(h_k\), and \(j_k\) denote a station's dry, hot, and joint frequencies in period \(k\), with \(k=0\) for 1991–2007 and \(k=1\) for 2008–2024. Write \(c_k=j_k-d_kh_k\), the excess joint frequency relative to the product of the two marginal frequencies. The exact identity

$$j_1-j_0=\underbrace{(d_1-d_0)\frac{h_1+h_0}{2}}_{A_D}+\underbrace{(h_1-h_0)\frac{d_1+d_0}{2}}_{A_H}+\underbrace{(c_1-c_0)}_{A_C}$$

allocates the change in the marginal-frequency product symmetrically between dry and hot frequencies. The three terms sum exactly to the joint-frequency difference at every station and in their equally weighted network mean. Multiplication by 100 expresses the terms in percentage points.

This is an accounting identity, not a causal model. In particular, \(c\) is covariance between binary threshold indicators; its change does not isolate copula change. Indicator covariance may vary with marginal frequencies even under an unchanged copula. We therefore call \(A_C\) the *excess-joint term*, and interpret \(A_D\) and \(A_H\) as frequency components rather than fractions of externally forced change.

### 2.6. Resampling, multiplicity, and reproducibility

For the new partition we generated 4,999 circular moving-block bootstrap replicates, separately within each 17-year period. Four-year blocks were sampled with replacement and truncated to the original period length. Within each period, the same sampled years were used for every station and both variables, retaining observed within-year spatial and cross-variable dependence. The two periods were resampled independently. Baseline thresholds were re-estimated in every replicate before evaluating both period samples, so the intervals include variation in the estimated thresholds. A conditional-threshold sensitivity holds the original cutoffs fixed; two- and six-year blocks assess the resampling choice.

We report 95% percentile intervals and wider, nominal Bonferroni intervals using quantiles 0.003125 and 0.996875 for eight primary quantities: four partition terms under two event definitions. These are approximate intervals from short, discretized samples, not exact simultaneous coverage guarantees. Sensitivity scenarios are robustness checks rather than additional discoveries. No local significance symbols are placed on the new station maps. Resampling assumes that blocks represent variability within each period and does not model secular nonstationarity inside the periods or dependence across their boundary.

The extension is controlled by `publication_config.yaml` and reproduced with `python run_publication.py`. Input hashes, package versions, station membership, thresholds, bootstrap network draws, summary tables, and vector figures are archived under `outputs/publication_v2`. `python validate_publication.py` checks the probability identity and calendar coverage and rebuilds both historical thermal index sets from daily observations. Historical bootstrap and clustering outputs retain their original provenance; their complete pipeline was not rerun for this extension.

### 2.7. Climate classification and stratified analysis

We assigned station coordinates to the 1-km present-climate Köppen–Geiger raster of Beck et al. (2018), representing 1980–2016. The six reporting groups are BWh hot desert (36 stations), BWk cold desert (18), BSh hot steppe (10), BSk cold steppe (35), temperate Csa/Cfa (15), and Dsa cold dry-summer (10). The temperate group combines 13 Csa and two Cfa stations for reporting; this aggregation does not imply climatic equivalence. One station required a nearest-valid raster cell, as recorded with classification confidence in the station table. These are fixed climate labels, not an analysis of climate-type transitions. The newer map collection of Beck et al. (2023) provides a relevant alternative baseline, but its classes were not substituted for the archived 2018 assignments.

For each group we summarized station-specific thermal slopes and fixed-baseline changes. Historical label-permutation tests compared the maximum minus minimum group mean across six regimes, using 999 permutations and FDR across 32 prespecified index–metric combinations. Their exchangeability assumption is limited by spatial dependence; we therefore use these results as supporting diagnostics. Of the 32 combinations, 16 concern time trends, four fixed-baseline changes, and 12 responses to a regional temperature anomaly. These are omnibus tests of between-group variation, not 32 pairwise comparisons.

For the compound analysis, climate groups are intersected with the balanced annual or summer station set. We calculate group means of each probability-partition term and obtain exploratory 95% within-group intervals from the same synchronized 4,999 bootstrap fields used for the network analysis. Group-specific replication therefore preserves cross-station and cross-group covariance. These intervals are not multiplicity-adjusted between-regime tests, and their overlap or separation is not used to declare one climate regime more affected than another. Counts of zero dry thresholds are reported alongside the summer results.

## 3. Results

### 3.1. Thermal change differs across count quantiles

Warm indices increased and cool indices declined, with the largest lower-to-upper slope contrasts in the daytime indices (Fig. 2; Table 1). In the annual network-mean warm-day series, slopes rose from 9.45 days decade⁻¹ at the 0.10 quantile to 15.55 at the median and 22.50 at the 0.90 quantile. Warm-night slopes were also positive, but their median and upper-quantile values were similar, at 12.14 and 12.21 days decade⁻¹. Thus, the strongest change in warm-day counts occurred toward the upper part of the conditional annual distribution, whereas warm nights exhibited a smaller slope contrast.

![Quantile profiles](../outputs/publication_v2/figures/fig02_quantile_profiles.png)

*Figure 2. Quantile slopes for annual network-mean counts (solid lines) and network-mean OLS slopes (dashed lines). Shading is the interquartile range of station-specific slopes at each quantile, representing spatial heterogeneity rather than a confidence interval for the solid line. The response is an annual count of threshold-exceeding days, not daily temperature intensity. Network focal slopes are independently recomputed from the archived annual counts.*

**Table 1. Network-mean thermal trends and station-level median FDR screening.**

| Index | OLS | q0.10 | q0.50 | q0.90 | Δ₁ | Median FDR-retained stations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Warm days | 14.02 | 9.45 | 15.55 | 22.50 | 13.06 | 115/124 |
| Warm nights | 10.96 | 8.74 | 12.14 | 12.21 | 3.47 | 104/124 |
| Cool days | −12.60 | −5.60 | −12.00 | −15.37 | −9.77 | 97/124 |
| Cool nights | −10.86 | −5.91 | −9.83 | −12.06 | −6.16 | 86/124 |

*All trend coefficients and Δ₁ are in days decade⁻¹. Δ₁ is computed before rounding. FDR applies separately across stations for each index at the median; these counts do not represent a joint field-significance test. Network slopes and station inferential counts describe different estimands.*

Cool-day slopes became more negative from −5.60 to −15.37 days decade⁻¹ between the lower and upper quantiles; cool-night slopes changed from −5.91 to −12.06. Negative cool-day slopes occurred at every station at the median and upper quantile. Nevertheless, direction and precision were different: at the upper quantile, historical pointwise bootstrap intervals excluded zero for 73 warm-day, 88 warm-night, 57 cool-day, and 61 cool-night stations. Widespread signs therefore did not imply precise local tail estimates everywhere.

### 3.2. Spatial heterogeneity and fixed-baseline consistency

Station asymmetry maps show broad warming-consistent directions together with substantial variation in magnitude (Fig. 3). Positive warm-day Δ₁ occurred at 118 stations and negative cool-day Δ₁ at all 124. The maps retain individual-station values; they do not imply a continuous trend surface. Historical Moran diagnostics gave the most consistent spatial organization for cool-day slopes, while clustering was sensitive enough to remain a supplementary description.

![Thermal asymmetry maps](../outputs/publication_v2/figures/fig03_thermal_asymmetry_maps.png)

*Figure 3. Station-specific Δ₁ for the four annual thermal indices. All panels use the same zero-centered color scale in days decade⁻¹. Positive values indicate a larger upper-quantile slope; negative values indicate a more negative upper-quantile slope. Colors encode point estimates without local significance claims or spatial interpolation.*

Applying fixed 1991–2007 daily thresholds reproduced the main directions (Fig. 4). Network mean annual counts changed by +36.78 warm days, +28.41 warm nights, −10.69 cool days, and −8.19 cool nights in 2008–2024. The expected direction occurred at 124, 117, 124, and 110 stations, respectively. Climate-regime differences are examined in Section 3.3 and Figure 5. Supplementary Figures S1–S5 document quality, robustness and spatial diagnostics, while Figures S6–S8 describe exploratory clustering and representative station profiles.

![Fixed-baseline thermal changes](../outputs/publication_v2/figures/fig04_fixed_baseline_changes.png)

*Figure 4. Station mean annual-count differences between 2008–2024 and 1991–2007 using fixed daily thresholds estimated in 1991–2007. Dots are station differences, diamonds are equal-station means, and horizontal segments and central marks show station interquartile ranges and medians. These segments describe between-station spread, not uncertainty in a network mean.*

Historical checks support these broad thermal patterns. Excluding 39 detrended-homogeneity-flagged stations changed network Δ₁ from 13.06 to 13.05 days decade⁻¹ for warm days and from 3.47 to 3.45 for warm nights; cool-day and cool-night contrasts changed from −9.77 to −10.39 and from −6.16 to −4.87. Masking inconsistent temperatures changed the focal network slopes by at most 0.026 days decade⁻¹. Bootstrap means were stable when increasing replication from 200 to 400, although upper-quantile results were more sensitive to the resampling method. Stability of signs therefore provides firmer evidence than fine precision in every tail coefficient.

### 3.3. Climate regimes organize the spatial contrasts

All six climate groups had positive mean upper-quantile warm-day and warm-night slopes and negative mean cool-day and cool-night slopes (Fig. 5; Table 2). The strongest mean upper-quantile warm-day trends occurred in BSh hot steppe (23.85 days decade⁻¹) and Dsa cold dry-summer settings (22.23), whereas BWh hot desert had a smaller warm-day slope (15.84) but the largest warm-night slope (16.64). Thus, the network-wide daytime contrast did not describe every climate setting: hot-desert stations exhibited comparatively strong nighttime change.

![Climate classification and thermal response](../outputs/publication_v2/figures/fig05_climate_regimes.png)

*Figure 5. (a) Station membership in six fixed Köppen–Geiger reporting groups, using the Beck et al. (2018) 1980–2016 climate raster. The 15-station temperate group pools Csa and Cfa for reporting. (b) Group means of station-specific 0.90-quantile slopes for annual thermal-event counts, with a common zero-centered scale. These are descriptive means, not slopes fitted to group-average series, classification-change estimates, or proof of between-group significance.*

**Table 2. Climate-regime composition and mean station upper-quantile thermal slopes.**

<!-- REGIME_THERMAL_TABLE -->
| Regime | N | Mean elevation (m) | Warm days | Warm nights | Cool days | Cool nights |
| --- | --- | --- | --- | --- | --- | --- |
| BWh hot desert | 36 | 575.05 | 15.84 | 16.64 | -20.59 | -18.36 |
| BWk cold desert | 18 | 1401.71 | 16.75 | 11.62 | -15.28 | -13.48 |
| BSh hot steppe | 10 | 495.39 | 23.85 | 12.36 | -16.85 | -15.34 |
| BSk cold steppe | 35 | 1442.13 | 21.62 | 12.72 | -17.56 | -12.82 |
| Csa/Cfa temperate | 15 | 671.57 | 20.44 | 13.64 | -18.37 | -17.09 |
| Dsa cold dry-summer | 10 | 1917.11 | 22.23 | 11.99 | -19.62 | -14.09 |
<!-- END_REGIME_THERMAL_TABLE -->

*Elevation is the group mean in meters; all four slopes are in days decade⁻¹. Station counts refer to the complete thermal network. Compound analyses use smaller balanced subsets (Table 4). These summaries do not identify the climatic cause of a trend.*

BWh stations also had the most negative mean upper-quantile cool-day slope (−20.59 days decade⁻¹), compared with −15.28 in BWk cold desert. These contrasts describe the data but are less decisive inferentially: no time-trend or fixed-baseline regime statistic survived FDR across the historical 32-metric family. The two retained statistics (q = 0.048 each) concerned thermal-index responses to the internally derived regional temperature anomaly, not trends per decade. Consequently, the climate classification improves geographical interpretation without establishing six statistically distinct trend populations.

### 3.4. Explicit dry–hot concurrence increased on stable station sets

With both marginal thresholds required, annual dry–hot frequency increased from 9.16% to 29.81% between periods, a change of 20.64 percentage points across 104 stations (95% interval 9.22–40.84). June–September frequency increased from 5.54% to 18.39%, or 12.85 percentage points across 103 stations (7.20–27.18; Fig. 6). The wider primary-family intervals also remained above zero. These frequencies describe equal-weight station networks, not percentages of Iranian land area.

![Compound network extent](../outputs/publication_v2/figures/fig06_fixed_compound_extent.png)

*Figure 6. Percentage of balanced-network stations classified as dry, hot, or simultaneously dry AND hot in each year. Annual and June–September definitions use separate station sets and their own fixed 1991–2007 precipitation and temperature quartiles. Shading identifies the baseline period. Strict inequalities are used; joint-event extent cannot exceed either marginal extent. Counts are based on screened daily observations and calendar-day coverage.*

Summer dryness classification was constrained by precipitation zeros: 28 of the 103 stations had a baseline lower quartile of zero. Such stations cannot satisfy the strict dry inequality and contribute zero joint frequency in the primary classification. No annual station had a zero lower quartile. Sensitivities changing the tie convention and station inclusion are therefore necessary for interpreting the seasonal difference; the primary annual and summer magnitudes should not be read as a controlled comparison of identical event populations.

### 3.5. The hot-frequency term accounts for most of the summer point change

The exact partition attributes the summer point difference to +2.89 percentage points in the dry-frequency term, +9.91 in the hot-frequency term, and +0.06 in the excess-joint term (Fig. 7; Table 3). The corresponding annual components were +8.18, +12.48, and −0.02 percentage points. These components sum to the observed joint changes before rounding. The hot-frequency term had the larger point estimate under both definitions, but the overlapping component intervals do not by themselves constitute a paired test that it exceeds the dry-frequency term.

![Compound frequency partition](../outputs/publication_v2/figures/fig07_compound_partition.png)

*Figure 7. Partition of the late-minus-early joint-frequency difference. Points are equal-station estimates; thick intervals are 95% percentile intervals and thin intervals are nominal Bonferroni intervals across eight primary quantities. The 4,999 synchronized circular year-block replicates re-estimate baseline thresholds. The three colored component estimates sum exactly to the total change. The excess-joint term is a covariance difference between binary indicators, not a causal or copula-isolated contribution.*

**Table 3. Compound-frequency differences and partition terms, in percentage points.**

| Definition | Stations | Joint change [95% interval] | Dry-frequency term [95% interval] | Hot-frequency term [95% interval] | Excess-joint term [95% interval] |
| --- | ---: | --- | --- | --- | --- |
| Annual | 104 | 20.64 [9.22, 40.84] | 8.18 [4.18, 18.61] | 12.48 [6.16, 23.44] | −0.02 [−4.35, 3.42] |
| June–September | 103 | 12.85 [7.20, 27.18] | 2.89 [0.53, 10.43] | 9.91 [5.79, 16.29] | 0.06 [−1.79, 2.37] |

*Intervals include baseline-threshold re-estimation and synchronized resampling within each period. The wider nominal primary-family intervals are archived in `compound_partition_primary.csv` and drawn in Fig. 7. In particular, the summer dry-frequency term does not exclude zero under that wider interval. No component is interpreted as a physical attribution fraction.*

The excess-joint intervals included zero under both definitions, providing no resolved network-wide change in this term. This differs from claiming that temperature and precipitation are independent, or that their dependence structure is unchanged. Station components varied spatially (Fig. 8), while the network average could conceal compensating local patterns. These short station records do not justify local significance claims from inspection of the mapped colors.

![Spatial partition](../outputs/publication_v2/figures/fig08_compound_partition_maps.png)

*Figure 8. June–September station estimates of the joint-frequency difference and its three partition terms for the 103-station balanced network. All panels share a symmetric color scale in percentage points. Gray crosses identify the 28 stations with zero baseline precipitation thresholds, which cannot satisfy the strict dry definition; their zero-valued contributions are retained in network means. Colored circles show the remaining station estimates without interpolation or significance symbols. Threshold values and flags are supplied in the station table.*

### 3.6. Compound-frequency changes within climate regimes

The annual joint-frequency point difference was positive in every climate group, ranging from 16.86 percentage points in the temperate group to 36.13 in BSh hot steppe (Fig. 9). Summer point differences ranged from 5.20 in BWh hot desert to 21.32 in Dsa cold dry-summer settings; the temperate and BSk cold-steppe groups had changes of 20.78 and 16.22 percentage points, respectively (Table 4). These patterns extend the climate-regime interpretation from individual thermal indices to explicit compound concurrence.

![Compound frequency changes by climate regime](../outputs/publication_v2/figures/fig09_compound_climate_regimes.png)

*Figure 9. Annual and June–September joint-frequency changes by climate regime. Points are equally weighted within-group station estimates; intervals are exploratory 95% percentile intervals from synchronized primary bootstrap fields with threshold re-estimation. The displayed n is the balanced sample for that definition and group. Intervals describe each group separately and are not adjusted tests comparing groups. Summer estimates include stations with zero lower precipitation quartiles.*

**Table 4. Summer compound-frequency change and its components by climate regime.**

<!-- REGIME_COMPOUND_TABLE -->
| Regime | N | N₀ | Joint change [95% interval] | Dry-frequency term | Hot-frequency term | Excess-joint term |
| --- | --- | --- | --- | --- | --- | --- |
| BWh hot desert | 26 | 17 | 5.20 [2.71, 14.03] | 1.87 | 3.85 | -0.52 |
| BWk cold desert | 13 | 2 | 9.50 [2.69, 23.53] | 1.06 | 8.54 | -0.11 |
| BSh hot steppe | 8 | 5 | 5.88 [2.21, 16.91] | 1.08 | 5.06 | -0.26 |
| BSk cold steppe | 33 | 3 | 16.22 [7.13, 34.22] | 2.78 | 12.46 | 0.99 |
| Csa/Cfa temperate | 15 | 1 | 20.78 [6.27, 44.71] | 6.32 | 14.74 | -0.28 |
| Dsa cold dry-summer | 8 | 0 | 21.32 [10.29, 45.59] | 4.95 | 17.06 | -0.69 |
<!-- END_REGIME_COMPOUND_TABLE -->

*All changes and components are percentage points. N₀ denotes stations with zero baseline precipitation quartiles. The joint-change interval is an exploratory 95% within-group interval. Component point estimates sum to the joint change before rounding; complete intervals for every component and both seasons are in the supplementary tables.*

The hot-frequency point term exceeded the dry-frequency term in all six summer groups, while every group-specific excess-joint interval included zero. However, lower summer changes in hot-arid groups partly reflect the strict precipitation definition: zero dry thresholds occurred at 17 of 26 balanced BWh stations and five of eight BSh stations, compared with none of the eight Dsa stations. The smaller BWh estimate therefore cannot be interpreted as evidence that its stations face less absolute heat exposure or drought hazard. Within-group samples of eight stations and broad intervals further limit ranking the regimes.

### 3.7. Definition and resampling sensitivity

All examined scenarios retained positive network joint-frequency point changes under both definitions (Fig. 10). Using the 20th/80th or 30th/70th percentile thresholds changed the annual estimate to 17.76 or 20.36 percentage points and the summer estimate to 10.97 or 14.45. Raising daily coverage to 90% retained 86 annual and 101 summer stations and yielded changes of 22.09 and 12.93 percentage points. Two- and six-year blocks altered interval endpoints without changing the observed point estimates.

Threshold uncertainty materially widened some intervals. For example, holding thresholds fixed gave a summer interval of 6.17–19.59 percentage points, compared with 7.20–27.18 when thresholds were re-estimated. Inclusive treatment of threshold ties changed the summer point estimate to 17.36 percentage points. Excluding zero-threshold sites retained 75 summer stations and yielded a change of 17.65 percentage points (9.73–33.80). On the common 101-station network, annual and summer changes were 20.62 and 13.10 percentage points. These alternatives demonstrate sensitivity of absolute magnitude to definition and station population; they are not independent replications of one estimand.

![Sensitivity of joint-frequency changes](../outputs/publication_v2/figures/fig10_partition_sensitivity.png)

*Figure 10. Joint-frequency differences and 95% synchronized-block intervals under ten analytical scenarios. Point estimates coincide for alternative block lengths because the event definition is unchanged. Alternative thresholds, coverage, tie rules and station selection change the estimand. These are sensitivity checks rather than independent confirmatory tests.*

## 4. Discussion

### 4.1. Thermal changes in relation to earlier observational studies

The thermal and compound analyses describe related but distinct aspects of recent change. The thermal indices show that annual warm-day counts shifted unevenly across their conditional distribution, with a larger positive slope at the upper than lower quantile. This pattern is consistent with the motivation for quantile analysis in Barbosa et al. (2011) and Fan (2014): reporting only OLS or median trends would conceal changes among high- and low-count years. It does not mean that the annual quantile itself has become more probable; a fixed conditional quantile retains its probability definition while its count value changes with time.

Fan (2014) offers a particularly useful comparison because both studies model annual extreme-temperature counts. Across 549 Chinese stations during 1960–2008, the reported national warm-night slopes were 4.63, 8.15, and 8.83 days decade⁻¹ at the 0.10, 0.50, and 0.90 quantiles. Our corresponding estimates are 8.74, 12.14, and 12.21 days decade⁻¹ (Table 1). Both profiles indicate a larger shift among high-count than low-count years, but the larger Iranian coefficients cannot be assigned to regional climate sensitivity without harmonizing periods, percentile construction, station weighting, and preprocessing. Fan used homogenized records, whereas our daily archive has not undergone metadata-supported homogenization. Barbosa et al. (2011), in contrast, analyzed daily mean temperature itself; their quantile slopes in temperature units provide a conceptual precedent, not a numerical benchmark for count trends.

Within Iran, Soltani et al. (2016) examined 50 stations over 1975–2010 and reported warm-day and warm-night increases of 12 and 14 days decade⁻¹ for the 1995–2010 subperiod, alongside decreases of 4 and 3 days decade⁻¹ in cold days and nights. Our 1991–2024 OLS estimates are 14.02, 10.96, −12.60, and −10.86 days decade⁻¹, respectively. The common directions support consistency in the broad thermal signal, while the day–night ordering and cooling-index magnitudes differ. These are comparisons across different periods and index implementations, not formal evidence of acceleration or a reversal of physical controls.

The more recent Iranian analysis of Ghasemi (2026) used 78 stations with monthly records spanning 1955–2024 and reported warming of 0.54 °C decade⁻¹ at the cold minimum-temperature tail and 0.34 °C decade⁻¹ at the hot maximum-temperature tail. That stronger cold-tail temperature change does not contradict our larger warm-day count contrast: monthly temperature quantiles and quantiles of annual threshold-exceedance counts measure different responses. A small temperature shift can produce a large count change where many daily observations lie near the threshold. This precedent also narrows our contribution: applying quantile regression in Iran is already established; linking annual count distributions, explicit compound occurrence, and climate-stratified component estimates is the focus here.

**Table 5. Comparison with relevant observational and compound-event studies.**

| Study | Evidence or principal finding | Relationship to this study and comparison limit |
| --- | --- | --- |
| Barbosa et al. (2011) | Central European daily-temperature quantile trends and spatial clustering | Supports distributional analysis; temperature slopes cannot be compared numerically with annual-count slopes. |
| Fan (2014) | Chinese annual warm-night counts show larger upper- than lower-quantile increases | Similar profile shape; record length, homogenization and baseline differ. Numerical comparison is given in Section 4.1. |
| Soltani et al. (2016) | Earlier Iranian warm-event increases and cold-event decreases | Signs agree; published subperiod trends are not estimates for our 1991–2024 window. |
| Ghasemi (2026) | Iranian monthly-temperature tails change asymmetrically | Direct national QR precedent; temperature intensity and annual event frequency remain distinct estimands. |
| Alizadeh et al. (2020) | Increasing US compound dry–hot extent and changing event composition | Directionally consistent; joint-rarity classes, seasonal windows and much longer records prevent pooling magnitudes. |
| Bevacqua et al. (2022) | Precipitation trends shape projected compound hot–dry occurrence | A future-model result is compatible with a large observed hot-frequency term; thresholds and climate state matter. |
| Schmutz et al. (2026) | Heat-index changes dominate increasing compound occurrence in much of Europe, with regional exceptions | Qualitatively consistent with our larger hot-frequency point term; our remainder is not a full dependence decomposition. |
| Jamali et al. (2026) | Iranian station climate classifications show transitions toward drier types | Motivates climate context; our fixed classes do not estimate climate-zone migration. |

*Comparisons concern published estimands and reported findings. They are not a meta-analysis or a statistical test of differences between studies. Sources and bibliographic verification are documented in the accompanying reference audit.*

Warm-night trends were also widespread, but their smaller quantile contrast indicates a different distributional expression from warm days. Temperature minima and maxima respond to radiation, moisture availability, boundary-layer processes, and circulation in different ways. Those processes provide hypotheses for subsequent work, not mechanisms identified by these regressions. The absence of humidity, soil moisture, radiation, or circulation observations prevents a direct explanation of the daytime–nighttime difference. Annual event counts also do not resolve the duration, intensity, or persistence of individual heatwaves.

The fixed-threshold compound analysis adds a hydroclimatic condition: high seasonal or annual temperature must coincide with precipitation below a baseline cutoff. Its increase on balanced station sets indicates that changing data availability alone does not account for the estimated expansion. Yet the magnitudes remain conditional on the reference climate and observing locations. Relative quartile events are useful for comparing occurrence through time at each site; they are not equivalent to an absolute heat-stress threshold, a soil-moisture drought, or a societal impact.

### 4.2. More joint events do not require a resolved excess-joint increase

The summer partition identifies a large positive hot-frequency term and an excess-joint term close to zero at network scale. The interpretation is statistical: an increase in the marginal frequency of hot seasons can contribute substantially to more frequent dry–hot seasons without a detectable increase in excess concurrence beyond the marginal product. This is compatible with established compound-event theory, which distinguishes marginal changes from the dependence that joins them (Zscheischler and Seneviratne, 2017). The present binary partition is deliberately simpler than the full distributional decomposition in Schmutz et al. (2026); it does not estimate a copula or separate dependence changes from all marginal effects on indicator covariance.

The annual dry-frequency term was also positive, and the summer estimate should not lead to the general claim that precipitation is unimportant. Bevacqua et al. (2022) showed why precipitation change and internal variability can be central to compound-event occurrence. Differences in time scale, baseline, station population, and threshold position matter when comparing such findings. Here, the annual analysis integrates precipitation through the year, whereas the summer classification encounters frequent zero totals at the lower threshold. A smaller summer dry term can therefore reflect the statistical definition as well as the behavior of precipitation.

Alizadeh et al. (2020) provides the closest supplied precedent for examining changes in compound-event extent and composition. However, their much longer archive supports a different temporal perspective. The present 34-year record cannot supply stable estimates of rare return periods or distinguish recent change from all modes of multidecadal variability. The historical joint-rarity outputs are retained as complementary descriptive evidence, while the revised main analysis asks the more explicit question of how often both marginal thresholds were crossed. An apparent difference in relative annual and summer expansion between the two definitions is therefore informative about the estimand, rather than evidence that one set of calculations must be numerically wrong.

Specifically, Alizadeh et al. used United States climate divisions during 1896–2017, with October–September annual and March–August warm-season windows. Their increasing compound extent and greater heat contribution are qualitatively consistent with our growing joint occurrence and large summer hot-frequency term. However, their empirical joint-rarity classes are not our fixed marginal AND events, and their warm season is not June–September. A numerical ratio between the reported changes would combine differences in temporal aggregation, spatial support, record length, and event definition. The comparison supports the scientific relevance of separating event composition, while leaving the Iranian magnitudes conditional on this study's design.

### 4.3. What climate classification adds to the interpretation

The Köppen–Geiger grouping connects the station results to recognizable temperature–moisture environments. It is useful even when a formal difference between groups is not resolved: the joint-change estimates, sample sizes, and zero-precipitation thresholds can be inspected together rather than hidden inside one national mean. In summer, the relatively small BWh and BSh joint changes coexist with many structurally zero dry-event classifications. Conversely, larger point changes in temperate and Dsa groups describe departures from their own reference climates, not larger absolute aridity. Climate classification thus exposes a threshold-related interpretation problem as well as spatial variation.

Jamali et al. (2026) analyzed 279 Iranian stations grouped by record availability and reported shifts toward arid types at 27% of long-record and 18% of mid-record stations. Their time-varying classifications address climate transitions, whereas our Beck et al. (2018) assignments hold the geographical stratification fixed. The results are complementary, but our maps cannot establish that stations crossed climate boundaries or attribute the observed compound-frequency components to aridification. Beck et al. (2023) provides updated historical classification periods that could support a future assignment-sensitivity analysis; that newer product was not substituted for the raster actually used here.

The current classification also has statistical limits. None of the 16 historical time-trend group comparisons survived the specified FDR family, and the new compound intervals are within-group uncertainty summaries rather than adjusted contrasts between groups. Sparse classes and spatial dependence limit the power and interpretation of regime rankings. The defensible result is heterogeneous point estimates within a common warming-consistent pattern, not a set of statistically distinct climate responses. Keeping the classification in the main analysis makes those limits visible and provides a reproducible basis for testing them with longer records.

### 4.4. Uncertainty and spatial interpretation

Resampling whole annual station fields is important because treating stations as independent replicates would overstate the information in a spatially coherent climate signal. Our synchronized blocks preserve the observed spatial arrangement within sampled years, and re-estimation of baseline thresholds acknowledges that those cutoffs were estimated from only 17 observations. The resulting intervals are appreciably wider than some conditional-threshold intervals. Their asymmetry also reflects discrete event counts and threshold variability; reporting only a bootstrap standard deviation would obscure this behavior.

The bootstrap nevertheless has limits. Four-year blocks leave few effective blocks in each period, and circular resampling treats the period boundary as adjacent for sampling purposes. Within-period trends challenge the approximate stationarity assumption. The Bonferroni intervals address multiplicity in the primary family only to the extent that the marginal bootstrap intervals have adequate coverage; they cannot correct the short record or make sensitivity analyses confirmatory. Failure to resolve an excess-joint change is consequently absence of clear evidence at this scale, rather than proof of stable dependence.

Point maps and equal-station summaries also describe a monitoring network rather than the whole land surface. Dense station clusters have more influence than sparsely sampled terrain of equal area. Spatially opposing components can cancel in network averages, while a small set of local colors cannot establish coherent regional change. The maps therefore use common scales and avoid interpolated surfaces or significance markers. Area-weighted products, independent gridded observations, and explicit spatial models would be needed to estimate physical affected area and formally evaluate regional contrasts.

### 4.5. Scope, transferability, and next steps

The framework is transferable where paired temperature and precipitation observations exist: define the marginal thresholds, retain a stable sampling network, calculate station-first probabilities, and partition their difference with synchronized uncertainty estimation. Its scientific value lies in preventing a broad statement about increasing joint occurrence from being mistaken for evidence of increasing coupling. The reproducible outputs also make alternative threshold choices inspectable, which is particularly important in arid regions where zero precipitation affects event classification.

Several limitations remain before the work can support stronger conclusions. Daily observations have not undergone metadata-supported homogenization, precipitation totals remain sensitive to nonrandom missing days, and the baseline is shorter than a conventional 30-year climatology. The balanced subsets improve comparability but exclude stations with poorer coverage, potentially changing geographical representation. Existing thermal median probabilities are analytic and their FDR screening does not fully resolve temporal and spatial dependence concerns. The present extension strengthens compound-frequency interpretation without retroactively solving every inferential limitation of the historical thermal analysis.

The most useful next extensions would add independent evidence: longer homogenized records, gridded data for area-weighted comparison, or soil-moisture and evaporative-demand measurements for process-oriented definitions. A copula-based counterfactual decomposition would address a different and more demanding dependence question, requiring adequate sample size and model diagnostics. Daily heatwave duration and simultaneous warm-day/warm-night persistence would also complement the annual-frequency estimands. These analyses should be evaluated against explicit questions rather than added simply to increase methodological complexity.

## 5. Conclusions

This study connected distributional changes in annual thermal-event counts with explicitly defined dry–hot concurrence across the Iranian observing network during 1991–2024. Its contribution is the joint interpretation of thermal quantile profiles, fixed-baseline compound probabilities, and climate-stratified component estimates within a reproducible workflow. These complementary diagnostics distinguish changes among high-count thermal years from changes in the frequency with which low precipitation and high temperature occur together.

The thermal response was warming-consistent but uneven across the count distribution. Warm-day slopes increased from 9.45 days decade⁻¹ at the lower quantile to 22.50 at the upper quantile, whereas warm-night slopes ranged from 8.74 to 12.21. A mean-only summary would therefore miss much of the daytime distributional contrast. Fixed-baseline comparisons also showed 36.78 additional warm days and 28.41 additional warm nights per year in the later period. These are changes in annual event frequency, with separate implications from the intensity or duration of individual heatwaves.

Explicit dry–hot occurrence increased on balanced station networks. Annual frequency rose from 9.16% to 29.81%, and June–September frequency from 5.54% to 18.39%. The summer increase of 12.85 percentage points had a 95% synchronized-block interval of 7.20–27.18 points. Its statistical partition comprised a 9.91-point hot-frequency term, a 2.89-point dry-frequency term, and a 0.06-point excess-joint term whose interval included zero. Increasing compound occurrence can therefore be described without claiming a detected strengthening of dependence or physical coupling. The positive annual dry-frequency term also prevents generalizing the summer balance into a claim that precipitation change is unimportant.

Climate classification adds spatial and methodological context. All six groups had positive joint-frequency point differences, while the summer changes and prevalence of zero dry thresholds varied substantially among groups. The larger temperate and Dsa estimates describe relative changes against local baselines, not greater absolute hazard than in desert environments. Group contrasts remain exploratory: historical thermal differences did not survive the specified time-trend multiplicity screening, and the new intervals do not formally test differences between climate groups. The classification is consequently useful for identifying where event definitions and sample composition affect interpretation.

The direction of network compound-frequency change persisted across the examined threshold, coverage, tie-treatment, station-selection, and resampling alternatives, although its magnitude and interval width varied. This supports the directional finding within the evaluated designs while emphasizing that definitions are part of the result. Inference remains restricted by the 34-year record, 17-year reference period, unequal spatial sampling, precipitation missingness, and absence of metadata-supported daily homogenization. Longer homogenized observations, independent spatial products, and soil-moisture or evaporative-demand measurements would be the most informative next steps. They would enable stronger tests of regional contrasts and process explanations while retaining the central distinction between more frequent joint events and changes in their dependence structure.

## Data and code availability

Daily observations were supplied by the Iran Meteorological Organization (IRIMO). **[AUTHOR TO CONFIRM ACCESS CONDITIONS AND REDISTRIBUTION LICENSE.]** The reproducible implementation comprises the original thermal pipeline and the separate publication extension. **[INSERT PERMANENT CODE/DATA ARCHIVE DOI OR PUBLIC REPOSITORY RELEASE.]** Input hashes and package versions are included with the extension outputs.

## Author declarations

**[INSERT AUTHORS AND AFFILIATIONS, CONTRIBUTIONS, FUNDING, COMPETING INTERESTS, AND ACKNOWLEDGEMENTS.]**

## Supplementary material

The [supplementary index](Supplementary_Q1_2026.md) provides complete component intervals, all sensitivity scenarios, station thresholds and membership, bootstrap draws, eleven supporting figures, eight tables, validation scope, and a catalog of retained numerical evidence. Supplementary Figures S9–S11 provide explicitly exploratory warming associations, historical joint-rarity sensitivity and geographical associations.

## References

Alexander LV, Zhang X, Peterson TC, Caesar J, Gleason B, Klein Tank AMG, Haylock M, Collins D, Trewin B, Rahimzadeh F, Tagipour A, Rupa Kumar K, Revadekar J, Griffiths G, Vincent L, Stephenson DB, Burn J, Aguilar E, Brunet M, Taylor M, New M, Zhai P, Rusticucci M, Vazquez-Aguirre JL (2006). Global observed changes in daily climate extremes of temperature and precipitation. *Journal of Geophysical Research: Atmospheres* 111:D05109. https://doi.org/10.1029/2005JD006290

Alexandersson H (1986). A homogeneity test applied to precipitation data. *Journal of Climatology* 6:661–675. https://doi.org/10.1002/joc.3370060607

Alizadeh MR, Adamowski J, Nikoo MR, AghaKouchak A, Dennison P, Sadegh M (2020). A century of observations reveals increasing likelihood of continental-scale compound dry-hot extremes. *Science Advances* 6:eaaz4571. https://doi.org/10.1126/sciadv.aaz4571

Barbosa SM, Scotto MG, Alonso AM (2011). Summarising changes in air temperature over Central Europe by quantile regression and clustering. *Natural Hazards and Earth System Sciences* 11:3227–3233. https://doi.org/10.5194/nhess-11-3227-2011

Beck HE, Zimmermann NE, McVicar TR, Vergopolan N, Berg A, Wood EF (2018). Present and future Köppen–Geiger climate classification maps at 1-km resolution. *Scientific Data* 5:180214. https://doi.org/10.1038/sdata.2018.214

Beck HE, McVicar TR, Vergopolan N, Berg A, Lutsko NJ, Dufour A, Zeng Z, Jiang X, van Dijk AIJM, Miralles DG (2023). High-resolution (1 km) Köppen-Geiger maps for 1901–2099 based on constrained CMIP6 projections. *Scientific Data* 10:724. https://doi.org/10.1038/s41597-023-02549-6

Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B* 57:289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

Bevacqua E, Zappa G, Lehner F, Zscheischler J (2022). Precipitation trends determine future occurrences of compound hot–dry events. *Nature Climate Change* 12:350–355. https://doi.org/10.1038/s41558-022-01309-5

Buishand TA (1982). Some methods for testing the homogeneity of rainfall records. *Journal of Hydrology* 58:11–27. https://doi.org/10.1016/0022-1694(82)90066-X

Dunn RJH, Alexander LV, Donat MG, Zhang X, Bador M, Herold N, Lippmann T, Allan R, Aguilar E, Barry AA, Brunet M, Caesar J, Chagnaud G, Cheng V, Cinco T, Durre I, de Guzman R, Htay TM, Wan Ibadullah WM, Bin Ibrahim MKI, Khoshkam M, Kruger A, Kubota H, Leng TW, Lim G, Li‐Sha L, Marengo J, Mbatha S, McGree S, Menne M, de los Milagros Skansi M, Ngwenya S, Nkrumah F, Oonariya C, Pabon‐Caicedo JD, Panthou G, Pham C, Rahimzadeh F, Ramos A, Salgado E, Salinger J, Sané Y, Sopaheluwakan A, Srivastava A, Sun Y, Timbal B, Trachow N, Trewin B, van der Schrier G, Vazquez‐Aguirre J, Vasquez R, Villarroel C, Vincent L, Vischel T, Vose R, Bin Hj Yussof MN (2020). Development of an updated global land in situ-based dataset of temperature and precipitation extremes: HadEX3. *Journal of Geophysical Research: Atmospheres* 125:e2019JD032263. https://doi.org/10.1029/2019JD032263

Fan LJ (2014). Quantile trends in temperature extremes in China. *Atmospheric and Oceanic Science Letters* 7:304–308. https://doi.org/10.3878/j.issn.1674-2834.13.0102

Ghasemi AR (2026). Analyzing spatiotemporal patterns of extreme temperatures in Iran using principal component analysis and quantile regression. *Earth and Space Science* 13:e2025EA004860. https://doi.org/10.1029/2025EA004860

Jamali M, Eslamian S, Shayannejad M, Gohari A (2026). Observed warming–driven aridification and climate-type transitions across Iran. *Journal of Arid Environments* 235:105606. https://doi.org/10.1016/j.jaridenv.2026.105606

Katz RW, Brown BG (1992). Extreme events in a changing climate: variability is more important than averages. *Climatic Change* 21:289–302. https://doi.org/10.1007/BF00139728

Koenker R, Bassett G Jr (1978). Regression quantiles. *Econometrica* 46:33–50. https://doi.org/10.2307/1913643

Koenker R (2005). *Quantile Regression*. Cambridge University Press, Cambridge. https://doi.org/10.1017/CBO9780511754098

Künsch HR (1989). The jackknife and the bootstrap for general stationary observations. *The Annals of Statistics* 17:1217–1241. https://doi.org/10.1214/aos/1176347265

Pettitt AN (1979). A non-parametric approach to the change-point problem. *Applied Statistics* 28:126–135. https://doi.org/10.2307/2346729

Schmutz J, Vrac M, François B, Bulut B (2026). Spatial structures of emerging hot and dry compound events over Europe from 1950 to 2023. *Natural Hazards and Earth System Sciences* 26:881–900. https://doi.org/10.5194/nhess-26-881-2026

Soltani M, Laux P, Kunstmann H, Stan K, Sohrabi MM, Molanejad M, Sabziparvar AA, Ranjbar SaadatAbadi A, Ranjbar F, Rousta I, Zawar-Reza P, Khoshakhlagh F, Soltanzadeh I, Babu CA, Azizi GH, Martin MV (2016). Assessment of climate variations in temperature and precipitation extreme events over Iran. *Theoretical and Applied Climatology* 126:775–795. https://doi.org/10.1007/s00704-015-1609-5

Vinod HD (2006). Maximum entropy ensembles for time series inference in economics. *Journal of Asian Economics* 17:955–978. https://doi.org/10.1016/j.asieco.2006.09.001

Zhang X, Aguilar E, Sensoy S, Melkonyan H, Tagiyeva U, Ahmed N, Kutaladze N, Rahimzadeh F, Taghipour A, Hantosh TH, Albert P, Semawi M, Karam Ali M, Al-Shabibi MHS, Al-Oulan Z, Zatari T, Al Dean Khelet I, Hamoud S, Sagir R, Demircan M, Eken M, Adiguzel M, Alexander LV, Peterson TC, Wallis T (2005a). Trends in Middle East climate extreme indices from 1950 to 2003. *Journal of Geophysical Research: Atmospheres* 110:D22104. https://doi.org/10.1029/2005JD006181

Zhang X, Hegerl G, Zwiers FW, Kenyon J (2005b). Avoiding inhomogeneity in percentile-based indices of temperature extremes. *Journal of Climate* 18:1641–1651. https://doi.org/10.1175/JCLI3366.1

Zhang X, Alexander L, Hegerl GC, Jones P, Klein Tank AMG, Peterson TC, Trewin B, Zwiers FW (2011). Indices for monitoring changes in extremes based on daily temperature and precipitation data. *WIREs Climate Change* 2:851–870. https://doi.org/10.1002/wcc.147

Zscheischler J, Seneviratne SI (2017). Dependence of drivers affects risks associated with compound events. *Science Advances* 3:e1700263. https://doi.org/10.1126/sciadv.1700263
