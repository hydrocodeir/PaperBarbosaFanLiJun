# Quantile-dependent changes in thermal extremes and compound dry–hot conditions across Iran, 1991–2024

**Running title:** Thermal and dry–hot changes across Iran

## Abstract

Changes in mean temperature do not fully characterize thermal extremes, particularly when trends differ between ordinary and extreme years or when heat occurs concurrently with precipitation deficit. Using daily observations from 124 Iranian meteorological stations during 1991–2024, we evaluated annual warm-day, warm-night, cool-day, and cool-night counts with station-level quantile regression, moving-block bootstrap intervals, false-discovery-rate control, spatial diagnostics, and Köppen–Geiger stratification. An empirical joint-probability analysis quantified annual and June–September dry–hot conditions. Our results show that warm days underwent the strongest distributional change, with the annual network-mean slope increasing from 9.45 days decade⁻¹ at the 0.10 quantile to 22.50 days decade⁻¹ at the 0.90 quantile. By contrast, warm-night slopes were positive but less asymmetric (8.74 and 12.21 days decade⁻¹, respectively), while cool-day and cool-night slopes reached −15.37 and −12.06 days decade⁻¹ at the 0.90 quantile. At the median, false-discovery-rate control retained 115 warm-day, 104 warm-night, 97 cool-day, and 86 cool-night station trends. Thresholds derived from 1991–2007 produced late-period changes of +36.78 warm days, +28.41 warm nights, −10.69 cool days, and −8.19 cool nights yr⁻¹. Warm-season dry–hot events with empirical joint return periods of at least 10 years affected an average of 10.35 stations yr⁻¹ in 1991–2007 and 41.53 stations yr⁻¹ in 2008–2024. The affected-station trend remained positive under a 4-year moving-block null (Kendall τ = 0.583, p = 0.0002), while event composition shifted from predominantly dry-dominant to predominantly hot-dominant. The findings indicate widespread but distributionally and climatically heterogeneous warming of thermal extremes across Iran, accompanied by an expanding warm-season dry–hot signal; however, the short record, absence of metadata-supported daily homogenization, and empirical definition of rarity limit inference about long-return-period hazards and external forcing.

**Keywords:** temperature extremes; quantile regression; compound events; empirical joint probability; climate regimes; Iran

## 1. Introduction

Changes in climate risk depend on the behavior of extremes as well as averages. A shift in mean temperature can coincide with changes in variability, skewness, or tail behavior, so the frequency of unusually warm or cool conditions need not change uniformly across the distribution (Katz and Brown, 1992; Katz et al., 2002). Observational studies have documented increasing warm extremes and declining cold extremes across many regions, while event-attribution and projection studies show that continued warming increases the likelihood of high-temperature records (Alexander et al., 2006; Seneviratne et al., 2014; Fischer and Knutti, 2015; Fischer et al., 2021). Quantile regression is suited to this problem because it estimates changes in different parts of a response distribution rather than reducing the temporal signal to a change in its mean (Koenker and Bassett, 1978; Koenker, 2005; Reich, 2012).

Percentile-based counts of warm days, warm nights, cool days, and cool nights provide an interpretable measure of changes in daily temperature occurrence (Frich et al., 2002; Zhang et al., 2011; Dunn et al., 2020). Their annual distributions can nevertheless change unevenly. A trend at the upper quantile of annual warm-day counts, for example, describes years in which warm days are already frequent and may differ substantially from the trend in years with relatively few warm days. Analyses confined to means or medians can therefore miss amplification or contraction in the tails (Barbosa et al., 2011; Fan, 2014). This distinction is relevant in arid and topographically complex regions, where land–atmosphere coupling, circulation, elevation, and maritime influence can produce marked spatial differences in hot extremes (Seneviratne et al., 2010; McKinnon et al., 2021; Rousi et al., 2022; Vautard et al., 2023).

Heat is also more consequential when it occurs with limited moisture supply. Dry and hot conditions can reinforce one another through reduced evaporative cooling and increased sensible-heat flux, although the strength and direction of this coupling vary among regions and events (Leonard et al., 2014; Zscheischler and Seneviratne, 2017; Zscheischler et al., 2018; Maraun et al., 2025). Compound-event studies consequently treat precipitation and temperature as jointly distributed drivers rather than independent hazards (Hao et al., 2013; Sarhadi et al., 2018; Alizadeh et al., 2020). The spatial extent of joint events is particularly relevant because simultaneous impacts at many locations can strain regional response capacity. Alizadeh et al. (2020), for example, found that compound dry–hot extremes in the contiguous United States became more frequent and spatially connected over a century of observations, with heat becoming a more prominent event component in recent decades.

Iran spans humid Caspian lowlands, the Alborz and Zagros mountain systems, interior plateaus and deserts, and hot coasts bordering the Persian Gulf and Gulf of Oman. Previous studies have reported increasing temperature extremes in Iran and the wider Middle East and North Africa, but estimates vary with index definition, record length, spatial scale, and data source (Zhang et al., 2005a; Soltani et al., 2016; Vaghefi et al., 2019; Zittis et al., 2016; Francis and Fonseca, 2024; Naderi et al., 2024). Regional studies also identify growing heat and humid-heat hazards in southwest Asia (Pal and Eltahir, 2016; Raymond et al., 2020, 2024). Less is known about whether recent station trends across Iran differ systematically between lower and upper quantiles, whether those differences are consistent across climate regimes, and whether the thermal signal is accompanied by a change in the network-wide occurrence of dry–hot years.

Here, we address these questions using daily records from 124 stations for 1991–2024, with three objectives: to quantify trends across the annual distributions of four percentile-based temperature indices, to assess their spatial organization and sensitivity to uncertainty, threshold, and data-quality choices, and to determine whether annual and warm-season dry–hot conditions have become more widespread across the station network. Because the study is observational, agreement among the consistency checks strengthens the interpretation of recent change but does not constitute detection and attribution of anthropogenic forcing.

## 2. Data and methods

### 2.1. Station observations and quality screening

We analyzed daily minimum, maximum, and mean temperature together with daily precipitation at 124 Iranian meteorological stations from 1991 through 2024 (Fig. 1), providing 34 years of observations across coastal, lowland, plateau, and mountain environments. Station coordinates and elevations supported mapping and physiographic comparisons; however, the formal data provider, acquisition procedure, original units, and access conditions must be specified before submission: **[DATA SOURCE AND ACCESS STATEMENT REQUIRED]**.

![Figure 1. Study region and station network.](../outputs/figures/ijoc_study_area_regional_context_new.png)

*Figure 1. Iranian station network used in the analysis. Station color denotes elevation; neighboring countries and adjacent seas provide regional context.*

Median daily completeness was 98.96% for minimum temperature, 99.42% for maximum temperature, 99.35% for mean temperature, and 98.94% for precipitation, with no duplicate station dates. Across the 1,539,956 station-day records, missing values numbered 23,135 for minimum temperature, 13,944 for maximum temperature, 20,750 for mean temperature, and 16,347 for precipitation. Quality screening identified 45 records at 32 stations for which minimum temperature exceeded maximum temperature and, after these overlapping conflicts were excluded, 2,165 records at 61 stations for which mean temperature fell outside the interval defined by the daily minimum and maximum. Because the source archive lacks the metadata required to correct individual observations, we retained the supplied values in the primary analysis but repeated the regional thermal and compound analyses after masking inconsistent values, thereby testing their influence without modifying any source record.

We applied Pettitt, standard normal homogeneity, and Buishand-type tests to annual mean temperature as diagnostic screens (Pettitt, 1979; Buishand, 1982; Alexandersson, 1986). At least one test flagged 119 raw annual series, whereas only 39 detrended series were flagged; because breakpoint tests can respond to genuine trend as well as non-climatic discontinuity, and because the archive lacks station histories and a reference-network homogenization design, we did not impose automatic daily adjustments. We instead repeated regional summaries after excluding every station flagged by a detrended test and omitted leap days before calculating day-of-year thresholds.

### 2.2. Annual temperature-extreme indices

For each station, we calculated day-of-year thresholds from the 1991–2024 record using observations within ±5 days of each calendar day, which produced an 11-day moving window with circular treatment at the beginning and end of the year. The 10th and 90th percentiles of maximum temperature defined cool- and warm-day thresholds, while the corresponding minimum-temperature percentiles defined cool- and warm-night thresholds; a station-wide reference sample replaced a local window only when fewer than 15 valid observations were available. We then summed daily exceedance indicators by year and set an annual value to missing when fewer than 80% of the relevant daily observations were valid.

Because the primary thresholds use the complete analysis period, the resulting indices describe redistribution relative to the 1991–2024 empirical climate rather than exceedance of an independent historical normal. We therefore estimated a second set of thresholds from 1991–2007, applied them unchanged to all years, and compared mean annual counts in 2008–2024 with those in the baseline period, which tests whether the direction of change persists when late-period observations do not contribute to threshold estimation.

### 2.3. Quantile trends and uncertainty

For annual index value \(Y\) and time \(T\), the conditional quantile model was

$$
Q_Y(\tau \mid T)=\beta_0(\tau)+\beta_1(\tau)T,
$$

where time was expressed in decades. We fitted station-level models from \(\tau=0.10\) to \(0.90\) in increments of 0.01 and focused inference on \(\tau=0.10\), 0.50, and 0.90, while ordinary least-squares slopes provided a mean-trend comparison. Network-mean annual series summarized national-scale behavior, although station-level models remained the basis for inferential counts.

Distributional asymmetry was summarized as

$$
\Delta_1=\hat\beta_1(0.90)-\hat\beta_1(0.10).
$$

For a warm index, positive \(\Delta_1\) indicates that high-count years increased faster than low-count years, whereas negative \(\Delta_1\) for a declining cool index indicates stronger contraction at the upper quantile.

We estimated station-level uncertainty with a moving-block bootstrap using 200 primary replicates and block length \(\lceil n^{1/3}\rceil\), bounded between 2 and 8 years, which yielded 4-year blocks for the 34-year series. A 400-replicate rerun assessed Monte Carlo stability, while maximum-entropy bootstrap results provided a method sensitivity check (Kunsch, 1989; Hall et al., 1995; Lahiri, 2003; Vinod, 2006). Because analytic tail intervals were incomplete, we used bootstrap percentile intervals at the tail quantiles; at \(\tau=0.50\), analytic station-level probabilities were adjusted separately for each index using the Benjamini–Hochberg false discovery rate (FDR) at 0.05 (Benjamini and Hochberg, 1995).

The quantile estimates minimize the asymmetric absolute-loss function

$$
\hat{\boldsymbol{\beta}}(\tau)=\arg\min_{\boldsymbol{\beta}}\sum_{i=1}^{n}\rho_{\tau}\!\left(y_i-\mathbf{x}_i^{\mathsf{T}}\boldsymbol{\beta}\right),\qquad
\rho_{\tau}(u)=u\left[\tau-\mathbf{1}(u<0)\right],
$$

where \(\mathbf{x}_i=(1,T_i)^{\mathsf{T}}\). This makes explicit that lower and upper quantiles are estimated from asymmetric losses rather than from separate transformations of a mean model.

### 2.4. Spatial and climate-regime analyses

We measured spatial autocorrelation in station slopes with Moran's \(I\), using five nearest neighbors and 499 label permutations (Moran, 1950), but treated the resulting maps as displays of station values rather than spatially continuous inference. Hierarchical clustering provided an exploratory description of similar quantile profiles, with stability evaluated through reduced features, alternative linkage and distance choices, and the adjusted Rand index (Hubert and Arabie, 1985). Because several partitions were sensitive to feature selection or lacked geographic compactness, detailed cluster assignments and representative-station panels are reported in the Supplementary Material and do not constitute primary evidence.

We assigned station locations to the 1-km present-day Köppen–Geiger classification of Beck et al. (2018), retaining six groups with sufficient stations for summary: hot desert (BWh; n = 36), cold desert (BWk; n = 18), hot steppe (BSh; n = 10), cold steppe (BSk; n = 35), temperate Csa/Cfa (n = 15), and cold dry-summer Dsa (n = 10). The two Cfa stations were combined with Csa, after which differences between regime means were tested with 999 station-label permutations and FDR correction across 32 contrasts; these comparisons describe association with climate setting rather than a causal effect of Köppen–Geiger class.

### 2.5. Regional warming and fixed-threshold consistency checks

We calculated annual mean-temperature anomalies relative to each station's 1991–2007 mean, averaged them across available stations for each year, and related annual index counts to the resulting regional station-temperature anomaly with quantile models. The coefficient provides an internally scaled association in days yr⁻¹ °C⁻¹, but it is not independent evidence of forcing because the anomaly and indices derive from the same station archive.

We also defined station-level signal-to-noise ratios as the absolute bootstrap mean slope divided by its bootstrap standard deviation, using a ratio of at least 2 together with the expected warming-consistent sign as a descriptive screening threshold. By contrast, composite “fingerprint scores” produced during exploratory analysis were excluded from primary evidence because their components are dependent and their equal weighting lacks external calibration.

Additional diagnostics were retained for transparency and are reported in the Supplementary Material: representative-station extraction from cluster centroids, alternative linkage and distance specifications, interpolation-method comparisons for display maps, simple associations of trend metrics with latitude, longitude, and elevation, and signal-emergence summaries based on the same sign and signal-to-noise rule. These diagnostics do not expand the causal scope of the station analysis.

### 2.6. Compound dry–hot analysis

Following the empirical joint-probability logic used in studies of concurrent dryness and heat (Corbella and Stretch, 2012; Sarhadi et al., 2018; Alizadeh et al., 2020), we examined annual precipitation with annual mean temperature and June–September precipitation with mean daily maximum temperature over the same months. A station-year was retained when both variables had at least 80% daily coverage, and a station required at least 25 valid years. We use the term *dry–hot* because precipitation deficit alone does not represent the meteorological, soil-moisture, hydrological, and human dimensions encompassed by drought (Van Loon et al., 2016).

For precipitation \(P_y\) and temperature \(T_y\) in year \(y\), the empirical probability of conditions at least as dry and hot was

$$
\hat p_{AND,y}=\frac{1}{n}\sum_{i=1}^{n}\mathbf{1}(P_i\leq P_y,\;T_i\geq T_y),
$$

with empirical joint return period \(\widehat{RP}_{AND,y}=1/\hat p_{AND,y}\). We used the complete station record to define ranks and report thresholds of 5, 10, and 20 years, focusing on the latter two; however, with only 34 years, these values represent empirical rarity classes rather than stable estimates of rare long-return-period hazards.

For each year and threshold, we counted affected stations, calculated their percentage among valid stations, and summarized temporal change with Kendall's \(\tau\) and linear slopes. Because annual extent can be serially dependent, we added a residual moving-block sensitivity with 4-year blocks and 4,999 replicates: detrended residual blocks generated a no-trend null for two-sided probabilities, whereas adding resampled blocks to the fitted trend produced 95% intervals. Kolmogorov–Smirnov, Cramér–von Mises, and Anderson–Darling tests assessed differences between 1991–2007 and 2008–2024, while Moran's \(I\), five-nearest-neighbor weights, and 499 permutations evaluated annual binary event fields.

For descriptive comparison of event components, we also calculated marginal empirical return periods for precipitation deficit and heat excess, labeling a joint event dry-dominant when its marginal dry return period exceeded its marginal hot return period, hot-dominant under the reverse ordering, and co-dominant when the two were equal. Because the label compares marginal rarity, it identifies an event component rather than a causal driver.

## 3. Results

### 3.1. Quantile-dependent thermal changes

Our analysis shows that warm indices increased and cool indices declined across most of the station network, although their slopes differed among quantiles (Fig. 2; Table 1). For the annual network-mean warm-day series, the slope increased from 9.45 days decade⁻¹ at \(\tau=0.10\) to 15.55 at \(\tau=0.50\) and 22.50 at \(\tau=0.90\), yielding \(\Delta_1=13.06\) days decade⁻¹ and indicating greater change in high-count years. By contrast, warm nights increased with nearly equal median and upper-quantile slopes (12.14 and 12.21 days decade⁻¹), which produced a smaller \(\Delta_1\) of 3.47 days decade⁻¹.

![Figure 2. Regional quantile trends.](../outputs/figures/ijoc_regional_quantile_panels.png)

*Figure 2. Quantile-regression slope profiles fitted to annual network-mean counts of warm days, warm nights, cool days, and cool nights. Slopes are expressed in days per decade. These profiles summarize regional behavior; station-level models provide the inferential basis.*

**Table 1.** Network-mean trends and station-level median field significance.

| Index | OLS slope | \(\beta_1(0.10)\) | \(\beta_1(0.50)\) | \(\beta_1(0.90)\) | \(\Delta_1\) | FDR-retained stations at \(\tau=0.50\) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Warm days | 14.02 | 9.45 | 15.55 | 22.50 | 13.06 | 115/124 |
| Warm nights | 10.96 | 8.74 | 12.14 | 12.21 | 3.47 | 104/124 |
| Cool days | −12.60 | −5.60 | −12.00 | −15.37 | −9.77 | 97/124 |
| Cool nights | −10.86 | −5.91 | −9.83 | −12.06 | −6.16 | 86/124 |

*Slopes and \(\Delta_1\) are in days decade⁻¹. OLS, ordinary least squares; FDR, Benjamini–Hochberg false discovery rate at 0.05. FDR was applied separately within each index at the median quantile.*

The cool indices changed in the opposite direction, with cool-day slopes becoming more negative from the lower to the upper quantile (−5.60 to −15.37 days decade⁻¹) and cool-night slopes changing from −5.91 to −12.06 days decade⁻¹. These gradients yielded \(\Delta_1\) values of −9.77 and −6.16 days decade⁻¹, respectively, so the strongest distributional asymmetry occurred in warm-day and cool-day counts, whereas nighttime changes were more evenly distributed among quantiles.

Station-level results show that the regional patterns were not produced by a small subset of sites. All 124 warm-day OLS slopes were positive, 123 stations had positive slopes at \(\tau=0.10\) and 0.50, and 122 remained positive at \(\tau=0.90\), while warm-night slopes were positive at 122, 117, and 116 stations across the three focal quantiles. Every station had a negative cool-day slope at \(\tau=0.50\) and 0.90, compared with 121 at \(\tau=0.10\), and cool-night slopes were also predominantly negative. The sign of \(\Delta_1\) was consequently warming-consistent at 118 warm-day stations, 107 warm-night stations, all 124 cool-day stations, and 116 cool-night stations.

Bootstrap intervals excluded zero at \(\tau=0.90\) for 73 warm-day, 88 warm-night, 57 cool-day, and 61 cool-night stations, whereas the corresponding counts at \(\tau=0.10\) were 47, 69, 74, and 68. Because fewer stations had tail intervals excluding zero than had directionally consistent slopes, network-wide sign consistency was stronger than station-specific tail precision.

### 3.2. Spatial organization and climate regimes

Tail asymmetry varied geographically (Fig. 3). Warm-day \(\Delta_1\) was generally positive and reached its largest values at Esfahan (34.11 days decade⁻¹), Zanjan (28.44), Daran (28.15), Abumusa Island (27.42), and Bandar-e Anzali (26.15), whereas every cool-day \(\Delta_1\) estimate was negative. The strongest upper-tail cool-day contractions occurred at Bandar-e Lengeh (−33.70), Chabahar (−32.87), Bushehr (coastal; −30.53), Orumiyeh (−26.39), and Jask (−25.54); however, these station estimates do not define continuous regional surfaces.

![Figure 3. Station-level tail asymmetry.](../outputs/figures/ijoc_main_delta1_maps.png)

*Figure 3. Station values of \(\Delta_1=\beta_1(0.90)-\beta_1(0.10)\) for the four temperature indices. Positive values denote stronger change at the upper quantile; negative values denote upper-quantile contraction relative to the lower quantile.*

Moran's \(I\) was significant in 6 of the 12 index–quantile fields and was most consistent for cool days, with \(I=0.290\), 0.192, and 0.348 at \(\tau=0.10\), 0.50, and 0.90 (permutation p = 0.002–0.004). Warm-day slopes were spatially autocorrelated at \(\tau=0.10\) (\(I=0.258\), p = 0.002) and 0.90 (\(I=0.166\), p = 0.008), but not at the median. In contrast, warm nights showed significant structure only at \(\tau=0.10\), and cool nights did not reach p < 0.05 at any focal quantile. Although exploratory clusters of day indices were more geographically compact than random station labels, night-index clusters were not; we therefore treat cluster membership as a secondary description rather than evidence of fixed climatic regions.

All six Köppen–Geiger groups exhibited warming-consistent changes, but their magnitudes differed (Fig. 4; Table 2). The mean warm-day slope at \(\tau=0.90\) was largest at hot-steppe stations (23.85 days decade⁻¹), followed by cold dry-summer (22.23) and cold-steppe stations (21.62). Hot-desert stations had the strongest upper-quantile decline in cool days (−20.59 days decade⁻¹) together with the largest upper-quantile increase in warm nights and decline in cool nights, whereas cold-desert stations generally showed the weakest changes.

![Figure 4. Climate-regime quantile trends.](../outputs/figures/advanced_climate_regimes/climate_regime_quantile_profiles.png)

*Figure 4. Mean station-level quantile slopes within six Köppen–Geiger climate groups. BWh, hot desert; BWk, cold desert; BSh, hot steppe; BSk, cold steppe; C, temperate Csa/Cfa; Dsa, cold dry-summer.*

**Table 2.** Selected climate-regime results.

| Climate regime | Stations | Mean elevation (m) | Warm-day \(\beta_1(0.90)\) | Cool-day \(\beta_1(0.90)\) | Fixed-threshold warm-day change | Fixed-threshold cool-day change |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| BWh hot desert | 36 | 575 | 15.84 | −20.59 | 32.01 | −11.12 |
| BWk cold desert | 18 | 1,402 | 16.75 | −15.28 | 28.73 | −8.04 |
| BSh hot steppe | 10 | 495 | 23.85 | −16.85 | 46.62 | −11.98 |
| BSk cold steppe | 35 | 1,442 | 21.62 | −17.56 | 38.36 | −10.59 |
| C temperate | 15 | 672 | 20.44 | −18.37 | 41.12 | −10.13 |
| Dsa cold dry-summer | 10 | 1,917 | 22.23 | −19.62 | 46.63 | −13.81 |

*Quantile slopes are in days decade⁻¹. Fixed-threshold changes compare mean annual counts in 2008–2024 with 1991–2007 and are in days yr⁻¹.*

After adjustment for 32 regime contrasts, only the median warm-day response to regional temperature anomaly and the upper-quantile cool-day response had FDR q < 0.05 (q = 0.048 for each); the remaining differences are therefore descriptive, particularly for groups containing only 10–18 stations.

### 3.3. Fixed thresholds, regional warming, and robustness

The network-mean station temperature anomaly increased by 0.446 °C decade⁻¹ (95% confidence interval: 0.272–0.621 °C decade⁻¹). Thresholds derived from 1991–2007 produced late-period changes of +36.78 warm days, +28.41 warm nights, −10.69 cool days, and −8.19 cool nights yr⁻¹ in 2008–2024 (Fig. 5). Because the corresponding direction was consistent at 124, 117, 124, and 110 stations, the primary signs were not an artifact of using the complete 1991–2024 record to define percentile thresholds.

![Figure 5. Fixed-threshold early–late changes.](../outputs/figures/advanced_climate_change_signal/fixed_baseline_period_change_summary.png)

*Figure 5. Changes in annual warm-day, warm-night, cool-day, and cool-night counts in 2008–2024 relative to thresholds and mean counts derived from 1991–2007. Bars summarize station-level changes across the 124-station network.*

The index response per degree of regional station warming also had the expected direction. At \(\tau=0.90\), slopes were +29.62 warm days, +22.60 warm nights, −21.03 cool days, and −19.66 cool nights yr⁻¹ °C⁻¹, while at the median the expected sign occurred at every warm-day, warm-night, and cool-day station and at 119 of 124 cool-night stations. These associations agree with the temporal trends but do not provide independent attribution estimates.

The main results were stable across several analysis choices. Excluding the 39 stations flagged by any detrended homogeneity test changed network \(\Delta_1\) from 13.06 to 13.05 days decade⁻¹ for warm days, 3.47 to 3.45 for warm nights, −9.77 to −10.39 for cool days, and −6.16 to −4.87 for cool nights. Increasing bootstrap depth from 200 to 400 replicates produced station-level correlations above 0.996 for the primary bootstrap means; in contrast, moving-block and maximum-entropy summaries differed more at \(\tau=0.90\), with mean absolute differences of 1.68–3.34 days decade⁻¹, indicating method-dependent tail uncertainty despite stable principal signs.

Masking logically inconsistent daily temperatures had still smaller effects, as the largest absolute change across the four indices and three focal quantiles was 0.026 days decade⁻¹ for a network slope and 0.024 days decade⁻¹ for \(\Delta_1\). In the compound analysis, screening changed early- or late-period affected-station means by no more than 0.53 station and linear trend slopes by no more than 0.38 station decade⁻¹; thus, contradictory records do not explain the reported network patterns, although they remain part of the broader data-quality limitation.

The split-period slopes did not indicate uniform acceleration at every quantile. Warm-day OLS and median slopes increased from 10.94 and 10.43 days decade⁻¹ in 1991–2007 to 20.01 and 24.50 in 2008–2024, whereas the upper-quantile slope decreased from 23.65 to 14.97; warm nights showed a similar contrast. Negative cool-day and cool-night OLS slopes were larger in magnitude during the first half of the record, but given the 17-year subperiods, these estimates indicate descriptive temporal nonlinearity rather than tested breakpoints.

![Figure 8. Robustness and split-period diagnostics.](../outputs/figures/ijoc_robustness_synthesis.png)

*Figure 8. Synthesis of bootstrap-method, bootstrap-depth, homogeneity-exclusion, and related robustness diagnostics. The figure is used to compare sensitivity of estimates, not to create an additional significance test.*

![Figure 9. Split-period quantile comparison.](../outputs/figures/ijoc_split_period_comparison.png)

*Figure 9. Network-mean OLS and focal-quantile slopes for 1991–2007 and 2008–2024. Differences among quantiles are interpreted as descriptive temporal nonlinearity because each subperiod contains 17 years.*

### 3.4. Expansion of compound dry–hot conditions

The number of stations meeting the empirical dry–hot thresholds increased under both event definitions (Fig. 6; Table 3). Under the annual definition, \(RP\geq10\) events affected an average of 12.71 stations yr⁻¹ in 1991–2007 and 39.82 in 2008–2024, with a linear trend of 15.63 stations decade⁻¹. The raw Kendall result (\(\tau=0.318\), p = 0.0086) remained positive under the moving-block null (p = 0.0152). Annual \(RP\geq20\) events likewise increased from 3.12 to 17.12 affected stations yr⁻¹, for which the dependence-aware Kendall p-value was 0.0258.

![Figure 6. Extent of compound dry–hot conditions.](../outputs/compound_dry_hot/figures/compound_dry_hot_extent_timeseries.png)

*Figure 6. Annual percentage of valid stations meeting empirical dry–hot joint-return-period thresholds. The annual definition combines annual precipitation and mean temperature; the warm-season definition combines June–September precipitation and mean maximum temperature.*

Warm-season changes were larger. Events with \(RP\geq10\) affected an average of 10.35 stations yr⁻¹ in the early period and 41.53 in the late period, with \(\tau=0.583\) and a slope of 19.24 stations decade⁻¹. Under the 4-year block null, p = 0.0002 for both statistics and the 95% slope interval was 13.94–24.53 stations decade⁻¹. Warm-season \(RP\geq20\) events increased from 2.59 to 16.12 stations yr⁻¹ (\(\tau=0.539\); slope = 8.72 stations decade⁻¹), with dependence-aware p-values of 0.0002 and a slope interval of 5.66–11.68 stations decade⁻¹. The largest extent occurred in 2021, when 99 stations met the \(RP\geq10\) threshold.

**Table 3.** Trends in the number of stations meeting compound dry–hot thresholds.

| Definition | Threshold | 1991–2007 mean | 2008–2024 mean | Kendall \(\tau\) | Block-null p for \(\tau\) | Slope (stations decade⁻¹) | 95% block-bootstrap slope interval |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Annual precipitation + mean temperature | \(RP\geq10\) | 12.71 | 39.82 | 0.318 | 0.0152 | 15.63 | 7.11–24.05 |
| Annual precipitation + mean temperature | \(RP\geq20\) | 3.12 | 17.12 | 0.320 | 0.0258 | 7.95 | 2.71–13.23 |
| June–September precipitation + maximum temperature | \(RP\geq10\) | 10.35 | 41.53 | 0.583 | 0.0002 | 19.24 | 13.94–24.53 |
| June–September precipitation + maximum temperature | \(RP\geq20\) | 2.59 | 16.12 | 0.539 | 0.0002 | 8.72 | 5.66–11.68 |

*Return periods are empirical rarity classes within the 34-year record. Dependence-aware probabilities and intervals use 4-year residual moving blocks and 4,999 replicates.*

Early–late distribution tests supported the warm-season contrast: at \(RP\geq10\), mean affected-station percentage increased from 8.50% to 34.04%, with p-values of 0.00024 for Kolmogorov–Smirnov, 0.000057 for Cramér–von Mises, and at most 0.001 for Anderson–Darling. At \(RP\geq20\), mean extent increased from 2.12% to 13.22%, with corresponding p-values of 0.0046, 0.00057, and at most 0.001.

Event-component classification also changed between periods (Fig. 7). Among warm-season \(RP\geq10\) station-events, the early period contained 147 dry-dominant, 27 hot-dominant, and 2 co-dominant cases, whereas the late period contained 141, 535, and 30, respectively. At \(RP\geq20\), the corresponding counts changed from 39 dry-dominant and 4 hot-dominant cases to 51 and 205. Spatial connectedness of the annual warm-season event fields increased in parallel, with Moran's \(I\) slopes of 0.075 decade⁻¹ at \(RP\geq10\) (p = 0.0047) and 0.087 decade⁻¹ at \(RP\geq20\) (p = 0.0017), while the annual definition showed weaker and less stable changes.

![Figure 7. Change in the relative rarity of dry and hot event components.](../outputs/compound_dry_hot/figures/compound_dry_hot_driver_shift.png)

*Figure 7. Fractions of compound station-events classified as dry-dominant, hot-dominant, or co-dominant in 1991–2007 and 2008–2024. Dominance denotes the larger marginal empirical return period and is not a causal attribution.*

![Figure 10. Spatial frequency of compound dry–hot station-years.](../outputs/compound_dry_hot/figures/compound_dry_hot_station_frequency_maps.png)

*Figure 10. Station-level frequency of annual and June–September compound dry–hot events at empirical (RP\geq10) and (RP\geq20) thresholds. Values are percentages of valid station-years; they are not area-weighted estimates.*

![Figure 11. Spatial connectedness of compound events.](../outputs/compound_dry_hot/figures/compound_dry_hot_moran_connectedness.png)

*Figure 11. Annual Moran's (I) diagnostics for the spatial connectedness of compound dry–hot event fields. Positive trends indicate increasing similarity among neighboring station classifications, subject to the network geometry and permutation design.*

## 4. Discussion

### 4.1. Distributional and regional expression of thermal change

Our results show a widespread increase in warm-event frequency and decline in cool-event frequency, consistent with previous assessments of Iran, the Middle East, and global land observations (Zhang et al., 2005a; Alexander et al., 2006; Soltani et al., 2016; Dunn et al., 2020), but the quantile analysis reveals that this regional warming signal was not uniform across annual count distributions. Warm-day trends increased from 9.45 days decade⁻¹ at the lower quantile to 22.50 days decade⁻¹ at the upper quantile, whereas cool-day trends changed from −5.60 to −15.37 days decade⁻¹. Consequently, years that already contained many warm days changed faster than low-warm-frequency years, while years with many cool-day occurrences underwent the largest reductions.

These quantiles describe annual event counts rather than daily temperature magnitude: a high quantile represents a year near the upper part of the annual frequency distribution, not the hottest observed day or the intensity of a heatwave. The coefficients therefore indicate that high-warm-frequency years became more common at a faster rate and high-cool-frequency years diminished more rapidly, whereas trends in event intensity, duration, persistence, or within-event temperature severity require separate indices. This distinction is important because an annual frequency coefficient can otherwise be misinterpreted as a trend in daily temperature magnitude.

Although OLS slopes captured the common direction, they compressed the contrast between low- and high-count years into a single coefficient, and median slopes provided no information about the two tails. The large warm-day and cool-day values of \(\Delta_1\) therefore show why mean-only reporting is incomplete, consistent with quantile-based studies in which tail changes differed from central tendencies (Barbosa et al., 2011; Fan, 2014). Two stations can have similar mean slopes but markedly different changes in high-frequency years; as a result, monitoring based only on means or medians would not distinguish their distributional response.

Daytime and nighttime indices also followed different distributional trajectories. Warm-night and cool-night trends were widespread, but their upper-to-lower quantile contrasts were smaller than those of the day indices, and the upper-quantile warm-day slope exceeded the warm-night slope by 10.29 days decade⁻¹. We interpret this difference as evidence that the clearest asymmetry in the archive occurs in daytime event counts, whereas warm-night change more closely resembles a distribution-wide shift; it does not imply that nighttime warming is negligible, particularly at low-elevation hot-desert and coastal stations where high nighttime temperature and humid heat can persist after sunset (Raymond et al., 2024).

Several processes could generate day-night differences, because soil-moisture limitation, reduced evaporative cooling, cloud and radiation changes, boundary-layer behavior, and circulation persistence can influence maximum and minimum temperatures differently (Seneviratne et al., 2010; McKinnon et al., 2021; Vautard et al., 2023; Maraun et al., 2025). Iran's aridity and relief make these mechanisms plausible; however, the archive contains no direct measurements of soil moisture, radiation, cloud, circulation, irrigation, or surface-energy partitioning, so the present analysis identifies the asymmetry and its geography without assigning a dominant physical cause.

The direction of change was shared across all six Köppen–Geiger groups, yet effect magnitudes differed among climate settings: hot-steppe and cold dry-summer stations had the largest mean upper-quantile warm-day trends, whereas hot-desert stations had the largest upper-quantile cool-day reduction and the strongest night-index changes. Spatial diagnostics support part of this heterogeneity, as cool-day slopes were autocorrelated at all focal quantiles and warm-day slopes at both tails; in contrast, the absence of significant Moran's \(I\) in other fields indicates that spatial organization was neither uniform nor equally strong across indices. National means therefore describe the common direction, but maps should remain displays of station estimates rather than continuous representations of regional climate.

Climate-regime summaries provide an interpretable context for these differences without implying that a single geographical gradient explains the network. Recent evidence of warming-related climate-type transitions in Iran makes this comparison relevant, although that study addressed changing classifications rather than the quantile slopes examined here (Jamali et al., 2026). Our inference is narrower: stations assigned to different present-day climates share a warming-consistent direction but differ in average magnitude, while the classification itself neither explains individual-station change nor constitutes a causal exposure.

Only two of 32 between-regime contrasts remained below the FDR threshold, partly because some groups contain only 10–18 stations, neighboring sites are dependent, and a 1-km classification raster cannot resolve every site-scale influence. The remaining differences are therefore descriptive rather than evidence for six distinct trend populations, and the merger of two Cfa stations with Csa was a sample-size decision rather than a claim of climatic equivalence. Cluster assignments require similar restraint because several partitions changed with feature set, linkage, or distance measure and night-index clusters lacked geographic compactness; we consequently retain clustering as an exploratory supplement while basing the principal spatial interpretation on station estimates, regime summaries, and Moran's \(I\).

### 4.2. Expansion and changing composition of compound dry–hot conditions

The thermal-index results acquire a broader hydroclimatic context when evaluated alongside precipitation. Affected-station extent increased under both compound definitions, but the increase was larger and more consistent in June–September, when mean extent at the 10-year empirical rarity threshold rose from 10.35 to 41.53 stations between the two 17-year periods. Because the trend remained positive when detrended residuals were resampled in four-year blocks, this expansion cannot be explained solely by applying an independence-based test to a serially ordered annual series.

Warm-season concentration is physically credible because precipitation deficit, soil-moisture availability, and high temperature can interact during the period of greatest evaporative demand (Zscheischler and Seneviratne, 2017; Alizadeh et al., 2020; Lv et al., 2026); nevertheless, concurrence does not quantify feedback strength. The analysis includes neither soil moisture, evapotranspiration, vapor-pressure deficit, nor surface fluxes, and precipitation deficit represents only one dimension of drought. We therefore use *dry–hot* rather than *drought–heat*, consistent with the view that drought definitions depend on the process and impact domain under investigation (Van Loon et al., 2016).

Changes in event composition provide information beyond affected-station counts. Early warm-season events were predominantly dry-dominant, meaning that precipitation deficit was marginally rarer than high temperature within each station record, whereas late-period events were predominantly hot-dominant; heat thus became the more unusual component in a much larger fraction of joint events. This observation agrees with the independently estimated increase in warm-event frequencies, although it neither establishes that heat caused the precipitation deficit nor separates thermodynamic change from circulation variability.

Our comparison with Alizadeh et al. (2020) is informative because both studies use empirical joint probabilities to examine compound-event extent and relative component rarity, but the difference in record length precludes direct comparison of hazard magnitude. Their 122-year US climate-division archive represents a much wider range of temporal variability and supports longer return-period analysis, whereas the 34-year Iranian station archive resolves only within-record rarity classes. The common result is therefore qualitative: heat became a more prominent component of recent dry–hot events, while the present study further shows that this transition coincided with quantile-dependent changes in four temperature-frequency indices and was strongest during the warm season.

Increasing Moran's \(I\) suggests that warm-season events also became more connected within the observing network, which is relevant because simultaneous conditions at many sites differ from isolated local occurrence. Station connectedness is not, however, equivalent to affected land area: spatial sampling is uneven, every station receives equal weight, and a dense group of sites can influence the extent metric more than a sparsely monitored region of equal size. Gridded observations, area weighting, or spatial event reconstruction would be required to estimate physical extent, whereas the present metric supports the narrower conclusion that recent joint events involved a larger and more organized share of the station network.

The thermal indices and compound analysis answer related but distinct questions. Percentile indices describe how annual warm- and cool-event frequencies changed across their distributions, while empirical joint probabilities identify station-years that were simultaneously dry and hot relative to the local record. Agreement between the two lines of evidence strengthens the interpretation that exceptional heat has become more prominent in recent Iranian climate conditions, but it does not constitute attribution or establish future risk without assumptions about stationarity and continued change.

### 4.3. Robustness and methodological implications

The sensitivity analyses address different vulnerabilities rather than repeating the same confirmation. Fixed thresholds test whether full-period percentile estimation created the trend direction; homogeneity exclusion assesses whether potentially discontinuous stations dominate regional summaries; bootstrap comparisons evaluate Monte Carlo stability and dependence on the resampling model; internal-consistency masking tests contradictory daily temperatures; and residual blocks address serial dependence in compound extent. Because these checks target separate parts of the analytical chain, their agreement provides a stronger basis for the main direction than would repeated application of a single model.

Thresholds estimated from 1991–2007 reproduced the direction of the full-period indices, reducing concern that later observations influenced the percentiles in a way that generated the reported signs. This sensitivity does not establish a conventional climatological normal, because 17 years provide fewer day-of-year observations than a standard 30-year reference period and early-period sampling variability affects the thresholds. Similarly, excluding the 39 stations flagged by a detrended homogeneity test left regional \(\Delta_1\) estimates close to their full-network values, but annual breakpoint diagnostics cannot replace metadata-supported daily homogenization or detect every seasonal and tail-specific discontinuity. The combined evidence shows that neither full-period thresholds nor flagged stations control the national sign pattern, while local tail magnitudes remain more uncertain.

Increasing moving-block bootstrap depth from 200 to 400 replicates changed the bootstrap means little, whereas maximum-entropy and moving-block results differed more for several upper-tail summaries. Monte Carlo depth was therefore less influential than the resampling model, and directional agreement is more secure than fine comparison of interval endpoints. The internal-consistency screen led to similarly small changes, no greater than 0.026 days decade⁻¹ in focal network slopes or 0.024 days decade⁻¹ in \(\Delta_1\); these records should still be resolved if corrected source data become available, but they do not account for the regional findings.

Split-period slopes did not support uniform acceleration at every quantile. Although some OLS and median warm-index slopes were larger in 2008–2024, upper-quantile warm-day slopes and several cool-index slopes followed a different contrast, which indicates that a single linear coefficient does not capture every temporal feature of the record. Because each subperiod contains only 17 years and is sensitive to individual events and internal variability, we regard these estimates as descriptive temporal nonlinearity rather than evidence for a breakpoint or regime shift.

More broadly, the analysis shows why distribution-aware methods can add information to national assessments. Mean trends estimate average change, whereas quantile regression identifies whether low-, typical-, and high-frequency years evolve at comparable rates; in Iran, upper-to-lower differences for warm and cool days were of similar order to the mean slopes themselves. Reporting a small set of pre-specified quantiles alongside OLS therefore retains both average and distributional information rather than substituting one incomplete summary for another.

Applying quantile models across stations, indices, and probability levels also creates multiplicity and dependence. We limited the primary interpretation to three focal quantiles, applied FDR at the median, and used bootstrap intervals at the tails, but this design does not provide simultaneous field significance over the full 0.10–0.90 grid. Spatial quantile models, hierarchical shrinkage, or simultaneous confidence bands could extend the analysis by representing dependence among stations and quantiles explicitly (Reich, 2012), especially where station density is uneven.

Empirical joint probabilities offer a complementary transferable method because they make few distributional assumptions and place stations with different climates on a relative scale, although their largest resolvable rarity is constrained by sample length and full-period ranks mix early and late climate states. Parametric copulas could extrapolate further but would introduce assumptions that are difficult to validate with 34 annual observations; consequently, the empirical approach is appropriate for within-record comparison, whereas longer records or nonstationary multivariate models are necessary for design-level hazard estimates.

We excluded the exploratory composite fingerprint score from the main evidence because equal weighting of dependent indicators can create an apparently precise scale without external calibration. The revised framework instead keeps effect size, sign consistency, multiplicity-adjusted station counts, uncertainty, spatial structure, climate-regime summaries, and compound extent visible, allowing each element to be evaluated separately. This structure is transferable to other dryland or topographically complex regions, although climate classes, seasons, bootstrap blocks, and rarity thresholds must be adapted to the local record and research question.

### 4.4. Scope, relevance, and research priorities

The findings are relevant to climate monitoring because they identify which parts of the annual frequency distribution are changing most rapidly. For Iran, reporting only mean temperature or a single index trend would miss the pronounced upper-quantile response of warm days and contraction of cool days, whereas a compact set of lower-, median-, and upper-quantile estimates could test whether this structure persists as the record is updated. Station-level results remain necessary because national and regime averages conceal local departures along coasts, within hot deserts, and across high-elevation terrain.

The meteorological analysis does not quantify impacts because it contains no observations of population exposure, health outcomes, crops, ecosystems, electricity demand, reservoir operation, or adaptive capacity. It can identify locations and event types for subsequent impact studies, but it cannot rank societal vulnerability or establish adaptation priorities. Coastal nighttime change, for example, may contribute to heat exposure, yet humid-heat assessment requires humidity-based indices and exposure information rather than temperature counts alone (Raymond et al., 2020, 2024).

Record length remains the central limit on inference: 34 years can characterize recent observational change but provide a restricted basis for tail regression, subperiod comparison, and empirical joint return periods. The 10- and 20-year labels are within-record rarity categories rather than stationary century-scale risk estimates, while full-period thresholds describe redistribution relative to the analyzed climate and the 17-year fixed baseline serves only as a sensitivity check. Neither construction supports extrapolation beyond the observed period without additional assumptions.

Spatial and physical interpretation is constrained for related reasons. The regional temperature anomaly and indices derive from the same archive, station counts are not area weighted, the Köppen–Geiger raster simplifies microclimates and siting, and no metadata-supported daily homogenization was available. Because circulation, humidity, radiation, soil moisture, evapotranspiration, irrigation, urbanization, and land-use observations are also absent, associations with geography, elevation, climate regime, or regional warming remain descriptive rather than causal.

Future work should first assemble station histories and neighboring reference series for daily homogenization, then compare the station estimates with reanalysis and spatial quantile models to evaluate coverage gaps and covariance. Process observations and model experiments are also needed to test alternative circulation and land–atmosphere mechanisms, while longer records, soil-moisture and evaporative-demand measures, event duration, and impact data would determine whether the expanding station footprint translates into hydrological, agricultural, ecological, or health consequences. Formal anthropogenic attribution would additionally require counterfactual climate-model ensembles.

These constraints narrow, rather than negate, the central conclusion. Across several definitions and sensitivity tests, recent Iranian thermal change cannot be represented as a uniform shift in mean conditions: the largest changes occurred in particular parts of the annual event-count distribution and varied among climate settings, while concurrent heat and precipitation deficit expanded across the station network during the warm season.

## 5. Conclusions

This study provides a station-based assessment of how the distribution of annual thermal-extreme counts changed across Iran between 1991 and 2024. The principal result is not a single national warming coefficient but a systematic change in the shape of the frequency distributions. Network warm-day slopes increased from 9.45 days decade⁻¹ at the 0.10 quantile to 22.50 days decade⁻¹ at the 0.90 quantile, while cool-day slopes changed from −5.60 to −15.37 days decade⁻¹. Thus, years already characterized by frequent warm days intensified faster than low-frequency years, and years with many cool days experienced the strongest contraction. Warm nights and cool nights changed in the same respective directions, but their smaller (Delta_1) contrasts indicate a more distribution-wide nighttime response.

The direction of change was shared across all six Köppen–Geiger groups, but its magnitude and quantile structure varied. Hot-steppe and cold dry-summer stations had the largest mean upper-quantile warm-day slopes, whereas hot-desert stations showed the strongest upper-quantile cool-day reduction and pronounced nighttime changes. Only two of 32 regime contrasts remained after FDR correction, so these differences are best treated as climatic context rather than evidence that the classes define distinct causal responses. Moran's (I) nevertheless showed geographic organization for several daytime fields, particularly cool-day slopes. National means should therefore be read with station-level estimates, uncertainty intervals, and the spatial diagnostics rather than as a continuous uniform surface.

The compound analysis identified a parallel expansion of concurrent dry–hot conditions. Under the warm-season (RP\geq10) empirical class, the mean number of affected stations increased from 10.35 yr⁻¹ in 1991–2007 to 41.53 yr⁻¹ in 2008–2024; the trend remained positive under four-year residual moving-block resampling, with (	au=0.583), a slope of 19.24 stations decade⁻¹, and a 95% interval of 13.94–24.53. The (RP\geq20) class showed the same direction. Event composition also changed from predominantly dry-dominant to predominantly hot-dominant station-events. These are within-record rarity results, not stationary estimates of 10- or 20-year hazards, and the component shift does not establish a causal heat–precipitation pathway.

Confidence in the principal direction rests on convergence across distinct checks. FDR screening retained 115 warm-day, 104 warm-night, 97 cool-day, and 86 cool-night median station trends; fixed 1991–2007 thresholds reproduced the signs of the full-period analysis; excluding 39 detrended-homogeneity-flagged stations changed regional (Delta_1) by at most 1.29 days decade⁻¹ in absolute value; and masking contradictory daily temperature records changed focal network slopes by no more than 0.026 days decade⁻¹. Increasing bootstrap depth had little effect on mean estimates, although differences among resampling models at the upper quantile show that tail interval endpoints remain method-dependent. These checks support the direction and broad spatial pattern, not arbitrary precision in every local tail estimate.

The scope of inference remains deliberately narrow. A 34-year record limits tail precision, subperiod comparisons, and any extrapolation to century-scale variability; the empirical return-period classes cannot substitute for stationary design estimates. Daily observations were screened but not metadata-homogenized, the network is not area weighted, and the warming covariate is derived from the same archive as the indices. The analysis also lacks direct measurements of humidity, soil moisture, radiation, circulation, evaporative demand, irrigation, exposure, or impacts. Consequently, the evidence supports an observational conclusion that recent Iranian thermal extremes changed in a warming-consistent but distributionally and climatically differentiated manner, accompanied by a larger warm-season station footprint of concurrent dry–hot conditions. It does not identify a dominant physical mechanism or quantify societal risk.

The most direct scientific extension is to combine longer, homogenized station records with spatial quantile models and independent gridded products, then test whether the observed day–night and regime contrasts persist after accounting for spatial covariance and measurement changes. Process observations and model experiments are needed to separate thermodynamic warming from circulation and land–atmosphere feedbacks, while humidity, soil-moisture, evaporative-demand, and impact data are required to determine whether the station-network signal translates into health, agricultural, ecological, or water-management consequences. Until those data are available, the present results offer a reproducible baseline for monitoring distributional thermal change and compound dry–hot concurrence across Iran.

## Data availability

The station-data provider, persistent access route, license, and any redistribution restrictions must be supplied before submission. Derived tables and figure-generation code are organized in the project repository. **[INSERT DATA-PROVIDER AND REPOSITORY STATEMENTS]**

## Code availability

The analysis pipeline, configuration file, dependency specification, and scripts used to generate the derived outputs accompany this manuscript. A permanent public repository and archived release should be cited at submission. **[INSERT REPOSITORY URL AND ARCHIVE DOI]**

## Author contributions

**[INSERT CRediT AUTHOR-CONTRIBUTION STATEMENT]**

## Funding

**[INSERT FUNDING INFORMATION OR “This research received no external funding.”]**

## Competing interests

**[INSERT JOURNAL-COMPLIANT COMPETING-INTERESTS DECLARATION]**

## Acknowledgements

**[INSERT ACKNOWLEDGEMENTS, INCLUDING THE DATA PROVIDER IF REQUIRED]**

## Supplementary material

The Supplementary Material contains station-level quantile results, robustness summaries, climate-regime contrasts, cluster-stability diagnostics, compound-event tests at all rarity thresholds, and internal-consistency sensitivity tables. A structured index of these artifacts is provided in [Supplementary_Material_Q1.md](Supplementary_Material_Q1.md), while the corresponding scientific audit is reported in [Scientific_Audit_Q1.md](Scientific_Audit_Q1.md).

## References

Alexander LV, Zhang X, Peterson TC, Caesar J, Gleason B, Klein Tank AMG, Haylock M, Collins D, Trewin B, Rahimzadeh F, Tagipour A, Kumar KR, Revadekar J, Griffiths G, Vincent L, Stephenson DB, Burn J, Aguilar E, Brunet M, Taylor M, New M, Zhai P, Rusticucci M, Vazquez-Aguirre JL (2006). Global observed changes in daily climate extremes of temperature and precipitation. *Journal of Geophysical Research: Atmospheres* 111:D05109. https://doi.org/10.1029/2005JD006290

Alexandersson H (1986). A homogeneity test applied to precipitation data. *Journal of Climatology* 6:661–675. https://doi.org/10.1002/joc.3370060607

Alizadeh MR, Adamowski J, Nikoo MR, AghaKouchak A, Dennison P, Sadegh M (2020). A century of observations reveals increasing likelihood of continental-scale compound dry-hot extremes. *Science Advances* 6:eaaz4571. https://doi.org/10.1126/sciadv.aaz4571

Barbosa SM, Scotto MG, Alonso AM (2011). Summarising changes in air temperature over Central Europe by quantile regression and clustering. *Natural Hazards and Earth System Sciences* 11:3227–3233. https://doi.org/10.5194/nhess-11-3227-2011

Beck HE, Zimmermann NE, McVicar TR, Vergopolan N, Berg A, Wood EF (2018). Present and future Köppen–Geiger climate classification maps at 1-km resolution. *Scientific Data* 5:180214. https://doi.org/10.1038/sdata.2018.214

Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B* 57:289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

Buishand TA (1982). Some methods for testing the homogeneity of rainfall records. *Journal of Hydrology* 58:11–27. https://doi.org/10.1016/0022-1694(82)90066-X

Corbella S, Stretch DD (2012). Multivariate return periods of sea storms for coastal erosion risk assessment. *Natural Hazards and Earth System Sciences* 12:2699–2708. https://doi.org/10.5194/nhess-12-2699-2012

Dunn RJH, Alexander LV, Donat MG, Zhang X, Bador M, Herold N, Lippmann T, Allan R, Aguilar E, Brunet M, Caesar J, Chagnaud G, Cheng V, Cinco T, Durre I, Htay MM, Hoang L, Hung NQ, Johnson F, Kruger A, Lau K, Leng TW, Loikith PC, Marengo J, Mbatha S, McGree S, Menne M, Skansi M, Trewin B, Villarroel C, Vincent LA, Vose RS, Yeo R, Zhang P (2020). Development of an updated global land in situ-based dataset of temperature and precipitation extremes: HadEX3. *Journal of Geophysical Research: Atmospheres* 125:e2019JD032263. https://doi.org/10.1029/2019JD032263

Fan LJ (2014). Quantile trends in temperature extremes in China. *Atmospheric and Oceanic Science Letters* 7:304–308. https://doi.org/10.3878/j.issn.1674-2834.13.0102

Fischer EM, Knutti R (2015). Anthropogenic contribution to global occurrence of heavy-precipitation and high-temperature extremes. *Nature Climate Change* 5:560–564. https://doi.org/10.1038/nclimate2617

Fischer EM, Sippel S, Knutti R (2021). Increasing probability of record-shattering climate extremes. *Nature Climate Change* 11:689–695. https://doi.org/10.1038/s41558-021-01092-9

Francis D, Fonseca R (2024). Recent and projected changes in climate patterns in the Middle East and North Africa region. *Scientific Reports* 14:10279. https://doi.org/10.1038/s41598-024-60976-w

Frich P, Alexander LV, Della-Marta P, Gleason B, Haylock M, Klein Tank AMG, Peterson T (2002). Observed coherent changes in climatic extremes during the second half of the twentieth century. *Climate Research* 19:193–212. https://doi.org/10.3354/cr019193

Hall P, Horowitz JL, Jing B-Y (1995). On blocking rules for the bootstrap with dependent data. *Biometrika* 82:561–574. https://doi.org/10.1093/biomet/82.3.561

Hao Z, AghaKouchak A, Phillips TJ (2013). Changes in concurrent monthly precipitation and temperature extremes. *Environmental Research Letters* 8:034014. https://doi.org/10.1088/1748-9326/8/3/034014

Hubert L, Arabie P (1985). Comparing partitions. *Journal of Classification* 2:193–218. https://doi.org/10.1007/BF01908075

Jamali M, Eslamian S, Shayannejad M, Gohari A (2026). Observed warming-driven aridification and climate-type transitions across Iran. *Journal of Arid Environments* 235:105606. https://doi.org/10.1016/j.jaridenv.2026.105606

Katz RW, Brown BG (1992). Extreme events in a changing climate: variability is more important than averages. *Climatic Change* 21:289–302. https://doi.org/10.1007/BF00139728

Katz RW, Parlange MB, Naveau P (2002). Statistics of extremes in hydrology. *Advances in Water Resources* 25:1287–1304. https://doi.org/10.1016/S0309-1708(02)00056-8

Koenker R (2005). *Quantile Regression*. Cambridge University Press, Cambridge. https://doi.org/10.1017/CBO9780511754098

Koenker R, Bassett G Jr (1978). Regression quantiles. *Econometrica* 46:33–50. https://doi.org/10.2307/1913643

Kunsch HR (1989). The jackknife and the bootstrap for general stationary observations. *The Annals of Statistics* 17:1217–1241. https://doi.org/10.1214/aos/1176347265

Lahiri SN (2003). *Resampling Methods for Dependent Data*. Springer, New York. https://doi.org/10.1007/978-1-4757-3803-2

Leonard M, Westra S, Phatak A, Lambert M, van den Hurk B, McInnes K, Risbey J, Schuster S, Jakob D, Stafford-Smith M (2014). A compound event framework for understanding extreme impacts. *WIREs Climate Change* 5:113–128. https://doi.org/10.1002/wcc.252

Lv B, Wang S, Chen G, Xiang B (2026). Precipitation and soil moisture coupling constrains subseasonal predictability of a prolonged extreme heatwave. *Communications Earth & Environment* 7:323. https://doi.org/10.1038/s43247-026-03341-1

Maraun D, Schiemann R, Ossó A, Jury M (2025). Changes in event soil moisture–temperature coupling can intensify very extreme heat beyond expectations. *Nature Communications* 16:734. https://doi.org/10.1038/s41467-025-56109-0

McKinnon KA, Poppick A, Simpson IR (2021). Hot extremes have become drier in the United States Southwest. *Nature Climate Change* 11:598–604. https://doi.org/10.1038/s41558-021-01076-9

Moran PAP (1950). Notes on continuous stochastic phenomena. *Biometrika* 37:17–23. https://doi.org/10.1093/biomet/37.1-2.17

Naderi M, Saatsaz M, Behrouj Peely A (2024). Extreme climate events under global warming in Iran. *Hydrological Sciences Journal* 69:337–364. https://doi.org/10.1080/02626667.2024.2317269

Pal JS, Eltahir EAB (2016). Future temperature in southwest Asia projected to exceed a threshold for human adaptability. *Nature Climate Change* 6:197–200. https://doi.org/10.1038/nclimate2833

Pettitt AN (1979). A non-parametric approach to the change-point problem. *Applied Statistics* 28:126–135. https://doi.org/10.2307/2346729

Raymond C, Matthews T, Horton RM (2020). The emergence of heat and humidity too severe for human tolerance. *Science Advances* 6:eaaw1838. https://doi.org/10.1126/sciadv.aaw1838

Raymond C, Matthews T, Tuholske C (2024). Evening humid-heat maxima near the southern Persian/Arabian Gulf. *Communications Earth & Environment* 5:591. https://doi.org/10.1038/s43247-024-01763-3

Reich BJ (2012). Spatiotemporal quantile regression for detecting distributional changes in environmental processes. *Journal of the Royal Statistical Society: Series C (Applied Statistics)* 61:535–553. https://doi.org/10.1111/j.1467-9876.2011.01025.x

Rousi E, Kornhuber K, Beobide-Arsuaga G, Luo F, Coumou D (2022). Accelerated western European heatwave trends linked to more-persistent double jets over Eurasia. *Nature Communications* 13:3851. https://doi.org/10.1038/s41467-022-31432-y

Sarhadi A, Ausín MC, Wiper MP, Touma D, Diffenbaugh NS (2018). Multidimensional risk in a nonstationary climate: joint probability of increasingly severe warm and dry conditions. *Science Advances* 4:eaau3487. https://doi.org/10.1126/sciadv.aau3487

Seneviratne SI, Corti T, Davin EL, Hirschi M, Jaeger EB, Lehner I, Orlowsky B, Teuling AJ (2010). Investigating soil moisture–climate interactions in a changing climate: a review. *Earth-Science Reviews* 99:125–161. https://doi.org/10.1016/j.earscirev.2010.02.004

Seneviratne SI, Donat MG, Mueller B, Alexander LV (2014). No pause in the increase of hot temperature extremes. *Nature Climate Change* 4:161–163. https://doi.org/10.1038/nclimate2145

Soltani M, Laux P, Kunstmann H, Stan K, Sohrabi MM, Molanejad M, Sabziparvar AA, Ranjbar SaadatAbadi A, Ranjbar F, Rousta I, Zawar-Reza P, Khoshakhlagh F, Soltanzadeh I, Babu CA, Azizi GH, Martin MV (2016). Assessment of climate variations in temperature and precipitation extreme events over Iran. *Theoretical and Applied Climatology* 126:775–795. https://doi.org/10.1007/s00704-015-1609-5

Vaghefi SA, Keykhai M, Jahanbakhshi F, Sheikholeslami J, Ahmadi A, Yang H, Abbaspour KC (2019). The future of extreme climate in Iran. *Scientific Reports* 9:1464. https://doi.org/10.1038/s41598-018-38071-8

Van Loon AF, Stahl K, Di Baldassarre G, Clark J, Rangecroft S, Wanders N, Gleeson T, Van Dijk AIJM, Tallaksen LM, Hannaford J, Uijlenhoet R, Teuling AJ, Hannah DM, Sheffield J, Svoboda M, Verbeiren B, Wagener T, Van Lanen HAJ (2016). Drought in a human-modified world: reframing drought definitions, understanding, and analysis approaches. *Hydrology and Earth System Sciences* 20:3631–3650. https://doi.org/10.5194/hess-20-3631-2016

Vautard R, Cattiaux J, Happé T, Singh J, Bonnet R, Cassou C, Coumou D, D'Andrea F, Faranda D, Fischer EM, Ribes A, Sippel S, Yiou P (2023). Heat extremes in Western Europe increasing faster than simulated due to atmospheric circulation trends. *Nature Communications* 14:6803. https://doi.org/10.1038/s41467-023-42143-3

Vinod HD (2006). Maximum entropy ensembles for time series inference in economics. *Journal of Asian Economics* 17:955–978. https://doi.org/10.1016/j.asieco.2006.09.001

Zhang X, Aguilar E, Sensoy S, Melkonyan H, Tagiyeva U, Ahmed N, Kutaladze N, Rahimzadeh F, Taghipour A, Hantosh TH, Albert P, Semawi M, Karam Ali M, Al-Shabibi MHS, Al-Oulan Z, Zatari T, Al Dean Khelet I, Hamoud S, Sagir R, Demircan M, Eken M, Adiguzel M, Alexander LV, Peterson TC, Wallis T (2005a). Trends in Middle East climate extreme indices from 1950 to 2003. *Journal of Geophysical Research: Atmospheres* 110:D22104. https://doi.org/10.1029/2005JD006181

Zhang X, Alexander L, Hegerl GC, Jones P, Klein Tank AMG, Peterson TC, Trewin B, Zwiers FW (2011). Indices for monitoring changes in extremes based on daily temperature and precipitation data. *WIREs Climate Change* 2:851–870. https://doi.org/10.1002/wcc.147

Zittis G, Hadjinicolaou P, Lelieveld J (2016). Strongly increasing heat extremes in the Middle East and North Africa in the 21st century. *Climatic Change* 137:245–260. https://doi.org/10.1007/s10584-016-1665-6

Zscheischler J, Seneviratne SI (2017). Dependence of drivers affects risks associated with compound events. *Science Advances* 3:e1700263. https://doi.org/10.1126/sciadv.1700263

Zscheischler J, Westra S, van den Hurk BJJM, Seneviratne SI, Ward PJ, Pitman A, AghaKouchak A, Bresch DN, Leonard M, Wahl T, Zhang X (2018). Future climate risk from compound events. *Nature Climate Change* 8:469–477. https://doi.org/10.1038/s41558-018-0156-3
