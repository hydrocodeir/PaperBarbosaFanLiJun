# Quantile-Dependent Changes in Thermal Extremes Across Iran: A Robust Station-Based Assessment of Distributional Shifts, Uncertainty, and Regionalization
### or
# Beyond Mean Warming: Quantile-Dependent Thermal-Extreme Change, Uncertainty, and Regionalization Across Iran

**Running title:** Quantile-based thermal extremes across Iran

## Abstract

Mean change alone is often insufficient for understanding thermal extremes because the most societally relevant impacts emerge through changes in different parts of the distribution rather than through a uniform shift in the center alone. This study analyzes quantile-dependent trends in four annual thermal-extreme indices across Iran during 1991-2024 using observations from 124 meteorological stations. Warm days, warm nights, cool days, and cool nights were examined with station-wise quantile regression, moving-block bootstrap uncertainty estimation, false-discovery-rate field-significance screening, spatial diagnostics, and cluster-based regionalization. The results show a clear but non-uniform restructuring of thermal extremes. Warm days exhibit the strongest regional increase and the clearest upper-tail amplification, with regional network-mean slopes of `+9.45`, `+15.55`, and `+22.50` days per decade at `q0.10`, `q0.50`, and `q0.90`, respectively, and `Delta1 = +13.06`. Warm nights also increase, but more uniformly across the distribution, with corresponding slopes of `+8.74`, `+12.14`, and `+12.21` days per decade and `Delta1 = +3.47`. Cool days and cool nights decline throughout the country, with stronger contraction in upper quantiles, especially for cool days (`-5.60`, `-12.00`, `-15.37`; `Delta1 = -9.77`) and, to a lesser degree, cool nights (`-5.91`, `-9.83`, `-12.06`; `Delta1 = -6.16`). At `q0.50`, field-significance screening retained `115/124` warm-day, `104/124` warm-night, `97/124` cool-day, and `86/124` cool-night station results after FDR control. Moran's I was significant in `6/12` index-quantile fields, with the strongest and most consistent spatial organization in cool-day trends. Regionalization further shows dominant large clusters for warm nights and cool nights, but greater internal diversity for warm days and cool days. Sensitivity analyses indicate that the signs and broad ordering of the main results are robust to homogeneity-based station exclusion and deeper bootstrap resampling, whereas alternative bootstrap methods and some clustering specifications affect magnitude more than direction. Split-period comparisons point to temporal nonlinearity rather than a simple monotonic amplification across all quantiles. Overall, the study shows that recent thermal-extreme change across Iran is distributionally asymmetric, spatially heterogeneous, and driven most strongly by warm-day amplification and widespread contraction of cool extremes.

**Keywords:** thermal extremes; quantile regression; climate change; warm days; cool days; regionalization; Iran

## 1. Introduction

Changes in climate extremes are among the most consequential manifestations of contemporary climate change because their impacts propagate through agriculture, water resources, ecosystems, infrastructure, and human health `(Easterling et al., 2000; IPCC, 2021)`. Heat-related risk is especially sensitive to shifts in the tails of the thermal distribution, not only to changes in the mean `(Katz and Brown, 1992; Katz et al., 2002; Perkins, 2015)`. A region can therefore experience moderate mean change while simultaneously undergoing pronounced upper-tail intensification or accelerated loss of historically cool conditions. Distribution-aware methods are consequently essential for diagnosing the structure of thermal-extreme change `(Koenker and Bassett, 1978; Koenker, 2005; Barbosa et al., 2011; Fan, 2014; Reich, 2012)`.

This issue is especially important for annual warm and cool day-night indices. These indices are widely used in climatology because they convert daily temperature variability into interpretable annual counts of unusually warm or cool conditions `(Frich et al., 2002; Alexander et al., 2006; Zhang et al., 2011; Dunn et al., 2020)`. Yet many studies still emphasize mean or median behavior, even though the rates of change at lower, central, and upper quantiles may differ substantially. Such differences matter because they distinguish broad distribution-wide shifts from more selective amplification in the most extreme years `(Moberg and Jones, 2005; Kostopoulou and Jones, 2005; Fischer and Schär, 2010; Russo and Sterl, 2011)`.

Iran provides a strong test bed for this question. The country spans humid Caspian lowlands, the Zagros and Alborz mountain systems, central desert basins, and hot southern coasts influenced by the Persian Gulf and Gulf of Oman. These gradients in elevation, continentality, aridity, and marine influence favor spatially heterogeneous responses rather than a single national thermal signal. Previous studies have documented intensifying temperature extremes across Iran and the broader Middle East and North Africa `(Zhang et al., 2005a; Soltani et al., 2016; Vaghefi et al., 2019; Zittis et al., 2016; Francis and Fonseca, 2024)`, while regional work has further emphasized the growing relevance of heat stress and compound thermal hazards in southwest Asia `(Pal and Eltahir, 2016; Raymond et al., 2020; Zittis et al., 2021; Raymond et al., 2024)`. Even so, comparatively fewer station-based studies for Iran have combined distributional trend analysis, uncertainty quantification, field significance, and quantitative regionalization within a single workflow.

The present study addresses that gap using a much denser contemporary station network than is typical in long-record national syntheses. Using 124 meteorological stations for 1991-2024, we analyze annual warm days, warm nights, cool days, and cool nights with quantile regression, bootstrap uncertainty estimation, field-significance screening, spatial diagnostics, and cluster-based regionalization. The aim is not only to identify whether thermal extremes are changing, but to determine how those changes are structured across the distribution, how consistently they organize in space, and how stable the resulting interpretation remains under alternative methodological choices.

Five contributions are emphasized. First, the analysis resolves changes across the conditional distribution rather than relying on mean change alone. Second, it measures distributional asymmetry with explicitly defined tail-contrast metrics. Third, it quantifies uncertainty through dependent bootstrap resampling and complementary significance diagnostics. Fourth, it regionalizes station behavior through reproducible clustering and representative-station selection. Fifth, it evaluates the stability of the main findings with homogeneity-exclusion, bootstrap-depth, bootstrap-method, interpolation, clustering, and split-period sensitivity analyses. Together, these elements support a distribution-aware climatological synthesis of recent thermal-extreme change across Iran.

## 2. Data and methods

### 2.1. Data and study region

The analysis is based on daily temperature observations from 124 meteorological stations distributed across Iran for the common period 1991-2024. After preprocessing, each station-index series therefore contributes 34 annual values. The network spans coastal, lowland, plateau, and high-elevation settings and is dense enough to capture substantial physiographic heterogeneity while remaining fully observation-based.

![Figure 1. Study area and station network.](../outputs/figures/ijoc_study_area.png)

*Figure 1. Distribution of the 124 meteorological stations used in the analysis across Iran. Stations are colored by elevation to highlight the combined spatial and topographic coverage of the network. Source note: [ijoc_study_area.md](../outputs/output_docs/figures/ijoc_study_area.md).*

The station-based design is deliberate. All primary inference is anchored in observed station records rather than interpolated grid cells, which is especially important in a region where elevation contrasts and marine influence can generate sharp local thermal differences over short distances. The daily archive includes minimum, maximum, and, where available, mean temperature. Leap-day observations were removed before percentile-threshold estimation in order to maintain a consistent 365-day annual cycle.

Data-quality screening showed generally high completeness across the network, with median station completeness of `98.96%` for `tmin`, `99.42%` for `tmax`, and `99.35%` for `tmean`. No station contained duplicate dates. Internal consistency diagnostics nevertheless identified at least one `tmin > tmax` conflict at 32 stations and at least one `tmean` value outside the daily temperature range at 78 stations, indicating that screening remained necessary even in a high-coverage archive. Homogeneity diagnostics applied to annual mean temperature were highly sensitive in raw form, flagging `119/124` stations, whereas detrended versions still flagged `39/124` stations. The contrast between raw and detrended tests indicates that breakpoint-sensitive homogeneity statistics respond strongly in a region with pervasive warming and should not be interpreted mechanically as evidence for direct daily-series correction.

![Figure 2. Data completeness and detrended homogeneity diagnostics.](../outputs/figures/ijoc_data_quality_homogeneity.png)

*Figure 2. Combined overview of station-level completeness and detrended homogeneity diagnostics for the 1991-2024 station archive. The left panel summarizes completeness, and the right panel summarizes the number of stations flagged by detrended homogeneity tests. Source note: [ijoc_data_quality_homogeneity.md](../outputs/output_docs/figures/ijoc_data_quality_homogeneity.md).*

No automatic homogenization adjustment was imposed on the daily temperature records. That choice was deliberate. The available archive does not include the station-history metadata or neighbor-based reference framework needed for a defensible daily homogenization procedure, and the main breakpoint diagnostics were applied to annual mean temperature rather than directly to the percentile-based daily series used for the indices `(Zhang et al., 2005b)`. The main manuscript therefore retains the original daily observations, treats the diagnostics as screening information, and uses flagged-station exclusion as a conservative sensitivity bound rather than as a substitute for full homogenization. Supporting summaries are provided in Supplementary Figure S1 and Supplementary Tables S1-S3, while the exclusion sensitivity is reported in Supplementary Figure S8 and Supplementary Table S4.

### 2.2. Thermal-extreme indices

Four annual thermal-extreme indices were constructed from daily minimum (`tmin`) and maximum (`tmax`) temperature: warm days, warm nights, cool days, and cool nights. For each station, day-of-year-specific thresholds were estimated with a 5-day moving window and the 10th and 90th percentiles of the station's available 1991-2024 record. A minimum of 15 reference samples per day-of-year window was required; when local window support fell below that threshold, the station-wide reference sample was used as the fallback.

Warm days and cool days were defined from daily maximum temperature relative to the station-specific 90th and 10th percentile thresholds, respectively. Warm nights and cool nights were defined analogously from daily minimum temperature. Binary exceedance indicators were then summed by year to obtain annual counts `(Frich et al., 2002; Alexander et al., 2006; Zhang et al., 2011)`. Years with less than `80%` valid daily coverage for the corresponding temperature variable were set to missing at the annual-index stage.

Because the current run contains only 1991-2024 data, the threshold climatology is based on the full available analysis period rather than on a separate earlier reference interval. The resulting indices therefore describe how annual frequencies of warm and cool events evolved relative to the empirical distribution of the analyzed period itself. This produces internally consistent station-wise indices for the present dataset, but it also means the annual counts should be interpreted as distribution-relative diagnostics rather than as departures from a fixed pre-change climatological baseline.

### 2.3. Quantile regression and summary metrics

Quantile regression was used as the primary inferential tool because the goal is to resolve distribution-dependent change rather than only mean behavior `(Koenker and Bassett, 1978; Koenker, 2005)`. For an annual response `Y` and time covariate `T`, the `tau`-level conditional quantile is written as

$$
Q_Y(\tau \mid T) = \beta_0(\tau) + \beta_1(\tau) T, \qquad \tau \in (0,1).
$$

The estimator is obtained from the asymmetric absolute-loss objective

$$
\hat{\beta}(\tau) = \arg\min_{\beta} \sum_{i=1}^{n} \rho_{\tau}\!\left(y_i - x_i^\top\beta\right),
$$

where `x_i = (1, T_i)^\top` and

$$
\rho_{\tau}(u) = u\big(\tau - \mathbf{1}\{u < 0\}\big).
$$

This framework permits slope heterogeneity across the distribution, so `\beta_1(\tau)` can differ meaningfully between lower, central, and upper quantiles. Time was expressed in decades, making all slopes interpretable as changes in annual event counts per decade.

The full quantile profile was estimated from `q = 0.10` to `q = 0.90` at increments of `0.01`. Main-text interpretation focuses on `q0.10`, `q0.50`, and `q0.90`, together with an OLS benchmark slope. In addition to station-wise fits, descriptive regional quantile profiles were obtained by averaging each annual index across all stations and fitting quantile regression to that network-mean time series. These regional profiles are used for visual synthesis, whereas the main inferential backbone remains the station-level analysis.

To summarize distributional asymmetry, three contrasts were computed:

`Delta1 = slope(q0.90) - slope(q0.10)`

`Delta2 = slope(q0.90) - slope(q0.50)`

`Delta3 = slope(q0.50) - slope(q0.10)`

Positive `Delta1` indicates stronger upper-tail amplification, whereas negative `Delta1` indicates relative upper-tail weakening or, for declining indices, stronger upper-tail contraction.

### 2.4. Uncertainty, spatial significance, and regionalization

Bootstrap resampling was used to characterize uncertainty in the station-level quantile slopes. The baseline workflow uses a moving-block bootstrap with 200 replicates. Because annual extreme-index series may retain short-range temporal dependence, simple iid resampling was not treated as the default approach. Block length was selected automatically with an `n^(1/3)` rule and then bounded to the interval `[2, 8]`, following the general dependent-bootstrap logic discussed by `(Kunsch, 1989; Hall et al., 1995; Lahiri, 2003)`. For each focal quantile and for `Delta1`, the pipeline retained bootstrap means, standard deviations, and confidence intervals.

Inference in the tails relied primarily on the bootstrap summaries. This choice matters because the analytic tail outputs were insufficient for a stable national FDR screening layer in the current run, whereas bootstrap intervals remained available for the focal tail summaries. Tail-significance counts are therefore discussed from the bootstrap results, while field-significance screening is reported for the median quantile where the analytic station-level outputs were complete.

Field significance was assessed with false-discovery-rate control `(Benjamini and Hochberg, 1995)`, and spatial dependence was examined with Moran's I `(Moran, 1950)` using 499 permutations and 5-nearest-neighbor weights. In this context, FDR controls the expected proportion of false positives among the station-level results declared significant after multiple testing. This is important because the study evaluates many station-wise trend tests in parallel; without multiplicity control, some apparently significant results would be expected to arise by chance alone. FDR was preferred because it provides a practical balance between error control and statistical power, making it well suited to a dense multi-station setting where an overly strict family-wise correction could suppress spatially meaningful signals. Regionalization was performed with hierarchical clustering using average linkage and Euclidean distance. The clustering stage was treated as exploratory but reproducible: it was based on a screened feature set derived from focal quantile slopes, delta contrasts, and selected bootstrap summaries, and it was accompanied by representative-station extraction, reduced-feature reruns, alternative clustering specifications, and a permutation-based spatial compactness test.

### 2.5. Additional synthesis analyses

Several additional analyses were used to evaluate the stability and interpretability of the main results. First, a homogeneity-exclusion sensitivity reran the regional summaries after excluding all stations flagged by any detrended homogeneity diagnostic. Second, a bootstrap-depth sensitivity reran the bootstrap stage with 400 rather than 200 replicates. Third, a bootstrap-method comparison contrasted the baseline moving-block results against maximum-entropy bootstrap (`meboot`) summaries `(Efron and Tibshirani, 1993; Vinod, 2006)`. Fourth, interpolation sensitivity compared alternative spatial-display methods so that mapped surfaces could be interpreted as visualization products rather than hidden inferential layers. Fifth, alternative clustering solutions examined sensitivity to linkage, distance metric, and algorithmic family. Sixth, split-period comparisons contrasted the first and second halves of the 1991-2024 record, namely `1991-2007` and `2008-2024`, using regional network-mean annual series. Finally, a simple driver analysis related slope and `Delta1` metrics to latitude, longitude, and elevation.

### 2.6. Analytical workflow

All analyses were generated through a single pipeline that integrates data screening, annual-index construction, station-wise quantile regression, bootstrap uncertainty estimation, feature engineering, clustering, spatial diagnostics, and publication-oriented sensitivity modules. This integrated workflow reduces manual transfer error between stages and ensures that all tables, figures, and regional summaries used in the manuscript are traceable to the same underlying run configuration.

## 3. Results and discussion

### 3.1. Central claim and overview

The central result is that thermal-extreme change across Iran during 1991-2024 is strongly directional but not distributionally uniform. Warm indices rise across most of the network, cool indices decline across most of the network, and the magnitude of change depends systematically on quantile, index, and region. The most prominent signal is the amplification of warm days at upper quantiles, while the most spatially persistent cooling-related signal is the contraction of cool days. Warm nights also increase broadly, but with flatter tail asymmetry and a more dominant single regional cluster. Cool nights decline widely, although their spatial organization is weaker than that of cool days.

At the station level, the directional consistency is striking. All 124 warm-day OLS slopes are positive, and `123/124` warm-day stations remain positive at both `q0.10` and `q0.50`, with `122/124` still positive at `q0.90`. Warm nights are also broadly positive, but with more exceptions: `117/124` OLS slopes are positive, with `122/124`, `117/124`, and `116/124` positive at `q0.10`, `q0.50`, and `q0.90`, respectively. Cool days are almost uniformly negative, with all OLS, `q0.50`, and `q0.90` slopes negative and only three stations showing positive `q0.10` slopes. Cool nights are likewise predominantly negative, with `116/124` negative OLS slopes and negative `Delta1` at `116/124` stations. These sign structures already show that the manuscript's main interpretation depends on distributional shape rather than on isolated outliers.

### 3.2. Regional synthesis of quantile-dependent change

Figure 3 and Table 1 summarize the regional network-mean quantile profiles together with station-level field-significance and dominant-cluster information. Four points stand out. First, warm days have the strongest regional increase overall, with the slope rising from `+9.45` days per decade at `q0.10` to `+22.50` at `q0.90`. Second, warm nights also rise across the distribution, but their upper-tail slope (`+12.21`) is close to their median slope (`+12.14`), yielding much weaker tail amplification than for warm days. Third, both cool indices decline regionally, with cool days showing the stronger contraction, especially at `q0.90` (`-15.37`). Fourth, the sign of `Delta1` is index-specific and highly informative: it is positive for the warm indices and negative for the cool indices, which indicates upper-tail amplification of warming-related frequencies but upper-tail contraction of cooling-related frequencies.

![Figure 3. Regional quantile-coefficient panels.](../outputs/figures/ijoc_regional_quantile_panels.png)

*Figure 3. Regional network-mean quantile-regression profiles for warm days, warm nights, cool days, and cool nights. Each panel is fitted to the annual network-mean series and is descriptive rather than a substitute for station-level inference. Source note: [ijoc_regional_quantile_panels.md](../outputs/output_docs/figures/ijoc_regional_quantile_panels.md).*

**Table 1.** Regional synthesis of network-mean slopes, median-quantile field significance, and dominant-cluster behavior.

| Index | OLS | `q0.10` | `q0.50` | `q0.90` | `Delta1` | FDR-retained stations at `q0.50` | Dominant cluster |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Warm days | `+14.02` | `+9.45` | `+15.55` | `+22.50` | `+13.06` | `115/124` | Cluster 3 (`n = 69`, median `Delta1 = +14.48`) |
| Warm nights | `+10.96` | `+8.74` | `+12.14` | `+12.21` | `+3.47` | `104/124` | Cluster 4 (`n = 95`, median `Delta1 = +7.74`) |
| Cool days | `-12.60` | `-5.60` | `-12.00` | `-15.37` | `-9.77` | `97/124` | Cluster 2 (`n = 74`, median `Delta1 = -10.47`) |
| Cool nights | `-10.86` | `-5.91` | `-9.83` | `-12.06` | `-6.16` | `86/124` | Cluster 2 (`n = 114`, median `Delta1 = -8.75`) |

*Source note: regional slopes from [homogeneity_flag_exclusion_sensitivity.md](../outputs/output_docs/tables/homogeneity_flag_exclusion_sensitivity.md), median-quantile FDR counts from [station_significance_fdr.md](../outputs/output_docs/tables/station_significance_fdr.md), and cluster medians from [regional_cluster_composites.md](../outputs/output_docs/tables/regional_cluster_composites.md).*

The station-level synthesis reinforces the regional profiles. `Delta1` is positive at `118/124` warm-day stations and `107/124` warm-night stations, but negative at all 124 cool-day stations and at `116/124` cool-night stations. Bootstrap tail inference also shows substantial support for these directional changes: significant bootstrap intervals occur at `q0.90` for 73 warm-day, 88 warm-night, 57 cool-day, and 61 cool-night stations, while `q0.10` significance occurs at 47, 69, 74, and 68 stations, respectively. Field-significance control at `q0.50` further indicates that the regional signal is not an artifact of a few stations: the median slope remains significant after FDR control at most stations for all four indices.

### 3.3. Day-night asymmetry

The paired day-night comparison changes the physical interpretation in an important way. In this 1991-2024 network, the dominant warming signal is daytime rather than nighttime. Warm nights do increase, but warm-day slopes are larger at every regional summary metric. The warm-night minus warm-day contrast is `-3.06` for OLS, `-0.70` at `q0.10`, `-3.41` at `q0.50`, and `-10.29` at `q0.90`. The corresponding `Delta1` contrast is `-9.59`, showing that upper-tail amplification is much more strongly expressed in warm days than in warm nights.

For the cooling-related pair, the pattern is reversed. Cool nights also decline, but cool days decline more strongly, especially above the median. Cool-night minus cool-day contrasts are `+1.74` for OLS, `-0.30` at `q0.10`, `+2.17` at `q0.50`, and `+3.31` at `q0.90`, with a `Delta1` contrast of `+3.61`. In other words, the strongest loss of cool conditions is concentrated in cool days rather than cool nights, and the divergence becomes larger toward the upper part of the distribution.

**Table 2.** Regional day-night asymmetry, expressed as warm-night minus warm-day and cool-night minus cool-day.

| Metric | Warm-night minus warm-day | Cool-night minus cool-day |
| --- | ---: | ---: |
| OLS | `-3.06` | `+1.74` |
| `q0.10` | `-0.70` | `-0.30` |
| `q0.50` | `-3.41` | `+2.17` |
| `q0.90` | `-10.29` | `+3.31` |
| `Delta1` | `-9.59` | `+3.61` |

This asymmetry matters because it suggests that daytime warm-extreme frequency is the clearest leading edge of recent thermal change in Iran, whereas nighttime warming, although widespread, is more distributionally even. Likewise, the decline of cool conditions is more aggressive in the daytime indices, implying that the disappearance of unusually cool years is not simply a nocturnal phenomenon.

### 3.4. Representative stations and cluster-defined quantile geometries

Regionalization shows that the four indices do not organize in the same way. Warm days are distributed across four clusters with sizes `7`, `47`, `69`, and `1`, indicating strong internal diversity. Warm nights are much more concentrated, with a dominant cluster of 95 stations. Cool days are mainly divided between two larger groups (`45` and `74` stations), while cool nights are overwhelmingly concentrated in one cluster of 114 stations, with only a few highly atypical stations defining the remaining groups. This difference in cluster balance is already informative: daytime and cool-day responses are regionally differentiated, whereas cool-night and especially warm-night responses are more strongly dominated by a broad common pattern.

**Table 3.** Reproducibly selected representative stations, defined as the observed stations closest to the baseline cluster centroids.

| Index | Cluster | Representative station | Distance | `q0.10` | `q0.50` | `q0.90` | `Delta1` |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: |
| Warm days | 1 | Garmsar | `0.407` | `+4.44` | `+3.57` | `-0.00` | `-4.44` |
| Warm days | 2 | Fasa | `0.345` | `+7.14` | `+12.73` | `+13.94` | `+6.80` |
| Warm days | 3 | Omidiyeh (Aghajari) | `0.303` | `+10.00` | `+16.00` | `+27.14` | `+17.14` |
| Warm days | 4 | Esfahan | `0.000` | `+7.89` | `+13.33` | `+42.00` | `+34.11` |
| Warm nights | 1 | Dorudzan | `0.504` | `+1.18` | `-4.50` | `-4.50` | `-5.68` |
| Warm nights | 2 | Abadan | `0.423` | `+16.19` | `+16.67` | `+17.41` | `+1.22` |
| Warm nights | 3 | Semnan | `0.275` | `+5.36` | `+7.06` | `+7.00` | `+1.64` |
| Warm nights | 4 | Bushehr (Coastal) | `0.190` | `+8.18` | `+12.31` | `+15.20` | `+7.02` |
| Cool days | 1 | Arak | `0.333` | `-7.08` | `-12.11` | `-20.00` | `-12.92` |
| Cool days | 2 | Gorgan | `0.272` | `-5.00` | `-7.86` | `-15.00` | `-10.00` |
| Cool days | 3 | Abumusa Island | `0.495` | `-7.50` | `-20.00` | `-31.56` | `-24.06` |
| Cool days | 4 | Chahbahar | `0.834` | `-4.17` | `-9.50` | `-37.04` | `-32.87` |
| Cool nights | 1 | Dorudzan | `0.138` | `+6.32` | `+9.23` | `+2.86` | `-3.46` |
| Cool nights | 2 | Biyarjomand | `0.102` | `-7.33` | `-9.64` | `-16.79` | `-9.45` |
| Cool nights | 3 | Kish Island | `0.000` | `-0.91` | `-14.44` | `-49.69` | `-48.78` |
| Cool nights | 4 | Garmsar | `0.000` | `-26.36` | `-19.57` | `-24.64` | `+1.72` |

*Source note: [representative_station_selection.md](../outputs/output_docs/tables/representative_station_selection.md).*

![Figure 4. Representative-station comparisons.](../outputs/figures/ijoc_station_comparisons/main_figure_representative_stations.png)

*Figure 4. Representative-station panels for the four indices, with one station shown for each cluster according to the closest-observed-member rule. Source note: [main_figure_representative_stations.md](../outputs/output_docs/figures/ijoc_station_comparisons/main_figure_representative_stations.md).*

The representative stations clarify how the clusters differ physically. For warm days, Esfahan defines a distinct single-station cluster with extreme upper-tail amplification (`Delta1 = +34.11`), while Garmsar captures the only representative warm-day geometry with a negative `Delta1`, indicating little or no upper-tail intensification. Warm-night behavior is more consolidated: Bushehr (Coastal) typifies the dominant positive pattern, whereas Dorudzan represents a minority cluster with negative median and upper-tail slopes. For cool days, Chahbahar and Abumusa Island illustrate the most severe upper-tail contraction, while Arak and Gorgan represent more moderate but still clearly negative continental profiles. For cool nights, Biyarjomand captures the network's dominant declining pattern, Kish Island represents an extreme upper-tail-collapse outlier, and Garmsar stands out as a rare case with positive `Delta1` because its lower-quantile decline is even stronger than its upper-quantile decline.

### 3.5. Spatial expression and regionalization

The spatial pattern of `Delta1` confirms that distributional asymmetry is geographically organized rather than random. Warm-day `Delta1` is positive at `118/124` stations, and the largest values occur at Esfahan (`+34.11`), Zanjan (`+28.44`), Daran (`+28.15`), Abumusa Island (`+27.42`), and Bandar-E-Anzali (`+26.15`). Warm-night `Delta1` is also positive at most stations (`107/124`), but its geography is more mixed: large positive values occur at Zahak (`+21.05`), Khorramdareh (`+17.23`), Zanjan (`+14.90`), Nowshahr (`+14.48`), and Sarakhs (`+14.23`), while several stations show clearly negative values, including Ilam (`-18.25`) and Shiraz (`-16.45`). This mixed structure explains why warm nights rise strongly overall yet remain less sharply amplified in the upper tail than warm days.

For the cooling-related indices, the spatial pattern is even clearer. All 124 cool-day stations have negative `Delta1`, with the most extreme contractions at Bandar-E-Lengeh (`-33.70`), Chahbahar (`-32.87`), Bushehr (Coastal) (`-30.53`), Orumiyeh (`-26.39`), and Jask (`-25.54`). Cool-night `Delta1` is negative at `116/124` stations, and the most severe values occur at Kish Island (`-48.78`), Esfahan (`-25.79`), Sabzevar (`-25.64`), Shiraz (`-24.75`), and Abumusa Island (`-24.27`). In combination, these maps show that the most rapid erosion of cool extremes is concentrated at a subset of low-latitude, coastal, and inland-hot stations rather than being spatially uniform.

![Figure 5. Main Delta1 map panels.](../outputs/figures/ijoc_main_delta1_maps.png)

*Figure 5. Station maps of `Delta1 = slope(q0.90) - slope(q0.10)` for warm days, warm nights, cool days, and cool nights. Positive values indicate upper-tail amplification; negative values indicate relative upper-tail weakening or contraction. Source note: [ijoc_main_delta1_maps.md](../outputs/output_docs/figures/ijoc_main_delta1_maps.md).*

The cluster sizes provide a useful second view of this geography. The dominant warm-night cluster contains 95 stations and the dominant cool-night cluster contains 114 stations, implying that the night indices are characterized by large common regional modes with only a few atypical departures. Warm days and cool days, by contrast, partition more evenly among multiple groups, which is consistent with stronger regional differentiation and with the more distinct day-side spatial patterns visible in Figure 5.

### 3.5.1. Additional structural and spatial diagnostics

The spatial diagnostics sharpen this picture. Moran's I is significant in `6/12` index-quantile fields. Cool days show the strongest and most consistent spatial organization, with Moran's I of `0.290` at `q0.10`, `0.192` at `q0.50`, and `0.348` at `q0.90` (all `p = 0.002-0.004`). Warm days also show significant spatial structure at the tails (`0.258` at `q0.10`, `p = 0.002`; `0.166` at `q0.90`, `p = 0.008`) but not at the median. Warm nights show significant lower-tail structure only (`0.200`, `p = 0.004` at `q0.10`), whereas cool nights do not reach significance at `q0.10` or `q0.50` and are only marginal at `q0.90` (`p = 0.058`). This pattern means that spatial dependence is clearly real, but it is quantile- and index-specific rather than uniform across the entire analysis.

The permutation-based spatial compactness test yields a complementary result. Warm-day and cool-day clusters are more geographically compact than random labelings (`p = 0.002` and `0.008`, respectively), whereas warm-night and cool-night clusters are not (`p = 0.882` and `0.894`). This indicates that day-index regionalization corresponds more closely to physically contiguous spatial groupings, while night-index clustering often captures functionally similar quantile geometry across geographically dispersed stations. Supplementary dendrograms and quantile-map panels reinforce the same point: station similarity changes across quantiles, and the apparent spatial field is more coherent for some indices and quantiles than for others. Supporting visual summaries are provided in Supplementary Figures S2-S7 and S12-S13, and tabular diagnostics are provided in Supplementary Tables S10-S11.

### 3.6. Stronger explanatory context from geography and elevation

The driver analysis links these spatial contrasts to broad physiographic gradients. Warm-day `Delta1` decreases significantly with longitude (`std beta = -0.395`, `p < 0.001`; `rho = -0.426`), indicating that the strongest upper-tail amplification tends to occur farther west in the network, although several southern stations also display large positive values. Warm-night `Delta1` decreases with elevation (`std beta = -0.230`, `p = 0.012`; `rho = -0.244`), which means that upper-tail warm-night amplification is strongest at lower-elevation sites and weaker at highland stations.

Cool-day `Delta1` shows the clearest geographical structure of all four indices. Its association with latitude (`std beta = +0.450`, `p < 0.001`), longitude (`+0.410`, `p < 0.001`), and elevation (`+0.263`, `p = 0.001`) implies that the most negative `Delta1` values, and thus the strongest upper-tail loss of cool days, occur preferentially at lower, more southerly, and more westerly stations. Cool-night `Delta1` also becomes less negative with elevation (`std beta = +0.269`, `p = 0.003`; `rho = +0.295`), pointing again to especially strong contraction of cool-night extremes in low-elevation settings.

These relationships support a physically interpretable pattern. Low-elevation stations, many of them coastal or interior-hot, tend to show either stronger amplification of warm extremes or stronger contraction of cool extremes. Highland and northern stations still warm, but their distributional geometry is often flatter and less dominated by upper-tail extremes. The geographical controls are not strong enough to explain all station-level behavior by themselves, but they do indicate that the distributional trend structure is not arbitrary; it is meaningfully conditioned by Iran's large-scale physiography.

### 3.7. Robustness, methodological sensitivity, and temporal nonlinearity

The main directional conclusions remain stable under multiple sensitivity checks, but the magnitude of some summaries is method-dependent. Excluding all 39 stations flagged by any detrended homogeneity diagnostic leaves 85 stations and changes the regional network-mean `Delta1` values only modestly: warm days shift from `+13.06` to `+13.05`, warm nights from `+3.47` to `+3.45`, cool days from `-9.77` to `-10.39`, and cool nights from `-6.16` to `-4.87`. Thus, the manuscript's sign structure and overall ordering are not driven by the flagged subset alone, even though the exact magnitude of the cool-index asymmetry remains somewhat sensitive to station exclusion.

Bootstrap-depth sensitivity is likewise reassuring. Increasing the moving-block bootstrap from 200 to 400 replicates yields mean regional differences of only `-0.12` to `+0.09` across the four primary bootstrap mean metrics, with median absolute station-level differences mostly between `0.10` and `0.30` and station-to-station correlations above `0.996` in every case. By contrast, bootstrap-method sensitivity is more substantial in the tails: the mean absolute difference between moving-block and `meboot` summaries is about `0.97-1.16` days per decade at `q0.10` and `1.68-3.34` at `q0.90`. The implication is not that the results are unstable in sign, but that tail uncertainty magnitudes depend materially on bootstrap design. Interpolation comparisons also show that mapped surfaces are more reproducible between some methods than others. For example, linear RBF and linear interpolation agree strongly, with surface correlations of `0.84-0.93`, whereas the thin-plate spline solutions diverge much more from those alternatives. This supports interpreting the maps as spatial displays rather than as independent inferential products.

Regionalization robustness is uneven and therefore informative. In the reduced-feature rerun, agreement with the baseline clustering is very high for cool nights (`ARI = 0.94`), moderate for warm nights (`0.49`), and low for warm days (`0.15`) and cool days (`0.06`). Yet the alternative-method comparison shows that warm-day regionalization remains reasonably reproducible across several clustering algorithms (`ARI = 0.63-0.98`), warm nights are highly stable under the average-cityblock alternative (`0.94`) but less so under complete linkage and k-means, cool days remain only moderately reproducible (`0.36-0.46`), and cool nights range from exact reproduction under the average-cityblock variant (`1.00`) to low agreement under other alternatives (`0.12-0.21`). Taken together, these checks show that clustering is most robust for the dominant cool-night pattern and least uniquely determined for cool days, while warm-day clustering is sensitive to feature selection but still retains substantial structure across alternative hierarchical specifications.

In the current run, no alternative reference-period comparison was generated because the configured baseline and alternative settings collapse to the same effective period in the available dataset. Figure 6 therefore summarizes only the robustness components that produced distinct comparison outputs, namely bootstrap-method sensitivity, interpolation sensitivity, and clustering robustness.

![Figure 6. Robustness synthesis.](../outputs/figures/ijoc_robustness_synthesis.png)

*Figure 6. Multi-panel synthesis of the available methodological-sensitivity and clustering-robustness outputs generated in the current run. Depending on which comparison products exist, panels are shown dynamically for the corresponding sensitivity components. Source note: [ijoc_robustness_synthesis.md](../outputs/output_docs/figures/ijoc_robustness_synthesis.md).*

The split-period analysis shows temporal nonlinearity rather than a simple later-period intensification at every quantile. For warm days, the regional OLS slope increases from `+10.94` in `1991-2007` to `+20.01` in `2008-2024`, and the median quantile also strengthens (`+10.43` to `+24.50`), but the upper-tail slope weakens (`+23.65` to `+14.97`). Warm nights show a similar contrast: the OLS slope rises from `+11.15` to `+14.22`, while the upper-tail slope falls from `+15.60` to `+6.36`. For the cool indices, the strongest negative OLS slopes occur in the earlier half of the record, decreasing from `-27.62` to `-10.91` for cool days and from `-21.83` to `-10.47` for cool nights. These results indicate that the 34-year record is not well described by a single monotonic acceleration narrative. Instead, central and mean warming strengthen in the later half, but the tail geometry becomes more irregular, and the strongest cooling-related contractions are concentrated in the earlier subperiod.

![Figure 7. Split-period comparison.](../outputs/figures/ijoc_split_period_comparison.png)

*Figure 7. Comparison of regional network-mean `q0.10`, `q0.50`, `q0.90`, and OLS slopes between `1991-2007` and `2008-2024`. The figure documents temporal nonlinearity rather than uniform amplification across every quantile. Source note: [ijoc_split_period_comparison.md](../outputs/output_docs/figures/ijoc_split_period_comparison.md).*

Overall, the robustness results support the central conclusions while narrowing their proper scope. The signs and rank ordering of the main regional summaries are stable, but the exact magnitudes of tail contrasts and some clustering partitions should be interpreted as conditional on methodological choices rather than as uniquely fixed climatological constants. Supporting robustness outputs are compiled in Supplementary Figures S8-S11 and Supplementary Tables S4-S9.

### 3.8. Comparison with the wider climatological literature

The new results align with the broader literature in several ways while also refining it for Iran. First, the widespread increase in warm-event frequencies and decline in cool-event frequencies is consistent with prior evidence of intensifying heat extremes across Iran and the wider MENA domain `(Zhang et al., 2005a; Soltani et al., 2016; Vaghefi et al., 2019; Zittis et al., 2016; Francis and Fonseca, 2024)`. Second, the strong contraction of cool days and cool nights is consistent with the global reduction of cold extremes documented in large observational syntheses `(Alexander et al., 2006; Donat et al., 2013; Dunn et al., 2020)`.

The more distinctive contribution of the present analysis is that it shows recent Iranian warming to be strongly distributionally asymmetric and, in this dataset, led more clearly by daytime warm extremes than by nighttime ones. That result is physically plausible in an arid and topographically complex setting where land-surface feedbacks, cloud variability, and radiative conditions can favor strong daytime amplification `(Seneviratne et al., 2010; Sanchez-Lorenzo et al., 2017; Perkins-Kirkpatrick and Lewis, 2020)`. At the same time, the low-elevation concentration of strong warm-night amplification and strong cool-night contraction is consistent with the emerging importance of nocturnal and humid-heat risk in southern and coastal southwest Asia `(Pal and Eltahir, 2016; Raymond et al., 2020; Raymond et al., 2024)`.

The results also fit the broader quantile-regression literature. Earlier studies showed that quantile-specific temperature trends can reveal structures hidden by OLS or mean-only analysis `(Barbosa et al., 2011; Fan, 2014)`. The present study extends that logic to a dense national station network in Iran and shows that the same principle applies strongly in a dryland setting: warm-day amplification, warm-night flattening, and uneven cool-extreme contraction would all be partially obscured by mean trends alone. The manuscript therefore adds regional specificity to the wider literature by demonstrating that recent thermal-extreme change across Iran is not only warming in sign but also structurally differentiated in shape.

### 3.9. Limitations

Several limitations define the appropriate scope of inference. First, the record length is 34 years, which is sufficient for consistent recent-trend analysis but shorter than ideal for very stable tail estimation or for strong claims about multidecadal nonlinearity. Second, the thermal indices are built with thresholds derived from the available 1991-2024 record rather than from an independent earlier baseline, so the annual counts are distribution-relative within the analyzed period. Third, the daily records were screened for quality and homogeneity but were not formally homogenized with a metadata-supported daily adjustment procedure, which means some local tail magnitudes may still reflect unresolved non-climatic shifts.

Fourth, the main trend framework is station-wise rather than fully spatial, so covariance enters the study through diagnostics and regionalization rather than through a hierarchical spatial quantile model `(Reich, 2012)`. Fifth, the analytic tail outputs were insufficient for national FDR screening at `q0.10` and `q0.90`, so field-significance interpretation is anchored mainly in the median quantile while the tails rely on bootstrap intervals. Sixth, several clustering results are exploratory by design. Their stability varies by index, and the night-index clusters are not especially geographically compact relative to random labelings. Seventh, the driver analysis is intentionally simple and does not directly model circulation, soil moisture, humidity, radiation, or land-use effects. Finally, interpolated surfaces are used for display rather than as primary evidence, and the split-period comparison is descriptive rather than a formal breakpoint-detection framework.

These limitations do not overturn the main conclusions, but they do bound them. The manuscript is best interpreted as a robust station-based synthesis of recent distributional change across Iran rather than as a complete attribution or fully spatially explicit quantile model. The most useful next steps would combine metadata-supported daily homogenization, explicit spatial quantile modeling, and physically richer diagnostics of circulation, radiation, humidity, and land-surface processes.

## 4. Conclusions

This study demonstrates that recent thermal-extreme change across Iran is fundamentally distributional, not merely mean-based. Using 124 meteorological stations for 1991-2024, the analysis shows that warm days constitute the strongest warming-related signal, with regional network-mean slopes increasing from `+9.45` days per decade at `q0.10` to `+22.50` at `q0.90`, yielding `Delta1 = +13.06`. Warm nights also increase across the country, but more evenly across the distribution (`+8.74`, `+12.14`, and `+12.21` days per decade at `q0.10`, `q0.50`, and `q0.90`; `Delta1 = +3.47`). In contrast, cool days and cool nights both decline, with stronger upper-quantile contraction for cool days (`-5.60`, `-12.00`, `-15.37`; `Delta1 = -9.77`) than for cool nights (`-5.91`, `-9.83`, `-12.06`; `Delta1 = -6.16`). The central empirical message is therefore not only that Iran has warmed, but that the internal shape of annual thermal-extreme distributions has changed in a systematic and index-specific manner.

Several contributions of the study are especially important. First, the manuscript provides a dense station-based national assessment in which inference remains anchored in observations rather than in interpolated fields. Second, it combines quantile regression, bootstrap uncertainty estimation, field-significance screening, spatial autocorrelation diagnostics, and cluster-based regionalization within a single reproducible workflow. Third, it shows that the main results are not confined to a small subset of stations: at `q0.50`, FDR control retained `115/124` warm-day, `104/124` warm-night, `97/124` cool-day, and `86/124` cool-night station results. Fourth, it reveals that spatial organization differs by index. Cool days show the most persistent spatial autocorrelation, whereas warm days exhibit the clearest upper-tail amplification and the strongest evidence of geographically compact cluster structure. Warm nights and cool nights, by contrast, are characterized by larger dominant clusters, indicating broad common modes with fewer highly atypical departures.

These findings matter beyond descriptive trend reporting. Scientifically, they show that recent thermal change across Iran cannot be reduced to a single smooth national field or a single representative mean trend. The strongest changes emerge in the relative behavior of lower, central, and upper quantiles. Practically, this implies that heat-related risk assessment, agricultural planning, water-resource management, and climate adaptation strategies should pay particular attention to the increasing frequency of exceptionally warm-day years and the simultaneous erosion of cool-day and cool-night conditions. The results also suggest that low-elevation and several coastal or interior-hot stations are especially exposed to strong distributional shifts, making quantile-aware monitoring relevant for regional early warning and climate-risk prioritization.

More broadly, the study supports a methodological argument with wider relevance for climatology. In arid, topographically complex, and spatially heterogeneous environments, mean-only analyses can obscure the most consequential parts of the thermal signal. A quantile-aware, station-based framework is better suited to reveal whether change is distribution-wide, upper-tail dominated, or regionally uneven. In that sense, the present workflow offers a transferable template for other dryland and mountainous regions where national climate-change assessments may otherwise smooth over precisely the contrasts that matter most for extremes.

At the same time, the conclusions should be interpreted within the limits of the current dataset. The record length is 34 years, the thermal thresholds are derived from the available 1991-2024 period, the daily series were screened but not fully homogenized, and some regionalization outcomes remain exploratory rather than uniquely resolved. Future work should therefore extend this analysis with metadata-supported daily homogenization, explicit spatial quantile models, and physically richer diagnostics involving circulation, humidity, cloudiness, soil moisture, and surface energy balance. It would also be valuable to examine compound heat metrics and event-based heatwave characteristics so that the distributional changes identified here can be linked more directly to human and environmental impact pathways. Even with these caveats, the present evidence is clear: recent thermal-extreme change across Iran is asymmetric, regionally differentiated, and most fully understood when tail behavior is treated as a primary component of climate inference rather than as a secondary descriptive detail.

## Acknowledgements

Funding, acknowledgement, and data-availability statements should be inserted here in the journal's required format.

## Supplementary material

The Supplementary Material should contain at minimum the following items cited in the manuscript:

1. `Supplementary Figure S1`: [ijoc_data_quality_homogeneity.md](../outputs/output_docs/figures/ijoc_data_quality_homogeneity.md)
2. `Supplementary Figures S2-S4`: [figure_5_tau_0.10.md](../outputs/output_docs/figures/paper1_quantile_dendrograms/figure_5_tau_0.10.md), [figure_6_tau_0.50.md](../outputs/output_docs/figures/paper1_quantile_dendrograms/figure_6_tau_0.50.md), [figure_7_tau_0.90.md](../outputs/output_docs/figures/paper1_quantile_dendrograms/figure_7_tau_0.90.md)
3. `Supplementary Figures S5-S7`: [figure3_tau_0.10.md](../outputs/output_docs/figures/paper2_figure3_quantile_maps/figure3_tau_0.10.md), [figure3_tau_0.50.md](../outputs/output_docs/figures/paper2_figure3_quantile_maps/figure3_tau_0.50.md), [figure3_tau_0.90.md](../outputs/output_docs/figures/paper2_figure3_quantile_maps/figure3_tau_0.90.md)
4. `Supplementary Figure S8`: [ijoc_homogeneity_sensitivity.md](../outputs/output_docs/figures/ijoc_homogeneity_sensitivity.md)
5. `Supplementary Figure S9`: [ijoc_bootstrap_depth_sensitivity.md](../outputs/output_docs/figures/ijoc_bootstrap_depth_sensitivity.md)
6. `Supplementary Figure S10`: [ijoc_alternative_clustering_sensitivity.md](../outputs/output_docs/figures/ijoc_alternative_clustering_sensitivity.md)
7. `Supplementary Figure S11`: [cluster_spatial_validation.md](../outputs/output_docs/figures/advanced_regionalization/cluster_spatial_validation.md)
8. `Supplementary Figure S12`: [raw_vs_fdr_counts.md](../outputs/output_docs/figures/advanced_spatial_inference/raw_vs_fdr_counts.md)
9. `Supplementary Figure S13`: [moran_i_heatmap.md](../outputs/output_docs/figures/advanced_spatial_inference/moran_i_heatmap.md)
10. `Supplementary Tables S1-S3`: [data_quality_homogeneity_overview.md](../outputs/output_docs/tables/data_quality_homogeneity_overview.md), [data_quality_station_summary.md](../outputs/output_docs/tables/data_quality_station_summary.md), [data_homogeneity_tests_station_summary.md](../outputs/output_docs/tables/data_homogeneity_tests_station_summary.md)
11. `Supplementary Table S4`: [homogeneity_flag_exclusion_sensitivity.md](../outputs/output_docs/tables/homogeneity_flag_exclusion_sensitivity.md)
12. `Supplementary Tables S5-S6`: [bootstrap_depth_sensitivity_summary.md](../outputs/output_docs/tables/bootstrap_depth_sensitivity_summary.md), [bootstrap_depth_sensitivity_station_comparison.md](../outputs/output_docs/tables/bootstrap_depth_sensitivity_station_comparison.md)
13. `Supplementary Tables S7-S8`: [alternative_clustering_sensitivity_summary.md](../outputs/output_docs/tables/alternative_clustering_sensitivity_summary.md), [alternative_clustering_assignments.md](../outputs/output_docs/tables/alternative_clustering_assignments.md)
14. `Supplementary Table S9`: [regional_cluster_spatial_validation.md](../outputs/output_docs/tables/regional_cluster_spatial_validation.md)
15. `Supplementary Table S10`: [station_significance_fdr.md](../outputs/output_docs/tables/station_significance_fdr.md)
16. `Supplementary Table S11`: [spatial_autocorrelation_moran.md](../outputs/output_docs/tables/spatial_autocorrelation_moran.md)

## References

Alexander LV, Zhang X, Peterson TC, Caesar J, Gleason B, Klein Tank AMG, Haylock M, Collins D, Trewin B, Rahimzadeh F, Tagipour A, Kumar KR, Revadekar J, Griffiths G, Vincent L, Stephenson DB, Burn J, Aguilar E, Brunet M, Taylor M, New M, Zhai P, Rusticucci M, Vazquez-Aguirre JL (2006) Global observed changes in daily climate extremes of temperature and precipitation. *Journal of Geophysical Research: Atmospheres* 111:D05109. https://doi.org/10.1029/2005JD006290

Barbosa SM, Scotto MG, Alonso AM (2011) Summarising changes in air temperature over Central Europe by quantile regression and clustering. *Natural Hazards and Earth System Sciences* 11:3227-3233. https://doi.org/10.5194/nhess-11-3227-2011

Benjamini Y, Hochberg Y (1995) Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B* 57(1):289-300.

Donat MG, Alexander LV, Yang H, Durre I, Vose R, Caesar J, Gleason BE, Klein Tank AMG, Manton MJ, Molanus C, Peterson TC, Renom M, Brunet M, Caesar J, Zhang X, Kitching S, Haylock M (2013) Updated analyses of temperature and precipitation extreme indices since the beginning of the twentieth century: the HadEX2 dataset. *Journal of Geophysical Research: Atmospheres* 118(5):2098-2118. https://doi.org/10.1002/jgrd.50150

Easterling DR, Meehl GA, Parmesan C, Changnon SA, Karl TR, Mearns LO (2000) Climate extremes: observations, modeling, and impacts. *Science* 289(5487):2068-2074. https://doi.org/10.1126/science.289.5487.2068

Efron B, Tibshirani RJ (1993) *An Introduction to the Bootstrap*. Chapman & Hall, New York.

Frich P, Alexander LV, Della-Marta P, Gleason B, Haylock M, Klein Tank AMG, Peterson T (2002) Observed coherent changes in climatic extremes during the second half of the twentieth century. *Climate Research* 19:193-212. https://doi.org/10.3354/cr019193

Hall P, Horowitz JL, Jing B-Y (1995) On blocking rules for the bootstrap with dependent data. *Biometrika* 82(3):561-574. https://doi.org/10.1093/biomet/82.3.561

Hubert L, Arabie P (1985) Comparing partitions. *Journal of Classification* 2:193-218. https://doi.org/10.1007/BF01908075

IPCC (2021) Weather and climate extreme events in a changing climate. In: Masson-Delmotte V, Zhai P, Pirani A, Connors SL, Péan C, Berger S, Caud N, Chen Y, Goldfarb L, Gomis MI, Huang M, Leitzell K, Lonnoy E, Matthews JBR, Maycock TK, Waterfield T, Yelekçi O, Yu R, Zhou B (eds) *Climate Change 2021: The Physical Science Basis*. Cambridge University Press, Cambridge, Chapter 11. https://www.ipcc.ch/report/ar6/wg1/chapter/chapter-11/

Katz RW, Brown BG (1992) Extreme events in a changing climate: variability is more important than averages. *Climatic Change* 21(3):289-302. https://doi.org/10.1007/BF00139728

Katz RW, Parlange MB, Naveau P (2002) Statistics of extremes in hydrology. *Advances in Water Resources* 25(8-12):1287-1304. https://doi.org/10.1016/S0309-1708(02)00056-8

Kunsch HR (1989) The jackknife and the bootstrap for general stationary observations. *The Annals of Statistics* 17(3):1217-1241. https://doi.org/10.1214/aos/1176347265

Koenker R (2005) *Quantile Regression*. Cambridge University Press, Cambridge.

Koenker R, Bassett G Jr (1978) Regression quantiles. *Econometrica* 46(1):33-50. https://doi.org/10.2307/1913643

Kostopoulou E, Jones PD (2005) Assessment of climate extremes in the eastern Mediterranean. *Meteorology and Atmospheric Physics* 89:69-85. https://doi.org/10.1007/s00703-005-0122-2

Moberg A, Jones PD (2005) Trends in indices for extremes in daily temperature and precipitation in central and western Europe, 1901-99. *International Journal of Climatology* 25(9):1149-1171. https://doi.org/10.1002/joc.1163

Moran PAP (1950) Notes on continuous stochastic phenomena. *Biometrika* 37(1-2):17-23. https://doi.org/10.1093/biomet/37.1-2.17

Lahiri SN (2003) *Resampling Methods for Dependent Data*. Springer, New York.

Nairn JR, Fawcett RJB (2015) The excess heat factor: a metric for heatwave intensity and its use in classifying heatwave severity. *International Journal of Environmental Research and Public Health* 12(1):227-253. https://doi.org/10.3390/ijerph120100227

Pal JS, Eltahir EAB (2016) Future temperature in southwest Asia projected to exceed a threshold for human adaptability. *Nature Climate Change* 6:197-200. https://doi.org/10.1038/nclimate2833

Perkins SE (2015) A review on the scientific understanding of heatwaves: their measurement, driving mechanisms, and changes at the global scale. *Atmospheric Research* 164-165:242-267. https://doi.org/10.1016/j.atmosres.2015.05.014

Perkins-Kirkpatrick SE, Lewis SC (2020) Increasing trends in regional heatwaves. *Nature Communications* 11:3357. https://doi.org/10.1038/s41467-020-16970-7

Raei E, Nikoo MR, AghaKouchak A, Mazdiyasni O, Sadegh M (2018) GHWR, a multi-method global heatwave and warm-spell record and toolbox. *Scientific Data* 5:180206. https://doi.org/10.1038/sdata.2018.206

Raymond C, Matthews T, Horton RM (2020) The emergence of heat and humidity too severe for human tolerance. *Science Advances* 6(19):eaaw1838. https://doi.org/10.1126/sciadv.aaw1838

Raymond C, Matthews T, Tuholske C (2024) Evening humid-heat maxima near the southern Persian/Arabian Gulf. *Communications Earth & Environment* 5:591. https://doi.org/10.1038/s43247-024-01763-3

Reich BJ (2012) Spatiotemporal quantile regression for detecting distributional changes in environmental processes. *Journal of the Royal Statistical Society: Series C (Applied Statistics)* 61(4):535-553. https://doi.org/10.1111/j.1467-9876.2011.01025.x

Russo S, Sterl A (2011) Global changes in indices describing moderate temperature extremes from the daily output of a climate model. *Journal of Geophysical Research: Atmospheres* 116:D03104. https://doi.org/10.1029/2010JD014727

Sanchez-Lorenzo A, Enriquez-Alonso A, Calbó J, González JA, Wild M, Folini D, Norris JR, Vicente-Serrano SM (2017) Fewer clouds in the Mediterranean: consistency of observations and climate simulations. *Scientific Reports* 7:41475. https://doi.org/10.1038/srep41475

Seneviratne SI, Corti T, Davin EL, Hirschi M, Jaeger EB, Lehner I, Orlowsky B, Teuling AJ (2010) Investigating soil moisture-climate interactions in a changing climate: a review. *Earth-Science Reviews* 99(3-4):125-161. https://doi.org/10.1016/j.earscirev.2010.02.004

Soltani M, Laux P, Kunstmann H, Stan K, Sohrabi MM, Molanejad M, Sabziparvar AA, Ranjbar SaadatAbadi A, Ranjbar F, Rousta I, Zawar-Reza P, Khoshakhlagh F, Soltanzadeh I, Babu CA, Azizi GH, Martin MV (2016) Assessment of climate variations in temperature and precipitation extreme events over Iran. *Theoretical and Applied Climatology* 126:775-795. https://doi.org/10.1007/s00704-015-1609-5

Vaghefi SA, Keykhai M, Jahanbakhshi F, Sheikholeslami J, Ahmadi A, Yang H, Abbaspour KC (2019) The future of extreme climate in Iran. *Scientific Reports* 9:1464. https://doi.org/10.1038/s41598-018-38071-8

Vinod HD (2006) Maximum entropy ensembles for time series inference in economics. *Journal of Asian Economics* 17(6):955-978. https://doi.org/10.1016/j.asieco.2006.09.003

Ahmadalipour A, Moradkhani H (2018) Escalating heat-stress mortality risk due to global warming in the Middle East and North Africa (MENA). *Environment International* 117:215-225. https://doi.org/10.1016/j.envint.2018.05.014

Dunn RJH, Alexander LV, Donat MG, Zhang X, Bador M, Herold N, Lippmann T, Allan R, Aguilar E, Brunet M, Caesar J, Chagnaud G, Cheng V, Cinco T, Durre I, Htay MM, Hoang L, Hung NQ, Johnson F, Kruger A, Lau K, Leng TW, Loikith PC, Marengo J, Mbatha S, McGree S, Menne M, Skansi M, Trewin B, Villarroel C, Vincent LA, Vose RS, Yeo R, Zhang P (2020) Development of an updated global land in situ-based dataset of temperature and precipitation extremes: HadEX3. *Journal of Geophysical Research: Atmospheres* 125:e2019JD032263. https://doi.org/10.1029/2019JD032263

Fischer EM, Schär C (2010) Consistent geographical patterns of changes in high-impact European heatwaves. *Nature Geoscience* 3:398-403. https://doi.org/10.1038/ngeo866

Fan LJ (2014) Quantile trends in temperature extremes in China. *Atmospheric and Oceanic Science Letters* 7(4):304-308. https://doi.org/10.3878/j.issn.1674-2834.13.0102

Francis D, Fonseca R (2024) Recent and projected changes in climate patterns in the Middle East and North Africa (MENA) region. *Scientific Reports* 14:10279. https://doi.org/10.1038/s41598-024-60976-w

Zhang X, Aguilar E, Sensoy S, Melkonyan H, Tagiyeva U, Ahmed N, Kutaladze N, Rahimzadeh F, Taghipour A, Hantosh TH, Albert P, Semawi M, Karam Ali M, Al-Shabibi MHS, Al-Oulan Z, Zatari T, Al Dean Khelet I, Hamoud S, Sagir R, Demircan M, Eken M, Adiguzel M, Alexander LV, Peterson TC, Wallis T (2005a) Trends in Middle East climate extreme indices from 1950 to 2003. *Journal of Geophysical Research: Atmospheres* 110:D22104. https://doi.org/10.1029/2005JD006181

Zhang X, Hegerl G, Zwiers FW, Kenyon J (2005b) Avoiding inhomogeneity in percentile-based indices of temperature extremes. *Journal of Climate* 18(11):1641-1651. https://doi.org/10.1175/JCLI3366.1

Zhang X, Alexander L, Hegerl GC, Jones P, Klein Tank AMG, Peterson TC, Trewin B, Zwiers FW (2011) Indices for monitoring changes in extremes based on daily temperature and precipitation data. *WIREs Climate Change* 2(6):851-870. https://doi.org/10.1002/wcc.147

Zittis G, Hadjinicolaou P, Lelieveld J (2016) Strongly increasing heat extremes in the Middle East and North Africa (MENA) in the 21st century. *Climatic Change* 137:245-260. https://doi.org/10.1007/s10584-016-1665-6

Zittis G, Almazroui M, Alpert P, Ciais P, Cramer W, Dahdal Y, Fnais M, Francis D, Hadjinicolaou P, Howari FM, Kucera P, Kvalevåg MM, Lin L, Liu Z, Mihalopoulos N, Mostamandi S, Niang I, Noone K, Odoulami RC, Panitz HJ, Ratti C, Said F, Saleh K, Sielecki LE, Stenchikov G, Tsiros IX, Zittis C, Lelieveld J (2021) Business-as-usual will lead to super and ultra-extreme heatwaves in the Middle East and North Africa. *npj Climate and Atmospheric Science* 4:20. https://doi.org/10.1038/s41612-021-00178-7
