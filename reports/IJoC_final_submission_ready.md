# Quantile-dependent reorganization of thermal extremes across Iran during 1961-2024

**Running title:** Quantile-dependent thermal extremes across Iran

## Abstract

Mean trends alone are often insufficient to describe changes in climate extremes because the most consequential impacts emerge through shifts in the tails of the distribution. This study investigates quantile-dependent trends in four annual thermal-extreme indices across Iran during 1961-2024 using observations from 30 meteorological stations. Warm days, warm nights, cool days, and cool nights were analysed with quantile regression, while ordinary least squares (OLS) provided a benchmark for mean change. The workflow combined bootstrap uncertainty estimation, false-discovery-rate field-significance testing, cluster-based regionalization, spatial visualization, and sensitivity analysis for reference period, bootstrap method, interpolation, and clustering configuration. The results show that thermal-extreme change across Iran is not well represented as a uniform mean shift. Warm days and warm nights increased across most of the network, but with substantial quantile dependence and marked station-to-station heterogeneity. Warm days exhibited clear upper-tail amplification, whereas warm nights showed stronger overall increases but a more uniform regional structure. Cool days and cool nights declined, with the strongest negative changes concentrated in upper quantiles, implying rapid loss of years historically characterized by many cool events. Regional OLS trends were +9.53 and +13.58 days year^-1 per decade for warm days and warm nights, and -2.56 and -2.98 days year^-1 per decade for cool days and cool nights, respectively. Field significance was strongest at intermediate quantiles, whereas Moran's I showed weak national-scale spatial autocorrelation, indicating that the observed changes are regionally interpretable but not reducible to a single smooth spatial field. Split-period analysis further showed substantial intensification of warm extremes and stronger collapse of cool extremes after 1990. Overall, the results indicate that thermal-extreme change across Iran is best understood as a quantile-dependent reorganization of warm and cool event distributions, with wider relevance for climatological interpretation in arid and topographically complex regions.

**Keywords:** thermal extremes; quantile regression; climate change; warm nights; cool nights; regionalization; Iran

## 1. Introduction

Changes in climate extremes are among the most consequential manifestations of contemporary climate change because their impacts propagate through agriculture, public health, water resources, ecosystems, and infrastructure `(Easterling et al., 2000; IPCC, 2021)`. Mean temperature trends provide a useful first-order summary of climate change, but they are often insufficient for understanding extremes, since the most climatically relevant impacts depend not only on shifts in the centre of the distribution but also on changes in its tails `(Katz and Brown, 1992; Katz et al., 2002)`. A region may therefore exhibit moderate mean warming while simultaneously undergoing strong upper-tail intensification or rapid erosion of historically cool conditions. Distribution-aware approaches are consequently becoming increasingly important in climatology `(Koenker and Bassett, 1978; Koenker, 2005)`.

This issue is especially important for thermal-extreme indices. Warm and cool day-night indices have been widely used to characterize climate variability and change `(Frich et al., 2002; Alexander et al., 2006; Donat et al., 2013; Zhang et al., 2011)`, yet their interpretation often remains anchored in mean or median trends. Such summaries can obscure upper-tail amplification, lower-tail stabilization, and asymmetric decline across the distribution `(Moberg and Jones, 2005; Kostopoulou and Jones, 2005)`. In climatically heterogeneous regions, these effects may also vary sharply across space because of topographic contrasts, continentality, land-atmosphere coupling, and differences in moisture availability.

Iran offers a particularly useful setting for testing these issues. The country includes strong gradients in elevation, aridity, latitude, and marine influence, spanning humid Caspian lowlands, interior continental environments, hot southern lowlands, and high-elevation terrain. This diversity makes Iran a natural laboratory for examining whether observed thermal-extreme change behaves as a coherent national signal or as a set of locally modified regional responses. Previous work has documented changes in climate extremes across Iran and across the wider dryland and Mediterranean domain `(Kostopoulou and Jones, 2005; Zhang et al., 2005a; Soltani et al., 2016; Zittis et al., 2016)`, but relatively few station-based studies for Iran have resolved change across the full distribution while simultaneously addressing uncertainty, field significance, and regionalization.

The broader relevance of this question extends beyond Iran. A regional study is most useful when it contributes to wider climatological understanding rather than remaining purely descriptive. In that spirit, the present study does not ask only whether thermal extremes in Iran have changed. It asks whether those changes are best interpreted as a uniform mean shift or as a **quantile-dependent reorganization of the distribution**, and whether that reorganization is spatially coherent, physically interpretable, and robust to methodological choices.

Using 30 meteorological stations for 1961-2024, this study analyses four annual indices: warm days, warm nights, cool days, and cool nights. The specific contributions are fivefold. First, the study resolves trends across the full conditional distribution rather than relying on mean change alone. Second, it quantifies upper-tail amplification and contraction using tail-contrast metrics. Third, it integrates uncertainty through bootstrap-based diagnostics. Fourth, it uses clustering to regionalize quantile signatures rather than only mapping individual slopes. Fifth, it evaluates the robustness of the principal conclusions to alternative analytical choices. Collectively, these steps move the manuscript from a descriptive account of warming to a more explicit characterization of how the thermal-extremes distribution has been reorganized across Iran, with relevance for other arid and topographically complex regions.

## 2. Data and methods

### 2.1. Data and study region

The analysis is based on daily observations from 30 meteorological stations distributed across Iran. Station metadata include latitude, longitude, and elevation, thereby allowing the thermal-extreme signals to be interpreted in a physiographic context. After preprocessing, the common analysis period spans 1961-2024, yielding 64 annual values for each station-index series.

![Figure 1. Study area and station network.](../outputs/figures/ijoc_study_area.png)

*Figure 1. Distribution of the 30 meteorological stations used in the analysis across Iran. Stations are coloured by elevation to highlight the physiographic diversity of the network and the coexistence of lowland, mid-elevation, and highland environments.*

The station-based design is deliberate. National and regional summaries are derived from station estimates rather than gridded products, thereby preserving local climatic heterogeneity and avoiding interpolation assumptions at the primary inference stage. This choice is especially important in Iran, where elevation contrasts, continentality, and varying marine influence can generate strong short-range spatial differences in thermal behaviour. The input variables were daily minimum and maximum temperature, organized by station and calendar date, with leap-day observations removed before percentile-threshold estimation in order to preserve a consistent 365-day annual cycle.

Before index construction, the station records were screened for completeness, duplicate dates, and simple physical consistency, including `tmin <= tmax` and agreement of archived `tmean` with the daily temperature range where all three variables were available. Completeness was high across the network (median values about 98.7% for `tmin`, 99.2% for `tmax`, and 98.5% for `tmean` or a `tmin`-`tmax` surrogate), and no station contained duplicate calendar dates. Potential inhomogeneity was then screened using Pettitt-, SNHT-, and Buishand-type breakpoint diagnostics applied to annual mean temperature series, with linearly detrended versions used as the main screening tests to reduce confounding by long-term climatic trends. These diagnostics were used to identify potentially sensitive records for cautious interpretation rather than as automatic exclusion criteria, and no station was removed solely on the basis of a homogeneity flag. To assess the influence of flagged records on the principal results, a separate sensitivity analysis was conducted in which all stations flagged by any detrended homogeneity test were excluded from the regional summaries. A summary of completeness and homogeneity screening is provided in Supplementary Figure S1 and Supplementary Tables S1-S2, and the exclusion sensitivity is summarized in Supplementary Figure S10 and Supplementary Table S11.
### 2.2. Thermal-extreme indices

Four annual indices were constructed from daily minimum (`tmin`) and maximum (`tmax`) temperature observations: warm days, warm nights, cool days, and cool nights. Daily thresholds were estimated separately for each station using the 1961-1990 reference period, a 5-day moving calendar window, and the 10th and 90th percentiles of station-specific daily temperature distributions. A minimum reference-sample requirement was enforced at the day-of-year level to reduce instability in threshold estimation near data-sparse dates.

Warm days and cool days were defined from daily maximum temperature relative to the station-specific 90th and 10th percentile thresholds, respectively. Warm nights and cool nights were defined analogously from daily minimum temperature. Daily binary exceedance indicators were then summed within year to obtain annual counts, following the logic of percentile-based temperature-extreme diagnostics `(Frich et al., 2002; Alexander et al., 2006; Zhang et al., 2011)`. The resulting annual series are therefore interpretable as frequency-based measures of how often thermally unusual days or nights occurred in a given year at each station.

### 2.3. Quantile regression and summary metrics

Quantile regression was used as the principal inferential framework because the objective is to resolve distribution-dependent change rather than only mean change `(Koenker and Bassett, 1978; Koenker, 2005)`. For an annual index response `Y` and time covariate `T`, the `tau`-level conditional quantile is modeled as

$$
Q_Y(\tau \mid T) = \beta_0(\tau) + \beta_1(\tau) T, \qquad \tau \in (0,1).
$$

In contrast to OLS, which estimates `E[Y|T]` by minimizing squared residuals, quantile regression estimates `Q_Y(\tau|T)` by minimizing an asymmetric absolute-loss objective:

$$
\hat{\beta}(\tau) = \arg\min_{\beta} \sum_{i=1}^{n} \rho_{\tau}\!\left(y_i - x_i^\top\beta\right),
$$

where `x_i = (1, T_i)^\top`, and the check (pinball) loss is

$$
\rho_{\tau}(u) = u\,\big(\tau - \mathbf{1}\{u < 0\}\big)
=
\begin{cases}
\tau u, & u \ge 0,\\
(\tau - 1)u, & u < 0.
\end{cases}
$$

This formulation permits slope heterogeneity across the response distribution. Specifically, `\beta_1(\tau)` quantifies the decadal trend at quantile level `tau`, so differences between lower-tail, central, and upper-tail slopes reveal whether change is distribution-wide or concentrated in specific parts of the distribution.

Operationally, for each station-index series, linear trends were estimated across the conditional quantiles of annual counts, with time expressed in decades so that slopes are interpretable in days year^-1 per decade. Estimation was performed over the full quantile grid from `q = 0.05` to `q = 0.95` in increments of `0.01`, while focal interpretation emphasized `q0.05`, `q0.10`, `q0.50`, `q0.90`, and `q0.95`. OLS trends were estimated in parallel as a benchmark for mean change.

The quantile grid serves two different purposes in the manuscript. First, the dense grid is used to characterize the geometry of change across the full distribution and to construct regional quantile-coefficient panels. Second, the focal quantiles are used to build interpretable summary metrics, representative-station comparisons, maps, and significance diagnostics. This dual strategy preserves distributional detail without overloading the main text with an excessive number of parameters.

To summarize tail asymmetry, the following contrast was computed:

`Delta1 = slope(q0.95) - slope(q0.05)`

In addition, two secondary contrasts were used in feature engineering: `Delta2 = slope(q0.95) - slope(q0.50)` and `Delta3 = slope(q0.50) - slope(q0.05)`. Positive `Delta1` indicates stronger upper-tail intensification, whereas negative `Delta1` indicates relative upper-tail weakening or, for declining indices, stronger upper-tail contraction. These contrasts allow the study to distinguish broad-based distributional shifts from specifically upper-tail or lower-tail reorganization.

### 2.4. Uncertainty, spatial significance, and regionalization

Bootstrap resampling was used to characterize uncertainty in station-level quantile slopes. The main manuscript workflow employs 200 moving-block bootstrap replicates, because annual climate-index series may retain short-range temporal dependence and simple iid resampling would not preserve that structure adequately. Block length was not fixed a priori; instead, it was chosen data-adaptively from series length using a bounded cube-root rule so that the resampling window remained responsive to the effective record length while avoiding implausibly short or long blocks. For each focal quantile, bootstrap means, standard deviations, and confidence intervals were retained so that both point estimates and uncertainty structure could enter the later clustering and robustness stages. Maximum-entropy bootstrap (`meboot`) was retained as a sensitivity comparison because it is often useful for short, non-Gaussian series `(Efron and Tibshirani, 1993; Vinod, 2006)`, but it was treated as a secondary robustness check rather than as the main inferential baseline. In this study, the bootstrap is used primarily to compare uncertainty structure across indices, quantiles, and methodological variants rather than to claim ultra-precise tail interval estimation at individual stations; for that reason, the emphasis is placed on comparative robustness rather than on exact nominal interval coverage in the far tails. Because tail-depth adequacy remained a plausible reviewer concern, a targeted bootstrap-depth sensitivity analysis was also performed for the two central warm-extreme indices by re-estimating their station-level bootstrap summaries with 400 replicates instead of the default 200; this check was designed to test the stability of regional bootstrap means and tail-contrast metrics rather than to rerun the full manuscript workflow, and its results are summarized in Supplementary Figure S11 and Supplementary Tables S12-S13.

Field significance was assessed through false-discovery-rate control `(Benjamini and Hochberg, 1995)` across stations for the focal quantiles, and spatial dependence was examined using Moran's I `(Moran, 1950)` with permutation testing. Regionalization was performed using hierarchical clustering with average linkage and Euclidean distance on standardized feature vectors derived from station-level slope and uncertainty metrics. Hierarchical clustering was preferred because it preserves relative station-similarity structure, supports nested exploratory regionalization, and yields climatologically interpretable groupings without requiring a priori assumptions about the geometry of the clusters. The main clustering configuration used an uncertainty-aware feature set that included focal slopes, Delta contrasts, bootstrap means, and bootstrap variability. Robustness was assessed through a reduced-feature rerun, adjusted Rand index `(Hubert and Arabie, 1985)`, and a targeted alternative-method comparison in which the baseline clustering was compared with complete-linkage hierarchical clustering, Ward hierarchical clustering, average-linkage clustering with cityblock distance, and k-means. This additional check was designed to distinguish genuinely stable station-group structure from patterns that depend strongly on a single clustering convention.

### 2.5. Additional synthesis analyses

Three additional synthesis steps were used to strengthen interpretation. First, day-night asymmetry was quantified by directly comparing the regional behaviour of day and night indices. Second, split-period analysis was used to test temporal nonlinearity by comparing trends in 1961-1990 and 1991-2024. This split was treated as an interpretable synthesis diagnostic rather than as a formal breakpoint attribution exercise. Third, station behaviour was examined by elevation class to complement the geographic regression analysis and provide a more interpretable physiographic context for the observed patterns. These analyses were designed as explanatory layers rather than stand-alone statistical modules; their role is to help interpret the quantile-regression results in climatological rather than purely algorithmic terms.

Spatial visualization was treated in the same spirit. Station maps of `Delta1` were used as the main spatial summaries, while interpolated quantile maps were used only as visualization tools. The dedicated `paper2_figure3` map set emphasized a thin-plate-spline representation for selected quantiles, whereas the broader sensitivity framework compared alternative interpolation choices to test whether the major qualitative spatial messages were method-dependent.

### 2.6. Analytical workflow

Figure 2 summarizes the full analytical chain implemented in the pipeline, from station observations and percentile-threshold estimation to quantile regression, uncertainty estimation, regionalization, and publication-oriented synthesis products.

![Figure 2. Workflow of the analytical framework.](../outputs/figures/ijoc_workflow_flowchart.png)

*Figure 2. Workflow of the analytical framework used in the study. The sequence links daily station observations, percentile-based index construction, full-grid quantile regression, bootstrap uncertainty estimation, significance diagnostics, hierarchical regionalization, spatial visualization, and synthesis analyses used to generate the manuscript outputs.*

This workflow is important methodologically because the manuscript is not built around one isolated model fit. Instead, it integrates several mutually reinforcing steps: dense quantile estimation to resolve distributional structure, bootstrap resampling to assess uncertainty, regionalization to summarize station heterogeneity, map-based visualization to interpret spatial expression, and robustness checks to evaluate sensitivity to analytical choices. The final results should therefore be read as the outcome of an internally connected inference framework rather than as a collection of independent plots.

## 3. Results and discussion

### 3.1. Central claim and overview

The results support a clear central claim: **thermal-extreme change across Iran is not a simple mean shift, but a quantile-dependent reorganization of warm and cool event distributions**. Warm extremes increase and cool extremes decline, but the strongest changes often occur away from the centre of the distribution and vary sharply among indices, stations, and regional regimes. The main value of the analysis is therefore not only to confirm directional warming, but to show where in the distribution the reorganization is most intense and how that structure differs across space.

### 3.2. Regional synthesis of quantile-dependent change

At the regional scale, warm days and warm nights both increase, whereas cool days and cool nights decline. However, the distributional structure of change differs markedly among indices.

![Figure 3. Regional quantile-coefficient panels.](../outputs/figures/ijoc_regional_quantile_panels.png)

*Figure 3. Regional quantile-regression coefficients for (a) warm days, (b) warm nights, (c) cool days, and (d) cool nights. Black points show quantile-specific slopes, grey shading indicates 95% confidence intervals, and red lines show the OLS mean trend and its confidence bounds. The display range is restricted to q = 0.05-0.95.*

**Table 1.** Regional synthesis of mean trends, focal quantile trends, field significance, and dominant-cluster behaviour. The column `Approx. 1961-2024 change` is a simple slope-based translation over the full record length and is included only as an approximate measure of climatic magnitude.

| Index | OLS slope | q0.05 | q0.50 | q0.95 | Approx. 1961-2024 change | FDR sig. q0.10 | FDR sig. q0.50 | FDR sig. q0.90 | Dominant cluster | Dominant-cluster median Delta1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Warm days | +9.53 | +7.81 | +8.15 | +12.05 | +60.0 | 27 | 27 | 17 | 2 | +2.26 |
| Warm nights | +13.58 | +13.96 | +13.65 | +12.55 | +85.5 | 29 | 28 | 25 | 2 | +1.70 |
| Cool days | -2.56 | -1.51 | -2.01 | -8.05 | -16.1 | 0 | 17 | 26 | 1 | -6.64 |
| Cool nights | -2.98 | -1.77 | -2.47 | -6.78 | -18.8 | 12 | 23 | 25 | 2 | -8.21 |

Three points are especially important. First, warm nights show the strongest mean increase, corresponding to an approximate slope-based rise of about 85.5 days per year across the 1961-2024 record when the OLS estimate is translated over the full record length. Second, warm days display a clear upper-tail ordering (`q0.95 > q0.50 > q0.05`), indicating classical upper-tail amplification. Third, cool days and cool nights show the opposite structure: their strongest declines occur in upper quantiles, implying a rapid loss of years historically characterized by many cool events.

This structure moves the interpretation beyond significance testing alone. For example, the warm-day OLS slope of +9.53 days year^-1 per decade is not merely statistically non-zero; under a simple linear translation across the record it corresponds to roughly two additional months of warm-day conditions per year relative to the beginning of the period. The regional climate signal is therefore not only detectable, but climatologically large.

### 3.3. Day-night asymmetry

One of the most informative features of the dataset is the contrast between day and night indices.

**Table 2.** Regional day-night asymmetry in OLS and focal-quantile slopes, expressed as night-index trend minus day-index trend for the warm pair and cool-night minus cool-day trend for the cool pair.

| Comparison | OLS difference | q0.05 difference | q0.50 difference | q0.95 difference |
|---|---:|---:|---:|---:|
| Warm nights minus warm days | +4.05 | +6.14 | +5.50 | +0.50 |
| Cool nights minus cool days | -0.42 | -0.26 | -0.47 | +1.27 |

Warm nights increase more strongly than warm days at the mean, low, and median quantiles. This implies that nocturnal warming is not merely a tail phenomenon but a broad shift in the warm-night regime. Such behaviour is consistent with wider dryland and Mediterranean evidence showing enhanced nocturnal heat accumulation under warming `(Kostopoulou and Jones, 2005; Zittis et al., 2016)`.

By contrast, warm days show stronger upper-tail amplification than warm nights. The regional warm-day difference between q0.95 and q0.05 is +4.24 days year^-1 per decade, whereas warm nights show a slightly negative q0.95-q0.05 contrast (-1.41). Daytime extremes are therefore more strongly tilted toward upper-tail amplification, while nighttime warming is stronger overall but more distributionally uniform at the regional scale. This asymmetry is one of the clearest signs that day and night warming should not be treated as interchangeable manifestations of the same thermal trend.

### 3.4. Representative stations and cluster-defined quantile geometries

Regional means necessarily conceal strong station-scale contrasts. To connect the clustering results to physically interpretable slope geometries, representative stations were selected from the cluster structures of each index.

**Table 3.** Representative stations selected to illustrate contrasting cluster membership and focal-quantile behaviour. Stations were chosen to span the main cluster-defined geometries rather than to maximize a single metric.

| Index | Station | Cluster | q0.05 | q0.50 | q0.95 | Delta1 | Interpretation |
|---|---|---:|---:|---:|---:|---:|---|
| Warm days | Esfahan | 1 | +11.92 | +16.29 | +28.50 | +16.58 | Extreme upper-tail amplification |
| Warm days | Bushehr (Airport) | 2 | +3.83 | +4.00 | -2.50 | -6.33 | Upper-tail weakening |
| Warm days | Khoy | 3 | +11.72 | +12.43 | +8.83 | -2.89 | High baseline with downward-tilting upper tail |
| Warm days | Zanjan | 4 | +2.50 | +9.44 | +16.61 | +14.11 | Strong upward quantile transition |
| Warm nights | Ramsar | 1 | +25.20 | +25.59 | +19.82 | -5.38 | High baseline, muted upper tail |
| Warm nights | Tabriz | 2 | +9.73 | +11.96 | +15.23 | +5.50 | Moderate positive amplification |
| Warm nights | Tehran (Mehrabad Airport) | 3 | +9.11 | +15.88 | +21.03 | +11.91 | Strong upper-tail intensification |
| Warm nights | Shiraz | 4 | +9.77 | +14.00 | +35.68 | +25.90 | Exceptional upper-tail amplification |
| Cool days | Sanandaj | 1 | -1.46 | -4.44 | -8.11 | -6.64 | Typical upper-tail contraction |
| Cool days | Zabol | 2 | +1.40 | -2.86 | -8.75 | -10.15 | Positive lower tail, negative upper tail |
| Cool days | Orumiyeh | 3 | +0.17 | -3.21 | -16.77 | -16.95 | Extreme high-quantile collapse |
| Cool days | Torbat-E Heydariyeh | 4 | -0.26 | +3.26 | -12.35 | -12.10 | Positive median, collapsing upper tail |
| Cool nights | Khorramabad | 1 | +1.48 | +2.70 | -19.52 | -21.01 | Strongest upper-tail contraction |
| Cool nights | Qazvin | 2 | +0.50 | -3.56 | -13.26 | -13.76 | Classical contraction profile |
| Cool nights | Gorgan | 3 | +4.00 | +6.43 | +1.52 | -2.48 | Positive low- and mid-quantile regime |
| Cool nights | Shahrekord | 4 | +5.60 | +14.62 | +17.56 | +11.96 | Positive-tail outlier |

![Figure 4. Representative-station comparisons.](../outputs/figures/ijoc_station_comparisons/main_figure_representative_stations.png)

*Figure 4. Multi-panel comparison of representative stations selected to contrast cluster membership and focal-quantile behaviour for (a) warm days, (b) warm nights, (c) cool days, and (d) cool nights. The selected stations show that the clusters correspond to distinct quantile-slope geometries rather than only differences in average trend magnitude.*

These representative stations clarify the logic of the cluster structure. Warm-day change ranges from extreme upper-tail amplification (Esfahan, Zanjan) to relative upper-tail weakening (Bushehr Airport, Khoy). Warm nights span high-baseline but downward-tilting profiles (Ramsar), strong monotonic amplification (Tehran), and exceptional upper-tail growth (Shiraz). The cool indices likewise span dominant contraction profiles, mixed low- versus high-quantile behaviour, and isolated positive outliers such as Shahrekord for cool nights. The clusters are therefore climatologically interpretable rather than merely statistically convenient, because they correspond to distinct geometries of distributional change.

### 3.5. Spatial expression and regionalization

The spatial pattern of `Delta1` is coherent enough to support regional interpretation, but not smooth enough to justify treating the country as a single spatial field.

![Figure 5. Main Delta1 map panels.](../outputs/figures/ijoc_main_delta1_maps.png)

*Figure 5. Multi-panel station maps of `Delta1 = slope(q0.95) - slope(q0.05)` for (a) warm days, (b) warm nights, (c) cool days, and (d) cool nights. Warm indices generally show positive Delta1 values, whereas cool indices show widespread negative Delta1 values, consistent with upper-tail contraction.*

Warm days show a broad tendency toward positive `Delta1`, but with strong contrast between central, western, and southern stations. Warm nights also show positive `Delta1` at several stations, though less uniformly, and with one exceptional outlier in Shiraz. Cool days exhibit widespread negative `Delta1`, indicating upper-tail collapse across most of the network. Cool nights are similarly negative at most stations, but with notable positive anomalies at Shahrekord and, more modestly, Esfahan.

This pattern aligns with the formal field-significance results: intermediate quantiles are the most consistently significant across stations, yet Moran's I is not significant for any of the 20 slope fields tested. Accordingly, the spatial pattern is best interpreted as a set of regionally meaningful but non-smooth station responses, consistent with topographic and climatic heterogeneity rather than a single national-scale spatial mode.

### 3.5.1. Additional structural and spatial diagnostics

The main-text regionalization results are further supported by two additional diagnostic layers presented in the Supplementary Material. First, quantile-specific dendrograms at `q0.05`, `q0.50`, and `q0.95` (Supplementary Figures S2-S4) show that station similarity is not fixed across the distribution. Lower-tail structures are comparatively compressed, whereas upper-tail structures become more strongly separated, especially for warm nights and, in a different direction, for cool nights. This supports the interpretation that the clustering is not merely a multivariate artefact; rather, meaningful station-group separation is already visible in one-dimensional quantile-specific slope space.

Second, thin-plate-spline quantile maps for `q0.05`, `q0.10`, `q0.50`, `q0.90`, and `q0.95` (Supplementary Figures S5-S9) clarify how the geography of change evolves from lower to upper quantiles. Lower quantiles show smoother and more regionally coherent behaviour, whereas upper quantiles reveal sharper hotspots for warm indices and more focused contraction zones for cool indices. Because Moran's I does not support a strong national-scale spatial field, these interpolated panels are interpreted as visualization aids rather than as primary inferential products. Even so, they provide an important consistency check: the main station-based conclusions remain visible when the slope fields are viewed as continuous spatial summaries.

### 3.6. Stronger explanatory context from geography and elevation

The driver analysis already indicated that warm-night slopes are negatively associated with elevation, with standardized betas of -0.507 at q0.05 and -0.479 at q0.50. An elevation-class summary reinforces that result. Mean warm-night q0.50 slopes are about +18.43 days year^-1 per decade in lowland stations, +14.25 in mid-elevation stations, and only +6.58 in highland stations. Lower-elevation environments therefore appear especially prone to broad-based nocturnal warming.

The elevation-class contrasts also refine the interpretation of day versus night behaviour. Warm-day upper-tail amplification is strongest in the highland class (`Delta1 ≈ +7.48`), whereas warm-night amplification is strongest in the mid-elevation class (`Delta1 ≈ +6.25`) but weak in both lowland and highland classes. For cool nights, the highland class shows a much weaker contraction (`Delta1 ≈ -3.44`) than the lowland and mid-elevation classes (both about `-8.5` to `-8.7`). These contrasts suggest that the geography of day-night asymmetry cannot be reduced to a single monotonic elevation effect. Different parts of the thermal-extreme distribution appear to be modulated by different local controls.

From a physical standpoint, the warm-night patterns are consistent with greater nocturnal heat retention at lower elevations and in settings where radiative cooling may be reduced, but the present analysis does not directly identify the responsible mechanisms. The warm-day results, by contrast, suggest that some elevated and interior stations remain especially sensitive to upper-tail daytime amplification, possibly because exceptionally hot years intensify more rapidly than the broader mean regime. The cool-day and cool-night indices show the opposite direction of change, but again with stronger contraction in upper quantiles, indicating that the most historically cool years are disappearing fastest. These interpretations should therefore be read as physically plausible explanations consistent with the statistical patterns, not as formal process attribution.

### 3.7. Robustness, methodological sensitivity, and temporal nonlinearity

One of the strongest assets of the study is that the principal conclusions do not depend on a single analytical choice.

![Figure 6. Robustness synthesis.](../outputs/figures/ijoc_robustness_synthesis.png)

*Figure 6. Synthesis of robustness diagnostics. Panel (a) summarizes reference-period sensitivity, panel (b) summarizes bootstrap-method sensitivity with moving-block resampling as the baseline and `meboot` as the comparison, panel (c) shows interpolation-method agreement for the configured spatial-visualization workflow, and panel (d) reports clustering robustness as adjusted Rand index. The principal structural results remain stable despite some metric-level sensitivity.*

The robustness synthesis shows that warm-night metrics are most sensitive to reference-period choice, whereas cool-day and cool-night bootstrap summaries are more sensitive to the choice between moving-block and `meboot` resampling, especially in upper-tail means and confidence bounds. Interpolation agreement is strongest for cool nights and warm nights under the `linear_rbf` versus `linear` comparison, and clustering robustness is particularly high for warm nights (ARI = 0.93) and cool nights (ARI = 0.82). These diagnostics support a confident qualitative interpretation while still motivating caution for the most extreme tails.

An additional sensitivity test excluded all stations flagged by any detrended homogeneity diagnostic (Supplementary Figure S10; Supplementary Table S11). The principal directional conclusions were retained: warm indices remained positive and cool indices remained negative at the regional scale, and warm-day upper-tail amplification remained evident. However, some magnitude changes were non-negligible, particularly for warm nights, where the regional tail-contrast metric became more positive after exclusion. This result is useful rather than problematic: it indicates that the manuscript's broad climatic conclusions are not driven by a small subset of flagged stations, while also showing that exact tail-asymmetry magnitudes should be interpreted with appropriate caution.

A second targeted sensitivity check addressed bootstrap depth directly by increasing the number of resamples from 200 to 400 for `warm_days` and `warm_nights` (Supplementary Figure S11; Supplementary Tables S12-S13). The result is reassuringly stable. Regional bootstrap means for `q0.05`, `q0.50`, `q0.95`, and `Delta1` changed by no more than about `0.06 days year^-1 per decade` in absolute mean value, and station-level correspondence between the 200- and 400-replicate summaries remained extremely high (`r = 0.998-1.000`). Confidence-interval widths changed somewhat more than the means, as expected, but their station-level correlations still remained high (`r = 0.975-0.992`). This additional check does not remove all tail uncertainty, but it does show that the manuscript's main warm-extreme messages are not an artefact of a shallow bootstrap configuration.

A third sensitivity check examined whether the regionalization itself depends strongly on the chosen clustering method (Supplementary Figure S12; Supplementary Tables S14-S15). The answer is mixed, and therefore informative. `warm_days` remained comparatively stable across alternative methods (`ARI = 0.63-0.92`), and `warm_nights` also retained a strong common structure for most alternatives (`ARI = 0.88-0.93`), although complete linkage produced a noticeably less similar partition (`ARI = 0.43`). `cool_nights` showed moderate method dependence, with one alternative reproducing the baseline almost exactly (`ARI = 1.00`) but others yielding only partial agreement (`ARI = 0.34-0.48`). `cool_days` were the least stable (`ARI = 0.05-0.35`), indicating that their regional grouping should be read as exploratory rather than as a uniquely resolved partition. This pattern does not invalidate the station-scale signal; rather, it shows that some regional structures, especially for cooling-related indices, are weaker and less crisply separated than those for the warm indices.

Temporal nonlinearity provides an additional layer of evidence.

![Figure 7. Split-period comparison.](../outputs/figures/ijoc_split_period_comparison.png)

*Figure 7. Comparison of regional trend estimates between 1961-1990 and 1991-2024 for the focal quantiles and OLS trend. The later period shows marked intensification of warm extremes and stronger collapse of cool extremes, especially in upper quantiles.*

The split-period analysis shows that all four indices behave differently after 1990. Warm-day OLS slopes increase from about +0.63 days year^-1 per decade in 1961-1990 to +19.47 in 1991-2024, while warm-night OLS slopes rise from +4.44 to +19.49. Cool-day OLS slopes shift from +1.78 in the early period to -7.86 in the later period, and cool nights shift from -0.44 to -3.33. The strongest acceleration appears in the warm-day upper tail and in the cool-day upper-tail decline. The distributional reorganization detected in the full-period analysis is therefore not an artefact of long averaging; it reflects a marked strengthening of recent change.

### 3.8. Comparison with the wider climatological literature

The present findings are broadly consistent with the wider literature in three respects. First, the strong rise in warm extremes agrees with previous assessments over Iran and neighbouring dryland regions `(Soltani et al., 2016; Zhang et al., 2005a)`. Second, the stronger and more widespread increase in warm nights is consistent with Mediterranean and Middle Eastern evidence suggesting strong nocturnal warming and increasing heat stress `(Kostopoulou and Jones, 2005; Zittis et al., 2016)`. Third, the collapse of cool extremes, especially in upper quantiles, is consistent with the broader shift away from historically cool thermal regimes documented in global and regional extreme-index studies `(Alexander et al., 2006; Donat et al., 2013)`.

The distinctive contribution of the present study is not simply that it confirms warming over Iran. Rather, it shows that the most informative part of the signal lies in the **shape** of change: warm indices do not merely rise, and cool indices do not merely fall; both are internally reorganized across quantiles. In this respect, the study extends the interpretation beyond standard trend reporting and contributes to the international literature on how extremes evolve in arid and topographically complex settings.

### 3.9. Limitations

The limitations of the study should be stated explicitly. The station network is finite, the strongest tails remain the most uncertain part of the distribution, and the driver analysis is intentionally simple. Interpolated surfaces are used as visualization tools rather than primary inferential products. The period split is informative but still descriptive and does not amount to formal breakpoint attribution. These limitations do not invalidate the main conclusions, but they do define the appropriate scope of interpretation.

## 4. Conclusions

This study shows that thermal-extreme change across Iran is best understood not as a uniform shift in the mean, but as a **quantile-dependent reorganization** of warm and cool event distributions. Warm days and warm nights both increase, but with different internal structures: warm nights are stronger overall, whereas warm days show clearer upper-tail amplification. Cool days and cool nights both decline, with the strongest contractions concentrated in upper quantiles, indicating rapid loss of years historically characterized by many cool events.

Three conclusions are especially important. First, the mean signal alone understates how strongly the distribution has been reorganized. Second, station-scale heterogeneity is climatologically meaningful and can be summarized effectively through cluster-defined quantile geometries. Third, the main conclusions are robust to alternative methodological choices and are strengthened, rather than weakened, by split-period evidence showing stronger recent change.

For climatology more broadly, the results suggest that observational studies of extremes in arid and topographically complex regions benefit substantially from quantile-aware analysis. In such settings, the most scientifically and societally relevant information may lie less in whether a mean trend is positive or negative than in how different parts of the distribution are changing relative to one another. The wider implication is that regional climate-change assessments are likely to be more informative when they treat tail behaviour as a primary component of inference rather than as a secondary descriptive detail.

## Acknowledgements

Funding, acknowledgement, and data-availability statements should be inserted here in the journal's required format.

## Supplementary material

The Supplementary Material should contain at minimum the following items cited in the manuscript:

1. `Supplementary Figure S1`: `outputs/figures/ijoc_data_quality_homogeneity.png`
2. `Supplementary Figures S2-S4`: `outputs/figures/paper1_quantile_dendrograms/figure_5_tau_0.05.png`, `outputs/figures/paper1_quantile_dendrograms/figure_6_tau_0.50.png`, `outputs/figures/paper1_quantile_dendrograms/figure_7_tau_0.95.png`
3. `Supplementary Figures S5-S9`: `outputs/figures/paper2_figure3_quantile_maps/figure3_tau_0.05.png`, `outputs/figures/paper2_figure3_quantile_maps/figure3_tau_0.10.png`, `outputs/figures/paper2_figure3_quantile_maps/figure3_tau_0.50.png`, `outputs/figures/paper2_figure3_quantile_maps/figure3_tau_0.90.png`, `outputs/figures/paper2_figure3_quantile_maps/figure3_tau_0.95.png`
4. `Supplementary Figure S10`: `outputs/figures/ijoc_homogeneity_sensitivity.png`
5. `Supplementary Tables S1-S3`: `outputs/tables/data_quality_homogeneity_overview.csv`, `outputs/tables/data_quality_station_summary.csv`, `outputs/tables/data_homogeneity_tests_station_summary.csv`
6. `Supplementary Table S11`: `outputs/tables/homogeneity_flag_exclusion_sensitivity.csv`
7. `Supplementary Figure S11`: `outputs/figures/ijoc_bootstrap_depth_sensitivity.png`
8. `Supplementary Tables S12-S13`: `outputs/tables/bootstrap_depth_sensitivity_summary.csv`, `outputs/tables/bootstrap_depth_sensitivity_station_comparison.csv`
9. `Supplementary Figure S12`: `outputs/figures/ijoc_alternative_clustering_sensitivity.png`
10. `Supplementary Tables S14-S15`: `outputs/tables/alternative_clustering_sensitivity_summary.csv`, `outputs/tables/alternative_clustering_assignments.csv`

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

Zhang X, Aguilar E, Sensoy S, Melkonyan H, Tagiyeva U, Ahmed N, Kutaladze N, Rahimzadeh F, Taghipour A, Hantosh TH, Albert P, Semawi M, Karam Ali M, Al-Shabibi MHS, Al-Oulan Z, Zatari T, Al Dean Khelet I, Hamoud S, Sagir R, Demircan M, Eken M, Adiguzel M, Alexander LV, Peterson TC, Wallis T (2005a) Trends in Middle East climate extreme indices from 1950 to 2003. *Journal of Geophysical Research: Atmospheres* 110:D22104. https://doi.org/10.1029/2005JD006181

Zhang X, Hegerl G, Zwiers FW, Kenyon J (2005b) Avoiding inhomogeneity in percentile-based indices of temperature extremes. *Journal of Climate* 18(11):1641-1651. https://doi.org/10.1175/JCLI3366.1

Zhang X, Alexander L, Hegerl GC, Jones P, Klein Tank AMG, Peterson TC, Trewin B, Zwiers FW (2011) Indices for monitoring changes in extremes based on daily temperature and precipitation data. *WIREs Climate Change* 2(6):851-870. https://doi.org/10.1002/wcc.147

Zittis G, Hadjinicolaou P, Lelieveld J (2016) Strongly increasing heat extremes in the Middle East and North Africa (MENA) in the 21st century. *Climatic Change* 137:245-260. https://doi.org/10.1007/s10584-016-1665-6



