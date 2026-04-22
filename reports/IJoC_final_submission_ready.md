# An integrated quantile-based climatological synthesis of thermal extremes across Iran during 1961-2024

**Running title:** Quantile-based thermal extremes across Iran

## Abstract

Mean trends alone are often insufficient to describe changes in climate extremes because the most consequential impacts emerge through shifts in the tails of the distribution. This study investigates quantile-dependent trends in four annual thermal-extreme indices across Iran during 1961-2024 using observations from 30 meteorological stations. Warm days, warm nights, cool days, and cool nights were analysed with quantile regression, while ordinary least squares (OLS) provided a benchmark for mean change. The workflow combined data-quality and homogeneity diagnostics, bootstrap uncertainty estimation, false-discovery-rate field-significance testing, cluster-based regionalization, spatial visualization, and pipeline-integrated sensitivity analyses covering homogeneity-based station exclusion, bootstrap depth, alternative clustering specifications, and advanced publication-oriented robustness checks for reference period, bootstrap method, interpolation, and regionalization. The results show that thermal-extreme change across Iran is not well represented as a uniform mean shift. Warm days and warm nights increased across most of the network, but with substantial quantile dependence and marked station-to-station heterogeneity. Warm days exhibited clear upper-tail amplification, whereas warm nights showed stronger overall increases but a more uniform regional structure. Cool days and cool nights declined, with the strongest negative changes concentrated in upper quantiles, implying rapid loss of years historically characterized by many cool events. Regional OLS trends were +9.53 and +13.58 days year^-1 per decade for warm days and warm nights, and -2.56 and -2.98 days year^-1 per decade for cool days and cool nights, respectively. Field significance was strongest at intermediate quantiles, whereas Moran's I showed weak national-scale spatial autocorrelation, indicating that the observed changes are regionally interpretable but not reducible to a single smooth spatial field. Split-period analysis further showed substantial intensification of warm extremes and stronger collapse of cool extremes after 1990. Overall, the study provides an integrated quantile-based climatological synthesis with robust evaluation, showing how warm and cool event distributions have changed across Iran and offering a reproducible regional template for arid and topographically complex settings.

**Keywords:** thermal extremes; quantile regression; climate change; warm nights; cool nights; regionalization; Iran

## 1. Introduction

Changes in climate extremes are among the most consequential manifestations of contemporary climate change because their impacts propagate through agriculture, public health, water resources, ecosystems, and infrastructure `(Easterling et al., 2000; IPCC, 2021)`. Heatwaves and warm extremes have also intensified in many regions in ways that are not fully described by mean warming alone `(Perkins, 2015; Perkins-Kirkpatrick and Lewis, 2020)`. Mean temperature trends provide a useful first-order summary of climate change, but they are often insufficient for understanding extremes, since the most climatically relevant impacts depend not only on shifts in the centre of the distribution but also on changes in its tails `(Katz and Brown, 1992; Katz et al., 2002)`. A region may therefore exhibit moderate mean warming while simultaneously undergoing strong upper-tail intensification or rapid erosion of historically cool conditions. Distribution-aware approaches are consequently becoming increasingly important in climatology `(Koenker and Bassett, 1978; Koenker, 2005; Barbosa et al., 2011; Fan, 2014; Reich, 2012)`.

This issue is especially important for thermal-extreme indices. Warm and cool day-night indices have been widely used to characterize climate variability and change `(Frich et al., 2002; Alexander et al., 2006; Donat et al., 2013; Zhang et al., 2011; Dunn et al., 2020)`, yet their interpretation often remains anchored in mean or median trends. Such summaries can obscure upper-tail amplification, lower-tail stabilization, and asymmetric decline across the distribution `(Moberg and Jones, 2005; Kostopoulou and Jones, 2005; Fischer and Schär, 2010; Russo and Sterl, 2011)`. In climatically heterogeneous regions, these effects may also vary sharply across space because of topographic contrasts, continentality, land-atmosphere coupling, and differences in moisture availability.

Iran offers a particularly useful setting for testing these issues. The country includes strong gradients in elevation, aridity, latitude, and marine influence, spanning humid Caspian lowlands, interior continental environments, hot southern lowlands, and high-elevation terrain. This diversity makes Iran a natural laboratory for examining whether observed thermal-extreme change behaves as a coherent national signal or as a set of locally modified regional responses. Previous work has documented changes in climate extremes across Iran and across the wider dryland and Mediterranean domain `(Kostopoulou and Jones, 2005; Zhang et al., 2005a; Soltani et al., 2016; Zittis et al., 2016; Vaghefi et al., 2019; Ahmadalipour and Moradkhani, 2018; Francis and Fonseca, 2024)`, while broader regional assessments indicate intensifying heat stress, increasingly severe heatwaves, and growing humid-heat risk around the Persian Gulf and the wider MENA region `(Pal and Eltahir, 2016; Raymond et al., 2020; Zittis et al., 2021; Raymond et al., 2024)`. Even so, relatively few station-based studies for Iran have resolved change across the full distribution while simultaneously addressing uncertainty, field significance, and regionalization.

The broader relevance of this question extends beyond Iran. A regional study is most useful when it contributes to wider climatological understanding rather than remaining purely descriptive. In that spirit, the present study does not ask only whether thermal extremes in Iran have changed. It asks whether those changes are better characterized through an integrated quantile-based synthesis rather than through a uniform mean-shift summary alone, and whether the resulting distribution-aware interpretation is spatially coherent, physically interpretable, and robust to methodological choices.

Using 30 meteorological stations for 1961-2024, this study analyses four annual indices: warm days, warm nights, cool days, and cool nights. The specific contributions are fivefold. First, the study resolves trends across the full conditional distribution rather than relying on mean change alone. Second, it quantifies upper-tail amplification and contraction using tail-contrast metrics. Third, it integrates uncertainty through bootstrap-based diagnostics. Fourth, it uses clustering to regionalize quantile signatures rather than only mapping individual slopes. Fifth, it evaluates the robustness of the principal conclusions to alternative analytical choices. Collectively, these steps position the manuscript as an integrated quantile-based climatological synthesis and robust regional evaluation of thermal-extreme change across Iran, rather than as the introduction of a fundamentally new quantile concept.

## 2. Data and methods

### 2.1. Data and study region

The analysis is based on daily observations from 30 meteorological stations distributed across Iran. Station metadata include latitude, longitude, and elevation, thereby allowing the thermal-extreme signals to be interpreted in a physiographic context. After preprocessing, the common analysis period spans 1961-2024, yielding 64 annual values for each station-index series.

![Figure 1. Study area and station network.](../outputs/figures/ijoc_study_area.png)

*Figure 1. Distribution of the 30 meteorological stations used in the analysis across Iran. Stations are coloured by elevation to highlight the physiographic diversity of the network and the coexistence of lowland, mid-elevation, and highland environments.*

The station-based design is deliberate. National and regional summaries are derived from station estimates rather than gridded products, thereby preserving local climatic heterogeneity and avoiding interpolation assumptions at the primary inference stage, even though gridded extreme-index archives such as HadEX3 remain invaluable for broader-scale benchmarking `(Dunn et al., 2020)`. This choice is especially important in Iran, where elevation contrasts, continentality, and varying marine influence can generate strong short-range spatial differences in thermal behaviour. The input variables were daily minimum and maximum temperature, organized by station and calendar date, with leap-day observations removed before percentile-threshold estimation in order to preserve a consistent 365-day annual cycle.

Before index construction, the station records were screened for completeness, duplicate dates, and simple physical consistency, including `tmin <= tmax` and agreement of archived `tmean` with the daily temperature range where all three variables were available. Completeness was high across the network (median values about 98.7% for `tmin`, 99.2% for `tmax`, and 98.5% for `tmean` or a `tmin`-`tmax` surrogate), and no station contained duplicate calendar dates. Potential inhomogeneity was then screened using Pettitt-, SNHT-, and Buishand-type breakpoint diagnostics applied to annual mean temperature series, with linearly detrended versions used as the main screening tests to reduce confounding by long-term climatic trends. The raw annual-mean tests flagged all 30 stations, whereas the detrended screening still flagged 19 of 30 stations, indicating that breakpoint-sensitive diagnostics are highly responsive in a strongly warming regional climate and should not be interpreted mechanically as a mandate for station-wise correction.

No automatic homogenization adjustment was applied to the daily temperature series. That decision was deliberate rather than incidental. The available archive does not include the station metadata history or a neighboring reference-series framework needed for a defensible reference-based homogenization, and the breakpoint diagnostics were applied to annual mean temperature rather than directly to the daily `tmin`/`tmax` series used for percentile-based extremes. Under those conditions, a simple SNHT-style shift correction would require translating annual-mean break estimates into daily tail adjustments, which could itself distort percentile thresholds and extreme-count behaviour `(Zhang et al., 2005b)`. The homogeneity diagnostics were therefore used as screening tools to identify potentially sensitive records, while the main manuscript retained the original daily observations and treated flagged-station exclusion as a conservative bound on the conclusions rather than as a substitute for full homogenization. A summary of completeness and homogeneity screening is provided in Supplementary Figure S1 and Supplementary Tables S1-S2, and the exclusion sensitivity is summarized in Supplementary Figure S10 and Supplementary Table S11.
### 2.2. Thermal-extreme indices

Four annual indices were constructed from daily minimum (`tmin`) and maximum (`tmax`) temperature observations: warm days, warm nights, cool days, and cool nights. Daily thresholds were estimated separately for each station using the 1961-1990 reference period, a 5-day moving calendar window, and the 10th and 90th percentiles of station-specific daily temperature distributions. A minimum reference-sample requirement was enforced at the day-of-year level to reduce instability in threshold estimation near data-sparse dates.

Warm days and cool days were defined from daily maximum temperature relative to the station-specific 90th and 10th percentile thresholds, respectively. Warm nights and cool nights were defined analogously from daily minimum temperature. Daily binary exceedance indicators were then summed within year to obtain annual counts, following the logic of percentile-based temperature-extreme diagnostics `(Frich et al., 2002; Alexander et al., 2006; Zhang et al., 2011)`. The resulting annual series are therefore interpretable as frequency-based measures of how often thermally unusual days or nights occurred in a given year at each station.

### 2.3. Quantile regression and summary metrics

Quantile regression was used as the principal inferential framework because the objective is to resolve distribution-dependent change rather than only mean change `(Koenker and Bassett, 1978; Koenker, 2005)`. This choice is also consistent with earlier climate applications showing that quantile-specific temperature trends can reveal structures that are hidden in ordinary least-squares summaries `(Barbosa et al., 2011; Fan, 2014)`. More elaborate spatiotemporal quantile models have also been proposed for environmental-process change detection `(Reich, 2012)`, but the present study intentionally retains a transparent station-wise linear specification so that distributional slope contrasts remain directly interpretable for each index and station. For an annual index response `Y` and time covariate `T`, the `tau`-level conditional quantile is modeled as

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

The quantile grid serves two different purposes in the manuscript. First, the dense grid is used to characterize the geometry of change across the full distribution and to construct regional quantile-coefficient panels. Those panels are treated as descriptive pooled summaries because they are derived from annual station means and therefore can smooth over spatial heterogeneity. Second, the focal quantiles are used to build interpretable station-level summaries, representative-station comparisons, maps, and significance diagnostics. This dual strategy preserves distributional detail without overloading the main text with an excessive number of parameters while keeping primary inference anchored in the station-specific fits.

To summarize tail asymmetry, the following contrast was computed:

`Delta1 = slope(q0.95) - slope(q0.05)`

In addition, two secondary contrasts were used in feature engineering: `Delta2 = slope(q0.95) - slope(q0.50)` and `Delta3 = slope(q0.50) - slope(q0.05)`. Positive `Delta1` indicates stronger upper-tail intensification, whereas negative `Delta1` indicates relative upper-tail weakening or, for declining indices, stronger upper-tail contraction. These contrasts allow the study to distinguish broad-based distributional shifts from specifically upper-tail or lower-tail differentiation within the observed trend structure.

### 2.4. Uncertainty, spatial significance, and regionalization

Bootstrap resampling was used to characterize uncertainty in station-level quantile slopes. The main manuscript workflow employs 200 moving-block bootstrap replicates, because annual climate-index series may retain short-range temporal dependence and simple iid resampling would not preserve that structure adequately. Block length was not fixed a priori; instead, it was chosen data-adaptively from series length using a bounded cube-root rule so that the resampling window remained responsive to the effective record length while avoiding implausibly short or long blocks. For each focal quantile, bootstrap means, standard deviations, and confidence intervals were retained so that both point estimates and uncertainty structure could enter the later clustering and robustness stages. Inference for the focal tail quantiles was based primarily on these bootstrap confidence intervals, especially for `q0.05` and `q0.95`, where asymptotic analytic intervals can be unstable for moderately short annual series. Analytic quantile-regression confidence intervals were retained as a secondary comparison and robustness diagnostic rather than as the main inferential criterion; this distinction was tracked explicitly by comparing analytic and bootstrap significance flags in the station-summary outputs. Maximum-entropy bootstrap (`meboot`) was retained as a sensitivity comparison because it is often useful for short, non-Gaussian series `(Efron and Tibshirani, 1993; Vinod, 2006)`, but it was treated as a secondary robustness check rather than as the main inferential baseline. In this study, the bootstrap is used primarily to compare uncertainty structure across indices, quantiles, and methodological variants rather than to claim ultra-precise tail interval estimation at individual stations; for that reason, the emphasis is placed on comparative robustness rather than on exact nominal interval coverage in the far tails. Because tail-depth adequacy remained a plausible reviewer concern, a targeted bootstrap-depth sensitivity analysis was also performed for the two central warm-extreme indices by re-estimating their station-level bootstrap summaries with 400 replicates instead of the default 200; this check was designed to test the stability of regional bootstrap means and tail-contrast metrics rather than to rerun the full manuscript workflow, and its results are summarized in Supplementary Figure S11 and Supplementary Tables S12-S13.

Field significance was assessed through false-discovery-rate control `(Benjamini and Hochberg, 1995)` across stations for the focal quantiles, and spatial dependence was examined using Moran's I `(Moran, 1950)` with permutation testing. Because bootstrap inference was prioritized only for the focal tail summaries, whereas field-significance screening was required across the broader focal-quantile set, the FDR workflow remained based on the analytic quantile-regression intervals and their derived p-values; those results are therefore interpreted as a complementary screening layer rather than as the primary tail-inference device. Regionalization was performed using hierarchical clustering with average linkage and Euclidean distance, but the clustering stage was treated explicitly as exploratory rather than fully inferential. Before clustering, candidate features were screened for near-collinearity using an absolute pairwise-correlation threshold so that almost-duplicate summaries would not dominate the partition simply because they encode the same slope geometry in multiple ways. The main clustering configuration therefore used a parsimonious slope-based feature set consisting of `q0.05`, `q0.50`, and `q0.95` station slopes, while richer uncertainty-aware feature sets were retained only for sensitivity reruns. Hierarchical clustering was preferred because it preserves relative station-similarity structure, supports nested exploratory regionalization, and yields climatologically interpretable groupings without requiring a priori assumptions about the geometry of the clusters. Robustness was assessed through an expanded uncertainty-aware rerun, adjusted Rand index `(Hubert and Arabie, 1985)`, and a targeted alternative-method comparison in which the baseline clustering was compared with complete-linkage hierarchical clustering, Ward hierarchical clustering, average-linkage clustering with cityblock distance, and k-means. This additional check was designed to distinguish genuinely stable station-group structure from patterns that depend strongly on a single clustering convention or on one particular feature encoding.

### 2.5. Additional synthesis analyses

Three additional synthesis steps were used to strengthen interpretation. First, day-night asymmetry was quantified by directly comparing the descriptive cross-station summaries of day and night indices rather than by replacing the underlying station-level evidence. Second, split-period analysis was used to test temporal nonlinearity by comparing trends in 1961-1990 and 1991-2024. This split was treated as an interpretable synthesis diagnostic rather than as a formal breakpoint attribution exercise. Third, station behaviour was examined by elevation class to complement the geographic regression analysis and provide a more interpretable physiographic context for the observed patterns. These analyses were designed as explanatory layers rather than stand-alone statistical modules; their role is to help interpret the station-level quantile-regression results in climatological rather than purely algorithmic terms.

Spatial visualization was treated in the same spirit. Station maps of `Delta1` were used as the main spatial summaries, while interpolated quantile maps were used only as visualization tools. The dedicated `paper2_figure3` map set emphasized a thin-plate-spline representation for selected quantiles, whereas the broader sensitivity framework compared alternative interpolation choices to test whether the major qualitative spatial messages were method-dependent. Importantly, these interpolated surfaces were not treated as smooth estimators of a latent national field; they were used only to provide visual orientation between sparse station locations, with scientific interpretation anchored primarily in the station points themselves.

### 2.6. Analytical workflow

Figure 2 summarizes the full analytical chain implemented in the pipeline, from station observations and percentile-threshold estimation to quantile regression, uncertainty estimation, regionalization, and publication-oriented synthesis products.

![Figure 2. Workflow of the analytical framework.](../outputs/figures/ijoc_workflow_flowchart.png)

*Figure 2. Workflow of the analytical framework used in the study. The sequence links daily station observations, percentile-based index construction, full-grid quantile regression, bootstrap uncertainty estimation, significance diagnostics, hierarchical regionalization, spatial visualization, and synthesis analyses used to generate the manuscript outputs.*

This workflow is important methodologically because the manuscript is not built around one isolated model fit. Instead, it integrates several mutually reinforcing steps: dense quantile estimation to resolve distributional structure, bootstrap resampling to assess uncertainty, regionalization to summarize station heterogeneity, map-based visualization to interpret spatial expression, and robustness checks to evaluate sensitivity to analytical choices. The final results should therefore be read as the outcome of an internally connected quantile-based synthesis and evaluation framework rather than as a collection of independent plots.

Just as importantly, these robustness layers were implemented as operational parts of the pipeline rather than as informal afterthoughts. The executable workflow runs data-quality and homogeneity screening before index construction, then carries forward dedicated homogeneity-exclusion sensitivity, targeted bootstrap-depth sensitivity, alternative-clustering sensitivity, and advanced publication analyses for spatial inference, method sensitivity, driver analysis, and regionalization. The corresponding outputs are written automatically to the supplementary figure and table paths listed below so that the manuscript's robustness claims are backed by reproducible artifacts rather than by narrative reassurance alone.

## 3. Results and discussion

### 3.1. Central claim and overview

The results support a clear central claim: **thermal-extreme change across Iran is not adequately described by a simple mean shift alone, and a quantile-based synthesis reveals important distributional asymmetries in warm and cool event trends**. Warm extremes increase and cool extremes decline, but the strongest changes often occur away from the centre of the distribution and vary sharply among indices, stations, and regional regimes. The main value of the analysis is therefore not only to confirm directional warming, but to show where in the distribution the most pronounced changes occur and how that structure differs across space.

### 3.2. Regional synthesis of quantile-dependent change

As a descriptive cross-station summary, warm days and warm nights both increase, whereas cool days and cool nights decline. However, the distributional structure of change differs markedly among indices. These pooled panels are useful for orienting the reader to the overall geometry of change, but they are not treated as the main inferential layer because averaging annual counts across stations can mute climatically meaningful spatial heterogeneity.

![Figure 3. Regional quantile-coefficient panels.](../outputs/figures/ijoc_regional_quantile_panels.png)

*Figure 3. Descriptive regional-average quantile-regression profiles for (a) warm days, (b) warm nights, (c) cool days, and (d) cool nights. For each year, station values were first averaged across the network and quantile regression was then fitted to that pooled annual series. Black points show quantile-specific slopes, grey shading indicates analytic 95% confidence intervals from the quantile-regression fit, and red lines show the OLS mean trend and its confidence bounds. The display range is restricted to q = 0.05-0.95. These panels are descriptive summaries, whereas primary inference in the manuscript rests on station-level quantile regression and cross-station synthesis diagnostics.*

**Table 1.** Regional synthesis of mean trends, focal quantile trends, field significance, and dominant-cluster behaviour. The column `Approx. 1961-2024 change` is a simple slope-based translation over the full record length and is included only as an approximate measure of climatic magnitude.

| Index | OLS slope | q0.05 | q0.50 | q0.95 | Approx. 1961-2024 change | FDR sig. q0.10 | FDR sig. q0.50 | FDR sig. q0.90 | Dominant cluster | Dominant-cluster median Delta1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Warm days | +9.53 | +7.81 | +8.15 | +12.05 | +60.0 | 27 | 27 | 17 | 2 | +2.26 |
| Warm nights | +13.58 | +13.96 | +13.65 | +12.55 | +85.5 | 29 | 28 | 25 | 2 | +1.70 |
| Cool days | -2.56 | -1.51 | -2.01 | -8.05 | -16.1 | 0 | 17 | 26 | 1 | -6.64 |
| Cool nights | -2.98 | -1.77 | -2.47 | -6.78 | -18.8 | 12 | 23 | 25 | 2 | -8.21 |

Three points are especially important at this descriptive pooled level. First, warm nights show the strongest mean increase, corresponding to an approximate slope-based rise of about 85.5 days per year across the 1961-2024 record when the OLS estimate is translated over the full record length. Second, warm days display a clear upper-tail ordering (`q0.95 > q0.50 > q0.05`), indicating classical upper-tail amplification. Third, cool days and cool nights show the opposite structure: their strongest declines occur in upper quantiles, implying a rapid loss of years historically characterized by many cool events.

This structure is useful as a compact orientation device, but it should be read together with the station-level results that follow. For example, the warm-day OLS slope of +9.53 days year^-1 per decade implies a climatologically large pooled increase when translated across the full record, yet the later station-based sections show that this increase is not spatially uniform and that its quantile geometry varies materially across the network.

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

Regional means necessarily conceal strong station-scale contrasts. To connect the clustering results to physically interpretable slope geometries without introducing manual selection bias, one representative station per cluster was selected algorithmically as the observed station closest to the corresponding cluster centroid in the standardized clustering feature space used for the baseline regionalization. Because the clustering itself already applies collinearity screening, the representative-station rule is tied to the same screened feature set rather than to an ad hoc post hoc criterion.

**Table 3.** Reproducibly selected representative stations, defined as the closest observed member to each cluster centroid in the screened baseline clustering feature space. The exact station list and associated distances are generated automatically in `outputs/tables/representative_station_selection.csv` and should be inserted from that file after each pipeline run.

![Figure 4. Representative-station comparisons.](../outputs/figures/ijoc_station_comparisons/main_figure_representative_stations.png)

*Figure 4. Multi-panel comparison of algorithmically selected representative stations for (a) warm days, (b) warm nights, (c) cool days, and (d) cool nights. Within each index, one station per cluster is shown, chosen as the closest observed member to the corresponding cluster centroid in the screened clustering feature space. The panels therefore visualize cluster-defined quantile-slope geometries using a reproducible medoid-like rule rather than hand-picked examples.*

These representative stations clarify the logic of the cluster structure while remaining fully reproducible. Because each panel shows centroid-nearest observed members rather than manually curated examples, the figure should be read as a visualization of typical cluster geometry. Warm-day clusters still separate upper-tail amplification from weakening or flatter profiles; warm-night clusters separate broad-based nocturnal intensification from stronger upper-tail growth; and the cool indices continue to distinguish contraction-dominated regimes from mixed or locally positive profiles. The clusters are therefore climatologically interpretable rather than merely statistically convenient, but that interpretation now rests on a transparent selection rule that can be regenerated exactly from the pipeline.

### 3.5. Spatial expression and regionalization

The spatial pattern of `Delta1` is coherent enough to support regional interpretation, but not smooth enough to justify treating the country as a single spatial field.

![Figure 5. Main Delta1 map panels.](../outputs/figures/ijoc_main_delta1_maps.png)

*Figure 5. Multi-panel station maps of `Delta1 = slope(q0.95) - slope(q0.05)` for (a) warm days, (b) warm nights, (c) cool days, and (d) cool nights. Warm indices generally show positive Delta1 values, whereas cool indices show widespread negative Delta1 values, consistent with upper-tail contraction.*

Warm days show a broad tendency toward positive `Delta1`, but with strong contrast between central, western, and southern stations. Warm nights also show positive `Delta1` at several stations, though less uniformly, and with one exceptional outlier in Shiraz. Cool days exhibit widespread negative `Delta1`, indicating upper-tail collapse across most of the network. Cool nights are similarly negative at most stations, but with notable positive anomalies at Shahrekord and, more modestly, Esfahan.

This pattern aligns with the formal field-significance results: intermediate quantiles are the most consistently significant across stations, yet Moran's I is not significant for any of the 20 slope fields tested. Those field-significance summaries should be read as analytic-screening diagnostics rather than as the main inferential basis for the tail quantiles, for which the manuscript gives priority to bootstrap confidence intervals. Accordingly, the spatial pattern is best interpreted as a set of regionally meaningful but non-smooth station responses, consistent with topographic and climatic heterogeneity rather than a single national-scale spatial mode.

### 3.5.1. Additional structural and spatial diagnostics

The main-text regionalization results are further supported by two additional diagnostic layers presented in the Supplementary Material. First, quantile-specific dendrograms at `q0.05`, `q0.50`, and `q0.95` (Supplementary Figures S2-S4) show that station similarity is not fixed across the distribution. Lower-tail structures are comparatively compressed, whereas upper-tail structures become more strongly separated, especially for warm nights and, in a different direction, for cool nights. This supports the interpretation that the clustering is not merely a multivariate artefact; rather, meaningful station-group separation is already visible in one-dimensional quantile-specific slope space.

Second, thin-plate-spline quantile maps for `q0.05`, `q0.10`, `q0.50`, `q0.90`, and `q0.95` (Supplementary Figures S5-S9) clarify how the geography of change evolves from lower to upper quantiles. Lower quantiles show smoother and more regionally coherent behaviour, whereas upper quantiles reveal sharper hotspots for warm indices and more focused contraction zones for cool indices. Because Moran's I does not support a strong national-scale spatial field, these interpolated panels are interpreted as visualization aids rather than as primary inferential products. Their scientific role is therefore limited: they help the reader see how station-level contrasts are arranged in space, but they should not be read as precise surface estimates between stations or as evidence for spatial significance away from the observation points. Even so, they provide an important consistency check: the main station-based conclusions remain visible when the slope fields are viewed as continuous spatial summaries.

### 3.6. Stronger explanatory context from geography and elevation

The driver analysis already indicated that warm-night slopes are negatively associated with elevation, with standardized betas of -0.507 at q0.05 and -0.479 at q0.50. An elevation-class summary reinforces that result. Mean warm-night q0.50 slopes are about +18.43 days year^-1 per decade in lowland stations, +14.25 in mid-elevation stations, and only +6.58 in highland stations. Lower-elevation environments therefore appear especially prone to broad-based nocturnal warming.

The elevation-class contrasts also refine the interpretation of day versus night behaviour. Warm-day upper-tail amplification is strongest in the highland class (`Delta1 ≈ +7.48`), whereas warm-night amplification is strongest in the mid-elevation class (`Delta1 ≈ +6.25`) but weak in both lowland and highland classes. For cool nights, the highland class shows a much weaker contraction (`Delta1 ≈ -3.44`) than the lowland and mid-elevation classes (both about `-8.5` to `-8.7`). These contrasts suggest that the geography of day-night asymmetry cannot be reduced to a single monotonic elevation effect. Different parts of the thermal-extreme distribution appear to be modulated by different local controls.

From a physical standpoint, the warm-night patterns are consistent with greater nocturnal heat retention at lower elevations and in settings where radiative cooling may be reduced, but the present analysis does not directly identify the responsible mechanisms. The warm-day results, by contrast, suggest that some elevated and interior stations remain especially sensitive to upper-tail daytime amplification, possibly because exceptionally hot years intensify more rapidly than the broader mean regime. The cool-day and cool-night indices show the opposite direction of change, but again with stronger contraction in upper quantiles, indicating that the most historically cool years are disappearing fastest. These interpretations should therefore be read as physically plausible explanations consistent with the statistical patterns, not as formal process attribution.

### 3.7. Robustness, methodological sensitivity, and temporal nonlinearity

One of the strongest assets of the study is that the principal conclusions do not depend on a single analytical choice.

![Figure 6. Robustness synthesis.](../outputs/figures/ijoc_robustness_synthesis.png)

*Figure 6. Synthesis of robustness diagnostics. Panel (a) summarizes reference-period sensitivity, panel (b) summarizes bootstrap-method sensitivity with moving-block resampling as the baseline and `meboot` as the comparison, panel (c) shows interpolation-method agreement for the configured spatial-visualization workflow, and panel (d) reports clustering robustness as adjusted Rand index. The principal structural results remain stable despite some metric-level sensitivity.*

This robustness section should be read together with the pipeline design itself. The executable workflow does not merely compute the main station-level quantile fits; it also runs data-quality and homogeneity diagnostics, a station-exclusion sensitivity based on those homogeneity flags, a higher-replication bootstrap-depth check, an alternative-clustering sensitivity module, and a suite of advanced publication analyses. Reporting those layers explicitly matters for publication quality because they demonstrate that the manuscript's core claims were stress-tested within the same reproducible analytical environment rather than through disconnected ad hoc reruns. It also clarifies the role of the homogeneity analysis: in this study it is a screening and bounding device, not a claim that daily homogenization was solved or unnecessary.

The robustness synthesis shows that warm-night metrics are most sensitive to reference-period choice, whereas cool-day and cool-night bootstrap summaries are more sensitive to the choice between moving-block and `meboot` resampling, especially in upper-tail means and confidence bounds. Interpolation agreement is strongest for cool nights and warm nights under the `linear_rbf` versus `linear` comparison, and clustering robustness is particularly high for warm nights (ARI = 0.93) and cool nights (ARI = 0.82) even when the parsimonious slope-only baseline is compared with a richer uncertainty-aware rerun. These diagnostics support a confident qualitative interpretation while still motivating caution for the most extreme tails.

An additional sensitivity test excluded all stations flagged by any detrended homogeneity diagnostic (19 of 30 stations; Supplementary Figure S10; Supplementary Table S11). The principal directional conclusions were retained: warm indices remained positive and cool indices remained negative at the regional scale, and warm-day upper-tail amplification remained evident. However, some magnitude changes were non-negligible, particularly for warm nights, where the regional tail-contrast metric became more positive after exclusion. This result should not be over-sold as a replacement for formal homogenization. Its role is narrower and more defensible: it shows that the manuscript's broad climatic conclusions are not driven solely by the flagged subset, while also making clear that exact tail-asymmetry magnitudes remain conditional on unresolved record homogeneity.

A second targeted sensitivity check addressed bootstrap depth directly by increasing the number of resamples from 200 to 400 for `warm_days` and `warm_nights` (Supplementary Figure S11; Supplementary Tables S12-S13). The result is reassuringly stable. Regional bootstrap means for `q0.05`, `q0.50`, `q0.95`, and `Delta1` changed by no more than about `0.06 days year^-1 per decade` in absolute mean value, and station-level correspondence between the 200- and 400-replicate summaries remained extremely high (`r = 0.998-1.000`). Confidence-interval widths changed somewhat more than the means, as expected, but their station-level correlations still remained high (`r = 0.975-0.992`). This additional check does not remove all tail uncertainty, but it does show that the manuscript's main warm-extreme messages are not an artefact of a shallow bootstrap configuration.

A third sensitivity check examined whether the regionalization itself depends strongly on the chosen clustering method and feature encoding (Supplementary Figure S12; Supplementary Tables S14-S15). The answer is mixed, and therefore informative. `warm_days` remained comparatively stable across alternative methods (`ARI = 0.63-0.92`), and `warm_nights` also retained a strong common structure for most alternatives (`ARI = 0.88-0.93`), although complete linkage produced a noticeably less similar partition (`ARI = 0.43`). `cool_nights` showed moderate method dependence, with one alternative reproducing the baseline almost exactly (`ARI = 1.00`) but others yielding only partial agreement (`ARI = 0.34-0.48`). `cool_days` were the least stable (`ARI = 0.05-0.35`), indicating that their regional grouping should be read as exploratory rather than as a uniquely resolved partition. This pattern does not invalidate the station-scale signal; rather, it shows that some regional structures, especially for cooling-related indices, are weaker and less crisply separated than those for the warm indices, and that the clusters should be interpreted as descriptive regional summaries rather than as sharply identified latent climate regimes.

Temporal nonlinearity provides an additional layer of evidence.

![Figure 7. Split-period comparison.](../outputs/figures/ijoc_split_period_comparison.png)

*Figure 7. Comparison of descriptive regional-average trend summaries between 1961-1990 and 1991-2024 for the focal quantiles and OLS trend. The later period shows marked intensification of warm extremes and stronger collapse of cool extremes, especially in upper quantiles.*

The split-period analysis shows that all four indices behave differently after 1990. Warm-day OLS slopes increase from about +0.63 days year^-1 per decade in 1961-1990 to +19.47 in 1991-2024, while warm-night OLS slopes rise from +4.44 to +19.49. Cool-day OLS slopes shift from +1.78 in the early period to -7.86 in the later period, and cool nights shift from -0.44 to -3.33. The strongest acceleration appears in the warm-day upper tail and in the cool-day upper-tail decline. The distribution-aware patterns detected in the full-period analysis are therefore not artefacts of long averaging; they reflect a marked strengthening of recent change.

### 3.8. Comparison with the wider climatological literature

The present findings are broadly consistent with the wider literature in four respects. First, the strong rise in warm extremes agrees with previous assessments over Iran and neighbouring dryland regions `(Soltani et al., 2016; Zhang et al., 2005a; Vaghefi et al., 2019)`. Second, the stronger and more widespread increase in warm nights is consistent with Mediterranean and Middle Eastern evidence suggesting strong nocturnal warming and increasing heat stress `(Kostopoulou and Jones, 2005; Zittis et al., 2016; Ahmadalipour and Moradkhani, 2018; Francis and Fonseca, 2024)`. Third, the marked upper-tail strengthening of warm days is aligned with studies showing that high-impact heat extremes and heatwave characteristics often intensify more rapidly than central-tendency measures `(Fischer and Schär, 2010; Perkins, 2015; Perkins-Kirkpatrick and Lewis, 2020)`. Fourth, the collapse of cool extremes, especially in upper quantiles, is consistent with the broader shift away from historically cool thermal regimes documented in global and regional extreme-index studies `(Alexander et al., 2006; Donat et al., 2013; Dunn et al., 2020)`.

The distinctive contribution of the present study is not simply that it confirms warming over Iran. Rather, it offers a carefully integrated regional synthesis showing that the most informative part of the signal lies in the **shape** of change: warm indices do not merely rise, and cool indices do not merely fall; both vary meaningfully across quantiles. This matters especially in a region where human heat exposure is already emerging as a critical climate risk and where hot-night and humid-heat hazards are becoming more prominent in the Persian Gulf sector `(Pal and Eltahir, 2016; Raymond et al., 2020; Raymond et al., 2024)`. It also builds directly on earlier quantile-based temperature studies from Europe and China, which demonstrated the added value of quantile regression for identifying non-uniform warming across the distribution `(Barbosa et al., 2011; Fan, 2014)`. In this respect, the study contributes a robust, reproducible regional evaluation for an arid and topographically complex setting rather than claiming a new conceptual quantile framework, and thereby extends the interpretation beyond standard trend reporting, including scenarios in which present-day MENA warming is projected to evolve toward much more severe heatwave regimes `(Zittis et al., 2021)`.

### 3.9. Limitations

The limitations of the study should be stated explicitly. The station network is finite, the strongest tails remain the most uncertain part of the distribution, and the driver analysis is intentionally simple. Most importantly, the daily temperature records were screened for possible inhomogeneity but were not formally homogenized with a reference-based daily adjustment procedure. That choice avoided imposing a poorly identified breakpoint correction on percentile-based extremes, but it also means that some tail magnitudes may still reflect unresolved non-climatic shifts. Bootstrap inference is more defensible than asymptotic analytic intervals for the focal tails, but even bootstrap intervals for short annual series should still be interpreted as uncertainty summaries rather than exact small-sample coverage guarantees. Interpolated surfaces are used as visualization tools rather than primary inferential products, and the quantile-map figures are now framed around station-level slopes rather than around any separate significance overlay on the interpolated surface. The clustering stage is exploratory even after collinearity screening, because regional partitions can still depend on metric choice, feature encoding, and the limited station sample. The period split is informative but still descriptive and does not amount to formal breakpoint attribution. More broadly, the study does not estimate compound heat-stress metrics such as wet-bulb temperature or explicit heatwave-event properties, both of which are increasingly important in MENA impact studies `(Nairn and Fawcett, 2015; Raei et al., 2018; Ahmadalipour and Moradkhani, 2018; Raymond et al., 2020)`. These limitations do not invalidate the main conclusions, but they do define the appropriate scope of interpretation. A logical next step would be a dedicated daily homogenization study using station metadata and neighbor-based reference construction before any attempt to sharpen the absolute tail magnitudes.

## 4. Conclusions

This study shows that thermal-extreme change across Iran is not adequately summarized by a uniform shift in the mean alone. An integrated quantile-based climatological synthesis instead shows that warm days and warm nights both increase, but with different internal structures: warm nights are stronger overall, whereas warm days show clearer upper-tail amplification. Cool days and cool nights both decline, with the strongest contractions concentrated in upper quantiles, indicating rapid loss of years historically characterized by many cool events.

Three conclusions are especially important. First, the mean signal alone understates the extent of distributional asymmetry in thermal-extreme change. Second, station-scale heterogeneity is climatologically meaningful and can be summarized effectively through exploratory cluster-defined quantile geometries, provided that those groupings are not over-interpreted as uniquely identified regimes. Third, the main conclusions are robust to alternative methodological choices and are strengthened, rather than weakened, by split-period evidence showing stronger recent change, although the exact magnitude of some tail summaries should remain conditional on future daily homogenization work.

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

Barbosa SM, Scotto MG, Alonso AM (2011) Summarising changes in air temperature over Central Europe by quantile regression and clustering. *Natural Hazards and Earth System Sciences* 11:3227-3233. https://doi.org/10.5194/nhess-11-3227-2011

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

Nairn JR, Fawcett RJB (2015) The excess heat factor: a metric for heatwave intensity and its use in classifying heatwave severity. *International Journal of Environmental Research and Public Health* 12(1):227-253. https://doi.org/10.3390/ijerph120100227

Pal JS, Eltahir EAB (2016) Future temperature in southwest Asia projected to exceed a threshold for human adaptability. *Nature Climate Change* 6:197-200. https://doi.org/10.1038/nclimate2833

Perkins SE (2015) A review on the scientific understanding of heatwaves: their measurement, driving mechanisms, and changes at the global scale. *Atmospheric Research* 164-165:242-267. https://doi.org/10.1016/j.atmosres.2015.05.014

Perkins-Kirkpatrick SE, Lewis SC (2020) Increasing trends in regional heatwaves. *Nature Communications* 11:3357. https://doi.org/10.1038/s41467-020-16970-7

Raei E, Nikoo MR, AghaKouchak A, Mazdiyasni O, Sadegh M (2018) GHWR, a multi-method global heatwave and warm-spell record and toolbox. *Scientific Data* 5:180206. https://doi.org/10.1038/sdata.2018.206

Raymond C, Matthews T, Horton RM (2020) The emergence of heat and humidity too severe for human tolerance. *Science Advances* 6(19):eaaw1838. https://doi.org/10.1126/sciadv.aaw1838

Raymond C, Matthews T, Tuholske C (2024) Evening humid-heat maxima near the southern Persian/Arabian Gulf. *Communications Earth & Environment* 5:591. https://doi.org/10.1038/s43247-024-01763-3

Reich BJ (2012) Spatiotemporal quantile regression for detecting distributional changes in environmental processes. *Journal of the Royal Statistical Society: Series C (Applied Statistics)* 61(4):535-553. https://doi.org/10.1111/j.1467-9876.2011.01025.x

Russo S, Sterl A (2011) Global changes in indices describing moderate temperature extremes from the daily output of a climate model. *Journal of Geophysical Research: Atmospheres* 116:D03104. https://doi.org/10.1029/2010JD014727

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




