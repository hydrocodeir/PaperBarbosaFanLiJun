# Quantile Regression + Bootstrap + Clustering Report

## Data audit
- Number of stations: **30**
- Number of years in uploaded data: **64**
- Year range: **1961-2024**
- Station-index series below QR minimum length (3 years): **0**
- Station-index series below publication recommendation (20+ years): **0**

## Highest Delta1 stations by index
### Warm Days
- Esfahan: Delta1 = 16.577
- Ahvaz: Delta1 = 15.633
- Zanjan: Delta1 = 14.107
- Abadan: Delta1 = 12.156
- Sanandaj: Delta1 = 11.390

### Warm Nights
- Shiraz: Delta1 = 25.903
- Bam: Delta1 = 12.996
- Tehran (Mehrabad Airport): Delta1 = 11.915
- Zabol: Delta1 = 11.379
- Sabzevar: Delta1 = 11.278

### Cool Days
- Birjand: Delta1 = -2.138
- Kerman: Delta1 = -2.986
- Zahedan: Delta1 = -3.298
- Arak: Delta1 = -4.474
- Gorgan: Delta1 = -4.537

### Cool Nights
- Shahrekord: Delta1 = 11.961
- Esfahan: Delta1 = 2.857
- Tabriz: Delta1 = -2.255
- Gorgan: Delta1 = -2.478
- Orumiyeh: Delta1 = -3.421

## Sensitivity Checks
### Quantile τ = 0.05
- insufficient: 120

### Quantile τ = 0.95
- insufficient: 120

## Cluster sizes
### Warm Days
- Cluster 1: 2 stations
- Cluster 2: 12 stations
- Cluster 3: 4 stations
- Cluster 4: 12 stations

### Warm Nights
- Cluster 1: 6 stations
- Cluster 2: 15 stations
- Cluster 3: 8 stations
- Cluster 4: 1 stations

### Cool Days
- Cluster 1: 27 stations
- Cluster 2: 1 stations
- Cluster 3: 1 stations
- Cluster 4: 1 stations

### Cool Nights
- Cluster 1: 2 stations
- Cluster 2: 26 stations
- Cluster 3: 1 stations
- Cluster 4: 1 stations

## Clustering Robustness
- Comparison between the configured clustering setup and a reduced-feature rerun.
### Warm Days
- Stations compared: 30
- Adjusted Rand Index: 0.632
- Same-label fraction: 0.833

### Warm Nights
- Stations compared: 30
- Adjusted Rand Index: 0.934
- Same-label fraction: 0.200

### Cool Days
- Stations compared: 30
- Adjusted Rand Index: 0.685
- Same-label fraction: 0.133

### Cool Nights
- Stations compared: 30
- Adjusted Rand Index: 0.821
- Same-label fraction: 0.100

## Spatial Inference
- τ = 0.05: raw significant station-results = 0, FDR-retained = 0
- τ = 0.10: raw significant station-results = 75, FDR-retained = 68
- τ = 0.50: raw significant station-results = 98, FDR-retained = 95
- τ = 0.90: raw significant station-results = 96, FDR-retained = 93
- τ = 0.95: raw significant station-results = 0, FDR-retained = 0

- Moran's I tests with p < 0.05: 0/20 station-slope fields

## Method Sensitivity
- Reference-period sensitivity:
- cool_days | slope_0.05: mean abs diff = 0.732
- cool_nights | slope_0.05: mean abs diff = 1.497
- cool_days | slope_0.50: mean abs diff = 1.790
- warm_days | Delta1: mean abs diff = 2.129

- Bootstrap-method sensitivity:
- cool_days | boot_ci_low_0.05: mean abs diff = 0.833
- cool_days | boot_mean_0.05: mean abs diff = 0.896
- cool_nights | boot_mean_0.05: mean abs diff = 1.659
- cool_nights | boot_ci_high_0.05: mean abs diff = 1.758

- Interpolation-method agreement (top correlations):
- cool_nights | τ=0.95 | linear_rbf vs linear: r = 0.883, RMSE = 1.977
- warm_nights | τ=0.05 | linear_rbf vs linear: r = 0.841, RMSE = 2.583
- cool_nights | τ=0.05 | linear_rbf vs linear: r = 0.832, RMSE = 0.877
- warm_nights | τ=0.95 | linear_rbf vs linear: r = 0.824, RMSE = 4.331

## Driver Analysis
- cool_days | Delta1 ~ longitude: std beta = +0.317, p = 0.130, rho = +0.246
- warm_days | Delta1 ~ latitude: std beta = -0.281, p = 0.185, rho = -0.378
- cool_nights | Delta1 ~ elevation: std beta = +0.217, p = 0.262, rho = +0.074
- warm_nights | Delta1 ~ longitude: std beta = +0.216, p = 0.294, rho = +0.310
- warm_nights | Delta1 ~ latitude: std beta = -0.176, p = 0.391, rho = -0.214
- cool_nights | Delta1 ~ latitude: std beta = -0.127, p = 0.545, rho = -0.133

## Regional Composites
- warm_nights | cluster 4: median Delta1 = +25.903 (n=1)
- warm_days | cluster 1: median Delta1 = +13.983 (n=2)
- cool_nights | cluster 4: median Delta1 = +11.961 (n=1)
- warm_nights | cluster 3: median Delta1 = +10.848 (n=8)
- warm_days | cluster 4: median Delta1 = +8.217 (n=12)
- warm_days | cluster 2: median Delta1 = +2.262 (n=12)

