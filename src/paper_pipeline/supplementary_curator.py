"""Curated supplementary evidence, rendered directly from audited numeric tables."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from .publication_figures import theme, save, map_base, INDEXES, NAMES, COLORS
from .compound_partition import file_hash

SPECS = [
('figS01_data_quality','Coverage and homogeneity diagnostics','Calendar-based valid-day fractions for temperature and precipitation, calculated from the screened aggregates for all available station-years; boxes summarize station medians. Detrended homogeneity flags use nominal p < 0.05 for annual mean temperature; flags are diagnostic, not proof of artificial breaks.'),
('figS02_thermal_robustness','Sensitivity of thermal estimates','(a) Mean station interval-width ratio for 400 versus 200 moving-block replicates. (b) Mean absolute difference between saved moving-block and maximum-entropy bootstrap means. (c) Change in network slopes after internal-consistency screening. (d) Change after excluding homogeneity-flagged stations. Units in (b–d) are days per decade. The two historical alternative bootstrap ensembles were checked through their saved station summaries, not regenerated.'),
('figS03_station_quantile_maps','Station quantile slopes','Station estimates at the 0.10, 0.50 and 0.90 quantiles of annual thermal counts, on one common symmetric scale. No interpolation, area weighting or local significance symbols are used. Units are days per decade.'),
('figS04_median_precision','Bootstrap precision of median trends','Absolute bootstrap mean divided by bootstrap standard deviation at the median. This descriptive precision ratio is dimensionless; it is neither an emergence date nor a calibrated detection probability. No threshold-based significance labels are shown.'),
('figS05_spatial_diagnostics','Spatial dependence and multiplicity','(a) Moran’s I with five-nearest-neighbor weights. (b) Numbers of locally retained tests after Benjamini–Hochberg adjustment separately within each index–quantile family. NA means analytic tail intervals and tests are unavailable, not that zero stations are significant. Historical analytic probabilities are approximate and do not adjust for serial dependence. These counts do not constitute a field-significance test; permutation probabilities are tabulated separately.'),
('figS06_cluster_composites','Exploratory cluster composites','Median station slopes and upper-minus-lower contrasts within four fitted clusters per index. Cluster labels are specific to each index and have no shared climatic ordering. Cluster separation is partly induced by the fitted features and is not independent validation; single-station groups remain visible through their sample sizes.'),
('figS07_cluster_maps','Exploratory cluster membership','Station assignments from standardized, screened quantile-slope features and average-linkage Euclidean clustering. Colors distinguish categories only within each panel; clusters are not fixed climate regions. Alternative-method agreement and within-cluster sample sizes are reported in the supplementary tables.'),
('figS08_representative_profiles','Representative station profiles','One observed station per index-specific cluster, selected by the archived nearest-centroid rule. Lines show the full 0.10–0.90 quantile profile. These selected illustrations are not additional independent discoveries; all station coefficients and bootstrap intervals are available in the data catalog.'),
('figS09_internal_warming_association','Association with the network temperature anomaly','(a) Internally derived annual network temperature anomaly relative to 1991–2007. (b–e) Fixed-baseline network thermal-count response coefficients at three quantiles, with available median pointwise analytic intervals, in days per degree Celsius. Tail analytic intervals are unavailable and must not be inferred from the unadorned points. The predictor and outcomes share observations and temporal trends; the associations do not attribute change to external forcing.'),
('figS10_historical_joint_rarity','Historical joint-rarity sensitivity','Percentage of the available station network exceeding three empirical inverse-joint-probability cutoffs. Labels denote within-record rarity classes, not stable design return periods. This historical calculation uses unscreened daily data, available-row coverage denominators and a varying station network. It is retained as a definition sensitivity and is distinct from the main paper’s fixed marginal AND event.'),
('figS11_geographical_associations','Geographical associations of thermal asymmetry','Standardized multiple-regression coefficients relating the upper-minus-lower slope contrast to latitude, longitude and elevation. These are descriptive associations; no causal driver attribution or spatially adjusted significance is implied.'),
]

def heat(ax,frame,title,limit=None,cmap='RdBu_r',fmt='.2f',bounds=None):
    z=frame.to_numpy(dtype=float)
    limit=float(np.nanmax(abs(z))) if limit is None else limit
    norm=TwoSlopeNorm(vmin=-limit,vcenter=0,vmax=limit) if cmap=='RdBu_r' else (plt.Normalize(*bounds) if bounds else None)
    im=ax.imshow(z,aspect='auto',cmap=cmap,norm=norm)
    ax.set(xticks=range(len(frame.columns)),xticklabels=frame.columns,yticks=range(len(frame)),yticklabels=frame.index,title=title)
    for y in range(z.shape[0]):
        for x in range(z.shape[1]):
            if not np.isfinite(z[y,x]):
                ax.text(x,y,'NA',ha='center',va='center',fontsize=6.5,color='#555555');continue
            rgb=im.cmap(im.norm(z[y,x]))[:3]; luminance=sum(a*b for a,b in zip(rgb,[.2126,.7152,.0722]))
            label=format(z[y,x],fmt)
            if label in ['-0.00','-0.0']:label=label[1:]
            ax.text(x,y,label,ha='center',va='center',fontsize=6.5,color='white' if luminance<.5 else '#222222')
    return im

def create_supplementary_figures(root,out):
    theme();records=[];inputs=set()
    def load(name,folder='tables'):
        p=root/'outputs'/folder/(name+'.csv');inputs.add(p);return pd.read_csv(p)
    stations=pd.read_csv(root/'data/stationsInfo.csv')
    countries=gpd.read_file(root/'data/Iran_Sea_Ne.geojson').to_crs(4326)
    qr=load('qr_focus_slopes_and_bootstrap_summary').merge(stations,on=['station_id','station_name'],validate='many_to_one')
    # S1: calendar coverage, explicitly avoiding the historical available-row denominator.
    aggregates=load('screened_annual_seasonal_aggregates','publication_v2/tables')
    annual=aggregates.loc[aggregates.definition=='annual'].copy()
    annual['Temperature']=100*annual.temperature_days/annual.expected_days
    annual['Precipitation']=100*annual.precip_days/annual.expected_days
    fig,axes=plt.subplots(1,2,figsize=(7.1,3),layout='constrained')
    selected=['Temperature','Precipitation']
    for j,col in enumerate(selected):
        values=annual.groupby('station_id')[col].median().dropna()
        axes[0].boxplot(values,positions=[j],widths=.45,patch_artist=True,boxprops={'facecolor':['#84ADB8','#D4B278'][j]},flierprops={'markersize':2})
    axes[0].set(xticks=range(len(selected)),xticklabels=[c.replace('_coverage','').replace('_pct','').replace('_',' ') for c in selected],ylabel='Median station-year coverage (%)',title='(a) Calendar-day coverage')
    h=load('data_homogeneity_tests_station_summary'); columns=[c for c in h if 'detrended_pvalue' in c]
    counts=[int((h[c]<.05).sum()) for c in columns];counts.append(int((h[columns]<.05).any(axis=1).sum()))
    axes[1].bar(range(4),counts,color=['#659DAE']*3+['#AE6357'])
    axes[1].set(xticks=range(4),xticklabels=['Pettitt','SNHT','Buishand','Any'],ylabel='Flagged stations',title='(b) Detrended screening')
    for i,n in enumerate(counts):axes[1].text(i,n+.5,str(n),ha='center',fontsize=7)
    save(fig,out,SPECS[0][0],records)
    # S2: four sensitivity families; avoid duplicating the main fixed-baseline plot.
    fig,axes=plt.subplots(2,2,figsize=(7.1,5.4),layout='constrained')
    d=load('bootstrap_depth_sensitivity_summary');d=d.loc[d.metric.str.startswith('boot_ci_width_')].copy();d['ratio']=d.mean_400/d.mean_200
    f=d.pivot(index='index_name',columns='metric',values='ratio').loc[INDEXES];f.columns=['q10','q50','q90','Δ']
    f.index=NAMES;heat(axes[0,0],f,'(a) Interval width: 400 / 200',cmap='YlGnBu',bounds=(.9,1.1))
    d=load('bootstrap_method_sensitivity_summary');f=d.loc[d.metric.str.startswith('boot_mean_')].pivot(index='index_name',columns='metric',values='mean_abs_diff').loc[INDEXES];f.columns=['q10','q90'];f.index=NAMES
    heat(axes[0,1],f,'(b) Bootstrap method difference',cmap='YlGnBu')
    d=load('temperature_internal_consistency_quantile_sensitivity').set_index('index_name').loc[INDEXES];f=d[['slope_0.10_difference','slope_0.50_difference','slope_0.90_difference']].copy();f.columns=['q10','q50','q90'];f.index=NAMES
    heat(axes[1,0],f,'(c) Screened minus observed',limit=2)
    d=load('homogeneity_flag_exclusion_sensitivity').set_index('index_name').loc[INDEXES];f=d[['slope_0.10_difference_exclude_minus_all','slope_0.50_difference_exclude_minus_all','slope_0.90_difference_exclude_minus_all']].copy();f.columns=['q10','q50','q90'];f.index=NAMES
    heat(axes[1,1],f,'(d) Excluded minus full network',limit=2)
    save(fig,out,SPECS[1][0],records)
    fig,axes=plt.subplots(4,3,figsize=(7.1,9.4),layout='constrained')
    limit=float(np.ceil(np.max(np.abs(qr[[f'slope_{q:.2f}' for q in [.1,.5,.9]]].to_numpy()))/5)*5)
    for i,idx in enumerate(INDEXES):
        d=qr.loc[qr.index_name==idx]
        for j,q in enumerate([.1,.5,.9]):
            ax=axes[i,j];map_base(ax,countries,f'{NAMES[i]} · q{q:.2f}')
            im=ax.scatter(d.longitude,d.latitude,c=d[f'slope_{q:.2f}'],s=9,edgecolor='white',linewidth=.15,cmap='RdBu_r',vmin=-limit,vmax=limit,zorder=4)
            ax.tick_params(labelsize=5)
    fig.colorbar(im,ax=axes,shrink=.6,pad=.01,label='Days per decade')
    save(fig,out,SPECS[2][0],records)
    precision=load('climate_signal_emergence_station_level');precision=precision.loc[precision.metric=='0.50'].merge(stations,on=['station_id','station_name'])
    fig,axes=plt.subplots(2,2,figsize=(7.1,6.1),layout='constrained'); vmax=float(np.ceil(precision.signal_to_noise.max()))
    for i,(idx,ax) in enumerate(zip(INDEXES,axes.flat)):
        d=precision.loc[precision.index_name==idx];map_base(ax,countries,f'({chr(97+i)}) {NAMES[i]}')
        im=ax.scatter(d.longitude,d.latitude,c=d.signal_to_noise,s=16,edgecolor='white',linewidth=.2,cmap='viridis',vmin=0,vmax=vmax,zorder=4)
    fig.colorbar(im,ax=axes,shrink=.7,pad=.01,label='|Bootstrap mean| / bootstrap SD')
    save(fig,out,SPECS[3][0],records)
    fig,axes=plt.subplots(1,2,figsize=(7.1,2.8),layout='constrained')
    moran=load('spatial_autocorrelation_moran');f=moran.pivot(index='index_name',columns='tau',values='moran_i').loc[INDEXES];f.index=NAMES;f.columns=['q10','q50','q90'];heat(axes[0],f,'(a) Moran’s I',limit=.5)
    fdr=load('station_significance_fdr');f=fdr.groupby(['index_name','tau']).fdr_reject.sum(min_count=1).unstack().loc[INDEXES];f.index=NAMES;f.columns=['q10','q50','q90'];heat(axes[1],f,'(b) Retained local tests',cmap='YlGnBu',fmt='.0f')
    save(fig,out,SPECS[4][0],records)
    composites=load('regional_cluster_composites');clusters=load('cluster_assignments').merge(stations,on=['station_id','station_name']);metrics=['slope_0.10','slope_0.50','slope_0.90','Delta1']
    fig,axes=plt.subplots(2,2,figsize=(7.1,5.1),layout='constrained');limit=float(np.ceil(composites['median'].abs().max()/5)*5)
    for i,(idx,ax) in enumerate(zip(INDEXES,axes.flat)):
        local=composites.loc[composites.index_name==idx];f=local.pivot(index='cluster',columns='metric',values='median')[metrics]
        counts=clusters.loc[clusters.index_name==idx].cluster.value_counts();f.index=[f'C{k} (n={counts[k]})' for k in f.index];f.columns=['q10','q50','q90','Δ']
        im=heat(ax,f,f'({chr(97+i)}) {NAMES[i]}',limit=limit,fmt='.1f')
    fig.colorbar(im,ax=axes,shrink=.7,pad=.02,label='Days per decade');save(fig,out,SPECS[5][0],records)
    colors=['#287E9B','#C47842','#7666A8','#468567'];fig,axes=plt.subplots(2,2,figsize=(7.1,6.2),layout='constrained')
    for i,(idx,ax) in enumerate(zip(INDEXES,axes.flat)):
        map_base(ax,countries,f'({chr(97+i)}) {NAMES[i]}')
        for cluster,color in enumerate(colors,1):
            d=clusters.loc[(clusters.index_name==idx)&(clusters.cluster==cluster)]
            ax.scatter(d.longitude,d.latitude,s=18,c=color,edgecolor='white',linewidth=.2,zorder=4,label=f'C{cluster}: n={len(d)}')
        ax.legend(loc='lower left',fontsize=5.5,framealpha=.85)
    save(fig,out,SPECS[6][0],records)
    representative=load('representative_station_selection');profiles=load('qr_all_quantiles_long')
    fig,axes=plt.subplots(2,2,figsize=(7.1,6.4),layout='constrained')
    for i,(idx,ax) in enumerate(zip(INDEXES,axes.flat)):
        for _,r in representative.loc[representative.index_name==idx].iterrows():
            d=profiles.loc[(profiles.index_name==idx)&(profiles.station_id==r.station_id)].sort_values('tau')
            ax.plot(d.tau,d.slope,lw=1.1,color=colors[int(r.cluster)-1],label=f'C{int(r.cluster)}: {r.station_name}')
        ax.axhline(0,color='gray',lw=.5);ax.set(title=f'({chr(97+i)}) {NAMES[i]}',xlabel='Quantile',ylabel='Days per decade',xticks=[.1,.3,.5,.7,.9]);ax.legend(loc='upper center',bbox_to_anchor=(.5,-.23),ncol=2,fontsize=5.3,frameon=False)
    save(fig,out,SPECS[7][0],records)
    anomaly=load('regional_temperature_anomaly');response=load('warming_link_network_quantile_response')
    fig=plt.figure(figsize=(7.1,6.5),layout='constrained');gs=fig.add_gridspec(3,2,height_ratios=[.8,1,1]);ax=fig.add_subplot(gs[0,:])
    ax.plot(anomaly.year,anomaly.regional_temperature_anomaly_c,color='#A64B3D',lw=1);ax.axhline(0,color='gray',lw=.5);ax.set(title='(a) Network temperature anomaly',ylabel='°C relative to baseline')
    for i,idx in enumerate(INDEXES):
        ax=fig.add_subplot(gs[1+i//2,i%2]);d=response.loc[(response.index_name==idx)&response.tau.notna()].sort_values('tau')
        for _,r in d.iterrows():ax.plot([r.tau,r.tau],[r.ci_low_per_c,r.ci_high_per_c],color=COLORS[i],lw=1.2)
        ax.plot(d.tau,d.slope_per_c,'o-',color=COLORS[i],ms=3,lw=1);ax.axhline(0,color='gray',lw=.5);ax.set(title=f'({chr(98+i)}) {NAMES[i]}',xlabel='Quantile',ylabel='Days per °C',xticks=[.1,.5,.9])
    save(fig,out,SPECS[8][0],records)
    rarity=load('compound_dry_hot_yearly_extent','compound_dry_hot/tables');fig,axes=plt.subplots(2,1,figsize=(7.1,4.8),layout='constrained',sharex=True)
    for i,(definition,ax) in enumerate(zip(['annual_tmean','warm_season_tmax'],axes)):
        for threshold,color in zip([5,10,20],colors):
            d=rarity.loc[(rarity.definition==definition)&(rarity.return_period_threshold_years==threshold)]
            ax.plot(d.year,d.affected_fraction_pct,color=color,lw=1,marker='o',ms=2,label=f'Inverse joint probability ≥ {threshold}')
        ax.set(title=f'({chr(97+i)}) '+['Annual','June–September'][i],ylabel='Available stations (%)',ylim=(0,105));ax.legend(fontsize=6,loc='upper left')
    axes[-1].set_xlabel('Year');save(fig,out,SPECS[9][0],records)
    driver=load('driver_analysis_summary');f=driver.loc[driver.metric=='Delta1'].pivot(index='index_name',columns='predictor',values='std_beta').loc[INDEXES,['latitude','longitude','elevation']];f.index=NAMES;f.columns=['Latitude','Longitude','Elevation']
    fig,ax=plt.subplots(figsize=(5.7,2.8),layout='constrained');im=heat(ax,f,'Geographical associations with thermal asymmetry',limit=max(.5,float(abs(f).max().max())));fig.colorbar(im,ax=ax,label='Standardized coefficient');save(fig,out,SPECS[10][0],records)
    pd.DataFrame(records).to_csv(out/'supplementary_figure_manifest.csv',index=False)
    inputs.update([root/'data/stationsInfo.csv',root/'data/Iran_Sea_Ne.geojson',Path(__file__)])
    (out/'supplementary_figure_sources.json').write_text(json.dumps({p.relative_to(root).as_posix():file_hash(p) for p in sorted(inputs)},indent=2),encoding='utf-8')
    return records

def markdown_table(frame):
    frame=frame.copy()
    labels={'definition':'Period','component':'Component','estimate_pp':'Estimate','ci_low_pp':'95% lower','ci_high_pp':'95% upper',
            'family_ci_low_pp':'Family lower','family_ci_high_pp':'Family upper','scenario':'Scenario','n_stations':'N',
            'screen':'Screen','rows_flagged':'Flagged rows','stations_affected':'Stations','action_in_sensitivity':'Sensitivity action',
            'index_name':'Index','tau':'Quantile','moran_i':'Moran’s I','p_perm_two_sided':'Permutation p','moran_p_perm':'Permutation p',
            'local_FDR_retained':'FDR retained','available_analytic_tests':'Available tests','n_clusters':'Clusters',
            'observed_mean_within_cluster_distance_km':'Observed distance (km)','permuted_mean_distance_km':'Permuted distance (km)',
            'p_perm_more_compact':'Compactness p','climate_regime':'Climate regime','n_zero_dry_threshold':'N₀',
            'average_cityblock':'Average / cityblock','complete_euclidean':'Complete / Euclidean','ward_euclidean':'Ward / Euclidean',
            'kmeans':'k-means','expanded_features':'Expanded features','0.10':'q10','0.50':'q50','0.90':'q90','Delta1':'Δ'}
    values=dict(zip(INDEXES,NAMES))
    values.update({'annual':'Annual','warm_season':'June–September','joint_change':'Joint frequency',
                   'dry_marginal':'Dry-frequency term','hot_marginal':'Hot-frequency term','excess_joint':'Excess-joint term',
                   'BWh_hot_desert':'BWh: hot desert','BWk_cold_desert':'BWk: cold desert','BSh_hot_steppe':'BSh: hot steppe',
                   'BSk_cold_steppe':'BSk: cold steppe','C_temperate':'C: temperate','Dsa_cold_dry_summer':'Dsa: cold, dry summer'})
    for c in frame:
        frame[c]=frame[c].map(lambda v:values.get(v,v.replace('_',' ') if isinstance(v,str) else v))
        if c=='local_FDR_retained':frame[c]=frame[c].map(lambda v:int(v) if isinstance(v,(float,int)) else v)
    frame=frame.rename(columns=labels)
    lines=['| '+' | '.join(str(c) for c in frame.columns)+' |','| '+' | '.join('---' for _ in frame.columns)+' |']
    for row in frame.itertuples(index=False,name=None):
        lines.append('| '+' | '.join(f'{v:.3f}' if isinstance(v,float) else str(v) for v in row)+' |')
    return '\n'.join(lines)

def build_supplementary_document(root,out,primary,scenarios):
    load=lambda name:pd.read_csv(root/'outputs/tables'/f'{name}.csv')
    text='''# Supplementary material: thermal extremes and explicit dry–hot concurrence

This supplement accompanies [the revised manuscript](Manuscript_Q1_2026.md). Its figures were selected after an inventory of every existing output and regenerated from machine-readable tables. The [output audit](Output_Audit_2026.md) records verification scope, corrections and removals. The [data catalog](Supplementary_Data_Catalog.md) links retained numerical evidence and explains its role. Original material removed from the active output tree is preserved in a recovery archive identified in the audit.

## S1. Additional methods and interpretation

Thermal slopes concern annual counts and use the historical day-of-year percentile construction described in the main paper. Annual station indices were independently rebuilt from daily observations. All saved 200-replicate thermal bootstrap summaries were recalculated from their individual draws. Alternative 400-replicate and maximum-entropy ensembles are represented by saved station summaries; their summary aggregation was checked, but the alternative draws were not archived and those ensembles were not rerun.

Homogeneity diagnostics use annual mean temperatures, with Pettitt, SNHT and Buishand tests after linear detrending. Their nominal flags do not identify or correct specific artificial breaks. The exclusion sensitivity describes a changed observing network, not a homogenized dataset. Calendar coverage in Figure S1 uses the new screened aggregates; historical available-row completeness is retained in the data catalog and is not substituted for calendar-day coverage.

Spatial diagnostics use five-nearest-neighbor weights and 499 label permutations. The reported Moran probabilities are exploratory across fields. Analytic station probabilities are adjusted within each index–quantile family using Benjamini–Hochberg at 0.05. Tail analytic intervals are unavailable in the archived estimates; NA is retained instead of converting missing tests to zero. The local retained counts do not establish field significance or repair serial-dependence limitations. The precision ratio in Figure S4 is an absolute bootstrap mean divided by bootstrap standard deviation; it must not be interpreted as a time of emergence.

Clustering uses standardized station quantile slopes, average linkage, Euclidean distance and four requested groups per index. Features with absolute correlation at least 0.95 are screened in the configured order. Alternative linkage, distance, k-means and expanded uncertainty-feature choices test sensitivity; adjusted Rand index is used because raw label agreement depends on arbitrary cluster numbering. Cluster composites describe fitted groups, including any singleton clusters, rather than independently validated climate regions. Representatives are observed stations nearest fitted cluster centroids. They illustrate heterogeneity without adding independent inferential evidence.

Geographical regressions and regressions on the internally derived temperature anomaly remain descriptive. The latter share observations and trends with their outcomes. Historical composite “fingerprint” scores combine overlapping diagnostics and their circular-shift null shifts indices independently. They are excluded from the active supplementary evidence because neither the aggregate score nor that null supplies an independent attribution test. Their original files are recorded in the recovery manifest.

Historical empirical joint-rarity analyses use an observation-specific inverse joint probability and available-row completeness. Their yearly networks vary; no stable design return period, fixed marginal AND event, physical affected area or causal “driver” is inferred. Figure S10 preserves this analysis only to demonstrate the effect of event definition.

## S2. Numerical evidence

### Table S1. Primary compound partition with pointwise and family intervals

All entries are percentage points. Pointwise intervals use the 2.5th and 97.5th percentiles of 4,999 synchronized circular-block replicates. Family intervals use the 0.3125th and 99.6875th percentiles for eight primary quantities. Thresholds are re-estimated in each replicate. Nominal multiplicity adjustment does not guarantee exact coverage in short discrete series.

'''
    text+=markdown_table(primary[['definition','component','estimate_pp','ci_low_pp','ci_high_pp','family_ci_low_pp','family_ci_high_pp']])
    text+='\n\n### Table S2. Joint-frequency sensitivity\n\nChanging coverage, thresholds or station selection changes the estimand. Scenario intervals are exploratory rather than independent confirmations.\n\n'
    joint=scenarios.loc[scenarios.component=='joint_change',['definition','scenario','n_stations','estimate_pp','ci_low_pp','ci_high_pp']]
    text+=markdown_table(joint)
    text+='\n\n### Table S3. Daily internal-consistency screening\n\nScreening actions apply to the new compound extension and the explicit historical sensitivity, not retroactively to all historical thermal results.\n\n'
    text+=markdown_table(load('temperature_internal_consistency_screening')[['screen','rows_flagged','stations_affected','action_in_sensitivity']])
    q=load('qr_focus_slopes_and_bootstrap_summary');rows=[]
    for idx,name in zip(INDEXES,NAMES):
        d=q.loc[q.index_name==idx];r={'Index':name,'N':len(d)}
        for metric in ['0.10','0.50','0.90','Delta1']:
            r[metric]=int(((d['boot_ci_low_'+metric]>0)|(d['boot_ci_high_'+metric]<0)).sum())
        rows.append(r)
    text+='\n\n### Table S4. Station bootstrap intervals excluding zero\n\nCounts use pointwise 95% percentile intervals, without multiplicity adjustment. They are not equivalent to the analytic FDR counts in Table S5.\n\n'+markdown_table(pd.DataFrame(rows))
    moran=load('spatial_autocorrelation_moran');fdr=load('station_significance_fdr')
    counts=fdr.groupby(['index_name','tau']).agg(local_FDR_retained=('fdr_reject',lambda s:s.sum(min_count=1)),available_analytic_tests=('analytic_p','count')).reset_index()
    text+='\n\n### Table S5. Spatial diagnostics and retained local tests\n\nMoran probabilities are nominal permutation results; dependence limitations are described in Section S1. Tail analytic intervals are unavailable: NA denotes no estimable test, not zero significant stations. Bootstrap tail uncertainty is provided separately in Table S4.\n\n'+markdown_table(moran.merge(counts,on=['index_name','tau']).fillna('NA'))
    alt=load('alternative_clustering_sensitivity_summary').pivot(index='index_name',columns='method_label',values='adjusted_rand_index')
    alt=alt.join(load('cluster_robustness_summary').set_index('index_name').adjusted_rand_index.rename('expanded_features'))
    alt.columns=[c.replace('hierarchical_','') for c in alt.columns]
    text+='\n\n### Table S6. Clustering sensitivity: adjusted Rand index\n\nAgreement with the baseline partition. Values near one indicate similar assignments; low agreement cautions against treating the partition as uniquely determined.\n\n'+markdown_table(alt.reset_index())
    spatial=load('regional_cluster_spatial_validation')
    text+='\n\n### Table S7. Exploratory spatial compactness of clusters\n\nThe statistic is mean within-cluster geographical distance, compared with 499 permutations of cluster labels. These nominal probabilities do not validate the fitted groups as physical climate regions.\n\n'+markdown_table(spatial[['index_name','n_clusters','observed_mean_within_cluster_distance_km','permuted_mean_distance_km','p_perm_more_compact']])
    regimes=pd.read_csv(out/'tables/climate_regime_compound_partition.csv')
    text+='\n\n### Table S8. Compound-frequency change in each climate regime\n\nEqual-station within-group means and exploratory 95% synchronized-block intervals; annual and summer sample sizes differ. The complete four-component results are in the linked data catalog.\n\n'+markdown_table(regimes.loc[regimes.component=='joint_change',['definition','climate_regime','n_stations','n_zero_dry_threshold','estimate_pp','ci_low_pp','ci_high_pp']])
    text+='\n\n## S3. Supplementary figures\n\n'
    for i,(stem,title,caption) in enumerate(SPECS,1):
        text+=f'### Figure S{i}. {title}\n\n![{title}](../outputs/publication_v2/figures/{stem}.png)\n\n*{caption}*\n\n'
    text+='''## S4. Reproducibility and complete machine-readable evidence

All retained research tables are indexed in [Supplementary_Data_Catalog.md](Supplementary_Data_Catalog.md), including station-level results, complete bootstrap draws, reference metadata and validation records. Main and supplementary graphics are available as vector PDF, editable SVG, 350-dpi PNG and 600-dpi TIFF. The [supplementary figure atlas](../outputs/publication_v2/Supplementary_Figure_Atlas.pdf) contains Figures S1–S11; the [main atlas](../outputs/publication_v2/Figure_Atlas.pdf) contains Figures 1–10. Map boundaries supply geographical context only; their external source and redistribution license still require author confirmation.

Build the extension with `python run_publication.py`, validate raw-derived thermal and compound quantities with `python validate_publication.py`, and build the curated supplementary graphics with `python build_supplementary.py`. Rebuild text and atlases with `python build_publication_docs.py`. The wider historical numerical audit is reproduced with `python audit_output_data.py recompute`; it writes to an isolated work directory and never overwrites research outputs. Cleanup follows a file-specific manifest, with a verified recovery archive before deletion.
'''
    return text
