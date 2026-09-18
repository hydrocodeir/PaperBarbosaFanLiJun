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
('figS02_thermal_robustness','Sensitivity of thermal estimates','(a) Mean station interval-width ratio for 400 versus 200 moving-block replicates. (b) Mean absolute difference between saved moving-block and maximum-entropy bootstrap means. (c) Change in network slopes after internal-consistency screening. (d) Change after excluding homogeneity-flagged stations. Units in (b–d) are days per decade. Panels (c–d) use the historical available network, not the primary fixed 108-station network. The two historical alternative bootstrap ensembles were checked through their saved station summaries, not regenerated.'),
('figS03_station_quantile_maps','Station quantile slopes','Station estimates at the 0.10, 0.50 and 0.90 quantiles of annual thermal counts, on one common symmetric scale. No interpolation, area weighting or local significance symbols are used. Units are days per decade.'),
('figS04_median_precision','Bootstrap precision of median trends','Absolute bootstrap mean divided by bootstrap standard deviation at the median. This descriptive precision ratio is dimensionless; it is neither an emergence date nor a calibrated detection probability. No threshold-based significance labels are shown.'),
('figS05_spatial_diagnostics','Spatial dependence and multiplicity','(a) Moran’s I with five-nearest-neighbor weights. (b) Numbers of locally retained tests after Benjamini–Hochberg adjustment separately within each index–quantile family. NA means analytic tail intervals and tests are unavailable, not that zero stations are significant. Historical analytic probabilities are approximate and do not adjust for serial dependence. These counts do not constitute a field-significance test; permutation probabilities are tabulated separately.'),
('figS06_cluster_composites','Exploratory cluster composites','Median station slopes and upper-minus-lower contrasts within four fitted clusters per index. Cluster labels are specific to each index and have no shared climatic ordering. Cluster separation is partly induced by the fitted features and is not independent validation; single-station groups remain visible through their sample sizes.'),
('figS07_cluster_maps','Exploratory cluster membership','Station assignments from standardized, screened quantile-slope features and average-linkage Euclidean clustering. Colors distinguish categories only within each panel; clusters are not fixed climate regions. Alternative-method agreement and within-cluster sample sizes are reported in the supplementary tables.'),
('figS08_representative_profiles','Representative station profiles','One observed station per index-specific cluster, selected by the archived nearest-centroid rule. Lines show the full 0.10–0.90 quantile profile. These selected illustrations are not additional independent discoveries; all station coefficients and bootstrap intervals are available in the data catalog.'),
('figS09_internal_warming_association','Association with the network temperature anomaly','(a) Internally derived annual network temperature anomaly relative to 1991–2007. (b–e) Fixed-baseline network thermal-count response coefficients at three quantiles, with available median pointwise analytic intervals, in days per degree Celsius. Tail analytic intervals are unavailable and must not be inferred from the unadorned points. The predictor and outcomes share observations and temporal trends; the associations do not attribute change to external forcing.'),
('figS10_historical_joint_rarity','Historical joint-rarity sensitivity','Percentage of the available station network exceeding three empirical inverse-joint-probability cutoffs. Labels denote within-record rarity classes, not stable design return periods. This historical calculation uses unscreened daily data, available-row coverage denominators and a varying station network. It is retained as a definition sensitivity and is distinct from the main paper’s fixed marginal AND event.'),
('figS11_geographical_associations','Geographical associations of thermal asymmetry','Standardized multiple-regression coefficients relating the upper-minus-lower slope contrast to latitude, longitude and elevation. These are descriptive associations; no causal driver attribution or spatially adjusted significance is implied.'),
('figS12_thermal_network_composition','Thermal network composition and coverage sensitivity','(a) Number of stations with valid annual counts in the available daytime and nighttime networks; the dashed line marks the fixed common set of 108. (b) Upper-minus-lower quantile slope point estimates for four network definitions. Index-specific fixed sets each contain 109 stations; their common intersection contains 108. The 365-day-equivalent diagnostic scales each common station-year count by 365/valid days and assumes representative missing days. These are sensitivity estimates, not additional significance tests. Units in (b) are days per decade.'),
('figS13_thermal_contrast_intervals','Paired thermal contrast uncertainty and block sensitivity','(a) Four upper-minus-lower quantile slope contrasts on the common 108-station network. (b) Paired daytime-minus-nighttime contrasts in that signed asymmetry, separately for warm and cool indices. Points and segments show estimates and pointwise 95% percentile intervals from 4,999 synchronized circular year-pairs block replicates for each block length. Original year covariates travel with responses; the same sampled years are used across indices. Four-year blocks define the primary analysis. Intervals are conditional on fixed daily thresholds and the observed station set, and all shown contrast intervals include zero. The wider nominal six-contrast family intervals for the primary analysis are in Table S10.'),
('figS14_index_definition_asymmetry','Asymmetry sensitivity to index construction','Upper-minus-lower quantile slope contrasts on the same 108 stations under six index constructions. Full denotes the 1991–2024 reference; early denotes 1991–2007. Numbers 11 and 5 are total calendar-day window widths. T7 and T8 denote Hyndman–Fan linear and median-unbiased sample quantiles. Corrected early indices average target-year event counts after excluding that year and duplicating each other baseline year in turn (16 replacements); out-of-base thresholds remain fixed. Segments are exploratory pointwise 95% intervals from the same 4,999 four-year pairs-block draws applied to each constructed annual series. Index construction is held fixed inside this uncertainty resampling. These correlated sensitivities are not independent confirmations or simultaneous tests.'),
('figS15_zero_threshold_tie_effects','Paired effects of dry and hot threshold ties','Summer joint-frequency-change differences on identical station sets and synchronized year samples. Panel (a) changes P < q25 to P ≤ q25 while retaining T > q75; panel (b) then changes T > q75 to T ≥ q75. Points and exploratory 95% percentile intervals use 4,999 synchronized four-year block replicates with threshold refitting. These are paired definition effects, not changes in the observed climate or formal between-regime tests. In zero-cutoff stations, including dry ties classifies zero-precipitation seasons as dry rather than identifying a departure below the cutoff.'),
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
    coverage=load('thermal_network_coverage','publication_v2/tables')
    comparison=load('thermal_network_comparison','publication_v2/tables')
    fig,axes=plt.subplots(1,2,figsize=(7.1,3.6),layout='constrained')
    for idx,label,color,ls in [('warm_days','Daytime','#B74335','-'),('warm_nights','Nighttime','#465AA3',':')]:
        d=coverage.loc[coverage.index_name==idx]
        axes[0].plot(d.year,d.n_available,label=label,color=color,lw=1.2,ls=ls)
    axes[0].axhline(108,color='#333333',ls='--',lw=1,label='Common fixed: 108')
    axes[0].set(title='(a) Annual station membership',xlabel='Year',ylabel='Valid stations',ylim=(106,126),yticks=[108,112,116,120,124])
    axes[0].legend(frameon=False,fontsize=6.5,loc='center left')
    for j,(network,label,marker,color) in enumerate([('available','Available','o','#555D63'),('index_fixed','Index fixed (109)','s','#A97A37'),('common_fixed','Common fixed (108)','D','#287E9B'),('common_365_equivalent','Common: 365-day equivalent','^','#7666A8')]):
        d=comparison.loc[comparison.network==network].set_index('index_name').loc[INDEXES]
        axes[1].scatter(d.Delta1,np.arange(4)+(j-1.5)*.14,label=label,marker=marker,color=color,s=22)
    axes[1].axvline(0,color='gray',lw=.6)
    axes[1].set(title='(b) Asymmetry point estimates',yticks=range(4),yticklabels=NAMES,xlabel='Δ₁ (days per decade)',ylim=(-.6,4.8))
    axes[1].invert_yaxis();axes[1].legend(frameon=False,fontsize=5.8,loc='lower right')
    save(fig,out,SPECS[11][0],records)
    intervals=load('thermal_network_intervals','publication_v2/tables')
    contrasts=load('thermal_network_contrasts','publication_v2/tables')
    fig,axes=plt.subplots(2,1,figsize=(7.1,5.6),layout='constrained',sharex=True,gridspec_kw={'height_ratios':[1.7,1]})
    for j,(block,color,marker) in enumerate([(2,'#88969E','s'),(4,'#287E9B','o'),(6,'#B17B3B','^')]):
        d=intervals.loc[(intervals.network=='common_fixed')&(intervals.metric=='Delta1')&(intervals.block_length==block)].set_index('index_name').loc[INDEXES]
        paired=contrasts.loc[(contrasts.metric=='Delta1')&(contrasts.block_length==block)].set_index('contrast').loc[['warm_day_minus_night','cool_day_minus_night']]
        for ax,frame in [(axes[0],d),(axes[1],paired)]:
            y=np.arange(len(frame))+(j-1)*.2
            ax.hlines(y,frame.ci_low,frame.ci_high,color=color,lw=1.5 if block==4 else .9)
            ax.scatter(frame.estimate,y,s=25,color=color,marker=marker,label=f'{block}-year blocks')
    for ax in axes:ax.axvline(0,color='#555D63',lw=.8,ls='--');ax.invert_yaxis()
    axes[0].set(title='(a) Upper-minus-lower slope contrasts',yticks=range(4),yticklabels=NAMES)
    axes[0].legend(frameon=False,ncol=3,loc='lower left',fontsize=7)
    axes[0].set_ylim(4,-.6)
    axes[1].set(title='(b) Daytime-minus-nighttime asymmetry',yticks=range(2),yticklabels=['Warm indices','Cool indices'],xlabel='Paired contrast (days per decade)')
    save(fig,out,SPECS[12][0],records)
    from .definition_figures import plot_definition_supplement
    plot_definition_supplement(out,records,save)
    inputs.update([out/'tables/index_definition_trends.csv',out/'tables/zero_threshold_paired_effects.csv',root/'src/paper_pipeline/definition_figures.py'])
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

The primary thermal network analysis is distinct from those historical station summaries. It uses the common 108 stations meeting annual-index coverage in all 34 years for all four indices. Each index-specific complete set has 109 stations. Available, index-specific fixed, common fixed and 365-day-equivalent series are archived separately. The latter scales observed counts by 365/valid days and assumes representative missing days; it is not an imputed record.

Network intervals use 4,999 synchronized circular moving-block pairs draws of original calendar year and response field. Original time covariates are retained rather than reassigning sampled responses to a new time axis. All indices and network definitions share the saved sampled year positions. Four-year blocks are primary; two- and six-year blocks are sensitivities. Daily thresholds and station selection are held fixed. OLS, three focal quantile slopes, within-index upper-minus-lower contrasts and paired day-minus-night contrasts are saved for every replicate. The full primary quantile grid has pointwise 95% bands. Nominal family intervals in Table S10 use a six-contrast Bonferroni adjustment, separate from the compound family in Table S1. No exact short-record coverage is claimed. Independent linear programming checks the weighted check-loss minimizer and saved replicate coefficients; membership, means, contrasts and interval arithmetic are also checked by `validate_thermal_network.py`.

Index-construction sensitivities use the same 108-station daily records and masks. Each early-reference target year is removed and replaced by each of the other 16 years in turn; events are counted separately and counts averaged. The late-period cutoffs use the unchanged original reference. This deterministic correction is distinct from the conditional uncertainty bootstrap applied to completed annual indices. The primary full-record construction remains the estimand of Figure 2/Table 1; corrected early-reference differences supersede uncorrected headline period changes in Figure 4. Definitions and quantitative comparisons are in Tables S12–S13.

Structural-zero analyses keep observed-zero/positive strata fixed while varying strict, dry-inclusive-only and both-inclusive rules on identical station sets. For observed strict cutoffs, every zero-cutoff station has zero dry/joint frequency and zero partition components. Full-network point estimates are therefore weighted positive-subset means, exactly. Fixed-threshold draws retain this identity. Refitted cutoffs can become positive at originally zero-cutoff stations, giving nonzero bootstrap contributions; interval endpoints then cannot be rescaled mechanically. Paired rule-effect intervals use the same year draws and retain cross-variable dependence. They remain exploratory: for example, the summer additional hot-tie estimate is −2.28 percentage points while its percentile interval is [−8.51, −2.34], showing bootstrap centering displacement in this discrete short-baseline statistic. This behavior is reported rather than treating every percentile interval as calibrated confirmation.

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
    comparison=pd.read_csv(out/'tables/thermal_network_comparison.csv')
    text+='\n\n### Table S9. Thermal trend sensitivity to network definition\n\nAll slopes and contrasts are in days per decade. Available networks contain 115–124 daytime or 114–124 nighttime stations per year; index-specific fixed sets contain 109 and the common set contains 108. Common 365-equivalent denotes the coverage-scaled diagnostic described in Section S1. These point estimates do not replace the primary intervals in Table 1.\n\n'+markdown_table(comparison[['index_name','network','OLS','q10','q50','q90','Delta1']])
    intervals=pd.read_csv(out/'tables/thermal_network_intervals.csv')
    primary_block=json.loads((out/'thermal_network_metadata.json').read_text())['settings']['primary_block_length']
    d=intervals.loc[(intervals.network=='common_fixed')&(intervals.metric=='Delta1')&(intervals.block_length==primary_block)].copy()
    d['Contrast']=d.index_name.map(dict(zip(INDEXES,NAMES)))+' Δ₁'
    paired=pd.read_csv(out/'tables/thermal_network_contrasts.csv')
    paired=paired.loc[(paired.block_length==primary_block)&paired.contrast.isin(['warm_day_minus_night','cool_day_minus_night'])].copy()
    paired['Contrast']=paired.contrast.map({'warm_day_minus_night':'Warm Δ₁: day minus night','cool_day_minus_night':'Cool Δ₁: day minus night'})
    columns=['Contrast','estimate','ci_low','ci_high','family_ci_low','family_ci_high']
    frame=pd.concat([d[columns],paired[columns]],ignore_index=True)
    frame.columns=['Contrast','Estimate','95% lower','95% upper','Family lower','Family upper']
    text+='\n\n### Table S10. Primary paired thermal contrasts with uncertainty\n\nDays per decade on the common 108-station network. Pointwise intervals use 2.5th/97.5th percentiles and nominal six-contrast family intervals use 0.4167th/99.5833rd percentiles of 4,999 synchronized four-year block replicates. Day-minus-night comparisons use signed Δ₁, including for cool indices. All six intervals include zero; this does not establish equal slopes or equal asymmetry. Family columns stored for other networks or metrics are not part of this primary inferential family.\n\n'+markdown_table(frame)
    influence=pd.read_csv(out/'tables/thermal_network_leave_one_year_out.csv')
    frame=influence.groupby('index_name').agg(q90_min=('q90','min'),q90_max=('q90','max'),Delta1_min=('Delta1','min'),Delta1_max=('Delta1','max')).reindex(INDEXES).reset_index()
    text+='\n\n### Table S11. Leave-one-year-out thermal point estimates\n\nMinimum and maximum estimates across 34 fits, each omitting one original year from the common-network series and retaining the other calendar-year covariates. Units are days per decade. These ranges describe influence and are not confidence intervals. Omitted-year coefficients for all focal quantiles are archived in the data catalog.\n\n'+markdown_table(frame)
    from .definition_figures import SCENARIOS, LABELS_INDEX
    from .publication_regimes import REGIMES
    settings=json.loads((out/'index_definition_metadata.json').read_text())['settings']
    design=[]
    for case in settings['scenarios']:
        design.append({'Scenario':case['name'],'Reference':f"1991–{case['reference_end']}",'Window (days)':case['window'],'Quantile type':'T7' if case['method']=='linear' else 'T8','In-base correction':'16 donor replacements' if case['correction'] else 'None'})
    text+='\n\n### Table S12. Controlled thermal-index constructions\n\nAll six scenarios use the same 108 stations, strict inequalities, no leap days and at least 80% valid annual days. They report observed annual counts and, separately, percentages of valid days. The 5-day/T8 corrected sensitivity follows the percentile and in-base replacement conventions used in climdex, but uses a 17-year baseline and the study’s annual coverage rule; it is not a fully standard ETCCDI implementation with standard baseline and monthly completeness requirements. Historical T7 uses linear interpolation between sample order statistics; T8 uses the median-unbiased convention. No sparse reference window triggered the minimum-15-observation guard.\n\n'+markdown_table(pd.DataFrame(design))
    period=pd.read_csv(out/'tables/index_definition_period_summary.csv')
    trend=pd.read_csv(out/'tables/index_definition_trends.csv')
    rows=[]
    for case in SCENARIOS:
        for idx in INDEXES:
            r=period.loc[(period.scenario==case)&(period.index_name==idx)].iloc[0]
            d=trend.loc[(trend.scenario==case)&(trend.index_name==idx)].set_index('metric')
            v=d.loc['Delta1']
            rows.append([case,idx,r.change_days,r.change_rate_pp,d.loc['q10','estimate'],d.loc['q90','estimate'],f'{v.estimate:.2f} [{v.ci_low:.2f}, {v.ci_high:.2f}]'])
    text+='\n\n### Table S13. Thermal period changes and trend sensitivity\n\nPeriod differences are 2008–2024 minus 1991–2007 means, averaging stations equally on the common 108-station network. Count changes are days per year; rate changes are percentage points of valid observed days. Slopes and Δ₁ intervals are days per decade. Corrected baseline counts may be fractional because event counts, not thresholds, are averaged over donor replacements. Intervals are conditional sensitivity intervals, not independent tests.\n\n'+markdown_table(pd.DataFrame(rows,columns=['Scenario','index_name','Count change','Valid-day rate change','q10 slope','q90 slope','Δ₁ [95% interval]']))
    zero=pd.read_csv(out/'tables/zero_threshold_summary.csv')
    z=zero.loc[(zero.threshold_mode=='refitted')&(zero.component=='joint_change')]
    rows=[]
    for definition in ['annual','warm_season']:
        for group in ['All',*REGIMES]:
            d=z.loc[(z.definition==definition)&(z.climate_regime==group)]
            allrow=d.loc[(d.rule=='strict')&(d.stratum=='all')].iloc[0]
            pos=d.loc[(d.rule=='strict')&(d.stratum=='positive')].iloc[0]
            dry=d.loc[(d.rule=='dry_inclusive_only')&(d.stratum=='all')].iloc[0]
            both=d.loc[(d.rule=='both_inclusive')&(d.stratum=='all')].iloc[0]
            fmt=lambda r:f'{r.estimate_pp:.2f} [{r.ci_low_pp:.2f}, {r.ci_high_pp:.2f}]'
            rows.append([definition,group,int(allrow.n_stations),int(allrow.n_zero),int(pos.n_stations),fmt(allrow),fmt(pos),fmt(dry),fmt(both)])
    text+='\n\n### Table S14. Structural-zero dilution and tie definitions within climate regimes\n\nJoint-frequency changes in percentage points with exploratory 95% refitted-threshold intervals. N+ is the fixed observed-positive-cutoff subset. All strict point estimates equal N+/N times the positive-subset point estimate: this is an exact reweighting identity, not independent robustness evidence. Strict and inclusive all-station rules retain identical stations. Dry-inclusive changes only precipitation equality; both-inclusive also changes temperature equality. Group intervals do not test between-group contrasts. The three-station summer BSh positive subset is especially imprecise. Fixed-threshold intervals, paired rule effects, observed-zero strata and bootstrap cutoff changes are archived separately.\n\n'+markdown_table(pd.DataFrame(rows,columns=['definition','climate_regime','N','N₀','N+','Strict all','Strict positive subset','Dry-inclusive all','Both-inclusive all']))
    text+='\n\n## S3. Supplementary figures\n\n'
    for i,(stem,title,caption) in enumerate(SPECS,1):
        text+=f'### Figure S{i}. {title}\n\n![{title}](../outputs/publication_v2/figures/{stem}.png)\n\n*{caption}*\n\n'
    text+='''## S4. Reproducibility and complete machine-readable evidence

All retained research tables are indexed in [Supplementary_Data_Catalog.md](Supplementary_Data_Catalog.md), including station-level results, complete bootstrap draws, reference metadata and validation records. Main and supplementary graphics are available as vector PDF, editable SVG, 350-dpi PNG and 600-dpi TIFF. The [supplementary figure atlas](../outputs/publication_v2/Supplementary_Figure_Atlas.pdf) contains Figures S1–S15; the [main atlas](../outputs/publication_v2/Figure_Atlas.pdf) contains Figures 1–10. Map boundaries supply geographical context only; their external source and redistribution license still require author confirmation.

Build the extension with `python run_publication.py`, validate raw-derived thermal and compound quantities with `python validate_publication.py`, and build the curated supplementary graphics with `python build_supplementary.py`. Reproduce the thermal network extension alone with `python run_thermal_network.py` and independently verify it with `python validate_thermal_network.py`. Run index-construction and structural-zero sensitivities with `python run_index_definition.py` and `python run_zero_threshold.py`; check both using `python validate_index_zero.py`. Rebuild text and atlases with `python build_publication_docs.py`. The wider historical numerical audit is reproduced with `python audit_output_data.py recompute`; it writes to an isolated work directory and never overwrites research outputs. Cleanup follows a file-specific manifest, with a verified recovery archive before deletion.
'''
    return text
