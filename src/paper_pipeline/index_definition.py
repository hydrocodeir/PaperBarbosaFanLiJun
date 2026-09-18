"""Controlled percentile-index sensitivities with exhaustive in-base replacement.

For each baseline target year, remove it and duplicate each other year in turn.
Count target-year events separately for each replacement, then average counts.
This is an index construction correction, not an uncertainty bootstrap.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from .compound_partition import file_hash
from .thermal_network import PairQuantileSolver, weighted_ols, weights_from_indices, INDEXES


def window_samples(base, window):
    days=(np.arange(365)[:,None]+np.arange(-(window//2),window//2+1))%365
    samples=base[:,days].transpose(1,0,2).reshape(365,-1)
    year=np.broadcast_to(np.repeat(np.arange(len(base)),window),samples.shape)
    order=np.argsort(samples,axis=1)
    return np.take_along_axis(samples,order,axis=1),np.take_along_axis(year,order,axis=1)


def repeated_quantiles(sorted_values, year_labels, multiplicities, method):
    """Quantiles of an explicitly repeated sample, without materializing repeats."""
    valid=np.isfinite(sorted_values)
    w=multiplicities[:,year_labels]*valid[None,:,:]
    cumulative=np.cumsum(w,axis=2);n=cumulative[:,:,-1]
    if np.any(n==0):raise ValueError('Empty reference window')
    answer=[]
    for q in [.1,.9]:
        h=(n-1)*q if method=='linear' else n*q+(1/3+q*(1-1/3-1/3))-1
        h=np.clip(h,0,n-1);low=np.floor(h).astype(int);high=np.ceil(h).astype(int)
        il=(cumulative>low[:,:,None]).argmax(axis=2)
        ih=(cumulative>high[:,:,None]).argmax(axis=2)
        row=np.arange(365)[None,:]
        vl=sorted_values[row,il];vh=sorted_values[row,ih]
        fraction=h-low
        # Match NumPy's symmetric interpolation exactly: last-bit differences
        # can change strict comparisons for observations tied to a cutoff.
        answer.append(np.where(fraction>=.5,vh-(vh-vl)*(1-fraction),vl+fraction*(vh-vl)))
    return np.stack(answer,axis=-1)


def construct_counts(values, nbase, window, method, corrected, minimum=15):
    base=values[:nbase]
    sorted_values,labels=window_samples(base,window)
    # Do not silently change the historical fallback on sparsely sampled windows.
    # All real windows must meet the explicit minimum for this fixed network audit.
    if np.any(np.isfinite(sorted_values).sum(axis=1)<minimum):
        raise ValueError('Sparse window needs an explicitly audited fallback')
    cut=repeated_quantiles(sorted_values,labels,np.ones((1,nbase),int),method)[0]
    valid=np.isfinite(values)
    counts=np.column_stack([((values<cut[:,0])&valid).sum(axis=1),((values>cut[:,1])&valid).sum(axis=1)]).astype(float)
    if corrected:
        for target in range(nbase):
            donors=np.delete(np.arange(nbase),target)
            mult=np.ones((nbase-1,nbase),int);mult[:,target]=0;mult[np.arange(nbase-1),donors]=2
            if np.any((mult[:,labels]*np.isfinite(sorted_values)[None,:,:]).sum(axis=2)<minimum):
                raise ValueError('Sparse replacement window requires explicit handling')
            cuts=repeated_quantiles(sorted_values,labels,mult,method)
            y=values[target]
            counts[target,0]=np.mean(((y[None,:]<cuts[:,:,0])&valid[target]).sum(axis=1))
            counts[target,1]=np.mean(((y[None,:]>cuts[:,:,1])&valid[target]).sum(axis=1))
    return counts,cut


def build_index_definition(root,out,cfg):
    years=np.arange(cfg['analysis_years'][0],cfg['analysis_years'][1]+1)
    raw=pd.read_csv(root/'data/data.csv')
    raw=raw.loc[raw.year.isin(years)&~((raw.month==2)&(raw.day==29))].copy()
    dates=pd.to_datetime(raw[['year','month','day']])
    raw['doy']=dates.dt.dayofyear-((dates.dt.is_leap_year)&(dates.dt.month>2)).astype(int)
    member=pd.read_csv(out/'tables/thermal_network_membership.csv')
    ids=member.loc[member.common_complete,'station_id'].tolist()
    rows=[];thresholds=[]
    for pos,station in enumerate(ids):
        sdf=raw.loc[raw.station_id==station]
        for variable,cold,warm in [('tmax','cool_days','warm_days'),('tmin','cool_nights','warm_nights')]:
            values=sdf.pivot(index='year',columns='doy',values=variable).reindex(index=years,columns=np.arange(1,366)).to_numpy()
            valid=np.isfinite(values).sum(axis=1)
            assert (valid/365>=cfg['annual_min_valid_fraction']).all()
            for scenario in cfg['scenarios']:
                nbase=scenario['reference_end']-years[0]+1
                counts,cut=construct_counts(values,nbase,scenario['window'],scenario['method'],scenario['correction'],cfg['minimum_window_samples'])
                for k,index in enumerate([cold,warm]):
                    for year,value,nvalid in zip(years,counts[:,k],valid):
                        rows.append(dict(station_id=station,year=year,index_name=index,scenario=scenario['name'],count=value,valid_days=nvalid,rate_pct=100*value/nvalid))
                # Out-of-base cutoffs only; corrected in-base indices average event counts.
                for day,c in enumerate(cut,1):thresholds.append(dict(station_id=station,variable=variable,scenario=scenario['name'],doy=day,q10=c[0],q90=c[1]))
        if (pos+1)%10==0 or pos==0:print(f'Index definitions: {pos+1}/{len(ids)} stations',flush=True)
    table=out/'tables';annual=pd.DataFrame(rows);annual.to_csv(table/'index_definition_annual_counts.csv',index=False)
    pd.DataFrame(thresholds).to_csv(table/'index_definition_thresholds.csv',index=False)
    station=annual.assign(period=np.where(annual.year<=cfg['baseline_years'][1],'early','late')).groupby(['scenario','index_name','station_id','period'])[['count','rate_pct']].mean().unstack('period')
    summary=[]
    for (scenario,index,sid),r in station.iterrows():
        summary.append(dict(scenario=scenario,index_name=index,station_id=sid,early_count=r['count','early'],late_count=r['count','late'],change_days=r['count','late']-r['count','early'],change_rate_pp=r['rate_pct','late']-r['rate_pct','early']))
    station=pd.DataFrame(summary);station.to_csv(table/'index_definition_station_changes.csv',index=False)
    station.groupby(['scenario','index_name']).agg(n_stations=('station_id','size'),early_count=('early_count','mean'),late_count=('late_count','mean'),change_days=('change_days','mean'),change_rate_pp=('change_rate_pp','mean')).reset_index().to_csv(table/'index_definition_period_summary.csv',index=False)
    sampled=pd.read_csv(table/'thermal_network_sampled_year_positions_L4.csv').to_numpy();weights=weights_from_indices(sampled,len(years));x=(years-years[0])/10
    network=annual.groupby(['scenario','index_name','year'])[['count','rate_pct']].mean().reset_index();network.to_csv(table/'index_definition_network_series.csv',index=False)
    trends=[];draw_rows=[]
    for (scenario,index),d in network.groupby(['scenario','index_name']):
        y=d.sort_values('year')['count'].to_numpy();solver=PairQuantileSolver(x,y)
        q=solver.fit(np.ones(len(x)),[.1,.5,.9])[0];b=solver.fit(weights,[.1,.5,.9]);ols=float(weighted_ols(x,y,np.ones(len(x)))[0])
        vals=np.column_stack([weighted_ols(x,y,weights),b,b[:,2]-b[:,0]])
        points=[ols,*q,q[2]-q[0]]
        for j,metric in enumerate(['OLS','q10','q50','q90','Delta1']):
            lo,hi=np.quantile(vals[:,j],[.025,.975]);trends.append(dict(scenario=scenario,index_name=index,metric=metric,estimate=points[j],ci_low=lo,ci_high=hi))
        draw_rows.append(pd.DataFrame(vals,columns=['OLS','q10','q50','q90','Delta1']).assign(scenario=scenario,index_name=index,replicate=np.arange(len(vals))))
    pd.DataFrame(trends).to_csv(table/'index_definition_trends.csv',index=False)
    pd.concat(draw_rows,ignore_index=True).to_csv(table/'index_definition_bootstrap.csv',index=False)
    sources=['data/data.csv','index_definition_config.yaml','src/paper_pipeline/index_definition.py','outputs/publication_v2/tables/thermal_network_membership.csv','outputs/publication_v2/tables/thermal_network_sampled_year_positions_L4.csv']
    meta={'settings':cfg,'n_common_stations':len(ids),'correction':'Exhaustive target-year exclusion and duplication of each other baseline year; average event counts, not thresholds. 17 x 16 replacements for early-reference scenarios.','quantiles':'linear = Hyndman-Fan type 7; median_unbiased = type 8; strict < q10 and > q90; circular calendar-day windows; no leap days.','uncertainty':'Conditional 4999 synchronized original-year pairs block draws (L4), held corrected indices and thresholds; sensitivity intervals, not extra confirmatory tests.','sources':{p:file_hash(root/p) for p in sources}}
    (out/'index_definition_metadata.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    print('Index-definition sensitivities complete.',flush=True)
