"""Balanced thermal networks and synchronized year-pairs block uncertainty.

Time labels travel with sampled response fields. This never permutes responses
onto a new time axis or treats stations as independent replicates. Intervals are
approximate, conditional on the archived annual indices and selected stations.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from .compound_partition import circular_blocks, file_hash

INDEXES=['warm_days','warm_nights','cool_days','cool_nights']
TAUS=np.round(np.arange(.1,.901,.01),2)
FOCAL=[.1,.5,.9]
METRICS=['OLS','q10','q50','q90','Delta1']

class PairQuantileSolver:
    """Exact two-parameter check-loss minimization over line intersections.

For a full-rank intercept/time design a finite optimum has a vertex at two
observations. If multiple vertices minimize loss, report the midpoint of the
minimum/maximum optimal slopes. Convexity makes that midpoint optimal too.
Bootstrap counts are nonnegative observation weights, so vertices can be reused.
    """
    def __init__(self,x,y):
        self.x=np.asarray(x,float);self.y=np.asarray(y,float)
        a,b=np.triu_indices(len(x),1)
        distinct=self.x[b]!=self.x[a];a,b=a[distinct],b[distinct]
        self.slopes=(self.y[b]-self.y[a])/(self.x[b]-self.x[a])
        intercept=self.y[a]-self.slopes*self.x[a]
        residual=self.y[None,:]-intercept[:,None]-self.slopes[:,None]*self.x[None,:]
        residual[np.abs(residual)<1e-11]=0
        self.residual=residual
        self.negative=np.maximum(-residual,0)

    def fit(self,weights,taus):
        w=np.atleast_2d(weights).astype(float)
        if np.any((w>0).sum(axis=1)<3):raise ValueError('At least three distinct sampled years are required.')
        out=np.empty((len(w),len(taus)))
        for start in range(0,len(w),400):
            stop=min(start+400,len(w));wt=w[start:stop].T
            neg=self.negative@wt;total=self.residual@wt
            for j,tau in enumerate(taus):
                loss=neg+float(tau)*total;best=loss.min(axis=0)
                optimal=np.abs(loss-best[None,:])<=1e-9*(1+np.abs(best[None,:]))
                low=np.where(optimal,self.slopes[:,None],np.inf).min(axis=0)
                high=np.where(optimal,self.slopes[:,None],-np.inf).max(axis=0)
                out[start:stop,j]=(low+high)/2
        return out

def weighted_ols(x,y,w):
    w=np.atleast_2d(w);s=w.sum(axis=1);mx=(w@x)/s;my=(w@y)/s
    return ((w@(x*y))-s*mx*my)/((w@(x*x))-s*mx*mx)

def weights_from_indices(indices,n):
    return np.stack([np.bincount(row,minlength=n) for row in indices])

def summarize(values,point,**keys):
    lo,hi=np.quantile(values,[.025,.975])
    flo,fhi=np.quantile(values,[.05/12,1-.05/12])
    return dict(**keys,estimate=float(point),ci_low=float(lo),ci_high=float(hi),
                family_ci_low=float(flo),family_ci_high=float(fhi),bootstrap_reps=len(values))

def build_thermal_network(root:Path,out:Path,cfg:dict):
    settings=cfg['thermal_network'];reps=int(settings['bootstrap_reps']);seed=int(settings['seed'])
    table=out/'tables';table.mkdir(parents=True,exist_ok=True)
    source=root/'outputs/tables/annual_extreme_indices.csv'
    annual=pd.read_csv(source);years=np.arange(*[cfg['analysis_years'][0],cfg['analysis_years'][1]+1]);x=(years-years[0])/10
    pivots={idx:annual.pivot(index='year',columns='station_id',values=idx).reindex(years) for idx in INDEXES}
    ids=pivots[INDEXES[0]].columns
    member=pd.DataFrame({'station_id':ids})
    for idx in INDEXES:
        member[idx+'_valid_years']=pivots[idx].notna().sum().reindex(ids).to_numpy()
        member[idx+'_complete']=member[idx+'_valid_years'].eq(len(years))
    member['common_complete']=member[[i+'_complete' for i in INDEXES]].all(axis=1)
    common=member.loc[member.common_complete,'station_id'].to_numpy()
    member.to_csv(table/'thermal_network_membership.csv',index=False)
    coverage=[];series=[];profiles=[];points={};solvers={};intervals=[];draw_frames=[];contrasts=[];influence=[]
    for idx in INDEXES:
        p=pivots[idx];complete=member.loc[member[idx+'_complete'],'station_id'].to_numpy()
        valid_col='valid_days_tmax' if 'days' in idx else 'valid_days_tmin'
        v=annual.pivot(index='year',columns='station_id',values=valid_col).reindex(years)
        networks={'available':p.mean(axis=1),'index_fixed':p[complete].mean(axis=1),
                  'common_fixed':p[common].mean(axis=1),
                  'common_365_equivalent':(365*p[common]/v[common]).mean(axis=1)}
        for year in years:
            coverage.append(dict(index_name=idx,year=year,n_available=int(p.loc[year].count()),
                                 n_index_fixed=len(complete),n_common_fixed=len(common),
                                 mean_valid_days_common=float(v.loc[year,common].mean())))
        for network,y in networks.items():
            n=len(common) if network.startswith('common') else len(complete) if network=='index_fixed' else np.nan
            for year,value in zip(years,y):series.append(dict(index_name=idx,network=network,year=year,annual_count=value,n_stations=n if np.isfinite(n) else int(p.loc[year].count())))
            solver=PairQuantileSolver(x,y.to_numpy());solvers[idx,network]=solver
            slopes=solver.fit(np.ones((1,len(x))),TAUS)[0];ols=float(weighted_ols(x,y.to_numpy(),np.ones((1,len(x))))[0])
            focal=[slopes[np.flatnonzero(TAUS==t)[0]] for t in FOCAL]
            points[idx,network]=dict(zip(METRICS,[ols,*focal,focal[2]-focal[0]]))
            for t,b in zip(TAUS,slopes):profiles.append(dict(index_name=idx,network=network,tau=t,slope=b,ols_slope=ols,n_stations=n,ci_low=np.nan,ci_high=np.nan))
        for omit,year in enumerate(years):
            w=np.ones((1,len(x)));w[0,omit]=0
            b=solvers[idx,'common_fixed'].fit(w,FOCAL)[0]
            influence.append(dict(index_name=idx,omitted_year=year,q10=b[0],q50=b[1],q90=b[2],Delta1=b[2]-b[0]))
    profiles=pd.DataFrame(profiles)
    for length in settings['block_lengths']:
        rng=np.random.default_rng(seed+int(length))
        sampled=circular_blocks(len(x),int(length),reps,rng)
        if ((np.apply_along_axis(lambda a:len(np.unique(a)),1,sampled))<3).any():raise ValueError('Degenerate bootstrap draw; review seed/design.')
        weights=weights_from_indices(sampled,len(x))
        pd.DataFrame(sampled,columns=[f'position_{i}' for i in range(len(x))]).to_csv(table/f'thermal_network_sampled_year_positions_L{length}.csv',index=False)
        values={}
        for network in ['common_fixed','available']:
            for idx in INDEXES:
                solver=solvers[idx,network]
                tau_grid=TAUS if network=='common_fixed' and length==settings['primary_block_length'] else FOCAL
                qr=solver.fit(weights,tau_grid)
                if len(tau_grid)>3:
                    mask=(profiles.index_name==idx)&(profiles.network==network)
                    profiles.loc[mask,'ci_low']=np.quantile(qr,.025,axis=0)
                    profiles.loc[mask,'ci_high']=np.quantile(qr,.975,axis=0)
                    qr=qr[:,[np.flatnonzero(TAUS==t)[0] for t in FOCAL]]
                ols=weighted_ols(x,solver.y,weights)
                array=np.column_stack([ols,qr,qr[:,2]-qr[:,0]])
                frame=pd.DataFrame(array,columns=METRICS);frame['replicate']=np.arange(reps);frame['index_name']=idx;frame['network']=network;frame['block_length']=length;draw_frames.append(frame)
                for k,metric in enumerate(METRICS):
                    values[idx,network,metric]=array[:,k]
                    intervals.append(summarize(array[:,k],points[idx,network][metric],index_name=idx,network=network,metric=metric,block_length=length))
                print(f'Thermal bootstrap L={length} {network} {idx}: {reps} draws',flush=True)
        for label,day,night in [('warm_day_minus_night','warm_days','warm_nights'),('cool_day_minus_night','cool_days','cool_nights')]:
            d=values[day,'common_fixed','Delta1']-values[night,'common_fixed','Delta1']
            point=points[day,'common_fixed']['Delta1']-points[night,'common_fixed']['Delta1']
            contrasts.append(summarize(d,point,contrast=label,index_name='paired_indices',metric='Delta1',block_length=length))
        for idx in INDEXES:
            for metric in METRICS:
                d=values[idx,'common_fixed',metric]-values[idx,'available',metric]
                contrasts.append(summarize(d,points[idx,'common_fixed'][metric]-points[idx,'available'][metric],contrast='fixed_minus_available',index_name=idx,metric=metric,block_length=length))
    pd.DataFrame(coverage).to_csv(table/'thermal_network_coverage.csv',index=False)
    pd.DataFrame(series).to_csv(table/'thermal_network_series.csv',index=False)
    profiles.to_csv(table/'thermal_network_profiles.csv',index=False)
    pd.concat(draw_frames,ignore_index=True).to_csv(table/'thermal_network_bootstrap_focal.csv',index=False)
    pd.DataFrame(intervals).to_csv(table/'thermal_network_intervals.csv',index=False)
    pd.DataFrame(contrasts).to_csv(table/'thermal_network_contrasts.csv',index=False)
    pd.DataFrame(influence).to_csv(table/'thermal_network_leave_one_year_out.csv',index=False)
    point_rows=[dict(index_name=idx,network=network,**metrics) for (idx,network),metrics in points.items()]
    pd.DataFrame(point_rows).to_csv(table/'thermal_network_comparison.csv',index=False)
    metadata={'seed':seed,'settings':settings,'common_station_count':len(common),'years':years.tolist(),
              'bootstrap':'Circular moving blocks of (original year, full network response vector) pairs; shared indices across networks and four indices. Observed station membership held fixed.',
              'solver':'Exact line-intersection check-loss minimization; midpoint of extreme optimal vertex slopes for nonunique solutions.',
              'units':'Days per decade; common_365_equivalent uses 365*observed count/valid days as a coverage diagnostic only.',
              'intervals':'Pointwise percentile 95%; nominal six-contrast Bonferroni endpoints for four common-network Delta1 values and two paired day-minus-night contrasts. Family columns for other quantities are not confirmatory.',
              'limitations':'Conditional on annual index thresholds and observed fixed station set; no daily threshold refitting, homogenization, or guarantee of finite-sample coverage. Time covariates retained; approximate block resampling in a trending short record.',
              'input_hashes':{source.relative_to(root).as_posix():file_hash(source),'publication_config.yaml':file_hash(root/'publication_config.yaml'),Path(__file__).relative_to(root).as_posix():file_hash(Path(__file__))}}
    (out/'thermal_network_metadata.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
    print(f'Thermal network completed: {len(common)} common stations.',flush=True)
