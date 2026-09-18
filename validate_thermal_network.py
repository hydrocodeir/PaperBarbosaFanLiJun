"""Independent LP, network-membership and synchronized-contrast checks."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from scipy.optimize import linprog
from src.paper_pipeline.thermal_network import PairQuantileSolver,weights_from_indices,INDEXES,METRICS,TAUS

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'outputs/publication_v2';T=OUT/'tables'

def lp_loss(x,y,w,t):
    keep=w>0;x=x[keep];y=y[keep];w=w[keep];n=len(x)
    a=np.column_stack([np.ones(n),x,np.eye(n),-np.eye(n)])
    fit=linprog(np.r_[0,0,t*w,(1-t)*w],A_eq=a,b_eq=y,bounds=[(None,None)]*2+[(0,None)]*(2*n),method='highs')
    assert fit.success
    return fit.fun

def loss_at_slope(x,y,w,t,b):
    residual=y-b*x;order=np.argsort(residual);v=residual[order];cum=np.cumsum(w[order]);a=v[np.searchsorted(cum,t*w.sum(),side='left')]
    r=y-a-b*x
    return np.sum(w*np.where(r>=0,t*r,(t-1)*r))

def main():
    checks={};rng=np.random.default_rng(891);largest=0
    x=np.arange(34)/10
    for trial in range(24):
        y=2+3*x+rng.normal(size=34)*(1+x) if trial%2 else rng.integers(-2,7,34).astype(float)
        w=rng.multinomial(34,np.full(34,1/34)).astype(float)
        solver=PairQuantileSolver(x,y)
        for t,b in zip([.1,.5,.9],solver.fit(w,[.1,.5,.9])[0]):
            error=abs(loss_at_slope(x,y,w,t,b)-lp_loss(x,y,w,t));largest=max(largest,error)
            assert error<1e-6
    linear=PairQuantileSolver(x,4+7*x).fit(np.ones(34),[.1,.5,.9])
    assert np.allclose(linear,7)
    checks['independent_linear_program_check_loss_72_weighted_fits']={'max_loss_difference':largest,'status':'PASS'}
    checks['known_linear_trend_all_quantiles']='PASS'
    raw=pd.read_csv(ROOT/'outputs/tables/annual_extreme_indices.csv')
    membership=pd.read_csv(T/'thermal_network_membership.csv');common=set(membership.loc[membership.common_complete,'station_id'])
    expected=[]
    for idx in INDEXES:
        p=raw.pivot(index='year',columns='station_id',values=idx)
        expected.append(set(p.columns[p.notna().all()]))
    assert common==set.intersection(*expected) and len(common)==108
    checks['common_membership_independently_reconstructed']=len(common)
    series=pd.read_csv(T/'thermal_network_series.csv');profiles=pd.read_csv(T/'thermal_network_profiles.csv')
    saved=pd.read_csv(T/'thermal_network_bootstrap_focal.csv');intervals=pd.read_csv(T/'thermal_network_intervals.csv');contrasts=pd.read_csv(T/'thermal_network_contrasts.csv')
    for idx in INDEXES:
        p=raw.pivot(index='year',columns='station_id',values=idx)
        s=series.loc[(series.index_name==idx)&(series.network=='common_fixed')].sort_values('year')
        assert np.allclose(s.annual_count,p[sorted(common)].mean(axis=1))
        assert s.n_stations.eq(108).all()
        for network in ['common_fixed','available']:
            y=series.loc[(series.index_name==idx)&(series.network==network)].sort_values('year').annual_count.to_numpy()
            solver=PairQuantileSolver(x,y)
            for length in [2,4,6]:
                sampled=pd.read_csv(T/f'thermal_network_sampled_year_positions_L{length}.csv').to_numpy()
                assert sampled.shape==(4999,34)
                assert np.all((np.diff(sampled[:,:length],axis=1)%34)==1)
                chosen=[0,11,233,4998];w=weights_from_indices(sampled[chosen],34)
                d=saved.loc[(saved.index_name==idx)&(saved.network==network)&(saved.block_length==length)].set_index('replicate')
                rebuilt=solver.fit(w,[.1,.5,.9]);assert np.allclose(rebuilt,d.loc[chosen,['q10','q50','q90']])
                for j in range(4):
                    for k,t in enumerate([.1,.5,.9]):assert abs(loss_at_slope(x,y,w[j],t,rebuilt[j,k])-lp_loss(x,y,w[j],t))<1e-6
                assert np.allclose(d.Delta1,d.q90-d.q10)
    checks['station_means_and_fixed_membership']='PASS'
    checks['288_real_bootstrap_fits_checked_against_LP']='PASS'
    checks['sampled_year_pairs_block_contiguity_and_reconstruction']='PASS'
    for r in intervals.itertuples():
        d=saved.loc[(saved.index_name==r.index_name)&(saved.network==r.network)&(saved.block_length==r.block_length),r.metric]
        assert np.allclose(np.quantile(d,[.025,.975]),[r.ci_low,r.ci_high])
        assert np.allclose(np.quantile(d,[.05/12,1-.05/12]),[r.family_ci_low,r.family_ci_high])
    for r in contrasts.itertuples():
        d=saved.loc[saved.block_length==r.block_length]
        if r.contrast=='fixed_minus_available':
            a=d.loc[(d.index_name==r.index_name)&(d.network=='common_fixed')].set_index('replicate')[r.metric]
            b=d.loc[(d.index_name==r.index_name)&(d.network=='available')].set_index('replicate')[r.metric]
        else:
            prefix='warm' if r.contrast.startswith('warm') else 'cool'
            a=d.loc[(d.index_name==prefix+'_days')&(d.network=='common_fixed')].set_index('replicate').Delta1
            b=d.loc[(d.index_name==prefix+'_nights')&(d.network=='common_fixed')].set_index('replicate').Delta1
        assert np.allclose(np.quantile(a-b,[.025,.975]),[r.ci_low,r.ci_high])
        assert np.allclose(np.quantile(a-b,[.05/12,1-.05/12]),[r.family_ci_low,r.family_ci_high])
    checks['all_reported_pointwise_intervals_and_paired_contrasts']='PASS'
    meta=json.loads((OUT/'thermal_network_metadata.json').read_text())
    sampled=pd.read_csv(T/f"thermal_network_sampled_year_positions_L{meta['settings']['primary_block_length']}.csv").to_numpy()
    weights=weights_from_indices(sampled,34)
    for idx in INDEXES:
        y=series.loc[(series.index_name==idx)&(series.network=='common_fixed')].sort_values('year').annual_count.to_numpy()
        solver=PairQuantileSolver(x,y)
        d=profiles.loc[(profiles.index_name==idx)&(profiles.network=='common_fixed')].sort_values('tau')
        assert np.allclose(d.tau,TAUS)
        assert np.allclose(solver.fit(np.ones(34),TAUS)[0],d.slope)
        draws=solver.fit(weights,TAUS)
        assert np.allclose(np.quantile(draws,[.025,.975],axis=0).T,d[['ci_low','ci_high']])
        for row in d.itertuples():
            assert abs(loss_at_slope(x,y,np.ones(34),row.tau,row.slope)-lp_loss(x,y,np.ones(34),row.tau))<1e-6
    checks['324_primary_profile_points_LP_checked_and_all_bands_reconstructed']='PASS'
    for rel,h in meta['input_hashes'].items():assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==h
    checks['source_hashes']='PASS'
    result={'checks':checks,'scope':'Independent optimization, arithmetic and provenance checks; not empirical confirmation of nominal interval coverage or raw-data homogenization.'}
    (OUT/'thermal_network_validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
