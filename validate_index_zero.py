"""Independent duplication, historical-index, mixture and interval checks."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from src.paper_pipeline.index_definition import window_samples,repeated_quantiles

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'outputs/publication_v2';T=OUT/'tables'

def main():
    checks={};rng=np.random.default_rng(437)
    base=rng.normal(size=(5,365));base[0,3:17]=np.nan
    for window in [5,11]:
        v,labels=window_samples(base,window)
        mult=np.array([[0,1,2,1,1],[1,2,1,0,1]])
        for method in ['linear','median_unbiased']:
            result=repeated_quantiles(v,labels,mult,method)
            for j,w in enumerate(mult):
                for d in [0,7,111,364]:
                    explicit=np.repeat(v[d],w[labels[d]])
                    assert np.allclose(result[j,d],np.nanquantile(explicit,[.1,.9],method=method))
    checks['32_repeated_sample_quantile_checks_T7_T8']='PASS'
    counts=pd.read_csv(T/'index_definition_annual_counts.csv');ids=counts.station_id.unique()
    for case,filename in [('full_w11_t7','annual_extreme_indices.csv'),('fixed_w11_t7_raw','fixed_baseline_annual_extreme_indices.csv')]:
        historical=pd.read_csv(ROOT/'outputs/tables'/filename).set_index(['station_id','year'])
        for idx,d in counts.loc[counts.scenario==case].groupby('index_name'):
            values=d.set_index(['station_id','year'])['count'].sort_index()
            assert np.array_equal(values,historical.loc[values.index,idx]),(case,idx)
    checks['all_29376_historical_station_year_count_matches']='PASS'
    # Literal year replacement and NumPy quantiles, independent of weighted-order implementation.
    raw=pd.read_csv(ROOT/'data/data.csv');raw=raw.loc[raw.station_id==ids[0]].copy()
    raw=raw.loc[~((raw.month==2)&(raw.day==29))]
    date=pd.to_datetime(raw[['year','month','day']]);raw['doy']=date.dt.dayofyear-((date.dt.is_leap_year)&(date.dt.month>2)).astype(int)
    for variable,cold,warm in [('tmax','cool_days','warm_days'),('tmin','cool_nights','warm_nights')]:
        values=raw.pivot(index='year',columns='doy',values=variable).reindex(index=range(1991,2025),columns=range(1,366)).to_numpy();base=values[:17]
        for window,method,case in [(11,'linear','fixed_w11_t7_corrected'),(5,'linear','fixed_w5_t7_corrected'),(5,'median_unbiased','fixed_w5_t8_corrected')]:
            days=(np.arange(365)[:,None]+np.arange(-(window//2),window//2+1))%365
            for target in [0,8,16]:
                totals=[]
                for donor in np.delete(np.arange(17),target):
                    replaced=base.copy();replaced[target]=base[donor]
                    sample=replaced[:,days].transpose(1,0,2).reshape(365,-1)
                    cut=np.nanquantile(sample,[.1,.9],axis=1,method=method)
                    totals.append([(values[target]<cut[0]).sum(),(values[target]>cut[1]).sum()])
                expected=np.mean(totals,axis=0)
                for j,idx in enumerate([cold,warm]):
                    actual=counts.loc[(counts.station_id==ids[0])&(counts.year==1991+target)&(counts.scenario==case)&(counts.index_name==idx),'count'].item()
                    assert abs(actual-expected[j])<1e-10,(case,idx,target,actual,expected[j])
    checks['literal_16_donor_replacement_counts_for_36_station_year_indices']='PASS'
    rawcase=counts.loc[(counts.scenario=='fixed_w11_t7_raw')&(counts.year>2007)].set_index(['station_id','year','index_name'])['count'].sort_index()
    corrected=counts.loc[(counts.scenario=='fixed_w11_t7_corrected')&(counts.year>2007)].set_index(['station_id','year','index_name'])['count'].sort_index()
    assert np.array_equal(rawcase,corrected)
    assert np.allclose(counts.rate_pct,100*counts['count']/counts.valid_days)
    checks['correction_leaves_out_of_base_counts_unchanged_and_rates_match']='PASS'
    b=pd.read_csv(T/'index_definition_bootstrap.csv');trends=pd.read_csv(T/'index_definition_trends.csv')
    assert np.allclose(b.Delta1,b.q90-b.q10)
    for r in trends.itertuples():
        d=b.loc[(b.scenario==r.scenario)&(b.index_name==r.index_name),r.metric]
        assert np.allclose(np.quantile(d,[.025,.975]),[r.ci_low,r.ci_high])
    checks['all_120_index_sensitivity_intervals_reconstructed']='PASS'
    identity=pd.read_csv(T/'zero_threshold_dilution_identity.csv')
    assert np.allclose(identity.all_estimate_pp,identity.positive_weight*identity.positive_estimate_pp)
    assert np.allclose(identity.positive_weight,identity.n_positive/identity.n_all)
    checks['all_56_network_regime_component_dilution_identities']='PASS'
    draws=pd.read_csv(T/'zero_threshold_joint_bootstrap.csv');summary=pd.read_csv(T/'zero_threshold_summary.csv')
    for r in summary.loc[summary.component=='joint_change'].itertuples():
        d=draws.loc[(draws.definition==r.definition)&(draws.rule==r.rule)&(draws.threshold_mode==r.threshold_mode)]
        values=d[r.climate_regime+'__'+r.stratum]
        assert np.allclose(np.quantile(values,[.025,.975]),[r.ci_low_pp,r.ci_high_pp])
    for r in identity.loc[identity.component=='joint_change'].itertuples():
        d=draws.loc[(draws.definition==r.definition)&(draws.rule=='strict')&(draws.threshold_mode=='fixed')]
        assert np.allclose(d[r.climate_regime+'__all'],r.positive_weight*d[r.climate_regime+'__positive'])
        if r.n_zero:assert np.allclose(d[r.climate_regime+'__zero'],0)
    effects=pd.read_csv(T/'zero_threshold_paired_effects.csv')
    for r in effects.itertuples():
        d=draws.loc[(draws.definition==r.definition)&(draws.threshold_mode==r.threshold_mode)]
        a,c=('dry_inclusive_only','strict') if r.contrast=='dry_tie_change' else ('both_inclusive','dry_inclusive_only')
        delta=d.loc[d.rule==a].sort_values('replicate')[r.climate_regime+'__all'].to_numpy()-d.loc[d.rule==c].sort_values('replicate')[r.climate_regime+'__all'].to_numpy()
        assert np.allclose(np.quantile(delta,[.025,.975]),[r.ci_low_pp,r.ci_high_pp])
    checks['joint_intervals_paired_tie_intervals_and_fixed_threshold_scaling']='PASS'
    stations=pd.read_csv(T/'zero_threshold_stations.csv');agg=pd.read_csv(T/'screened_annual_seasonal_aggregates.csv',float_precision='round_trip')
    for r in stations.itertuples():
        d=agg.loc[(agg.definition==r.definition)&(agg.station_id==r.station_id)]
        baseline=d.loc[d.year<=2007]
        pc=np.quantile(baseline.precip,.25);tc=np.quantile(baseline.temperature,.75)
        dry=d.precip<pc if r.rule=='strict' else d.precip<=pc
        hot=d.temperature>=tc if r.rule=='both_inclusive' else d.temperature>tc
        for period,mask in [('early',d.year<=2007),('late',d.year>=2008)]:
            for name,flag in [('dry',dry),('hot',hot),('joint',dry&hot)]:assert np.isclose(100*flag[mask].mean(),getattr(r,name+'_'+period+'_pct'))
    checks['all_621_station_rule_event_frequencies_independently_reclassified']='PASS'
    for name in ['index_definition_metadata.json','zero_threshold_metadata.json']:
        for rel,digest in json.loads((OUT/name).read_text())['sources'].items():assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==digest,rel
    checks['new_analysis_source_hashes']='PASS'
    result={'checks':checks,'scope':'Independent index-count/quantile/reclassification/arithmetic/provenance checks; not validation of nominal bootstrap coverage or homogenization.'}
    (OUT/'index_zero_validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
