"""Separate structural-zero dilution, sample selection and tie definitions."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .compound_partition import partition,circular_blocks,COMPONENTS,file_hash
from .publication_regimes import REGIMES

RULES=['strict','dry_inclusive_only','both_inclusive']

def rates(p,t,pc,tc,rule):
    dry=p<pc if rule=='strict' else p<=pc
    hot=t>=tc if rule=='both_inclusive' else t>tc
    return dry.mean(axis=-2),hot.mean(axis=-2),(dry&hot).mean(axis=-2)


def build_zero_threshold(root,out,cfg):
    table=out/'tables'
    aggregate=pd.read_csv(table/'screened_annual_seasonal_aggregates.csv',float_precision='round_trip')
    primary=pd.read_csv(table/'compound_partition_stations.csv',float_precision='round_trip')
    assignments=pd.read_csv(root/'outputs/tables/koppen_geiger_station_assignments.csv').set_index('station_id')
    reps=cfg['bootstrap_reps'];summaries=[];station_rows=[];draw_rows=[];effects=[];identities=[];zero_prob=[]
    for definition in ['annual','warm_season']:
        archived=primary.loc[primary.definition==definition].sort_values('station_id')
        ids=archived.station_id.to_numpy();climate=assignments.loc[ids,'climate_regime'].to_numpy()
        a=aggregate.loc[aggregate.definition==definition]
        years=np.arange(cfg['analysis_years'][0],cfg['analysis_years'][1]+1)
        p=a.pivot(index='year',columns='station_id',values='precip').reindex(index=years,columns=ids).to_numpy()
        t=a.pivot(index='year',columns='station_id',values='temperature').reindex(index=years,columns=ids).to_numpy()
        early=years<=cfg['baseline_years'][1];p0,p1,t0,t1=p[early],p[~early],t[early],t[~early]
        pc=np.quantile(p0,.25,axis=0);tc=np.quantile(t0,.75,axis=0);zero=pc==0
        assert np.allclose(pc,archived.precip_threshold)
        masks={}
        for group in ['All',*REGIMES]:
            mask=np.ones(len(ids),bool) if group=='All' else climate==group
            for stratum,select in [('all',mask),('positive',mask&~zero),('zero',mask&zero)]:masks[group,stratum]=select
        rng=np.random.default_rng(cfg['seed']+(0 if definition=='annual' else 10000))
        i0=circular_blocks(len(p0),4,reps,rng);i1=circular_blocks(len(p1),4,reps,rng)
        points={}
        for rule in RULES:
            a0=rates(p0,t0,pc,tc,rule);a1=rates(p1,t1,pc,tc,rule)
            points[rule]=100*partition(*a0,*a1)
            if rule=='strict':assert np.allclose(points[rule],archived[list(COMPONENTS)])
            for j,sid in enumerate(ids):
                record=dict(definition=definition,rule=rule,station_id=sid,climate_regime=climate[j],zero_threshold=bool(zero[j]),precip_threshold=pc[j],temperature_threshold=tc[j])
                record.update({c:points[rule][j,k] for k,c in enumerate(COMPONENTS)})
                for period,r in [('early',a0),('late',a1)]:
                    record.update({metric+'_'+period+'_pct':100*r[k][j] for k,metric in enumerate(['dry','hot','joint'])})
                station_rows.append(record)
        for group in ['All',*REGIMES]:
            m=masks[group,'all'];pos=masks[group,'positive'];z=masks[group,'zero'];weight=pos.sum()/m.sum()
            for k,c in enumerate(COMPONENTS):
                allpoint=points['strict'][m,k].mean();positive=points['strict'][pos,k].mean()
                assert np.allclose(points['strict'][z,k],0)
                residual=allpoint-weight*positive
                assert abs(residual)<1e-10
                identities.append(dict(definition=definition,climate_regime=group,component=c,n_all=int(m.sum()),n_positive=int(pos.sum()),n_zero=int(z.sum()),positive_weight=weight,all_estimate_pp=allpoint,positive_estimate_pp=positive,reweighted_positive_pp=weight*positive,residual_pp=residual))
        for threshold_mode in ['fixed','refitted']:
            ensembles={rule:np.empty((reps,len(ids),4)) for rule in RULES};changed=np.zeros(len(ids))
            for start in range(0,reps,100):
                stop=min(start+100,reps);bp0,bt0,bp1,bt1=p0[i0[start:stop]],t0[i0[start:stop]],p1[i1[start:stop]],t1[i1[start:stop]]
                if threshold_mode=='refitted':
                    bc=np.quantile(bp0,.25,axis=1)[:,None,:];hc=np.quantile(bt0,.75,axis=1)[:,None,:]
                    changed+=(bc[:,0,:]>0).sum(axis=0)
                else:bc,hc=pc,tc
                for rule in RULES:
                    ensembles[rule][start:stop]=100*partition(*rates(bp0,bt0,bc,hc,rule),*rates(bp1,bt1,bc,hc,rule))
            if threshold_mode=='refitted':
                for j in np.flatnonzero(zero):zero_prob.append(dict(definition=definition,station_id=ids[j],climate_regime=climate[j],refitted_positive_threshold_fraction=changed[j]/reps))
                saved=pd.read_csv(table/f'bootstrap_network_{definition}.csv').to_numpy()
                assert np.allclose(ensembles['strict'].mean(axis=1),saved)
            for rule,draws in ensembles.items():
                assert np.allclose(draws[:,:,0],draws[:,:,1:].sum(axis=2))
                frame=pd.DataFrame(dict(replicate=np.arange(reps),definition=definition,rule=rule,threshold_mode=threshold_mode))
                for (group,stratum),mask in masks.items():
                    if not mask.any():continue
                    values=draws[:,mask].mean(axis=1);point=points[rule][mask].mean(axis=0)
                    frame[group+'__'+stratum]=values[:,0]
                    for k,c in enumerate(COMPONENTS):
                        lo,hi=np.quantile(values[:,k],[.025,.975])
                        summaries.append(dict(definition=definition,rule=rule,threshold_mode=threshold_mode,climate_regime=group,stratum=stratum,n_stations=int(mask.sum()),n_zero=int((mask&zero).sum()),component=c,estimate_pp=point[k],ci_low_pp=lo,ci_high_pp=hi))
                draw_rows.append(frame)
            for group in ['All',*REGIMES]:
                mask=masks[group,'all']
                for label,a,b in [('dry_tie_change','dry_inclusive_only','strict'),('hot_tie_increment','both_inclusive','dry_inclusive_only')]:
                    diff=(ensembles[a][:,mask,0]-ensembles[b][:,mask,0]).mean(axis=1)
                    lo,hi=np.quantile(diff,[.025,.975])
                    effects.append(dict(definition=definition,threshold_mode=threshold_mode,climate_regime=group,contrast=label,n_stations=int(mask.sum()),estimate_pp=(points[a][mask,0]-points[b][mask,0]).mean(),ci_low_pp=lo,ci_high_pp=hi))
            print(f'Zero threshold: {definition}, {threshold_mode}, {reps} paired draws',flush=True)
    for name,rows in [('summary',summaries),('stations',station_rows),('paired_effects',effects),('dilution_identity',identities),('refit_status',zero_prob)]:pd.DataFrame(rows).to_csv(table/f'zero_threshold_{name}.csv',index=False)
    pd.concat(draw_rows,ignore_index=True).to_csv(table/'zero_threshold_joint_bootstrap.csv',index=False)
    sources=['outputs/publication_v2/tables/screened_annual_seasonal_aggregates.csv','outputs/publication_v2/tables/compound_partition_stations.csv','outputs/tables/koppen_geiger_station_assignments.csv','publication_config.yaml','src/paper_pipeline/zero_threshold.py']
    metadata={'seed':cfg['seed'],'bootstrap_reps':reps,'block_length':4,'rules':RULES,'selection':'Zero/positive strata are fixed from observed baseline precipitation quartiles, not reselected in bootstrap draws. All tie comparisons use identical stations and sampled years.','interpretation':'Strict all-network point estimate equals positive-subset estimate times subset fraction; not independent robustness evidence. Fixed-threshold draws obey the same scaling; refitted draws can give positive cutoffs to observed-zero stations, so interval scaling need not hold.','sources':{p:file_hash(root/p) for p in sources}}
    (out/'zero_threshold_metadata.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
