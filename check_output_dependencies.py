"""Cross-table checks complementing the isolated numerical rerun."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import adjusted_rand_score
from sklearn.preprocessing import StandardScaler
from src.paper_pipeline.clustering import screen_clustering_features
from audit_output_data import compare, record, RESULTS

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'outputs';PUB=OUT/'publication_v2';AUDIT=OUT/'audit_cleanup'

def main():
    RESULTS.extend(json.loads((AUDIT/'numerical_checks.json').read_text()))
    load=lambda name:pd.read_csv(OUT/'tables'/f'{name}.csv',float_precision='round_trip')
    cfg=yaml.safe_load((ROOT/'config.yaml').read_text())
    q=load('qr_focus_slopes_and_bootstrap_summary');f=load('clustering_feature_table')
    compare('clustering_feature_table_source_columns',q,f)
    compare('publication_summary_table_source_columns',f,load('publication_summary_table'))
    for idx,d in f.groupby('index_name'):
        assert np.allclose(d.Delta1,d['slope_0.90']-d['slope_0.10'])
    compare('clustering_feature_screening_summary',load('clustering_feature_screening_summary'),pd.concat([
        screen_clustering_features(f,cfg,feature_set_label='baseline'),
        screen_clustering_features(f,cfg,feature_cols=cfg['clustering']['robustness_check']['reduced_features'],feature_set_label='sensitivity_rerun')],ignore_index=True))
    a=load('cluster_assignments').merge(load('cluster_assignments_reduced_features'),on=['index_name','station_id','station_name'])
    for _,r in load('cluster_robustness_summary').iterrows():
        d=a.loc[a.index_name==r.index_name];assert np.isclose(adjusted_rand_score(d.cluster,d.cluster_reduced_features),r.adjusted_rand_index)
    record('cluster_robustness_summary','MATCH','ARI independently recomputed from assignments.')
    representatives=load('representative_station_selection')
    for idx,d in f.groupby('index_name'):
        cols=cfg['clustering']['simple_features'];screen=screen_clustering_features(d,cfg)
        cols=screen.loc[screen.status.str.startswith('kept'),'feature'].tolist()
        x=d[cols].fillna(d[cols].median());z=pd.DataFrame(StandardScaler().fit_transform(x),index=d.index,columns=cols)
        for cluster,g in d.groupby('cluster'):
            dist=np.linalg.norm(z.loc[g.index]-z.loc[g.index].mean(),axis=1)
            chosen=g.assign(distance=dist).sort_values(['distance','station_name']).iloc[0]
            r=representatives.loc[(representatives.index_name==idx)&(representatives.cluster==cluster)].iloc[0]
            assert chosen.station_id==r.station_id
            assert np.isclose(chosen.distance,r.distance_to_cluster_centroid)
    record('representative_station_selection','MATCH','Every representative and centroid distance independently reconstructed.')
    p=pd.read_csv(PUB/'tables/compound_partition_primary.csv');s=pd.read_csv(PUB/'tables/compound_partition_stations.csv')
    for definition in ['annual','warm_season']:
        draws=pd.read_csv(PUB/f'tables/bootstrap_network_{definition}.csv')
        for _,r in p.loc[p.definition==definition].iterrows():
            lo,hi=np.quantile(draws[r.component],[.025,.975]);flo,fhi=np.quantile(draws[r.component],[.003125,.996875])
            assert np.allclose([lo,hi,flo,fhi],[r.ci_low_pp,r.ci_high_pp,r.family_ci_low_pp,r.family_ci_high_pp])
            assert np.isclose(s.loc[s.definition==definition,r.component].mean(),r.estimate_pp)
    record('primary_compound_intervals_and_station_means','MATCH','All eight estimates and both interval families reproduced from saved fields.')
    regimes=pd.read_csv(PUB/'tables/climate_regime_compound_partition.csv');draws=pd.read_csv(PUB/'tables/climate_regime_compound_bootstrap.csv')
    for _,r in regimes.iterrows():
        d=draws.loc[(draws.definition==r.definition)&(draws.climate_regime==r.climate_regime),r.component]
        assert len(d)==4999
        assert np.allclose(np.quantile(d,[.025,.975]),[r.ci_low_pp,r.ci_high_pp])
    for (definition,component),d in regimes.groupby(['definition','component']):
        r=p.loc[(p.definition==definition)&(p.component==component)].iloc[0]
        assert d.n_stations.sum()==r.n_stations
        assert np.isclose(np.average(d.estimate_pp,weights=d.n_stations),r.estimate_pp)
    record('regime_components_intervals_and_network_closure','MATCH','48 group intervals and eight weighted network means checked.')
    scenarios=pd.read_csv(PUB/'tables/compound_partition_sensitivity.csv');stations=pd.read_csv(PUB/'tables/compound_partition_all_scenario_stations.csv')
    for _,r in scenarios.iterrows():
        d=stations.loc[(stations.definition==r.definition)&(stations.scenario==r.scenario)]
        assert len(d)==r.n_stations
        assert np.isclose(d[r.component].mean(),r.estimate_pp)
        assert r.ci_low_pp<=r.ci_high_pp
    record('all_sensitivity_station_means','MATCH','All 80 scenario-component means and sample sizes reproduced; interval ordering checked.')
    record('dependency_check_completed','DONE','No additional bootstrap realizations simulated in this check.')

if __name__=='__main__':main()
