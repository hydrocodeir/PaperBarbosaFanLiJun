"""Inventory every output and reproduce numerical dependencies in an isolated workspace.

No source output is overwritten or deleted by this script. Expensive historical
bootstrap ensembles are validated from saved draws, not silently claimed rerun.
"""
from pathlib import Path
from collections import defaultdict
import argparse, copy, hashlib, json, shutil, time, warnings
import numpy as np
import pandas as pd
import yaml
from PIL import Image
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'outputs'
AUDIT=OUT/'audit_cleanup'
WORK=ROOT/'.audit_work/recomputed'
RESULTS=[]

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
    return h.hexdigest()

def inventory():
    AUDIT.mkdir(parents=True,exist_ok=True)
    if (AUDIT/'inventory_before.csv').exists():
        raise FileExistsError('Preserve the original audit inventory. Use a separate release directory for a new inventory.')
    rows=[]; schemas={}
    for p in sorted(OUT.rglob('*')):
        if not p.is_file() or AUDIT in p.parents: continue
        row=dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=digest(p),extension=p.suffix,status='READABLE',details='')
        try:
            if p.suffix=='.csv':
                d=pd.read_csv(p,float_precision='round_trip')
                schemas[row['path']]={'rows':len(d),'columns':list(d.columns),'missing':d.isna().sum().to_dict(),'exact_duplicate_rows':int(d.duplicated().sum())}
                row['details']=f'{len(d)} rows; {len(d.columns)} columns; {int(d.duplicated().sum())} exact duplicates'
                if np.isinf(d.select_dtypes(include='number').to_numpy()).any(): row['status']='REVIEW_INFINITY'
            elif p.suffix.lower() in ['.png','.tiff']:
                with Image.open(p) as im: row['details']=str(im.size); im.verify()
            elif p.suffix=='.pdf': row['details']=f'{len(PdfReader(p).pages)} pages'
            elif p.suffix=='.json': json.loads(p.read_text(encoding='utf-8'))
            else: p.read_bytes()
        except Exception as exc: row['status']='READ_FAILED'; row['details']=str(exc)
        rows.append(row)
    pd.DataFrame(rows).to_csv(AUDIT/'inventory_before.csv',index=False)
    (AUDIT/'table_schemas.json').write_text(json.dumps(schemas,indent=2),encoding='utf-8')
    print(pd.DataFrame(rows).status.value_counts().to_dict(),flush=True)

def record(name,status,detail):
    RESULTS.append(dict(check=name,status=status,detail=detail))
    (AUDIT/'numerical_checks.json').write_text(json.dumps(RESULTS,indent=2,default=str),encoding='utf-8')
    print(name,status,str(detail)[:240],flush=True)

def compare(name,a,b):
    common=[c for c in a if c in b]
    issues=[]; max_diff=0.
    if len(a)!=len(b): issues.append(f'rows {len(a)} versus {len(b)}')
    else:
        a=a.reset_index(drop=True);b=b.reset_index(drop=True)
        for c in common:
            if pd.api.types.is_numeric_dtype(a[c]) and pd.api.types.is_numeric_dtype(b[c]):
                x=a[c].to_numpy(dtype=float);y=b[c].to_numpy(dtype=float)
                if np.any(np.isfinite(x)&np.isfinite(y)): max_diff=max(max_diff,float(np.nanmax(np.abs(x-y))))
                if not np.allclose(x,y,equal_nan=True,atol=1e-7,rtol=1e-7): issues.append(c)
            elif not a[c].fillna('').astype(str).equals(b[c].fillna('').astype(str)): issues.append(c)
    record(name,'MISMATCH' if issues else 'MATCH',dict(rows=len(a),columns_checked=len(common),max_abs_difference=max_diff,mismatched_columns=issues))

def runstep(label,fn):
    before={p:p.stat().st_mtime_ns for p in WORK.rglob('*.csv')}
    print('START',label,flush=True); start=time.monotonic()
    fn()
    for p in WORK.rglob('*.csv'):
        if p.stat().st_mtime_ns==before.get(p):continue
        source=OUT/p.relative_to(WORK)
        if source.exists(): compare(source.relative_to(OUT).as_posix(),pd.read_csv(source,float_precision='round_trip'),pd.read_csv(p,float_precision='round_trip'))
    print('END',label,round(time.monotonic()-start,1),'seconds',flush=True)

def recompute():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.figure
    import matplotlib.pyplot as plt
    # Numerical outputs are reproduced; legacy figure exports are deliberately suppressed.
    matplotlib.figure.Figure.savefig=lambda *a,**k:None
    from src.paper_pipeline import quantile as qr, clustering as cl, advanced_analysis as aa
    from src.paper_pipeline import climate_change_signal as cs, climate_regime_analysis as cr
    from src.paper_pipeline import compound_dry_hot as cd, qc_sensitivity as qc
    from src.paper_pipeline import data_quality as dq, homogeneity_sensitivity as hs, clustering_sensitivity as cls
    from src.paper_pipeline import bootstrap_depth_sensitivity as bd
    from src.paper_pipeline.year_config import filter_to_analysis_years
    cfg=yaml.safe_load((ROOT/'config.yaml').read_text())
    for folder in ['tables','compound_dry_hot/tables']:
        (WORK/folder).mkdir(parents=True,exist_ok=True)
        for p in (OUT/folder).glob('*.csv'): shutil.copy2(p,WORK/folder/p.name)
    raw=filter_to_analysis_years(pd.read_csv(ROOT/cfg['paths']['data_csv']),cfg)
    stations=pd.read_csv(ROOT/cfg['paths']['station_csv'])
    load=lambda n:pd.read_csv(OUT/'tables'/f'{n}.csv',float_precision='round_trip')
    annual=load('annual_extreme_indices'); summary=load('qr_focus_slopes_and_bootstrap_summary'); features=load('clustering_feature_table')
    # 99,200 saved replicates support exact recomputation of all stored bootstrap summaries.
    boot=load('bootstrap_distributions_long'); rows=[]
    for (idx,sid),group in boot.groupby(['index_name','station_id'],sort=False):
        row={'index_name':idx,'station_id':sid}; row.update(qr.summarize_bootstrap(group.drop(columns=['index_name','station_id','station_name']),cfg['bootstrap']['alpha'])); rows.append(row)
    rebuilt=pd.DataFrame(rows); expected=summary.set_index(['index_name','station_id']).loc[rebuilt.set_index(['index_name','station_id']).index].reset_index()
    compare('saved_bootstrap_all_summaries',expected,rebuilt)
    assert boot.groupby(['index_name','station_id']).size().eq(200).all()
    record('saved_bootstrap_replicate_counts','MATCH',len(boot))
    fast=copy.deepcopy(cfg);fast['bootstrap']['enabled']=False
    def rebuild_qr():
        full,focus,_=qr.run_station_qr(annual,fast)
        compare('all_40176_station_quantile_slopes',load('qr_all_quantiles_long'),full)
        compare('all_station_focal_slopes_and_analytic_intervals',summary,focus)
    runstep('station QR without bootstrap',rebuild_qr)
    runstep('quality and homogeneity',lambda:dq.run_data_quality_assessment(raw,cfg,WORK))
    runstep('homogeneity exclusion',lambda:hs.run_homogeneity_exclusion_sensitivity(annual,load('data_homogeneity_tests_station_summary'),cfg,WORK))
    assignments,_=cl.run_clustering(cl.build_feature_table(summary,cfg),cfg)
    compare('cluster_assignments',load('cluster_assignments'),assignments)
    reduced,_=cl.run_clustering(cl.build_feature_table(summary,cfg),cfg,feature_cols=cfg['clustering']['robustness_check']['reduced_features'],label_col='cluster_reduced_features')
    compare('cluster_assignments_reduced_features',load('cluster_assignments_reduced_features'),reduced)
    runstep('alternative clustering',lambda:cls.run_alternative_clustering_sensitivity(features,cfg,WORK))
    runstep('spatial inference',lambda:aa.run_spatial_inference(summary,stations,cfg,WORK))
    runstep('geographical associations',lambda:aa.run_driver_analysis(features,stations,cfg,WORK))
    runstep('cluster composites and spatial validation',lambda:aa.run_regionalization_analysis(features,stations,cfg,WORK))
    runstep('fixed baseline, warming associations and internal scores',lambda:cs.run_climate_change_signal_analysis(raw,annual,summary,features,stations,cfg,WORK))
    runstep('raster assignments and regime summaries',lambda:cr.run_climate_regime_analysis(raw,annual,summary,features,stations,cfg,WORK))
    runstep('historical joint rarity from raw data',lambda:cd.run_compound_dry_hot_analysis(raw,stations,cfg,WORK))
    runstep('screened thermal and compound sensitivity',lambda:qc.run_internal_consistency_sensitivity(raw,annual,cfg,WORK))
    # Reproduce sensitivity aggregation from cached station estimates; the alternative ensembles were not saved.
    depth=load('bootstrap_depth_sensitivity_station_comparison')
    for _,r in load('bootstrap_depth_sensitivity_summary').iterrows():
        d=depth.loc[depth.index_name==r.index_name];x=d[r.metric+'_200'];y=d[r.metric+'_400']
        vals=[x.mean(),y.mean(),(y-x).mean(),(y-x).abs().median(),(y-x).abs().max(),x.corr(y)]
        expected=[r.mean_200,r.mean_400,r.mean_difference,r.median_abs_station_difference,r.max_abs_station_difference,r.station_correlation]
        if not np.allclose(vals,expected,equal_nan=True,rtol=1e-7,atol=1e-7):record('depth_'+r.index_name+'_'+r.metric,'MISMATCH',[vals,expected])
    record('bootstrap_depth_scope','CHECKED_AGGREGATION','400-replicate station estimates aggregated; individual alternative draws unavailable, ensemble not rerun.')
    method=load('bootstrap_method_sensitivity_station_level'); result=load('bootstrap_method_sensitivity_summary')
    for _,r in result.iterrows():
        d=method.loc[(method.index_name==r.index_name)&(method.comparison==r.comparison)]
        x=d[r.metric+'_base'];y=d[r.metric+'_alt']; vals=[x.mean(),y.mean(),(x-y).abs().mean(),x.corr(y)]
        expected=[r.mean_base,r.mean_alt,r.mean_abs_diff,r.correlation]
        if not np.allclose(vals,expected,equal_nan=True,rtol=1e-7,atol=1e-7):record('method_'+r.index_name+'_'+r.metric,'MISMATCH',[vals,expected])
    record('bootstrap_method_scope','CHECKED_AGGREGATION','Saved meboot station results checked; alternative ensemble not rerun.')
    interpolation=aa._plot_interpolation_comparison(summary,stations,cfg,cfg['advanced_analyses']['method_sensitivity']['interpolation_methods'],WORK)
    compare('interpolation_method_sensitivity_summary',load('interpolation_method_sensitivity_summary'),interpolation)
    plt.close('all')
    record('completed','DONE','Review mismatches before changing or deleting outputs.')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('phase',choices=['inventory','recompute']);args=parser.parse_args()
    AUDIT.mkdir(parents=True,exist_ok=True)
    warnings.filterwarnings('ignore',category=RuntimeWarning)
    if args.phase=='inventory':inventory()
    else:recompute()
