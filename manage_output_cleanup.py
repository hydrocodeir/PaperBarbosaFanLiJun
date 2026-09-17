"""Plan and archive file-specific output removals; never deletes files itself.

Run before the separately reviewed PowerShell deletion, then run with --finalize.
The initial inventory is immutable evidence of the pre-cleanup state.
"""
from pathlib import Path
import argparse
import hashlib
import json
import zipfile
import pandas as pd

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'outputs'
AUDIT=OUT/'audit_cleanup'
ARCHIVE=ROOT/'archives/output_cleanup_20260917.zip'

def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()

def plan():
    original=pd.read_csv(AUDIT/'inventory_before.csv').fillna('')
    current=set(pd.read_csv(OUT/'publication_v2/figure_manifest.csv').figure)|set(pd.read_csv(OUT/'publication_v2/supplementary_figure_manifest.csv').figure)
    duplicates={'table02_compound_partition.csv','table03_compound_partition.csv','table04_climate_regime_compound.csv','tableS01_joint_sensitivity.csv'}
    rows=[]
    for row in original.to_dict('records'):
        p=Path(row['path']);s=p.as_posix();reason='Numerical source, provenance, or current publication export retained.';action='KEEP'
        if s.startswith('outputs/figures/') or s.startswith('outputs/compound_dry_hot/figures/'):
            action='REMOVE';reason='Legacy graphic replaced by curated main/supplementary panels; numerical evidence retained.'
            if 'station_figures/' in s:reason='Redundant station-by-station export; all coefficients/draws retained and representative profiles provided in Figure S8.'
            if 'fingerprint' in s:reason='Composite attribution interpretation unsupported; excluded with its numerical score tables.'
        elif s.startswith('outputs/output_docs/') and p.name!='README.md':
            action='REMOVE';reason='Obsolete per-output description replaced by the audited supplementary catalog.'
        elif s.startswith('outputs/publication_v2/figures/') and p.stem not in current:
            action='REMOVE';reason='Superseded publication export outside the final main/supplementary manifests.'
        elif s.startswith('outputs/publication_v2/tables/') and p.name in duplicates:
            action='REMOVE';reason='Duplicate or display-only subset; canonical compound partition/sensitivity/regime tables retained.'
        elif s.startswith('outputs/tables/') and 'fingerprint' in p.name:
            action='REMOVE';reason='Overlapping composite score and unsuitable attribution null; archived rather than presented as supplementary evidence.'
        elif p.name in ['run_status.txt','run_status_detail.txt','run_status_summary.txt','REPORT.md']:
            action='REMOVE';reason='Historical run narrative/status superseded by the current audit; run_metadata.json retained.'
        elif s=='outputs/output_docs/README.md':
            action='REPLACE';reason='Replace obsolete index with current manuscript, supplementary, catalog and audit links.'
        if row['status']=='READ_FAILED':reason+=' Original file failed PNG integrity verification.'
        row.update(action=action,reason=reason)
        rows.append(row)
    manifest=pd.DataFrame(rows)
    manifest.to_csv(AUDIT/'curation_manifest.csv',index=False)
    removals=manifest.loc[manifest.action=='REMOVE'].copy()
    removals.to_csv(AUDIT/'removal_manifest.csv',index=False)
    backup=manifest.loc[manifest.action.isin(['REMOVE','REPLACE'])]
    ARCHIVE.parent.mkdir(exist_ok=True)
    if not ARCHIVE.exists():
        with zipfile.ZipFile(ARCHIVE,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=3) as z:
            for r in backup.itertuples():
                p=(ROOT/r.path).resolve()
                assert p.is_relative_to(OUT.resolve())
                assert sha(p)==r.sha256, f'Original changed; review before archive: {p}'
                z.write(p,r.path)
            z.write(AUDIT/'curation_manifest.csv','curation_manifest.csv')
    with zipfile.ZipFile(ARCHIVE) as z:
        for r in backup.itertuples():
            assert hashlib.sha256(z.read(r.path)).hexdigest()==r.sha256,r.path
    summary={'original_files':len(original),'original_bytes':int(original.bytes.sum()),'remove_files':len(removals),'remove_bytes':int(removals.bytes.sum()),'archive':ARCHIVE.relative_to(ROOT).as_posix(),'archive_sha256':sha(ARCHIVE),'archived_files':len(backup),'archive_contents_sha256_verified':True,'deletion_completed':False}
    (AUDIT/'cleanup_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))
    catalog(manifest)
    report(summary)

def role(p):
    n=p.stem
    if 'publication_v2' in p.parts:
        if 'bootstrap' in n:return 'Reproducibility: synchronized bootstrap draws for primary/group intervals.'
        if 'aggregates' in n:return 'Supplementary data: screened annual/seasonal aggregates and calendar coverage (S1).'
        if 'regime' in n:return 'Main/S8 evidence: climate-regime summaries or four-component partitions.'
        if 'quantile' in n or 'thermal' in n:return 'Main evidence: independently recomputed network thermal profiles/trends.'
        return 'Main/S1–S2 evidence: fixed-threshold event definitions, station partitions or sensitivity scenarios.'
    if 'compound_dry_hot' in p.parts:return 'Historical definition sensitivity (S10); varying network and empirical joint rarity, not fixed AND events or causal drivers.'
    if 'warming' in n or 'temperature_anomaly' in n:return 'Exploratory association (S9); shared observations/trends, not external-forcing attribution. Station fits were reproduced with the historical iteration limit; no independent station-level convergence certification.'
    if 'emergence' in n:return 'Precision diagnostics (S4); legacy filename does not denote an emergence date or calibrated detection.'
    if 'interpolation' in n:return 'Reproducibility only: interpolation comparison; no interpolated surface is used as inferential evidence.'
    if 'bootstrap_depth' in n or 'bootstrap_method' in n:return 'Sensitivity (S2): saved alternative station estimates/aggregation checked; alternative random ensembles not regenerated.'
    if 'bootstrap_distributions' in n:return 'Reproducibility: all 99,200 saved station bootstrap draws; all summary columns recalculated.'
    if 'homogeneity' in n or 'quality' in n or 'consistency' in n:return 'Quality/sensitivity (S1–S3); diagnostics do not constitute station homogenization.'
    if 'significance' in n:return 'Spatial diagnostics (S5): local analytic FDR; missing tail tests remain NA, not zero.'
    if 'spatial' in n:return 'Exploratory spatial diagnostics (S5/S7); nominal permutation probabilities.'
    if 'koppen' in n or 'climate_regime' in n:return 'Climate classification and group diagnostics; descriptive partitions with multiplicity/size limitations.'
    if 'cluster' in n or 'representative' in n:return 'Exploratory regionalization (S6–S8): features, assignments, stability or representative selection.'
    if 'driver' in n:return 'Geographical association (S11); legacy driver terminology does not imply causality.'
    if 'fixed_baseline' in n:return 'Main fixed-baseline sensitivity: annual indices and station/network trends or period contrasts.'
    if 'qr_' in n or 'publication_summary' in n:return 'Thermal evidence (main/S3–S4/S8): station quantiles, slopes, intervals and feature summaries.'
    if 'annual_extreme' in n:return 'Reproducibility: full-record thermal annual station indices rebuilt from daily data.'
    return 'Reproducibility metadata; see audit for verification scope.'

def catalog(manifest):
    removed=set(manifest.loc[manifest.action=='REMOVE','path'])
    lines=['# Supplementary data catalog','', 'This catalog distinguishes numerical evidence from display exports. Files remain in their canonical locations so scripts can reuse them without duplicate copies. The supplement embeds eight selected tables and eleven figures; the CSV files below provide complete station-level and sensitivity evidence. Raw observations remain in `data/` and are not redistributed by this catalog.','', 'Verification means reproducibility from the stated inputs, not proof that observations are error-free or homogenized. See [the audit](Output_Audit_2026.md) for bootstrap scope and limitations.','']
    for folder,title in [('publication_v2/tables','Current publication evidence'),('tables','Historical thermal and diagnostic evidence'),('compound_dry_hot/tables','Historical event-definition sensitivity')]:
        lines+=['## '+title,'','| Dataset | Rows × columns | Role and interpretation |','| --- | --- | --- |']
        for p in sorted((OUT/folder).glob('*.csv')):
            if p.relative_to(ROOT).as_posix() in removed:continue
            d=pd.read_csv(p);lines.append(f'| [{p.name}](../{p.relative_to(ROOT).as_posix()}) | {len(d):,} × {len(d.columns)} | {role(p)} |')
        lines+=['']
    lines+=['## Provenance and validation','', '- [Numerical audit](../outputs/audit_cleanup/numerical_checks.json), [dependency checks](../check_output_dependencies.py), and [independent network warming solver check](../outputs/audit_cleanup/network_warming_solver_check.csv).','- [Raw-index and compound validation](../outputs/publication_v2/validation.json).','- [Source hashes for supplementary figures](../outputs/publication_v2/supplementary_figure_sources.json).','- [Original run metadata](../outputs/run_metadata.json) and [publication configuration](../publication_config.yaml).','- [Reference audit](Reference_Audit_2026.md) and [curated reference metadata](verified_references.json).','- [File-by-file curation decisions](../outputs/audit_cleanup/curation_manifest.csv) and [cleanup summary](../outputs/audit_cleanup/cleanup_summary.json).','', 'Archived historical reports may contain obsolete paths and claims. The current manuscript, supplement and this catalog define the active publication package. Running the full legacy pipeline can recreate superseded exports; rerun curation review before treating those exports as publication evidence.','']
    (ROOT/'reports/Supplementary_Data_Catalog.md').write_text('\n'.join(lines),encoding='utf-8')

def report(s):
    checks=json.loads((AUDIT/'numerical_checks.json').read_text())
    assert not any(c['status']=='MISMATCH' for c in checks)
    text=f'''# Output audit and curation — 17 September 2026

## Inventory and decision

All {s['original_files']} original files ({s['original_bytes']/1024**2:.1f} MiB) were inventoried, hashed and inspected for file readability. CSV schemas, missing values, exact row duplicates and numerical infinities were recorded. There were 704 readable files and one PNG with an invalid IDAT checksum. Every original file has a decision and reason in the [curation manifest](../outputs/audit_cleanup/curation_manifest.csv).

{s['remove_files']} obsolete, redundant or excluded files ({s['remove_bytes']/1024**2:.1f} MiB) were selected for removal from the active tree. An additional old output index was backed up before replacement. The archive is [output_cleanup_20260917.zip](../{s['archive']}); every archived source file was verified against its original SHA-256. The archive's SHA-256 is `{s['archive_sha256']}`. Deletion completed: **{s['deletion_completed']}**. The archive preserves recovery, including the originally corrupt image; it does not repair that image.

## Numerical verification

- All 40,176 station quantile estimates and 496 focal station–index summaries were rerun from the annual indices and matched the stored estimates.
- All 32 stored bootstrap-summary fields across 496 station–index groups were recalculated from 99,200 saved draws; maximum discrepancy was floating-point rounding. Saved 200-replicate draws were not generated again.
- Quality and detrended homogeneity diagnostics, exclusions, baseline/alternative clustering, representative selection, spatial tests, fixed-baseline indices, warming associations, climate-raster assignments and group summaries, historical joint-rarity analysis and screening sensitivities were reproduced in an isolated work directory. See the {len(checks)} [numerical records](../outputs/audit_cleanup/numerical_checks.json) for exact table coverage.
- The compound extension's eight primary estimates and pointwise/family intervals were independently recalculated from station values and saved synchronized bootstrap draws. All 48 climate-group intervals, eight group-weighted network closures and 80 sensitivity component means/sample sizes were checked.
- Both thermal annual index definitions were independently rebuilt from daily observations by [validate_publication.py](../validate_publication.py); its nine checks passed. This validation's historical bootstrap/clustering scope statement pertains to that script alone; the separate audit additionally reran clustering.
- The 12 network warming-response slopes displayed in S9 were also solved independently by linear programming; the largest absolute difference was below 0.000005 days per degree Celsius. Station-level warming fits were reproduced with the historical iteration limit; some iterative fits emitted convergence warnings and their convergence is not independently certified.

No unresolved numerical mismatch was found in these checks. This conclusion establishes consistency with recorded inputs and implementations, not observational truth or the validity of every historical statistical interpretation. Alternative 400-replicate and maximum-entropy bootstrap draws were unavailable: only their saved station results and aggregate summaries were checked, not new ensembles. The sensitivity-scenario intervals were checked for ordering; this audit did not independently rerun all ten 4,999-draw extension scenarios.

## Corrections and scientific selection

1. Tail analytic probabilities/intervals are absent in the stored station fits. The curated FDR display uses **NA**, replacing the misleading zero produced by summing missing values. Median retained counts are 115, 104, 97 and 86. Bootstrap tail intervals remain separately available.
2. “Signal emergence” maps were relabeled as descriptive bootstrap precision. The ratio is not an emergence date or calibrated detection probability.
3. Historical composite fingerprint outputs were removed from active evidence: overlapping components and independent shifts of related indices do not provide an external-forcing attribution test. Five numerical score tables and their graphics remain recoverable in the archive.
4. Interpolated surfaces and hundreds of repeated station exports were replaced by station-point panels and selected, reproducible representative profiles. All underlying retained station quantiles and bootstrap draws remain available.
5. Historical joint-rarity outputs are explicitly separated from fixed marginal AND events. Their varying station network, available-row denominator and inverse-probability definition are disclosed; “return period” and causal “driver” interpretations are not carried into the supplement.
6. Duplicate display tables were removed; canonical partition, sensitivity and regime tables supply the manuscript and supplement directly. Climate classification is retained in main Figures 5/9, main tables and supplementary Table S8.
7. One broken PNG (`40790_robat_e_poshtebadam_figure4.png`) was excluded with the superseded station exports.

All nonstation legacy figure families were visually inspected using contact sheets, and the complete curated supplementary set was visually reviewed. All original station images underwent integrity checks and their shared numerical sources were checked; this was not a separate full-resolution visual inspection of every station image. The final set has **10 main figures, 11 supplementary figures, 5 main tables and 8 supplementary tables**, plus the linked [machine-readable data catalog](Supplementary_Data_Catalog.md). Each figure is available as PDF, SVG, PNG and TIFF.

## Reproducibility and recovery

Use `python audit_output_data.py recompute`, then `python check_output_dependencies.py` and `python check_warming_solver.py` for the numerical audit. Use `python build_supplementary.py` and `python build_publication_docs.py` to rebuild the curated supplement and atlases. The initial inventory is immutable; do not rerun the inventory phase over this cleaned release. Historical score comparisons after deletion require restoring their archived inputs into a separate review copy.

The recovery ZIP stores project-relative paths. Extract selected files into a separate folder for inspection before choosing to restore them; a blind extraction over the current output tree would reintroduce superseded material. The raw `data/`, reference `assets/` and original manuscript were not modified by cleanup. The archive remains outside `outputs/` and therefore reduces clutter, not total project disk usage.

Data units, source-data and boundary redistribution permissions, station metadata/homogenization and submission-specific requirements remain author responsibilities as described in the manuscript. These limitations are not resolved by a successful numerical audit.
'''
    (ROOT/'reports/Output_Audit_2026.md').write_text(text,encoding='utf-8')

def finalize():
    summary=json.loads((AUDIT/'cleanup_summary.json').read_text())
    manifest=pd.read_csv(AUDIT/'curation_manifest.csv')
    removed=manifest.loc[manifest.action=='REMOVE']
    assert all(not (ROOT/p).exists() for p in removed.path)
    assert sha(ARCHIVE)==summary['archive_sha256']
    summary['deletion_completed']=True
    active=[p for p in OUT.rglob('*') if p.is_file() and AUDIT not in p.parents]
    summary['active_files_excluding_audit']=len(active)
    summary['active_bytes_excluding_audit']=sum(p.stat().st_size for p in active)
    (AUDIT/'cleanup_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    pd.DataFrame([{'path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(active)]).to_csv(AUDIT/'inventory_after.csv',index=False)
    catalog(manifest);report(summary)
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--finalize',action='store_true');args=parser.parse_args()
    finalize() if args.finalize else plan()
