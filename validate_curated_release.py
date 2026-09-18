"""Check publication manifests, file integrity, links and source hashes after cleanup."""
from pathlib import Path
import hashlib
import json
import re
import xml.etree.ElementTree as ET
import pandas as pd
from PIL import Image
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'outputs/publication_v2'
checks={}
stems=[]
for name,count in [('figure_manifest.csv',10),('supplementary_figure_manifest.csv',15)]:
    m=pd.read_csv(OUT/name)
    assert len(m)==count and m.figure.is_unique
    stems.extend(m.figure)
expected={f'{stem}.{ext}' for stem in stems for ext in ['png','pdf','svg','tiff']}
assert {p.name for p in (OUT/'figures').iterdir() if p.is_file()}==expected
for p in (OUT/'figures').iterdir():
    if p.suffix in ['.png','.tiff']:
        with Image.open(p) as im:
            required_dpi=350 if p.suffix=='.png' else 600
            assert min(im.info.get('dpi',(0,0)))>=required_dpi-.1,(p,im.info.get('dpi'))
            assert min(im.size)>500,(p,im.size)
            im.verify()
    elif p.suffix=='.pdf':assert len(PdfReader(p).pages)==1
    elif p.suffix=='.svg':ET.parse(p)
checks['all_100_curated_graphic_files_integrity']='PASS'
for name,count in [('Figure_Atlas.pdf',10),('Supplementary_Figure_Atlas.pdf',15)]:
    assert len(PdfReader(OUT/name).pages)==count
checks['atlas_page_counts']='PASS'
for name in ['figure_source_hashes.json','supplementary_figure_sources.json']:
    for rel,expected_hash in json.loads((OUT/name).read_text()).items():
        with (ROOT/rel).open('rb') as f:actual=hashlib.file_digest(f,'sha256').hexdigest()
        assert actual==expected_hash,(name,rel)
checks['all_recorded_figure_source_hashes']='PASS'
paths=['reports/Manuscript_Q1_2026.md','reports/Supplementary_Q1_2026.md','reports/Supplementary_Data_Catalog.md','reports/Output_Audit_2026.md','reports/Publication_Guide_FA.md','reports/Reference_Audit_2026.md','reports/Q1_Reviewer_Report_2026_FA.md','reports/Thermal_Network_Revision_2026_FA.md','reports/Index_Zero_Revision_2026_FA.md','outputs/README.md']
# The former per-output documentation tree was retired after curation.
if (ROOT/'outputs/output_docs/README.md').exists():
    paths.append('outputs/output_docs/README.md')
for rel in paths:
    p=ROOT/rel
    for target in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
        if not target.startswith(('http','#')):assert (p.parent/target).exists(),(rel,target)
checks['active_document_links']='PASS'
s=(ROOT/'reports/Supplementary_Q1_2026.md').read_text(encoding='utf-8')
assert [int(n) for n in re.findall(r'^### Table S(\d+)\.',s,re.M)]==list(range(1,15))
assert [int(n) for n in re.findall(r'^### Figure S(\d+)\.',s,re.M)]==list(range(1,16))
assert s.count('| NA | 0 |')==8
checks['supplement_numbering_and_missing_tail_tests']='PASS'
catalog=(ROOT/'reports/Supplementary_Data_Catalog.md').read_text(encoding='utf-8')
tables=[p for folder in ['outputs/tables','outputs/compound_dry_hot/tables','outputs/publication_v2/tables'] for p in (ROOT/folder).glob('*.csv')]
for p in tables:
    assert p.relative_to(ROOT).as_posix() in catalog,p
    d=pd.read_csv(p)
    assert not d.empty,p
checks['all_retained_research_tables_cataloged']=len(tables)
summary=json.loads((ROOT/'outputs/audit_cleanup/cleanup_summary.json').read_text())
assert summary['deletion_completed'] and summary['archive_contents_sha256_verified']
removals=pd.read_csv(ROOT/'outputs/audit_cleanup/removal_manifest.csv')
assert all(not (ROOT/p).exists() for p in removals.path)
checks['archived_removals_absent_from_active_tree']=len(removals)
checks['scope']='Integrity, numerical-source provenance, catalog completeness and document consistency; see numerical audit for scientific verification limits.'
(ROOT/'outputs/audit_cleanup/release_validation.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
print(json.dumps(checks,indent=2))
