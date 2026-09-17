"""Render the curated supplementary figure set from audited numerical outputs."""
from pathlib import Path
from src.paper_pipeline.supplementary_curator import create_supplementary_figures

if __name__=='__main__':
    root=Path(__file__).resolve().parent
    result=create_supplementary_figures(root,root/'outputs/publication_v2')
    print(f'Created {len(result)} supplementary figures in four formats.')
