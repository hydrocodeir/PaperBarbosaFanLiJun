from pathlib import Path
import yaml
from src.paper_pipeline.index_definition import build_index_definition

if __name__=='__main__':
    root=Path(__file__).resolve().parent
    build_index_definition(root,root/'outputs/publication_v2',yaml.safe_load((root/'index_definition_config.yaml').read_text()))
