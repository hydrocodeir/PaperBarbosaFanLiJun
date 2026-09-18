from pathlib import Path
import yaml
from src.paper_pipeline.zero_threshold import build_zero_threshold

if __name__=='__main__':
    root=Path(__file__).resolve().parent
    build_zero_threshold(root,root/'outputs/publication_v2',yaml.safe_load((root/'publication_config.yaml').read_text()))
