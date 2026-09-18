"""Rebuild balanced-network thermal inference without rerunning compound models."""
from pathlib import Path
import yaml
from src.paper_pipeline.thermal_network import build_thermal_network

if __name__=='__main__':
    root=Path(__file__).resolve().parent
    cfg=yaml.safe_load((root/'publication_config.yaml').read_text(encoding='utf-8'))
    build_thermal_network(root,root/cfg['output_dir'],cfg)
