"""Run the independent manuscript extension without altering historical outputs."""
import argparse
from pathlib import Path
import sys
import yaml

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="publication_config.yaml")
    parser.add_argument("--figures-only", action="store_true")
    args = parser.parse_args()
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = ROOT / config_path
    cfg = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    from paper_pipeline.compound_partition import run_partition
    out = ROOT / cfg["output_dir"]
    if not args.figures_only:
        run_partition(cfg, ROOT, config_path)
        from paper_pipeline.publication_regimes import build_regime_tables
        build_regime_tables(ROOT, out, cfg)
    from paper_pipeline.publication_figures import create_figures
    create_figures(ROOT, out)
    print(f"Publication outputs: {out}", flush=True)

if __name__ == "__main__":
    main()
