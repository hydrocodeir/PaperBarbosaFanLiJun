from __future__ import annotations

import argparse

from src.paper_pipeline import run_pipeline


def main():
    parser = argparse.ArgumentParser(description="Quantile regression + bootstrap + clustering pipeline for temperature extremes.")
    parser.add_argument("--config", type=str, default="config.yaml")
    args = parser.parse_args()
    outdir = run_pipeline(args.config)
    print(f"Done. Outputs written to: {outdir}")


if __name__ == "__main__":
    main()
