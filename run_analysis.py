from __future__ import annotations

import argparse

from src.paper_pipeline import run_pipeline


def main():
    parser = argparse.ArgumentParser(description="Quantile regression + bootstrap + clustering pipeline for temperature extremes.")
    parser.add_argument("--config", type=str, required=True)
    parser.add_argument("--start-phase", type=int, default=1, help="Resume the pipeline from a cached phase (supported: 1, 8, 9, 10).")
    args = parser.parse_args()
    outdir = run_pipeline(args.config, start_phase=args.start_phase)
    print(f"Done. Outputs written to: {outdir}")


if __name__ == "__main__":
    main()
