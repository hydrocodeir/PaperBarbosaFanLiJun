"""Assemble references, numerical supplements, evidence tables and figure atlas."""
from pathlib import Path
import re
import json
import pandas as pd
from pypdf import PdfWriter

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs/publication_v2"
REPORTS = ROOT / "reports"


def markdown_table(frame):
    # Avoid making the optional tabulate package a runtime dependency.
    cols = list(frame.columns)
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
    for row in frame.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(f"{v:.2f}" if isinstance(v, float) else str(v) for v in row) + " |")
    return "\n".join(lines)


def main():
    manuscript_path = REPORTS / "Manuscript_Q1_2026.md"
    manuscript = manuscript_path.read_text(encoding="utf-8").split("## References")[0]
    refs = json.loads((REPORTS / "verified_references.json").read_text(encoding="utf-8"))
    from src.paper_pipeline.publication_regimes import REGIMES, LABELS
    thermal = pd.read_csv(OUT / "tables/table02_climate_regime_thermal.csv").set_index("climate_regime").loc[REGIMES]
    t2 = thermal[["n_stations", "mean_elevation_m", "warm_days", "warm_nights", "cool_days", "cool_nights"]].copy()
    t2.insert(0,"Regime",LABELS)
    t2.columns = ["Regime", "N", "Mean elevation (m)", "Warm days", "Warm nights", "Cool days", "Cool nights"]
    compound = pd.read_csv(OUT / "tables/climate_regime_compound_partition.csv")
    compound = compound.loc[compound.definition == "warm_season"]
    rows = []
    for regime,label in zip(REGIMES,LABELS):
        data=compound.loc[compound.climate_regime==regime].set_index("component")
        joint=data.loc["joint_change"]
        rows.append([label,int(joint.n_stations),int(joint.n_zero_dry_threshold),
            f"{joint.estimate_pp:.2f} [{joint.ci_low_pp:.2f}, {joint.ci_high_pp:.2f}]",
            data.loc["dry_marginal","estimate_pp"],data.loc["hot_marginal","estimate_pp"],data.loc["excess_joint","estimate_pp"]])
    t4=pd.DataFrame(rows,columns=["Regime","N","N₀","Joint change [95% interval]","Dry-frequency term","Hot-frequency term","Excess-joint term"])
    for tag,frame in [("REGIME_THERMAL_TABLE",t2),("REGIME_COMPOUND_TABLE",t4)]:
        pattern=rf"<!-- {tag} -->.*?<!-- END_{tag} -->"
        manuscript,count=re.subn(pattern,f"<!-- {tag} -->\n{markdown_table(frame)}\n<!-- END_{tag} -->",manuscript,flags=re.S)
        assert count==1,tag
    intervals = pd.read_csv(OUT / "tables/thermal_network_intervals.csv")
    metadata = json.loads((OUT / "thermal_network_metadata.json").read_text())
    block = metadata["settings"]["primary_block_length"]
    intervals = intervals.loc[(intervals.network == "common_fixed") & (intervals.block_length == block)]
    labels = {"warm_days": "Warm days", "warm_nights": "Warm nights", "cool_days": "Cool days", "cool_nights": "Cool nights"}
    metrics = ["OLS", "q10", "q50", "q90", "Delta1"]
    records, display = [], []
    for idx, label in labels.items():
        d = intervals.loc[intervals.index_name == idx].set_index("metric")
        record = {"index_name": idx, "n_stations": metadata["common_station_count"]}
        cells = [label]
        for metric in metrics:
            r = d.loc[metric]
            for key in ["estimate", "ci_low", "ci_high"]:
                record[metric + "_" + key] = r[key]
            cells.append(f"{r.estimate:.2f} [{r.ci_low:.2f}, {r.ci_high:.2f}]")
        records.append(record)
        display.append(cells)
    t1 = pd.DataFrame(records)
    t1.to_csv(OUT / "tables/table01_thermal_trends.csv", index=False)
    frame = pd.DataFrame(display, columns=["Index", "OLS [95% interval]", "q10 [95% interval]", "q50 [95% interval]", "q90 [95% interval]", "Δ₁ [95% interval]"])
    manuscript, count = re.subn(r"<!-- THERMAL_NETWORK_TABLE -->.*?<!-- END_THERMAL_NETWORK_TABLE -->",
        "<!-- THERMAL_NETWORK_TABLE -->\n" + markdown_table(frame) + "\n<!-- END_THERMAL_NETWORK_TABLE -->", manuscript, flags=re.S)
    assert count == 1
    manuscript += "## References\n\n" + "\n\n".join(refs) + "\n"
    manuscript_path.write_text(manuscript, encoding="utf-8")
    primary = pd.read_csv(OUT / "tables/compound_partition_primary.csv")
    scenarios = pd.read_csv(OUT / "tables/compound_partition_sensitivity.csv")
    joint = scenarios.loc[scenarios.component == "joint_change", ["definition", "scenario", "n_stations", "estimate_pp", "ci_low_pp", "ci_high_pp"]]
    from src.paper_pipeline.supplementary_curator import build_supplementary_document
    report = build_supplementary_document(ROOT, OUT, primary, scenarios)
    (REPORTS / "Supplementary_Q1_2026.md").write_text(report, encoding="utf-8")
    manifest = pd.read_csv(OUT / "figure_manifest.csv")
    writer = PdfWriter()
    for name in manifest.figure:
        writer.append(OUT / f"figures/{name}.pdf")
    writer.write(OUT / "Figure_Atlas.pdf")
    writer.close()
    supplementary = pd.read_csv(OUT / "supplementary_figure_manifest.csv")
    writer = PdfWriter()
    for name in supplementary.figure:
        writer.append(OUT / f"figures/{name}.pdf")
    writer.write(OUT / "Supplementary_Figure_Atlas.pdf")
    writer.close()

    # Main figure links and numerical table entries are machine checked.
    checks = {}
    for path in [manuscript_path, REPORTS / "Supplementary_Q1_2026.md"]:
        for target in re.findall(r"\]\(([^)]+)\)", path.read_text(encoding="utf-8")):
            if not target.startswith("http"):
                assert (path.parent / target).exists(), (path, target)
        checks[f"{path.name}_links"] = "PASS"
    checks["reference_count"] = len(refs)
    checks["supplementary_figure_count"] = len(supplementary)
    figures = [int(x) for x in re.findall(r"\*Figure (\d+)\.", manuscript)]
    assert figures == list(range(1,11)), figures
    checks["main_figure_count"] = len(figures)
    tables = [int(x) for x in re.findall(r"\*\*Table (\d+)\.",manuscript)]
    assert tables == list(range(1,6)),tables
    checks["main_table_count"] = len(tables)
    checks["regime_tables_generated_from_CSV"] = "PASS"
    checks["conclusion_word_count"] = len(manuscript.split("## 5. Conclusions")[1].split("## Data and code")[0].split())
    body=manuscript.split("## References")[0]
    for ref in refs:
        author=ref.split()[0]
        year=re.search(r"\((\d{4}[ab]?)\)",ref).group(1)
        assert re.search(re.escape(author)+r"[^;\n]{0,65}"+year,body), (author,year)
    checks["all_references_cited_in_body"] = "PASS"
    for row in records:
        line = next(line for line in manuscript.splitlines() if line.startswith(f"| {labels[row['index_name']]} |"))
        cells = [c.strip().replace("−", "-") for c in line.split("|")[2:-1]]
        assert len(cells) == 5
        for metric, cell in zip(metrics, cells):
            numbers = [float(n) for n in re.findall(r"-?\d+\.\d+", cell)]
            expected = [row[metric + "_" + k] for k in ["estimate", "ci_low", "ci_high"]]
            assert len(numbers) == 3 and all(abs(a-b) <= .0051 for a,b in zip(numbers, expected))
    checks["table1_fixed_network_coefficients_and_intervals"] = "PASS"
    checks["thermal_network_station_count"] = metadata["common_station_count"]
    for definition, label in [("annual", "Annual"), ("warm_season", "June–September")]:
        line = next(line for line in manuscript.splitlines() if line.startswith(f"| {label} |"))
        cells = [c.strip().replace("−", "-") for c in line.split("|")[2:-1]]
        rows = primary.loc[primary.definition == definition]
        assert int(cells[0]) == rows.n_stations.iloc[0]
        for cell, (_, row) in zip(cells[1:], rows.iterrows()):
            numbers = [float(n) for n in re.findall(r"-?\d+\.\d+", cell)]
            assert all(abs(a-b) <= .0051 for a,b in zip(numbers, [row.estimate_pp,row.ci_low_pp,row.ci_high_pp]))
    checks["table3_all_components_and_intervals"] = "PASS"
    construction = pd.read_csv(OUT / "tables/index_definition_period_summary.csv")
    for case in ["fixed_w11_t7_raw", "fixed_w11_t7_corrected", "fixed_w5_t7_corrected", "fixed_w5_t8_corrected"]:
        for idx in ["warm_days", "warm_nights"]:
            value = construction.loc[(construction.scenario == case) & (construction.index_name == idx), "change_days"].item()
            assert f"{value:.2f}" in manuscript, (case, idx, value)
    checks["matched_network_uncorrected_and_corrected_headline_values"] = "PASS"
    assert "(75/103) × 17.647059" in manuscript and "no additional evidence" in manuscript
    checks["zero_cutoff_reweighting_explicitly_separated"] = "PASS"
    (OUT / "document_validation.json").write_text(json.dumps(checks, indent=2), encoding="utf-8")
    from manage_output_cleanup import catalog
    catalog(pd.read_csv(ROOT / "outputs/audit_cleanup/curation_manifest.csv"))
    print(json.dumps(checks, indent=2))

if __name__ == "__main__":
    main()
