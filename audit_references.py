"""Read DOI registration metadata; flag bibliographic mismatches for review.

Metadata verification is separate from verification of a paper's conclusions.
Raw Crossref records are retained; failures are reported, never treated as passes.
"""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from difflib import SequenceMatcher
from pathlib import Path
import json
import re
import unicodedata
from urllib.request import Request, urlopen
from urllib.parse import quote
import pandas as pd

ROOT = Path(__file__).resolve().parent
NEW_DOIS = ["10.1029/2025EA004860", "10.1016/j.jaridenv.2026.105606",
            "10.1038/s41597-023-02549-6", "10.1175/JCLI3366.1"]


def normalize(text):
    return re.sub(r"[^a-z0-9]", "", unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode().lower())


def retrieve(doi):
    url = "https://api.crossref.org/works/" + quote(doi, safe="/")
    try:
        request = Request(url, headers={"User-Agent": "ClimateManuscriptReferenceAudit/1.0"})
        with urlopen(request, timeout=45) as response:
            return doi, json.load(response)["message"], None
    except Exception as exc:
        return doi, {}, str(exc)


def main():
    text = (ROOT / "reports/Manuscript_Q1_2026.md").read_text(encoding="utf-8")
    refs = text.split("## References", 1)[1].strip().split("\n\n")
    mapping = {re.search(r"https://doi.org/(\S+)", r).group(1).rstrip("."): r for r in refs if "https://doi.org/" in r}
    records, rows = {}, []
    cache_path = ROOT / "outputs/publication_v2/references/crossref_records.json"
    cached = json.loads(cache_path.read_text(encoding="utf-8"))["records"] if cache_path.exists() else {}
    def fetch(doi):
        return (doi, cached[doi], None) if cached.get(doi) else retrieve(doi)
    with ThreadPoolExecutor(max_workers=2) as pool:
        for doi, record, error in pool.map(fetch, list(dict.fromkeys(list(mapping) + NEW_DOIS))):
            records[doi] = record
            title = record.get("title", [""])[0]
            reference = mapping.get(doi, "NEW CANDIDATE")
            original_title = reference.split("). ", 1)[-1].split(". *", 1)[0].strip("*")
            if reference.split("). ", 1)[-1].startswith("*"):
                original_title = reference.split("). ", 1)[-1].split("*")[1]
            years = {key: value.get("date-parts", []) for key, value in record.items() if key.startswith("published") or key == "issued"}
            rows.append(dict(doi=doi, status="FAILED" if error else "DOI_METADATA_FOUND",
                             registered_title=title, registered_authors="; ".join(x.get("family", "") for x in record.get("author", [])),
                             registered_journal="; ".join(record.get("container-title", [])),
                             dates=json.dumps(years), volume=record.get("volume", ""), page=record.get("page", record.get("article-number", "")),
                             title_similarity=SequenceMatcher(None,normalize(original_title),normalize(title)).ratio() if doi in mapping else None,
                             manuscript_reference=reference, error=error))
    out = ROOT / "outputs/publication_v2/references"
    out.mkdir(parents=True, exist_ok=True)
    (out / "crossref_records.json").write_text(json.dumps({"checked_utc":datetime.now(timezone.utc).isoformat(), "records": records}, indent=2, ensure_ascii=False), encoding="utf-8")
    pd.DataFrame(rows).to_csv(out / "reference_metadata_audit.csv", index=False)
    print(pd.DataFrame(rows)[["doi", "status", "registered_title", "registered_authors", "title_similarity"]].to_string(index=False))

if __name__ == "__main__":
    main()
