from pathlib import Path
import subprocess
import markdown

ROOT = Path(__file__).resolve().parents[1]
md_path = ROOT / "outputs" / "IJoC_final_submission_ready.md"
html_path = ROOT / "outputs" / "IJoC_final_submission_ready.html"
pdf_path = ROOT / "outputs" / "IJoC_final_submission_ready.pdf"
chrome_path = Path(r"C:\Program Files\Google\Chrome\Application\chrome.exe")

if not md_path.exists():
    raise FileNotFoundError(f"Markdown file not found: {md_path}")
if not chrome_path.exists():
    raise FileNotFoundError(f"Chrome not found: {chrome_path}")

md_text = md_path.read_text(encoding="utf-8")
body = markdown.markdown(md_text, extensions=["extra", "tables", "fenced_code", "toc"])

css = """
@page { size: A4; margin: 18mm 16mm 18mm 16mm; }
body { font-family: Georgia, 'Times New Roman', serif; line-height: 1.45; color: #111; }
h1, h2, h3, h4 { page-break-after: avoid; }
p, li { orphans: 3; widows: 3; }
img { max-width: 100%; height: auto; }
table { border-collapse: collapse; width: 100%; font-size: 10pt; }
th, td { border: 1px solid #777; padding: 4px 6px; vertical-align: top; }
code, pre { font-family: Consolas, monospace; }
pre { white-space: pre-wrap; }
"""

html = (
    "<!doctype html><html><head><meta charset=\"utf-8\">"
    "<title>IJoC Final Submission</title>"
    f"<style>{css}</style></head><body>{body}</body></html>"
)

html_path.write_text(html, encoding="utf-8")

subprocess.run(
    [
        str(chrome_path),
        "--headless=new",
        "--disable-gpu",
        "--allow-file-access-from-files",
        f"--print-to-pdf={pdf_path}",
        html_path.resolve().as_uri(),
    ],
    check=True,
)

print(f"Generated PDF: {pdf_path}")
