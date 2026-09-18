"""Build and validate the International Journal of Climatology submission files."""

from __future__ import annotations

import re
import sys
import zipfile
from pathlib import Path

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK, WD_COLOR_INDEX
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt


PACKAGE = Path(__file__).resolve().parents[1]
ROOT = PACKAGE.parent
sys.path.insert(0, str(ROOT))

from export_manuscript_word import (  # noqa: E402
    add_page_number,
    add_text_with_urls,
    build_document,
    set_run_font,
)


def set_page_break_before(paragraph) -> None:
    p_pr = paragraph._p.get_or_add_pPr()
    element = OxmlElement("w:pageBreakBefore")
    element.set(qn("w:val"), "true")
    p_pr.append(element)


def highlight_placeholders(document: Document) -> int:
    count = 0
    pattern = re.compile(r"\[[A-Z][A-Z0-9 /,;:()–—'\-]+\]")
    marker_terms = (
        "[AUTHOR",
        "[Department",
        "[CORRESPONDING",
        "[full postal",
        "[EMAIL]",
        "[ORCID]",
        "[Email]",
        "[UNITS",
        "[ACCESS",
        "[LICENSE",
        "[ACKNOWLEDGE",
        "[INSERT",
        "[AUTHORS TO CONFIRM",
        "[AUTHORS:",
        "[SUBMISSION DATE]",
        "[Affiliation]",
        "[Postal address]",
        "[NONE / DOI OR URL]",
    )
    for paragraph in list(document.paragraphs) + [
        p for table in document.tables for row in table.rows for cell in row.cells for p in cell.paragraphs
    ]:
        marked_paragraph = any(term in paragraph.text for term in marker_terms)
        for run in paragraph.runs:
            if marked_paragraph or pattern.search(run.text):
                run.font.highlight_color = WD_COLOR_INDEX.YELLOW
        if marked_paragraph:
            count += 1
        else:
            count += sum(len(pattern.findall(run.text)) for run in paragraph.runs)
    return count


def finish_export(path: Path, *, title: str, header_text: str, page_break_heading: str | None = None) -> int:
    document = Document(path)
    document.core_properties.title = title
    document.core_properties.subject = "International Journal of Climatology submission"
    section = document.sections[0]
    header = section.header.paragraphs[0]
    header.text = header_text
    header.alignment = WD_ALIGN_PARAGRAPH.CENTER
    for run in header.runs:
        set_run_font(run, "Times New Roman", 9)
        run.italic = True
    footer = section.footer.paragraphs[0]
    if not footer._p.xpath('.//w:fldChar'):
        add_page_number(footer)
    if page_break_heading:
        for paragraph in document.paragraphs:
            if paragraph.text.strip() == page_break_heading:
                set_page_break_before(paragraph)
                break
    count = highlight_placeholders(document)
    document.save(path)
    return count


def configure_simple_document(document: Document, title: str, header_text: str) -> None:
    section = document.sections[0]
    section.page_width = Cm(21)
    section.page_height = Cm(29.7)
    section.top_margin = Cm(2.5)
    section.bottom_margin = Cm(2.5)
    section.left_margin = Cm(2.5)
    section.right_margin = Cm(2.5)
    normal = document.styles["Normal"]
    normal.font.name = "Times New Roman"
    normal.font.size = Pt(11)
    normal._element.rPr.rFonts.set(qn("w:ascii"), "Times New Roman")
    normal._element.rPr.rFonts.set(qn("w:hAnsi"), "Times New Roman")
    normal.paragraph_format.space_after = Pt(6)
    normal.paragraph_format.line_spacing = 1.15
    for style_name, size in (("Title", 16), ("Heading 1", 13), ("Heading 2", 11)):
        style = document.styles[style_name]
        style.font.name = "Times New Roman"
        style.font.size = Pt(size)
        style.font.bold = True
        style.font.color.rgb = None
    header = section.header.paragraphs[0]
    header.text = header_text
    header.alignment = WD_ALIGN_PARAGRAPH.CENTER
    for run in header.runs:
        set_run_font(run, "Times New Roman", 9)
        run.italic = True
    add_page_number(section.footer.paragraphs[0])
    document.core_properties.title = title
    document.core_properties.subject = "International Journal of Climatology submission"


def strip_inline_markdown(text: str) -> str:
    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r"\1 (\2)", text)
    return text.replace("**", "").replace("*", "").replace("`", "")


def build_simple_markdown(source: Path, destination: Path, *, header_text: str) -> int:
    document = Document()
    title = source.stem.replace("_", " ")
    configure_simple_document(document, title, header_text)
    first_title = True
    for raw in source.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("# "):
            paragraph = document.add_paragraph(style="Title" if first_title else "Heading 1")
            first_title = False
            paragraph.add_run(strip_inline_markdown(line[2:]))
        elif line.startswith("## "):
            paragraph = document.add_paragraph(style="Heading 1")
            paragraph.add_run(strip_inline_markdown(line[3:]))
        elif line.startswith("### "):
            paragraph = document.add_paragraph(style="Heading 2")
            paragraph.add_run(strip_inline_markdown(line[4:]))
        elif re.match(r"^- \[[ xX]\] ", line):
            checked = line[3].lower() == "x"
            text = strip_inline_markdown(line[6:])
            paragraph = document.add_paragraph(style="List Bullet")
            paragraph.add_run(("☒ " if checked else "☐ ") + text)
        elif re.match(r"^\d+\. ", line):
            paragraph = document.add_paragraph(style="List Number")
            paragraph.add_run(strip_inline_markdown(re.sub(r"^\d+\. ", "", line)))
        elif line.startswith("- "):
            paragraph = document.add_paragraph(style="List Bullet")
            paragraph.add_run(strip_inline_markdown(line[2:]))
        else:
            paragraph = document.add_paragraph()
            add_text_with_urls(paragraph, strip_inline_markdown(line))
    placeholders = highlight_placeholders(document)
    document.save(destination)
    return placeholders


def word_count_for_journal(markdown: str) -> int:
    body = re.split(r"(?m)^## References\s*$", markdown)[0]
    body = re.sub(r"(?m)^\*(?:Figure|Table)\s+\d+\..*\*\s*$", " ", body)
    body = re.sub(r"(?m)^#+\s*", " ", body)
    body = re.sub(r"!\[[^\]]*\]\([^)]+\)", " ", body)
    body = re.sub(r"<!--.*?-->", " ", body, flags=re.S)
    return len(body.split())


def validate_docx(path: Path) -> tuple[int, int, int]:
    with zipfile.ZipFile(path) as archive:
        bad = archive.testzip()
        if bad:
            raise RuntimeError(f"Corrupt DOCX member in {path.name}: {bad}")
    document = Document(path)
    return len(document.paragraphs), len(document.tables), len(document.inline_shapes)


def main() -> None:
    manuscript_md = PACKAGE / "Manuscript_IJC_2026.md"
    supplement_md = PACKAGE / "Supporting_Information_IJC_2026.md"
    manuscript_docx = PACKAGE / "Manuscript_IJC_2026.docx"
    supplement_docx = PACKAGE / "Supporting_Information_IJC_2026.docx"
    cover_docx = PACKAGE / "Cover_Letter_IJC_2026.docx"
    checklist_docx = PACKAGE / "Submission_Checklist_IJC_2026.docx"

    build_document(manuscript_md, manuscript_docx)
    build_document(supplement_md, supplement_docx)
    placeholder_counts = {
        manuscript_docx.name: finish_export(
            manuscript_docx,
            title="Distributional thermal change and the components of increasing dry–hot concurrence across Iran, 1991–2024",
            header_text="Thermal extremes and dry–hot concurrence in Iran",
            page_break_heading="Abstract",
        ),
        supplement_docx.name: finish_export(
            supplement_docx,
            title="Supporting Information: thermal extremes and explicit dry–hot concurrence",
            header_text="Supporting Information | Thermal extremes and dry–hot concurrence",
        ),
        cover_docx.name: build_simple_markdown(
            PACKAGE / "Cover_Letter_IJC_2026.md",
            cover_docx,
            header_text="Cover letter | International Journal of Climatology",
        ),
        checklist_docx.name: build_simple_markdown(
            PACKAGE / "Submission_Checklist_IJC_2026.md",
            checklist_docx,
            header_text="Submission checklist | International Journal of Climatology",
        ),
    }

    manuscript = manuscript_md.read_text(encoding="utf-8")
    supplement = supplement_md.read_text(encoding="utf-8")
    abstract_match = re.search(r"(?ms)^## Abstract\s*(.*?)^\*\*Keywords:", manuscript)
    abstract_words = len(abstract_match.group(1).split()) if abstract_match else -1
    running = re.search(r"\*\*Running title:\*\*\s*(.+)", manuscript).group(1).strip()
    reference_count = len(re.findall(r"(?m)^[A-Z][^\n]+\(\d{4}[a-z]?\)\.", manuscript))
    main_figures = len(re.findall(r"(?m)^\*Figure \d+\.", manuscript))
    main_tables = len(re.findall(r"(?m)^\*\*Table \d+\.", manuscript))
    supp_figures = len(re.findall(r"(?m)^### Figure S\d+\.", supplement))
    supp_tables = len(re.findall(r"(?m)^### Table S\d+\.", supplement))
    new_refs = {
        "Donat et al. (2014)": "Donat et al., 2014",
        "Rahimzadeh et al. (2009)": "Rahimzadeh et al., 2009",
        "Rahimzadeh and Nassaji Zavareh (2014)": "Rahimzadeh and Nassaji Zavareh, 2014",
        "Yosef et al. (2021)": "Yosef et al., 2021",
        "Wu et al. (2021)": "Wu et al., 2021",
        "Najafi et al. (2025)": "Najafi et al., 2025",
    }
    body = re.split(r"(?m)^## References\s*$", manuscript)[0]
    new_ref_checks = {name: (token in body) for name, token in new_refs.items()}
    reference_text = re.split(r"(?m)^## References\s*$", manuscript)[1]
    reference_blocks = [block.strip() for block in re.split(r"\n\s*\n", reference_text) if block.strip()]
    missing_body_citations: list[str] = []
    for reference in reference_blocks:
        match = re.match(r"([^\s]+).*?\((\d{4}[a-z]?)\)", reference)
        if not match:
            continue
        author, year = match.groups()
        if not re.search(re.escape(author) + r".{0,100}?" + re.escape(year), body, flags=re.S):
            missing_body_citations.append(f"{author} ({year})")
    docx_stats = {path.name: validate_docx(path) for path in (manuscript_docx, supplement_docx, cover_docx, checklist_docx)}
    word_count = word_count_for_journal(manuscript)

    conditions = {
        "main_word_count_at_or_below_7500": word_count <= 7500,
        "abstract_at_or_below_300": 0 < abstract_words <= 300,
        "running_title_below_70_characters": len(running) < 70,
        "main_figures_1_to_10_present": main_figures == 10,
        "main_tables_1_to_4_present": main_tables == 4,
        "supplement_figures_S1_to_S15_present": supp_figures == 15,
        "supplement_tables_S1_to_S15_present": supp_tables == 15,
        "six_new_IJC_references_cited": all(new_ref_checks.values()),
        "all_reference_entries_cited_in_main_text": not missing_body_citations,
        "all_docx_files_valid": all(stats[0] > 0 for stats in docx_stats.values()),
    }
    if not all(conditions.values()):
        raise RuntimeError(f"Submission validation failed: {conditions}")

    lines = [
        "International Journal of Climatology submission validation",
        "==========================================================",
        f"Journal-counted manuscript words (references and figure/table captions excluded; table text retained): {word_count}",
        f"Abstract words: {abstract_words}",
        f"Running-title characters: {len(running)}",
        f"References: {reference_count}",
        f"Main figures / tables: {main_figures} / {main_tables}",
        f"Supporting figures / tables: {supp_figures} / {supp_tables}",
        "",
        "New International Journal of Climatology citations:",
    ]
    lines.extend(f"- {name}: {'PASS' if passed else 'FAIL'}" for name, passed in new_ref_checks.items())
    lines.append(f"- uncited reference entries: {', '.join(missing_body_citations) if missing_body_citations else 'none'}")
    lines.extend(["", "DOCX structure (paragraphs, tables, embedded figures):"])
    lines.extend(f"- {name}: {stats}" for name, stats in docx_stats.items())
    lines.extend(["", "Placeholder counts highlighted in yellow:"])
    lines.extend(f"- {name}: {count}" for name, count in placeholder_counts.items())
    lines.extend(["", "Automated checks:"])
    lines.extend(f"- {name}: {'PASS' if passed else 'FAIL'}" for name, passed in conditions.items())
    lines.extend(
        [
            "",
            "AUTHOR ACTION REQUIRED: highlighted placeholders contain factual information absent from the repository.",
            "Complete them and rerun this builder before submission.",
        ]
    )
    (PACKAGE / "Validation_Report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
