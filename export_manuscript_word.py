"""Export the Q1 manuscript Markdown file to a formatted Word document.

The converter is intentionally local and deterministic: it embeds every figure,
converts Markdown tables to native Word tables, preserves emphasis and links, and
renders the five displayed equations as editable Word math text.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import mistune
from docx import Document
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Inches, Pt, RGBColor


URL_RE = re.compile(r"(https?://[^\s<>]+)")


def set_run_font(run, name: str, size: float | None = None) -> None:
    run.font.name = name
    if size is not None:
        run.font.size = Pt(size)
    fonts = run._element.get_or_add_rPr().get_or_add_rFonts()
    fonts.set(qn("w:ascii"), name)
    fonts.set(qn("w:hAnsi"), name)
    fonts.set(qn("w:cs"), name)


def set_cell_shading(cell, fill: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = tc_pr.find(qn("w:shd"))
    if shd is None:
        shd = OxmlElement("w:shd")
        tc_pr.append(shd)
    shd.set(qn("w:fill"), fill)


def set_repeat_table_header(row) -> None:
    tr_pr = row._tr.get_or_add_trPr()
    header = OxmlElement("w:tblHeader")
    header.set(qn("w:val"), "true")
    tr_pr.append(header)


def prevent_row_split(row) -> None:
    tr_pr = row._tr.get_or_add_trPr()
    tr_pr.append(OxmlElement("w:cantSplit"))


def add_page_number(paragraph) -> None:
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = paragraph.add_run()
    begin = OxmlElement("w:fldChar")
    begin.set(qn("w:fldCharType"), "begin")
    instruction = OxmlElement("w:instrText")
    instruction.set(qn("xml:space"), "preserve")
    instruction.text = " PAGE "
    separate = OxmlElement("w:fldChar")
    separate.set(qn("w:fldCharType"), "separate")
    end = OxmlElement("w:fldChar")
    end.set(qn("w:fldCharType"), "end")
    run._r.extend((begin, instruction, separate, end))
    set_run_font(run, "Times New Roman", 9)


def add_hyperlink(paragraph, text: str, target: str) -> None:
    relation = paragraph.part.relate_to(
        target,
        "http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink",
        is_external=True,
    )
    hyperlink = OxmlElement("w:hyperlink")
    hyperlink.set(qn("r:id"), relation)
    run = OxmlElement("w:r")
    run_properties = OxmlElement("w:rPr")
    color = OxmlElement("w:color")
    color.set(qn("w:val"), "0563C1")
    underline = OxmlElement("w:u")
    underline.set(qn("w:val"), "single")
    run_properties.extend((color, underline))
    run.append(run_properties)
    value = OxmlElement("w:t")
    value.text = text
    run.append(value)
    hyperlink.append(run)
    paragraph._p.append(hyperlink)


def add_text_with_urls(paragraph, text: str, *, bold=False, italic=False, code=False) -> None:
    position = 0
    for match in URL_RE.finditer(text):
        if match.start() > position:
            add_run(paragraph, text[position : match.start()], bold=bold, italic=italic, code=code)
        url = match.group(1).rstrip(".,;)")
        suffix = match.group(1)[len(url) :]
        add_hyperlink(paragraph, url, url)
        if suffix:
            add_run(paragraph, suffix, bold=bold, italic=italic, code=code)
        position = match.end()
    if position < len(text):
        add_run(paragraph, text[position:], bold=bold, italic=italic, code=code)


def add_run(paragraph, text: str, *, bold=False, italic=False, code=False):
    run = paragraph.add_run(text)
    run.bold = bold
    run.italic = italic
    if code:
        set_run_font(run, "Courier New", 9)
        run.font.color.rgb = RGBColor(70, 70, 70)
    else:
        # Leave the size to the paragraph style so that Title, Heading and
        # Caption styles retain their configured hierarchy.
        set_run_font(run, "Times New Roman")
    return run


def plain_text(nodes: list[dict]) -> str:
    parts: list[str] = []
    for node in nodes:
        node_type = node.get("type")
        if node_type in {"text", "codespan", "inline_math"}:
            parts.append(node.get("raw", "").lstrip("$").rstrip("$"))
        elif node_type == "image":
            parts.append(plain_text(node.get("children", [])))
        else:
            parts.append(plain_text(node.get("children", [])))
    return "".join(parts)


def render_inline(paragraph, nodes: list[dict], *, bold=False, italic=False) -> None:
    for node in nodes:
        node_type = node.get("type")
        if node_type == "text":
            add_text_with_urls(paragraph, node.get("raw", ""), bold=bold, italic=italic)
        elif node_type == "strong":
            render_inline(paragraph, node.get("children", []), bold=True, italic=italic)
        elif node_type == "emphasis":
            render_inline(paragraph, node.get("children", []), bold=bold, italic=True)
        elif node_type == "codespan":
            add_text_with_urls(
                paragraph, node.get("raw", ""), bold=bold, italic=italic, code=True
            )
        elif node_type == "link":
            add_hyperlink(paragraph, plain_text(node.get("children", [])), node["attrs"]["url"])
        elif node_type == "inline_math":
            run = add_run(paragraph, equation_text(node.get("raw", "")), bold=bold, italic=italic)
            set_run_font(run, "Cambria Math", 11)
        elif node_type in {"softbreak", "linebreak"}:
            paragraph.add_run("\n")
        elif node.get("children"):
            render_inline(paragraph, node["children"], bold=bold, italic=italic)


def equation_text(raw: str) -> str:
    """Convert the manuscript's displayed LaTeX to readable Unicode math."""
    source = raw.strip().strip("$")
    equations = {
        r"Q_Y(\tau\mid t)=\beta_0(\tau)+\beta_1(\tau)t.": "Qᵧ(τ | t) = β₀(τ) + β₁(τ)t.",
        r"\Delta_1=\hat\beta_1(0.90)-\hat\beta_1(0.10).": "Δ₁ = β̂₁(0.90) − β̂₁(0.10).",
        r"D_{sy}=\mathbf{1}(P_{sy}<q^{P}_{s,0.25}),\qquad H_{sy}=\mathbf{1}(T_{sy}>q^{T}_{s,0.75}),\qquad J_{sy}=D_{sy}H_{sy},": (
            "Dₛᵧ = 𝟙(Pₛᵧ < qᴾₛ,₀.₂₅),    Hₛᵧ = 𝟙(Tₛᵧ > qᵀₛ,₀.₇₅),    Jₛᵧ = DₛᵧHₛᵧ,"
        ),
        r"\widehat{\Delta j}_{\mathrm{all}}=\frac{N_+}{N}\widehat{\Delta j}_{+},": (
            "Δĵₐₗₗ = (N₊/N)Δĵ₊,"
        ),
        r"j_1-j_0=\underbrace{(d_1-d_0)\frac{h_1+h_0}{2}}_{A_D}+\underbrace{(h_1-h_0)\frac{d_1+d_0}{2}}_{A_H}+\underbrace{(c_1-c_0)}_{A_C}": (
            "j₁ − j₀ = (d₁ − d₀)(h₁ + h₀)/2 [Aᴅ] + (h₁ − h₀)(d₁ + d₀)/2 [Aʜ] + (c₁ − c₀) [Aᴄ]"
        ),
    }
    if source in equations:
        return equations[source]
    replacements = {
        r"\tau": "τ",
        r"\beta": "β",
        r"\Delta": "Δ",
        r"\mid": "|",
        r"\qquad": "    ",
        "-": "−",
    }
    for old, new in replacements.items():
        source = source.replace(old, new)
    source = re.sub(r"\\mathrm\{([^{}]+)\}", r"\1", source)
    source = re.sub(r"\\mathbf\{([^{}]+)\}", r"\1", source)
    source = source.replace("{", "").replace("}", "")
    return source


def add_equation(document: Document, raw: str) -> None:
    paragraph = document.add_paragraph()
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    paragraph.paragraph_format.space_before = Pt(4)
    paragraph.paragraph_format.space_after = Pt(4)
    run = paragraph.add_run(equation_text(raw))
    set_run_font(run, "Cambria Math", 11)


def set_picture_alt_text(inline_shape, description: str) -> None:
    doc_pr = inline_shape._inline.docPr
    doc_pr.set("descr", description)
    doc_pr.set("title", description)


def add_figure(document: Document, manuscript_dir: Path, node: dict) -> None:
    source = (manuscript_dir / node["attrs"]["url"]).resolve()
    if not source.exists():
        raise FileNotFoundError(f"Figure not found: {source}")
    alt_text = plain_text(node.get("children", [])) or source.stem
    paragraph = document.add_paragraph()
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    paragraph.paragraph_format.space_before = Pt(6)
    paragraph.paragraph_format.space_after = Pt(3)
    shape = paragraph.add_run().add_picture(str(source), width=Inches(6.35))
    set_picture_alt_text(shape, alt_text)


def add_table(document: Document, node: dict) -> None:
    table_head = node["children"][0]
    table_body = node["children"][1]
    headers = table_head.get("children", [])
    body_rows = table_body.get("children", [])
    table = document.add_table(rows=1, cols=len(headers))
    table.style = "Table Grid"
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = True

    header_row = table.rows[0]
    for index, cell_node in enumerate(headers):
        cell = header_row.cells[index]
        cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
        set_cell_shading(cell, "D9EAF7")
        paragraph = cell.paragraphs[0]
        paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
        render_inline(paragraph, cell_node.get("children", []), bold=True)
    set_repeat_table_header(header_row)
    prevent_row_split(header_row)

    for row_index, row_node in enumerate(body_rows):
        row = table.add_row()
        prevent_row_split(row)
        for column_index, cell_node in enumerate(row_node.get("children", [])):
            cell = row.cells[column_index]
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            if row_index % 2:
                set_cell_shading(cell, "F5F7FA")
            paragraph = cell.paragraphs[0]
            paragraph.alignment = WD_ALIGN_PARAGRAPH.LEFT
            render_inline(paragraph, cell_node.get("children", []))

    for row in table.rows:
        for cell in row.cells:
            cell.margin_top = Cm(0.08)
            cell.margin_bottom = Cm(0.08)
            for paragraph in cell.paragraphs:
                paragraph.paragraph_format.space_before = Pt(0)
                paragraph.paragraph_format.space_after = Pt(0)
                paragraph.paragraph_format.line_spacing = 1.0
                for run in paragraph.runs:
                    set_run_font(run, "Times New Roman", 8)


def configure_document(document: Document, title: str) -> None:
    section = document.sections[0]
    section.page_width = Cm(21.0)
    section.page_height = Cm(29.7)
    section.top_margin = Cm(2.3)
    section.bottom_margin = Cm(2.3)
    section.left_margin = Cm(2.2)
    section.right_margin = Cm(2.2)
    section.header_distance = Cm(1.1)
    section.footer_distance = Cm(1.1)

    styles = document.styles
    normal = styles["Normal"]
    normal.font.name = "Times New Roman"
    normal.font.size = Pt(11)
    normal._element.rPr.rFonts.set(qn("w:ascii"), "Times New Roman")
    normal._element.rPr.rFonts.set(qn("w:hAnsi"), "Times New Roman")
    normal.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    normal.paragraph_format.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE
    normal.paragraph_format.space_after = Pt(6)

    title_style = styles["Title"]
    title_style.font.name = "Times New Roman"
    title_style.font.size = Pt(16)
    title_style.font.bold = True
    title_style.font.color.rgb = RGBColor(31, 78, 121)
    title_style.paragraph_format.space_after = Pt(10)

    for level, size in ((1, 14), (2, 12), (3, 11)):
        style = styles[f"Heading {level}"]
        style.font.name = "Times New Roman"
        style.font.size = Pt(size)
        style.font.bold = True
        style.font.color.rgb = RGBColor(31, 78, 121)
        style.paragraph_format.keep_with_next = True
        style.paragraph_format.space_before = Pt(12 if level == 1 else 8)
        style.paragraph_format.space_after = Pt(4)

    caption = styles["Caption"]
    caption.font.name = "Times New Roman"
    caption.font.size = Pt(9)
    caption.font.italic = False
    caption.font.color.rgb = RGBColor(45, 45, 45)
    caption.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    caption.paragraph_format.line_spacing = 1.0
    caption.paragraph_format.space_before = Pt(3)
    caption.paragraph_format.space_after = Pt(8)

    header = section.header.paragraphs[0]
    header.text = "Thermal extremes and dry–hot concurrence in Iran"
    header.alignment = WD_ALIGN_PARAGRAPH.CENTER
    for run in header.runs:
        set_run_font(run, "Times New Roman", 9)
        run.italic = True
        run.font.color.rgb = RGBColor(90, 90, 90)
    add_page_number(section.footer.paragraphs[0])

    properties = document.core_properties
    properties.title = title
    properties.subject = "Q1 journal manuscript"
    properties.keywords = "Iran; temperature extremes; dry-hot events; quantile regression"
    properties.comments = "Generated from reports/Manuscript_Q1_2026.md"


def build_document(source: Path, destination: Path) -> None:
    markdown = mistune.create_markdown(renderer="ast", plugins=["table", "math"])
    ast = markdown(source.read_text(encoding="utf-8"))
    title = plain_text(ast[0].get("children", []))
    document = Document()
    configure_document(document, title)
    in_references = False

    for block in ast:
        block_type = block.get("type")
        if block_type in {"blank_line", "block_html"}:
            continue

        if block_type == "heading":
            level = block.get("attrs", {}).get("level", 2)
            heading_text = plain_text(block.get("children", []))
            if level == 1:
                paragraph = document.add_paragraph(style="Title")
                paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
            else:
                paragraph = document.add_paragraph(style=f"Heading {min(level - 1, 3)}")
            render_inline(paragraph, block.get("children", []))
            in_references = heading_text.strip().lower() == "references"
            continue

        if block_type == "table":
            add_table(document, block)
            continue

        if block_type == "paragraph":
            children = block.get("children", [])
            if len(children) == 1 and children[0].get("type") == "image":
                add_figure(document, source.parent, children[0])
                continue
            if len(children) == 1 and children[0].get("type") == "inline_math":
                add_equation(document, children[0].get("raw", ""))
                continue

            text = plain_text(children).strip()
            is_caption = text.startswith("Figure ") or text.startswith("Table ")
            paragraph = document.add_paragraph(style="Caption" if is_caption else None)
            if text.startswith("Running title:"):
                paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
                paragraph.paragraph_format.space_after = Pt(10)
            elif in_references:
                paragraph.paragraph_format.left_indent = Cm(0.6)
                paragraph.paragraph_format.first_line_indent = Cm(-0.6)
                paragraph.paragraph_format.line_spacing = 1.0
                paragraph.paragraph_format.space_after = Pt(5)
            render_inline(paragraph, children)
            continue

        # Preserve uncommon block types as text rather than silently dropping them.
        fallback = plain_text(block.get("children", [])) or block.get("raw", "")
        if fallback.strip():
            paragraph = document.add_paragraph()
            add_run(paragraph, fallback)

    destination.parent.mkdir(parents=True, exist_ok=True)
    document.save(destination)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "source",
        nargs="?",
        type=Path,
        default=Path("reports/Manuscript_Q1_2026.md"),
    )
    parser.add_argument(
        "destination",
        nargs="?",
        type=Path,
        default=Path("reports/Manuscript_Q1_2026.docx"),
    )
    args = parser.parse_args()
    build_document(args.source.resolve(), args.destination.resolve())
    print(f"Created: {args.destination.resolve()}")


if __name__ == "__main__":
    main()
