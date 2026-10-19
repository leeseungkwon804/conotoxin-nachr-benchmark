#!/usr/bin/env /usr/bin/python3
"""Fix the Supplementary Materials docx: the wide tables (up to 31 columns) were
rendered with pandoc's default equal-width autofit-to-page-width layout, so header
and cell text wrapped one character per line. Fix: put the 'Supplementary Tables'
section in landscape orientation with narrower margins, and give every table
explicit fixed-width columns sized by each column's actual content length, at a
smaller font size.

Run ONCE, right after regenerating Supplementary_Materials.docx from
SUPP_COMBINED.md via pandoc (see build_supplementary.py in the session
scratchpad). Not idempotent -- re-running against an already-fixed file
inserts a second landscape section break.
"""
import docx
from docx.oxml.ns import qn
import copy

ROOT = "/Volumes/dkyoo_SSD1/review_partime/conotoxin_project"
SRC = ROOT + "/manuscript/supplementary/Supplementary_Materials.docx"
OUT = ROOT + "/manuscript/supplementary/Supplementary_Materials.docx"

d = docx.Document(SRC)
body = d.element.body

# ---- 1. locate the split point: last paragraph of the Figures section ----
target_text = "(Caption is rendered directly in the figure image above.)"
split_after = None
for p in d.paragraphs:
    if p.text.strip() == target_text:
        split_after = p._p
# split_after is now the LAST such paragraph (Figure S6's caption note)
assert split_after is not None

# ---- 2. build a portrait sectPr (ends section 1) and insert as its own paragraph ----
PORTRAIT_SECTPR = (
    '<w:sectPr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
    '<w:pgSz w:w="12240" w:h="15840"/>'
    '<w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" '
    'w:header="720" w:footer="720" w:gutter="0"/>'
    '</w:sectPr>'
)
from lxml import etree
sect_p = etree.SubElement(body, qn('w:p'))
ppr = etree.SubElement(sect_p, qn('w:pPr'))
ppr.append(etree.fromstring(PORTRAIT_SECTPR))
split_after.addnext(sect_p)

# ---- 3. change the body-level final sectPr to landscape, narrower margins ----
final_sectpr = body.find(qn('w:sectPr'))
for child in list(final_sectpr):
    final_sectpr.remove(child)
LANDSCAPE_INNER = (
    '<w:sectPr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
    '<w:pgSz w:w="15840" w:h="12240" w:orient="landscape"/>'
    '<w:pgMar w:top="900" w:right="720" w:bottom="900" w:left="720" '
    'w:header="720" w:footer="720" w:gutter="0"/>'
    '<w:footnotePr><w:numRestart w:val="eachSect"/></w:footnotePr>'
    '</w:sectPr>'
)
new_final = etree.fromstring(LANDSCAPE_INNER)
for child in list(new_final):
    final_sectpr.append(child)

AVAILABLE_TWIPS = 15840 - 720 - 720  # landscape width minus L/R margins = 14400

# ---- 4. resize every table's columns by content length, shrink font ----
W_NS = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'


def cell_text(tc):
    return "".join(t.text or "" for t in tc.iter(qn('w:t')))


def fix_table(tbl):
    tblPr = tbl.find(qn('w:tblPr'))
    tblGrid = tbl.find(qn('w:tblGrid'))
    rows = tbl.findall(qn('w:tr'))
    if not rows:
        return
    ncols = len(tblGrid.findall(qn('w:gridCol')))

    # sample column max content length (chars), across all rows
    col_max = [0] * ncols
    for tr in rows:
        tcs = tr.findall(qn('w:tc'))
        for i, tc in enumerate(tcs):
            if i >= ncols:
                break
            txt = cell_text(tc)
            # count each "word chunk" length capped, since wrapping happens on
            # whitespace normally -- but these are mostly single-token cells
            col_max[i] = max(col_max[i], len(txt))

    # weight = sqrt-damped content length so extreme outliers don't eat the row,
    # with a floor so short numeric columns stay legible
    weights = [max(3.2, min(cm, 22)) ** 0.75 for cm in col_max]
    total_w = sum(weights)
    MIN_TWIPS = 460
    widths = [max(MIN_TWIPS, round(AVAILABLE_TWIPS * w / total_w)) for w in weights]
    # renormalize to exactly fill AVAILABLE_TWIPS after applying the floor
    scale = AVAILABLE_TWIPS / sum(widths)
    widths = [round(w * scale) for w in widths]

    # rewrite tblGrid
    for old in list(tblGrid):
        tblGrid.remove(old)
    for w in widths:
        gc = etree.SubElement(tblGrid, qn('w:gridCol'))
        gc.set(qn('w:w'), str(w))

    # tblPr: fixed layout, absolute width in dxa (not pct)
    tblW = tblPr.find(qn('w:tblW'))
    tblW.set(qn('w:type'), 'dxa')
    tblW.set(qn('w:w'), str(AVAILABLE_TWIPS))
    layout = tblPr.find(qn('w:tblLayout'))
    if layout is None:
        layout = etree.SubElement(tblPr, qn('w:tblLayout'))
    layout.set(qn('w:type'), 'fixed')

    # per-cell explicit width (must match grid) + font size on every run
    n_data_rows = len(rows) - 1
    fsz = '15' if ncols <= 8 else ('14' if ncols <= 16 else '12')  # half-points: 7.5/7/6 pt
    for tr in rows:
        tcs = tr.findall(qn('w:tc'))
        for i, tc in enumerate(tcs):
            if i >= ncols:
                break
            tcPr = tc.find(qn('w:tcPr'))
            if tcPr is None:
                tcPr = etree.SubElement(tc, qn('w:tcPr'))
                tc.insert(0, tcPr)
            tcW = tcPr.find(qn('w:tcW'))
            if tcW is None:
                tcW = etree.SubElement(tcPr, qn('w:tcW'))
            tcW.set(qn('w:w'), str(widths[i]))
            tcW.set(qn('w:type'), 'dxa')
            for r in tc.iter(qn('w:r')):
                rpr = r.find(qn('w:rPr'))
                if rpr is None:
                    rpr = etree.SubElement(r, qn('w:rPr'))
                    r.insert(0, rpr)
                sz = rpr.find(qn('w:sz'))
                if sz is None:
                    sz = etree.SubElement(rpr, qn('w:sz'))
                sz.set(qn('w:val'), fsz)
                szCs = rpr.find(qn('w:szCs'))
                if szCs is None:
                    szCs = etree.SubElement(rpr, qn('w:szCs'))
                szCs.set(qn('w:val'), fsz)


for tbl in d.tables:
    fix_table(tbl._tbl)

d.save(OUT)
print("saved", OUT)
