#!/usr/bin/env python3
"""
Write the v3 captions into the figure PPTX, one per slide.

Slides 5-7 still carry the "Figure N. Title.  Text." placeholder and slides 1-4 and S1-S6 carry
the superseded long captions; both are replaced from figures_v3/CAPTIONS_v3.md, which is the
single source the manuscript also uses.

Run with /usr/bin/python3 — the homebrew interpreter on this machine has a broken pyexpat and
cannot parse the Office XML.

Usage: /usr/bin/python3 scripts/set_pptx_captions.py [--dry-run]
"""
import copy
import os
import re
import sys

from pptx import Presentation
from pptx.enum.text import MSO_AUTO_SIZE
from pptx.oxml.ns import qn
from pptx.util import Emu, Pt

ROOT = "/Volumes/dkyoo_SSD1/review_partime/conotoxin_project"
V3 = os.path.join(ROOT, "manuscript", "figures_v3")
PPTX = os.path.join(V3, "2026-08-14  Figures v3.pptx")
CAPS = os.path.join(V3, "CAPTIONS_v3.md")

SLIDE_ORDER = ["1", "2", "3", "4", "5", "6", "7", "S1", "S2", "S3", "S4", "S5", "S6"]
INLINE = re.compile(r"(\*\*.+?\*\*|\*.+?\*)")


def load_captions():
    """-> {'1': 'Figure 1. ...', 'S5': 'Figure S5. ...'} from the markdown source."""
    out = {}
    for line in open(CAPS, encoding="utf-8"):
        line = line.strip()
        m = re.match(r"\*\*Figure (S?\d+)\.\*\*\s*(.+)$", line)
        if m:
            out[m.group(1)] = line
    missing = [k for k in SLIDE_ORDER if k not in out]
    if missing:
        raise SystemExit("no caption found for: %s" % ", ".join(missing))
    return out


EMU_IN = 914400.0


def fitted_size(text, width_in, height_in, start_pt=12.0, floor_pt=7.0):
    """Largest point size (in 0.5 pt steps, down to floor_pt) whose wrapped text fits the box.

    python-pptx can set SHAPE_TO_FIT_TEXT but cannot reflow — only PowerPoint recomputes the
    height on open. A caption that overflows would sit on top of the figure in every other viewer,
    so the size is chosen here instead of trusting the autofit flag.

    Arial's average advance width is about 0.52 em over mixed-case prose; line spacing is 1.2 em.
    Both are approximations, hence the 10% headroom.
    """
    usable_w = width_in - 0.2          # 0.1 in text-frame inset each side
    usable_h = height_in - 0.1
    size = start_pt
    while size > floor_pt:
        chars_per_line = usable_w / (0.52 * size / 72.0)
        lines = -(-len(text) // int(chars_per_line))
        if lines * (1.2 * size / 72.0) <= usable_h * 0.9:
            return size
        size -= 0.5
    return floor_pt


def set_caption(tf, text, proto_run, size_pt=None):
    """Replace a text frame's content with text, honouring **bold**/*italic* markdown.

    Every new run is cloned from proto_run so the slide's own font, colour and typeface survive —
    the point is to change the words, not the typography DK set. Bold and italic are then set
    explicitly on every run: the prototype is the bold "Figure N." run, so a clone that is merely
    left alone comes out bold, which is how the whole caption went bold the first time.

    Runs are inserted BEFORE any a:endParaRPr. DrawingML fixes the child order of a:p as
    (a:pPr?, runs*, a:endParaRPr?); appending after that trailing element produces a file
    LibreOffice renders happily and PowerPoint rejects, dropping the caption. That is exactly what
    happened the first time this ran.
    """
    p = tf.paragraphs[0]
    for r in list(p.runs):
        r._r.getparent().remove(r._r)
    for extra in list(tf.paragraphs)[1:]:
        extra._p.getparent().remove(extra._p)
    tail = p._p.find(qn("a:endParaRPr"))

    for piece in INLINE.split(text):
        if not piece:
            continue
        bold = piece.startswith("**") and piece.endswith("**")
        ital = not bold and piece.startswith("*") and piece.endswith("*")
        body = piece[2:-2] if bold else piece[1:-1] if ital else piece
        r = copy.deepcopy(proto_run._r)
        if tail is not None:
            tail.addprevious(r)
        else:
            p._p.append(r)
        run = p.runs[-1]
        run.text = body
        run.font.bold = bold
        run.font.italic = ital
        if size_pt is not None:
            run.font.size = Pt(size_pt)


ORDER = ["pPr", "r", "br", "fld", "endParaRPr"]


def check_order(prs):
    """Fail before saving if any paragraph violates the a:p child order PowerPoint enforces."""
    bad = []
    for i, slide in enumerate(prs.slides, 1):
        for sh in slide.shapes:
            if not sh.has_text_frame:
                continue
            for para in sh.text_frame.paragraphs:
                seen = [c.tag.split("}")[1] for c in para._p]
                ranks = [ORDER.index(t) for t in seen if t in ORDER]
                if ranks != sorted(ranks):
                    bad.append("slide %d: %s" % (i, " ".join(seen)))
    if bad:
        raise SystemExit("a:p child order is wrong — PowerPoint would drop this text:\n  "
                         + "\n  ".join(bad))
    print("\na:p child order valid on every paragraph of every slide.")


def main():
    dry = "--dry-run" in sys.argv
    caps = load_captions()
    prs = Presentation(PPTX)
    if len(prs.slides) != len(SLIDE_ORDER):
        raise SystemExit("expected %d slides, found %d — slide order is no longer safe to assume"
                         % (len(SLIDE_ORDER), len(prs.slides)))

    slide_h = prs.slide_height / EMU_IN
    shrunk = []
    for slide, key in zip(prs.slides, SLIDE_ORDER):
        # the caption box is the widest text box on the slide; panel letters are ~0.34 in
        boxes = [sh for sh in slide.shapes
                 if sh.has_text_frame and sh.width > Emu(int(3 * 914400))]
        if not boxes:
            raise SystemExit("slide for Figure %s has no caption-width text box" % key)
        box = max(boxes, key=lambda sh: sh.width)
        old = box.text_frame.text.strip()
        proto = None
        for para in box.text_frame.paragraphs:
            if para.runs:
                proto = para.runs[0]
                break
        if proto is None:
            raise SystemExit("caption box on the Figure %s slide has no run to copy from" % key)

        plain = re.sub(r"\*+", "", caps[key])
        avail = slide_h - box.top / EMU_IN
        base = (proto.font.size.pt if proto.font.size else 12.0)
        size = fitted_size(plain, box.width / EMU_IN, avail, start_pt=base)
        print("Figure %-3s  %3d words  %.1f in of room  -> %.1f pt%s"
              % (key, len(plain.split()), avail, size,
                 "" if size == base else "  (was %.1f)" % base))
        if dry:
            continue

        set_caption(box.text_frame, caps[key], proto, size_pt=size)
        box.text_frame.word_wrap = True
        box.text_frame.auto_size = MSO_AUTO_SIZE.SHAPE_TO_FIT_TEXT
        # store a height matching the fitted text, so the box is right even in viewers that do
        # not recompute autofit; PowerPoint will refine it on open
        chars_per_line = (box.width / EMU_IN - 0.2) / (0.52 * size / 72.0)
        lines = -(-len(plain) // int(chars_per_line))
        box.height = Emu(int((lines * 1.2 * size / 72.0 + 0.1) * EMU_IN))
        def room_needed(pt):
            cpl = (box.width / EMU_IN - 0.2) / (0.52 * pt / 72.0)
            return (-(-len(plain) // int(cpl))) * 1.2 * pt / 72.0 / 0.9 + 0.1

        if room_needed(size) > avail:
            raise SystemExit("Figure %s does not fit even at the %.1f pt floor (%.2f in needed, "
                             "%.2f in available) — the caption would sit on the figure. Move the "
                             "panel block up before re-running." % (key, size, room_needed(size), avail))
        if size < base:
            shrunk.append("Figure %s: %.1f -> %.1f pt. %.2f in of room below the figure; "
                          "10 pt needs %.2f in (move panels up %.2f in), "
                          "12 pt needs %.2f in (move panels up %.2f in)."
                          % (key, base, size, avail,
                             room_needed(10), max(0, room_needed(10) - avail),
                             room_needed(12), max(0, room_needed(12) - avail)))

    if dry:
        return
    check_order(prs)
    prs.save(PPTX)
    print("\nsaved %s" % PPTX)
    if shrunk:
        print("\ncaptions set below the slide's own 12 pt to fit the space left under the figure —"
              "\nmove the panel block up if you would rather keep them at 12 pt:")
        for w in shrunk:
            print("  " + w)
    else:
        print("every caption fits at the slide's own point size.")


if __name__ == "__main__":
    main()
