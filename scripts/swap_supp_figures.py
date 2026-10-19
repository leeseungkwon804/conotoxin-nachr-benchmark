#!/usr/bin/env python3
"""
Replace the supplementary figure images in the PPTX with the re-rendered uppercase-letter versions.

S2, S4, S5 and S6 carried lowercase panel letters baked into the image while the main figures used
uppercase; the four scripts were changed to draw uppercase and re-run at 600 dpi. Only those four
slides are touched.

Each replacement keeps the picture's left edge and centre, and fits it inside the box the old image
occupied without distorting the new aspect ratio — the re-rendered figures are not the same shape as
the ones DK placed, so scaling to the old width and height would stretch them.

Run with /usr/bin/python3 (the homebrew interpreter cannot parse Office XML here).

Usage: /usr/bin/python3 scripts/swap_supp_figures.py
"""
import os

from PIL import Image
from pptx import Presentation
from pptx.util import Emu

ROOT = "/Volumes/dkyoo_SSD1/review_partime/conotoxin_project"
V3 = os.path.join(ROOT, "manuscript", "figures_v3")
PPTX = os.path.join(V3, "2026-08-14  Figures v3.pptx")
REGEN = os.path.join(V3, "_regen")

# slide number (1-based) -> re-rendered file
SWAPS = {9: "FigureS2.png", 11: "FigureS4.png", 12: "FigureS5.png", 13: "FigureS6.png"}

EMU_IN = 914400.0


def main():
    prs = Presentation(PPTX)
    for slide_no, fname in sorted(SWAPS.items()):
        slide = prs.slides[slide_no - 1]
        pics = [sh for sh in slide.shapes if sh.shape_type == 13]
        if len(pics) != 1:
            raise SystemExit("slide %d has %d pictures, expected exactly 1" % (slide_no, len(pics)))
        old = pics[0]
        box_l, box_t = old.left, old.top
        box_w, box_h = old.width, old.height

        # the caption box fixes both the usable width and the floor the image must stay above
        cap = max((sh for sh in slide.shapes
                   if sh.has_text_frame and sh.width > Emu(int(3 * EMU_IN))),
                  key=lambda sh: sh.width)
        avail_w = cap.width
        avail_h = cap.top - box_t - Emu(int(0.15 * EMU_IN))

        path = os.path.join(REGEN, fname)
        w, h = Image.open(path).size
        # the re-rendered figures are wider than the ones they replace, so they are refitted to the
        # usable width rather than left at the old box size, which would waste a third of the slide
        scale = min(avail_w / w, avail_h / h)
        new_w, new_h = int(w * scale), int(h * scale)
        new_l = cap.left + (avail_w - new_w) // 2
        new_t = box_t

        el = old._element
        el.getparent().remove(el)
        pic = slide.shapes.add_picture(path, Emu(new_l), Emu(new_t), Emu(new_w), Emu(new_h))
        # put it back at the bottom of the z-order so the caption stays on top
        slide.shapes._spTree.remove(pic._element)
        slide.shapes._spTree.insert(2, pic._element)

        print("slide %-2d %-13s  %.2f x %.2f in -> %.2f x %.2f in  (bottom %.2f in)"
              % (slide_no, fname, box_w / EMU_IN, box_h / EMU_IN,
                 new_w / EMU_IN, new_h / EMU_IN, (new_t + new_h) / EMU_IN))

    prs.save(PPTX)
    print("\nsaved %s" % PPTX)


if __name__ == "__main__":
    main()
