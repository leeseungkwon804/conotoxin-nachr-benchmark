#!/usr/bin/env python3
"""
Rough composites of Figures 5-7 (and Figure 2 with its new panel e) from the individual panels.

A starting layout for rearranging in PowerPoint, not the final figures.

Uses PIL rather than matplotlib: this machine's matplotlib font stack is broken (font_manager
imports plistlib, which imports a pyexpat built against a mismatched libexpat), and the pods that
would otherwise do the rendering are down. PIL only needs to paste and label, so it is enough.

Usage: python3 scripts/assemble_figures_v3_pil.py
"""
import os

from PIL import Image, ImageDraw, ImageFont

ROOT = "/Volumes/dkyoo_SSD1/review_partime/conotoxin_project"
V3 = os.path.join(ROOT, "manuscript", "figures_v3")
PAN = os.path.join(V3, "panels")

FONT_CANDIDATES = [
    "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
    "/System/Library/Fonts/Helvetica.ttc",
    "/Library/Fonts/Arial Bold.ttf",
]
GAP = 26          # px between panels
PADDING = 34      # px around the whole figure
LETTER_PAD = 8


def load_font(size):
    for p in FONT_CANDIDATES:
        if os.path.exists(p):
            try:
                return ImageFont.truetype(p, size)
            except Exception:
                pass
    return ImageFont.load_default()


def open_panel(name, width):
    """Panel scaled to a target width, or a placeholder if it does not exist."""
    p = os.path.join(PAN, name + ".png")
    if not os.path.exists(p):
        ph = Image.new("RGB", (width, int(width * 0.55)), "#FAFAFA")
        d = ImageDraw.Draw(ph)
        d.text((width // 2, int(width * 0.27)), "%s\n(not yet made)" % name,
               fill="#B00020", anchor="mm", align="center", font=load_font(28))
        return ph
    im = Image.open(p).convert("RGB")
    h = max(1, round(im.height * width / im.width))
    return im.resize((width, h), Image.LANCZOS)


def compose(out_name, rows, total_width=2400):
    """rows: list of lists of (panel_name, letter). Each row splits the width evenly."""
    font = load_font(52)
    imgs = []
    for row in rows:
        w = (total_width - GAP * (len(row) - 1)) // len(row)
        imgs.append([(open_panel(n, w), letter) for n, letter in row])

    height = PADDING * 2 + sum(max(im.height for im, _ in row) for row in imgs) + GAP * (len(imgs) - 1)
    canvas = Image.new("RGB", (total_width + PADDING * 2, height), "white")
    draw = ImageDraw.Draw(canvas)

    y = PADDING
    for row in imgs:
        x = PADDING
        rh = max(im.height for im, _ in row)
        for im, letter in row:
            canvas.paste(im, (x, y))
            if letter:
                draw.text((x + LETTER_PAD, y + LETTER_PAD), letter, fill="black", font=font)
            x += im.width + GAP
        y += rh + GAP

    canvas.save(os.path.join(V3, out_name + ".png"), dpi=(400, 400))
    canvas.save(os.path.join(V3, out_name + ".pdf"), "PDF", resolution=400)
    print("  wrote %s  (%d x %d)" % (out_name, canvas.width, canvas.height))


print("assembling into", V3)

# Figure 5 — contact map on its own row (it is tall), then the two supporting panels
compose("Figure5", [
    [("Figure5_a", "a")],
    [("Figure5_b", "b"), ("Figure5_c_structure", "c")],
])

# Figure 6
compose("Figure6", [[("Figure6_a", "a")], [("Figure6_b", "b")]])

# Figure 7
compose("Figure7", [[("Figure7_a", "a")], [("Figure7_b", "b")]])

# Figure 2 — existing a-d plus the new e
compose("Figure2_with_e", [
    [("Figure2_a", "a"), ("Figure2_b", "b")],
    [("Figure2_c", "c"), ("Figure2_d", "d")],
    [("Figure2_e", "e")],
])

# an alternative Figure 5 in which panel c is the quantitative version instead of the render
compose("Figure5_altc_quantitative", [
    [("Figure5_a", "a")],
    [("Figure5_b", "b"), ("Figure5_c", "c")],
])
print("done")
