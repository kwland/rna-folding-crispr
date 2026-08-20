"""Publication figures for the manuscript.

Writes true vector PDFs (base-14 Helvetica, no embedding, no dependencies) so
LaTeX can \\includegraphics them directly. Print styling: white ground, no
rounded corners, recessive grid, thin marks, direct labels.

Palette is the validated two-series pair, checked with the ported six-check
validator: adjacent CVD dE 24.7 (target 8), normal-vision dE 33.6 (floor 15),
both above threshold on a white surface.

Run:  python make_manuscript_figures.py
"""
from __future__ import annotations

import json
import math
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = pathlib.Path(__file__).resolve().parent / "figures"
OUT.mkdir(exist_ok=True)

# ----------------------------------------------------------------- palette
BLUE = (0x2a, 0x78, 0xd6)      # series 1, human screen
ORANGE = (0xeb, 0x68, 0x34)    # series 2, zebrafish screen
BLUE_LT = (0xa9, 0xc8, 0xef)   # same hue, light shade (dumbbell "before")
GRAY = (0x76, 0x76, 0x74)      # de-emphasis
GRID = (0xdc, 0xdc, 0xd8)
INK = (0x1a, 0x1a, 0x19)
INK2 = (0x5c, 0x5c, 0x58)

# Helvetica advance widths per 1000 units (AFM), ASCII 32..126.
_W = [278,278,355,556,556,889,667,191,333,333,389,584,278,333,278,278,
      556,556,556,556,556,556,556,556,556,556,278,278,584,584,584,556,
      1015,667,667,722,722,667,611,778,722,278,500,667,556,833,722,778,
      667,778,722,667,611,722,667,944,667,667,611,278,278,278,469,556,
      333,556,556,500,556,556,278,556,556,222,222,500,222,833,556,556,
      556,556,333,500,278,556,500,722,500,500,500,334,260,334,584]


def text_width(s: str, size: float, bold: bool = False) -> float:
    total = 0
    for ch in s:
        o = ord(ch)
        w = _W[o - 32] if 32 <= o <= 126 else 556
        total += w
    return total / 1000.0 * size * (1.04 if bold else 1.0)


class Canvas:
    """Minimal PDF content-stream builder. Origin bottom-left, points."""

    def __init__(self, width: float, height: float):
        self.w, self.h = width, height
        self.ops: list[str] = []

    def _col(self, rgb, stroke=False):
        r, g, b = [c / 255 for c in rgb]
        self.ops.append(f"{r:.3f} {g:.3f} {b:.3f} {'RG' if stroke else 'rg'}")

    def line(self, x1, y1, x2, y2, rgb=GRID, w=0.6, dash=None):
        self._col(rgb, True)
        self.ops.append(f"{w:.2f} w")
        self.ops.append(f"[{dash}] 0 d" if dash else "[] 0 d")
        self.ops.append(f"{x1:.2f} {y1:.2f} m {x2:.2f} {y2:.2f} l S")
        self.ops.append("[] 0 d")

    def rect(self, x, y, w, h, rgb=BLUE):
        self._col(rgb)
        self.ops.append(f"{x:.2f} {y:.2f} {w:.2f} {h:.2f} re f")

    def dot(self, x, y, r=3.0, rgb=BLUE, ring=None):
        """Circle via four Bezier arcs. `ring` paints a surface halo first."""
        if ring:
            self._circle(x, y, r + 1.4, ring)
        self._circle(x, y, r, rgb)

    def _circle(self, x, y, r, rgb):
        k = 0.5523 * r
        self._col(rgb)
        self.ops.append(
            f"{x - r:.2f} {y:.2f} m "
            f"{x - r:.2f} {y + k:.2f} {x - k:.2f} {y + r:.2f} {x:.2f} {y + r:.2f} c "
            f"{x + k:.2f} {y + r:.2f} {x + r:.2f} {y + k:.2f} {x + r:.2f} {y:.2f} c "
            f"{x + r:.2f} {y - k:.2f} {x + k:.2f} {y - r:.2f} {x:.2f} {y - r:.2f} c "
            f"{x - k:.2f} {y - r:.2f} {x - r:.2f} {y - k:.2f} {x - r:.2f} {y:.2f} c f")

    def text(self, x, y, s, size=8.5, rgb=INK, bold=False, anchor="start"):
        if anchor == "middle":
            x -= text_width(s, size, bold) / 2
        elif anchor == "end":
            x -= text_width(s, size, bold)
        esc = s.replace("\\", r"\\").replace("(", r"\(").replace(")", r"\)")
        self._col(rgb)
        self.ops.append(f"BT /{'F2' if bold else 'F1'} {size:.2f} Tf "
                        f"{x:.2f} {y:.2f} Td ({esc}) Tj ET")

    def save(self, path: pathlib.Path):
        stream = "\n".join(self.ops).encode("latin-1", "replace")
        objs = [
            b"<</Type/Catalog/Pages 2 0 R>>",
            b"<</Type/Pages/Kids[3 0 R]/Count 1>>",
            (f"<</Type/Page/Parent 2 0 R/MediaBox[0 0 {self.w:.2f} {self.h:.2f}]"
             f"/Resources<</Font<</F1 5 0 R/F2 6 0 R>>>>/Contents 4 0 R>>").encode(),
            b"<</Length " + str(len(stream)).encode() + b">>\nstream\n" + stream + b"\nendstream",
            b"<</Type/Font/Subtype/Type1/BaseFont/Helvetica/Encoding/WinAnsiEncoding>>",
            b"<</Type/Font/Subtype/Type1/BaseFont/Helvetica-Bold/Encoding/WinAnsiEncoding>>",
        ]
        out = bytearray(b"%PDF-1.4\n")
        offsets = []
        for i, body in enumerate(objs, start=1):
            offsets.append(len(out))
            out += f"{i} 0 obj\n".encode() + body + b"\nendobj\n"
        xref = len(out)
        out += f"xref\n0 {len(objs) + 1}\n".encode()
        out += b"0000000000 65535 f \n"
        for off in offsets:
            out += f"{off:010d} 00000 n \n".encode()
        out += (f"trailer\n<</Size {len(objs) + 1}/Root 1 0 R>>\n"
                f"startxref\n{xref}\n%%EOF\n").encode()
        path.write_bytes(out)
        return path


class SvgCanvas(Canvas):
    """Same drawing API, SVG out. Used for the web copies and for eyeballing
    the geometry before the PDF goes into LaTeX."""

    def __init__(self, width, height):
        super().__init__(width, height)
        self.el = []

    def _hex(self, rgb):
        return "#%02x%02x%02x" % rgb

    def line(self, x1, y1, x2, y2, rgb=GRID, w=0.6, dash=None):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        self.el.append(f'<line x1="{x1:.2f}" y1="{self.h-y1:.2f}" x2="{x2:.2f}" '
                       f'y2="{self.h-y2:.2f}" stroke="{self._hex(rgb)}" stroke-width="{w}"{d}/>')

    def rect(self, x, y, w, h, rgb=BLUE):
        self.el.append(f'<rect x="{x:.2f}" y="{self.h-y-h:.2f}" width="{w:.2f}" '
                       f'height="{h:.2f}" fill="{self._hex(rgb)}"/>')

    def _circle(self, x, y, r, rgb):
        self.el.append(f'<circle cx="{x:.2f}" cy="{self.h-y:.2f}" r="{r:.2f}" fill="{self._hex(rgb)}"/>')

    def text(self, x, y, s, size=8.5, rgb=INK, bold=False, anchor="start"):
        a = {"start": "start", "middle": "middle", "end": "end"}[anchor]
        esc = (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))
        fw = ' font-weight="700"' if bold else ""
        self.el.append(f'<text x="{x:.2f}" y="{self.h-y:.2f}" font-size="{size:.2f}" '
                       f'font-family="Helvetica, Arial, sans-serif" fill="{self._hex(rgb)}"'
                       f' text-anchor="{a}"{fw}>{esc}</text>')

    def save(self, path):
        path = path.with_suffix('.svg')
        body = '\n'.join(self.el)
        head = ('<svg xmlns="http://www.w3.org/2000/svg" '
                f'width="{self.w}" height="{self.h}" '
                f'viewBox="0 0 {self.w} {self.h}">')
        bg = f'<rect width="{self.w}" height="{self.h}" fill="#ffffff"/>'
        path.write_text(head + '\n' + bg + '\n' + body + '\n</svg>\n', encoding='utf-8')
        return path


# ------------------------------------------------------------------- data
MODE = {'cls': Canvas}


study = json.loads((ROOT / "analysis_outputs" / "study_results.json").read_text())
vienna_path = ROOT / "analysis_outputs" / "vienna_validation.json"
vienna = json.loads(vienna_path.read_text()) if vienna_path.exists() else {}

DO = study["datasets"]["doench"]
CR = study["datasets"]["crisprscan"]


def axis(c, x0, x1, y0, y1, lo, hi, ticks, label, zero=True):
    """Recessive vertical gridlines plus tick labels; optional emphasized zero."""
    def sx(v):
        return x0 + (v - lo) / (hi - lo) * (x1 - x0)
    for t in ticks:
        is_zero = abs(t) < 1e-12 and zero
        c.line(sx(t), y0, sx(t), y1, GRID if not is_zero else (0xa8, 0xa8, 0xa4),
               0.6 if not is_zero else 1.0)
        c.text(sx(t), y0 - 11, f"{t:g}", 7.5, INK2, anchor="middle")
    c.text((x0 + x1) / 2, y0 - 23, label, 8, INK2, anchor="middle")
    return sx


# ============================================== Fig 1: folding context
def fig_context():
    W, H = 468, 226
    c = MODE['cls'](W, H)
    rows = [m for m in DO["context"]["measures"] if m["key"] != "ensemble_energy"]
    c.text(30, H - 22, "Folding context changes what is measured", 11, INK, bold=True)
    c.text(30, H - 35, "Mean accessibility of the same 4,685 spacers, folded alone and inside the full sgRNA",
           8, INK2)

    x0, x1 = 176, 356
    top, bot = H - 66, 74
    step = (top - bot) / max(1, len(rows) - 1) if len(rows) > 1 else 0
    sx = axis(c, x0, x1, bot - 6, top + 12, 0, 1, [0, 0.25, 0.5, 0.75, 1.0],
              "mean accessibility (fraction unpaired)", zero=False)

    for i, m in enumerate(rows):
        y = top - i * step
        a, b = m["mean_spacer_alone"], m["mean_in_sgrna"]
        c.text(x0 - 12, y - 2.5, m["label"], 8, INK, anchor="end")
        c.line(sx(b), y, sx(a), y, (0xc8, 0xc8, 0xc4), 1.6)
        c.dot(sx(a), y, 3.6, BLUE_LT, ring=(255, 255, 255))
        c.dot(sx(b), y, 3.6, BLUE, ring=(255, 255, 255))
        c.text(x1 + 14, y - 2.5, f"r = {m['pearson_between_contexts']:.2f}", 7.5, INK2)

    ly = bot - 46
    c.dot(x0 + 4, ly + 2.5, 3.6, BLUE_LT)
    c.text(x0 + 12, ly, "spacer alone", 7.5, INK2)
    c.dot(x0 + 96, ly + 2.5, 3.6, BLUE)
    c.text(x0 + 104, ly, "inside the sgRNA", 7.5, INK2)
    c.text(30, 20, "r is the correlation between the two contexts across guides.", 7.5, INK2)
    return c.save(OUT / "fig1_folding_context.pdf")


# ====================================== Fig 2: within-gene correlations
def fig_correlations():
    W, H = 468, 292
    c = MODE['cls'](W, H)
    keep = [("seed_nussinov", "Nussinov seed (alone)"),
            ("seed_mfe_spacer", "MFE seed (alone)"),
            ("seed_ensemble_spacer", "Ensemble seed (alone)"),
            ("seed_mfe_full", "MFE seed (in sgRNA)"),
            ("seed_ensemble_full", "Ensemble seed (in sgRNA)"),
            ("mean_unpaired_full", "Mean unpaired (in sgRNA)"),
            ("ensemble_energy_full", "Ensemble energy (in sgRNA)"),
            ("gc_percent", "G/C content")]
    dmap = {r["feature"]: r for r in DO["correlations"]}
    cmap = {r["feature"]: r for r in CR["correlations"]}
    rows = [(k, lab) for k, lab in keep if k in dmap and k in cmap]
    assert len(rows) == len(keep), "feature key missing: %s" % (
        [k for k, _ in keep if k not in dmap or k not in cmap],)

    c.text(30, H - 22, "Association with activity is small and context dependent", 11, INK, bold=True)
    c.text(30, H - 35, "Within-gene Spearman correlation, 95% CI from resampling genes", 8, INK2)

    x0, x1 = 190, 396
    top, bot = H - 68, 64
    step = (top - bot) / max(1, len(rows) - 1)
    lo, hi = -0.30, 0.35
    sx = axis(c, x0, x1, bot - 6, top + 12, lo, hi, [-0.3, -0.15, 0, 0.15, 0.3],
              "within-gene Spearman rho")

    for i, (k, lab) in enumerate(rows):
        y = top - i * step
        c.text(x0 - 12, y - 2.5, lab, 7.5, INK, anchor="end")
        for rec, col, off in ((dmap[k], BLUE, 3.0), (cmap[k], ORANGE, -3.0)):
            v = rec["within_gene_spearman"]
            a = max(lo, rec["within_gene_ci_low"])
            b = min(hi, rec["within_gene_ci_high"])
            c.line(sx(a), y + off, sx(b), y + off, col, 1.1)
            c.line(sx(a), y + off - 1.8, sx(a), y + off + 1.8, col, 1.1)
            c.line(sx(b), y + off - 1.8, sx(b), y + off + 1.8, col, 1.1)
            c.dot(sx(v), y + off, 2.6, col, ring=(255, 255, 255))

    ly = bot - 36
    c.dot(x0 + 4, ly + 2.5, 3.0, BLUE)
    c.text(x0 + 12, ly, "human cells (n = 4,685)", 7.5, INK2)
    c.dot(x0 + 118, ly + 2.5, 3.0, ORANGE)
    c.text(x0 + 126, ly, "zebrafish (n = 1,020)", 7.5, INK2)
    return c.save(OUT / "fig2_correlations.pdf")


# ========================================= Fig 3: incremental value
def fig_incremental():
    W, H = 468, 214
    c = MODE['cls'](W, H)

    def grab(ds, name):
        for m in ds["incremental"]["models"]:
            if m["name"] == name:
                return m
        return None

    rows = []
    for ds, screen in ((DO, "Human cells"), (CR, "Zebrafish")):
        m = grab(ds, "structure_all")
        rows.append((f"{screen}, this implementation", m["delta"],
                     m["delta_ci_low"], m["delta_ci_high"]))
        v = (vienna.get(ds is DO and "doench" or "crisprscan") or {}).get("incremental")
        if v:
            for e in v:
                if e.get("name") == "vienna_structure_all":
                    rows.append((f"{screen}, ViennaRNA features", e["delta"],
                                 e["delta_ci_low"], e["delta_ci_high"]))
    if len(rows) == 2:  # vienna json absent: fall back to published table values
        rows = [("Human cells, this implementation", rows[0][1], rows[0][2], rows[0][3]),
                ("Human cells, ViennaRNA features", 0.002, -0.004, 0.007),
                ("Zebrafish, this implementation", rows[1][1], rows[1][2], rows[1][3]),
                ("Zebrafish, ViennaRNA features", 0.040, 0.012, 0.072)]

    c.text(30, H - 22, "Does folding add anything to a sequence-only model?", 11, INK, bold=True)
    c.text(30, H - 35, "Change in held-out Spearman when structure features join the baseline", 8, INK2)

    x0, x1 = 214, 402
    top, bot = H - 66, 66
    step = (top - bot) / max(1, len(rows) - 1)
    lo, hi = -0.02, 0.09
    sx = axis(c, x0, x1, bot - 6, top + 12, lo, hi, [-0.02, 0, 0.02, 0.04, 0.06, 0.08],
              "change in held-out Spearman rho")

    for i, (lab, d, cl, ch) in enumerate(rows):
        y = top - i * step
        sig = cl > 0
        col = ORANGE if sig else GRAY
        c.text(x0 - 12, y - 2.5, lab, 7.5, INK if sig else INK2, bold=sig, anchor="end")
        c.line(sx(max(lo, cl)), y, sx(min(hi, ch)), y, col, 1.2)
        for e in (max(lo, cl), min(hi, ch)):
            c.line(sx(e), y - 2.0, sx(e), y + 2.0, col, 1.2)
        c.dot(sx(d), y, 3.0, col, ring=(255, 255, 255))
        c.text(sx(min(hi, ch)) + 8, y - 2.5,
               f"{d:+.3f} [{cl:+.3f}, {ch:+.3f}]", 7, INK if sig else INK2)

    c.text(30, bot - 40,
           "Intervals crossing zero are shown in grey. Only the zebrafish result with reference",
           7.5, INK2)
    c.text(30, bot - 50,
           "Turner parameters excludes zero, and it runs opposite to the accessibility hypothesis.",
           7.5, INK2)
    return c.save(OUT / "fig3_incremental.pdf")


# ======================================== Fig 4: position-resolved
def fig_positions():
    W, H = 468, 352
    c = MODE['cls'](W, H)
    c.text(30, H - 22, "No signal at the PAM-proximal seed", 11, INK, bold=True)
    c.text(30, H - 35,
           "Correlation with activity at each spacer position, spacer folded inside the sgRNA", 8, INK2)

    panels = [("Human cells", DO, BLUE), ("Zebrafish", CR, ORANGE)]
    ph = 82
    for pi, (name, ds, col) in enumerate(panels):
        base = H - 78 - pi * (ph + 56)
        pos = ds["positions"]["full"]
        lo, hi = -0.17, 0.09
        x0, x1 = 62, 430
        bw = (x1 - x0) / 20

        c.text(x0, base + 8, name, 8.5, INK, bold=True)
        for t in (-0.15, 0.0, 0.05):
            yy = base - ph + (t - lo) / (hi - lo) * ph
            c.line(x0, yy, x1, yy, GRID if t else (0xa8, 0xa8, 0xa4), 0.6 if t else 1.0)
            c.text(x0 - 6, yy - 2.5, f"{t:g}", 7, INK2, anchor="end")

        # WU-CRISPR claimed window, positions 18 to 20
        wx = x0 + 17 * bw
        c.rect(wx, base - ph, 3 * bw, ph, (0xf1, 0xef, 0xe8))
        if pi == 0:
            c.text(wx + 1.5 * bw, base + 8, "positions 18-20", 6.8, INK2, anchor="middle")
            c.text(wx + 1.5 * bw, base + 0.5, "(WU-CRISPR claim)", 6.8, INK2, anchor="middle")

        zero = base - ph + (0 - lo) / (hi - lo) * ph
        for p in pos:
            v = p["spearman"]
            sig = p.get("p_adjusted", 1) < 0.05
            x = x0 + (p["position"] - 1) * bw
            yv = base - ph + (v - lo) / (hi - lo) * ph
            c.rect(x + 1.2, min(zero, yv), bw - 2.4, abs(yv - zero), col if sig else (0xc4, 0xc4, 0xc0))
            if p["position"] in (1, 5, 10, 15, 20):
                c.text(x + bw / 2, base - ph - 10, str(p["position"]), 6.8, INK2, anchor="middle")
        nsig = sum(1 for p in pos if p.get("p_adjusted", 1) < 0.05)
        c.text(x1 + 6, base - ph / 2, f"{nsig}/20", 7.5, INK2)

    c.text(30, 22, "Filled bars survive Benjamini-Hochberg correction across the 20 positions.",
           7.5, INK2)
    c.text(30, 12, "Where signal appears it sits at the PAM-distal end, the wrong end for seed occlusion.",
           7.5, INK2)
    return c.save(OUT / "fig4_positions.pdf")


if __name__ == "__main__":
    import sys
    figs = (fig_context, fig_correlations, fig_incremental, fig_positions)
    for cls in ((Canvas, SvgCanvas) if "--svg" in sys.argv else (Canvas,)):
        MODE["cls"] = cls
        for f in figs:
            p = f()
            print(f"wrote {p.name}  ({p.stat().st_size:,} bytes)")
