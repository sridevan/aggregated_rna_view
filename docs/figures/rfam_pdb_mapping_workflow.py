import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

out = sys.argv[1]
W, H = 24, 13.5
fig = plt.figure(figsize=(W, H), dpi=150)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
fig.patch.set_facecolor("white")

INK, MUTED = "#1f2937", "#57534e"
PDBE = ("#e6f2f1", "#2f7f79")      # fill, edge
RFAM = ("#fdf3e3", "#b7791f")
STEP = ("#eef2f7", "#475569")
OUTC = ("#eaf5ec", "#3b7a4a")
NOTE = ("#ffffff", "#a8a29e")

def box(x, y, w, h, title, body, col, title_size=15, body_size=11.5, dashed=False):
    fill, edge = col
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.18", fc=fill, ec=edge, lw=1.8,
                                ls=(0, (4, 3)) if dashed else "-"))
    ax.text(x + 0.22, y + h - 0.28, title, ha="left", va="top", fontsize=title_size, fontweight="bold", color=INK)
    if body:
        ax.text(x + 0.22, y + h - 0.28 - 0.42 * (title_size / 15) - 0.12, body, ha="left", va="top", fontsize=body_size, color=MUTED, linespacing=1.45)

def arrow(p, q, color=INK, dashed=False, rad=0.0, lw=1.8, label=None, lpos=0.5, loff=(0, 0)):
    a = FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=18, lw=lw, color=color, connectionstyle=f"arc3,rad={rad}",
                        ls=(0, (5, 3)) if dashed else "-", shrinkA=2, shrinkB=2)
    ax.add_patch(a)
    if label:
        mx, my = p[0] + (q[0] - p[0]) * lpos + loff[0], p[1] + (q[1] - p[1]) * lpos + loff[1]
        ax.text(mx, my, label, fontsize=10.5, color=color, ha="center", va="center", style="italic",
                bbox=dict(fc="white", ec="none", pad=1.5))

# ---------------------------------------------------------------- title
ax.text(0.6, H - 0.55, "Mapping PDB RNA chains to Rfam consensus positions", fontsize=24, fontweight="bold", color=INK, va="top")
ax.text(0.6, H - 1.2, "Chain regions are Infernal search hits (run by the Rfam team); we re-align each hit region so every nucleotide gets a consensus position or a label.",
        fontsize=13.5, color=MUTED, va="top")

# ---------------------------------------------------------------- layout: six step columns
tw, gap = 3.55, 0.37
xs = [0.6 + i * (tw + gap) for i in range(6)]
cx = [x + tw / 2 for x in xs]

# ---------------------------------------------------------------- sources (each above the step that reads it, so arrows drop straight down)
sy, sh = 9.55, 1.75
ax.text(0.6, sy + sh + 0.22, "SOURCES  (cached; retrieval date and URL recorded)", fontsize=11.5, color=MUTED, fontweight="bold")
box(xs[0], sy, tw, sh, "Infernal hits, via PDBe", "Rfam CSV of significant hits:\none family + region per\n(author) chain", RFAM, title_size=13.5, body_size=11)
box(xs[1], sy, tw, sh, "PDBe molecules, residues", "Entity sequence, author\nnumbering, insertion codes", PDBE, title_size=13.5, body_size=11)
box(xs[2], sy, tw, sh, "Rfam 15.1 family, CM", "Covariance model, CLEN,\nchecksum, clan", RFAM, title_size=13.5, body_size=11)
box(xs[3], sy, 2 * tw + gap, sh, "PDBe updated mmCIF  /  ModelServer per chain",
    "_pdbx_poly_seq_scheme and atom_site: identifiers and which residues\nhave coordinates (curated sets: one chain, model 1, never full mmCIF)", PDBE, title_size=13.5, body_size=11)
box(xs[5], sy, tw, sh, "Rfam family/structures", "The same hits as Rfam\npublishes them: region,\nCM span, bit score", RFAM, title_size=13.5, body_size=11)

# ---------------------------------------------------------------- steps
ty, th = 6.6, 1.95
ax.text(0.6, ty + th + 0.22, "WORKFLOW  (python3 -m rna_explorer build <ACC>)", fontsize=11.5, color=MUTED, fontweight="bold", zorder=10, bbox=dict(fc="white", ec="none", pad=2))
steps = [
    ("1  Select chains", "One instance per\nInfernal hit region; curated:\n(pdb, label_asym_id) list"),
    ("2  Extract sequence", "Hit region only;\nmodified residues as\nparent base"),
    ("3  Re-align region", "Unique sequences (md5);\ncmalign local mode,\nInfernal 1.1.5"),
    ("4  Assign states", "One state per residue and\nper consensus column;\njoin back to every chain"),
    ("5  Validate", "5 pass/fail checks;\nfailing chains kept,\nflagged ok: false"),
    ("6  Aggregate", "Per consensus position:\nunique sequences | chains\nwith coordinates"),
]
for x, (t, b) in zip(xs, steps):
    box(x, ty, tw, th, t, b, STEP, title_size=14.5, body_size=11.5)
for i in range(5):
    arrow((xs[i] + tw, ty + th / 2), (xs[i + 1], ty + th / 2), lw=2.2)
for i in (0, 1, 2):
    arrow((cx[i], sy), (cx[i], ty + th))
arrow((cx[3], sy), (cx[3], ty + th), label="coordinates", loff=(0, 0))
arrow((cx[4], sy), (cx[4], ty + th), label="identifiers", loff=(0, 0))
arrow((cx[5] - 0.3, sy), (cx[4] + 1.2, ty + th), color=RFAM[1], dashed=True, label="cross-check only", lpos=0.5, loff=(0.25, 0.05))

# ---------------------------------------------------------------- detail row
ny, nh = 3.55, 2.3
box(xs[0], ny, 3 * tw + 2 * gap, nh, "Example: 1Y26 chain X (adenine riboswitch, RF00167, CLEN 102)",
    "Infernal hit region label_seq_id 3–69  →  CM span 19–87, average PP 0.999\n"
    "71 nucleotides = 67 match + 4 flank (2 at each end)\n"
    "102 columns = 67 match + 2 deletion (64, 72) + 33 uncovered (1–18, 88–102)\n"
    "consensus 50  →  label_seq_id 34  →  author G46, coordinates present", NOTE, title_size=13, body_size=11.5)
box(xs[3], ny, tw, nh, "Residue states",
    "match → position 1..CLEN\ninsert → anchor + offset\nlocal_end → anchor + offset\nnot_submitted_flank\n+ coordinates flag (≥ 1 atom)", NOTE, title_size=12.5, body_size=11)
box(xs[4], ny, tw, nh, "Checks",
    "ids = _pdbx_poly_seq_scheme\ncoordinate flags agree\npositions monotonic\nall CLEN columns, once\ninsert offsets consistent", NOTE, title_size=12.5, body_size=11)
for i in (3, 4):
    ax.plot([cx[i], cx[i]], [ny + nh, ty], color=MUTED, lw=1.3, ls=(0, (2, 2)))

# ---------------------------------------------------------------- outputs band
oy, oh = 0.45, 2.3
ow = xs[5] + tw - 0.6
box(0.6, oy, ow, oh, "Outputs  data/families/<ACC>/", "", OUTC, title_size=13.5)
rows = [("mappings/<seq_name>.json", "residues[] and consensus_columns[] per chain"),
        ("validation.json", "per-chain checks; Rfam region and CM-span cross-checks"),
        ("consensus_stats.json", "per-position coverage (two denominators), mean PP"),
        ("family.json, instances.json", "chains and provenance: Rfam release, CM checksum, Infernal, sources, dates")]
for k, (f, d) in enumerate(rows):
    yy = oy + oh - 0.95 - k * 0.36
    ax.text(0.85, yy, f, fontsize=11.5, color=INK, family="monospace", va="center")
    ax.text(5.6, yy, d, fontsize=11.5, color=MUTED, va="center")
ax.text(15.3, oy + oh - 0.95, "Used by (all on consensus positions)", fontsize=11.5, color=INK, fontweight="bold", va="center")
ax.text(15.3, oy + oh - 1.3, "coverage tracks and 2D/3D selection sync;\nloops, ligand contacts, SASA, pockets;\ninteractions; structural clustering;\nresidue-number relabelling (display only)", fontsize=11.5, color=MUTED, va="top", linespacing=1.45)
arrow((cx[5], ty), (cx[5], oy + oh), lw=2.2)

fig.savefig(out, dpi=150, facecolor="white")
print("saved", out)
