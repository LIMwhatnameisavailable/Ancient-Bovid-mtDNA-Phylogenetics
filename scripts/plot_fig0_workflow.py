#!/usr/bin/env python3
"""Fig0: Workflow overview — Wenjiangduo bovine mitochondrial DNA analysis.

Output: 09_figures/fig0_workflow.png + .svg (300 dpi, 4:3 aspect ratio)
"""

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from pathlib import Path
import datetime

# ─── Paths ───
OUT_DIR = Path('/mnt/disk2/srtp2024/LIM/IVPP/analysis/09_figures')
LOG_DIR = Path('/mnt/disk2/srtp2024/LIM/IVPP/analysis/logs')
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ─── Style setup ───
plt.rcParams.update({
    'font.family':        'serif',
    'font.serif':         ['Times New Roman', 'Georgia', 'DejaVu Serif'],
    'font.size':          10,
    'figure.facecolor':   'white',
    'savefig.facecolor':  'white',
    'axes.unicode_minus': False,
})

# ─── Color constants (user-specified) ───
NAVY       = '#1B3A5C'
AMBER      = '#C97D2E'
AMBER_BG   = '#FDF3E3'
TEAL       = '#2A7F7F'
GRAY_BG    = '#E8E8E8'
GRAY_TEXT  = '#666666'
TEXT_PRI   = '#333333'
WHITE      = '#FFFFFF'
GRAY_BAR   = '#D0D0D0'
ARROW_C    = '#1B3A5C'
NAVY_LIGHT = '#3A6A9A'
BORDER_C   = '#CCCCCC'

VERT_DIV   = '#CCCCCC'  # vertical divider line in dual-column boxes


# ═══════════════════════════════════════════════════════════
#  Figure setup
# ═══════════════════════════════════════════════════════════
W, H = 12, 9
fig, ax = plt.subplots(1, 1, figsize=(W, H))
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.axis('off')
fig.patch.set_facecolor('white')

# ─── Layout parameters ───
BOX_W   = 8.0
BOX_CX  = 6.0
BOX_L   = BOX_CX - BOX_W/2
BOX_R   = BOX_CX + BOX_W/2

H_NORMAL = 0.95
H_STEP2  = 2.40
H_CONCL  = 1.00
GAP      = 0.22
MARGIN_B = 0.40

# Compute box bottom-y positions (bottom to top)
y = MARGIN_B
y6b = y; y6t = y6b + H_CONCL          # Box 6 — Conclusion
y5b = y6t + GAP; y5t = y5b + H_NORMAL # Box 5 — Haplotype
y4b = y5t + GAP; y4t = y4b + H_STEP2  # Box 4 — Phylo Tree (Step 2) – tallest
y3b = y4t + GAP; y3t = y3b + H_NORMAL # Box 3 — Species ID (Step 1)
y2b = y3t + GAP; y2t = y2b + H_NORMAL # Box 2 — QC
y1b = y2t + GAP; y1t = y1b + H_NORMAL # Box 1 — Data Input

boxes = [
    # (yb, yt, label, is_special)
    (y1b, y1t, 'Data Collection & Preparation', False),
    (y2b, y2t, 'Quality Control & Filtering',    False),
    (y3b, y3t, 'Species Identification  (Step 1)', False),
    (y4b, y4t, 'Phylogenetic Analysis  (Step 2)',  True),   # tallest, amber
    (y5b, y5t, 'Haplotype Analysis  (Step 3)',     False),
    (y6b, y6t, 'Key Findings',                     False),  # amber fill
]


# ═══════════════════════════════════════════════════════════
#  Helper functions
# ═══════════════════════════════════════════════════════════

def draw_box(ax, xl, xr, yb, yt, fill_color, edge_color, lw=1.2):
    """Draw a rectangle with optional rounded corners."""
    rect = FancyBboxPatch(
        (xl, yb), xr - xl, yt - yb,
        boxstyle="round,pad=0.04",
        facecolor=fill_color, edgecolor=edge_color, linewidth=lw,
    )
    ax.add_patch(rect)


def draw_step_capsule(ax, cx, y_top, label, fill=TEAL, text_color=WHITE):
    """Draw a small capsule badge above/at the top of a box."""
    fs = 7.5
    txt = ax.text(0, 0, label, fontsize=fs, ha='center', va='center',
                  fontweight='bold', color=text_color)
    bb = txt.get_window_extent().transformed(ax.transData.inverted())
    pad_x = 0.1
    pad_y = 0.04
    cw = bb.width + 2 * pad_x
    ch = bb.height + 2 * pad_y
    # Place centered at cx, slightly below y_top so it's inside the box
    cy = y_top - 0.01
    cap = FancyBboxPatch(
        (cx - cw / 2, cy - ch + 0.02),
        cw, ch,
        boxstyle="round,pad=0.02",
        facecolor=fill, edgecolor='none', zorder=5,
    )
    ax.add_patch(cap)
    txt.set_position((cx, cy - ch / 2 + 0.03))
    # Return the capsule for potential reuse
    return cap, txt


def draw_arrow_down(ax, cx, y_from, y_to, color=ARROW_C, lw=1.5):
    """Draw downward arrow between two boxes."""
    ax.annotate(
        '', xy=(cx, y_to), xytext=(cx, y_from),
        arrowprops=dict(arrowstyle='->', color=color, lw=lw),
    )


def draw_vertical_divider(ax, cx, yb, yt, color=VERT_DIV, lw=0.8):
    """Draw a thin vertical line inside a box to separate two columns."""
    ax.plot([cx, cx], [yb, yt], color=color, linewidth=lw, solid_capstyle='round')


def draw_multi_text(ax, lines, x, y, fontsize=9, color=TEXT_PRI, line_sp=0.16,
                    ha='left', va='top', weight='normal'):
    """Draw multiple lines of text starting from (x, y)."""
    for i, line in enumerate(lines):
        ax.text(x, y - i * line_sp, line, fontsize=fontsize, color=color,
                ha=ha, va=va, fontweight=weight)


def draw_centered_multi_text(ax, lines, x, y, fontsize=9, color=TEXT_PRI,
                             line_sp=0.16, va='top', weight='normal'):
    """Draw centered multi-line text."""
    for i, line in enumerate(lines):
        ax.text(x, y - i * line_sp, line, fontsize=fontsize, color=color,
                ha='center', va=va, fontweight=weight)


# ═══════════════════════════════════════════════════════════
#  Draw arrows between boxes
# ═══════════════════════════════════════════════════════════
draw_arrow_down(ax, BOX_CX, y1b, y2t)
draw_arrow_down(ax, BOX_CX, y2b, y3t)
draw_arrow_down(ax, BOX_CX, y3b, y4t)
draw_arrow_down(ax, BOX_CX, y4b, y5t)
draw_arrow_down(ax, BOX_CX, y5b, y6t)


# ═══════════════════════════════════════════════════════════
#  Box 1 — Data Input
# ═══════════════════════════════════════════════════════════
draw_box(ax, BOX_L, BOX_R, y1b, y1t, 'white', NAVY, 1.2)
draw_step_capsule(ax, BOX_CX, y1t, 'STEP 0', TEAL)

draw_multi_text(ax, [
    'Source: 18 ancient bovine mtDNA genomes  |  Wenjiangduo site (4th-5th cent.)',
    'Two sub-sites:  Tower context (n=2)  +  Main context (n=16)',
    'Sequence length: 16,322-16,341 bp  |  18 FASTA → merged alignment',
], BOX_L + 0.3, y1t - 0.30, fontsize=9.5, line_sp=0.22)


# ═══════════════════════════════════════════════════════════
#  Box 2 — QC
# ═══════════════════════════════════════════════════════════
draw_box(ax, BOX_L, BOX_R, y2b, y2t, 'white', NAVY, 1.2)
draw_step_capsule(ax, BOX_CX, y2t, 'QC', TEAL)

qc_left_x = BOX_L + 0.3
qc_center = BOX_CX
qc_right_x = BOX_R - 0.3

draw_multi_text(ax, [
    'Coverage stratification:',
    '  HIGH (≥100×):   8 samples  —  confident calls',
    '  MID  (20-99×):  8 samples  —  adequate for analysis',
    '  LOW  (<20×):    2 samples  —  excluded from pop. stats',
], qc_left_x, y2t - 0.28, fontsize=9, line_sp=0.19)

draw_multi_text(ax, [
    'N% warning (≥5%): 3 samples',
    '  Z9853 (30.0%), Z9859 (8.7%), Z9865 (7.4%)',
    'Removed: 6 lowQ from ref panel (sensitivity test)',
    '→ 16 samples for downstream (excl. LOW)',
], qc_center + 0.5, y2t - 0.28, fontsize=9, line_sp=0.19)


# ═══════════════════════════════════════════════════════════
#  Box 3 — Species Identification (Step 1)  — dual column
# ═══════════════════════════════════════════════════════════
draw_box(ax, BOX_L, BOX_R, y3b, y3t, 'white', NAVY, 1.2)
draw_step_capsule(ax, BOX_CX, y3t, 'STEP 1', TEAL)

# Left column — BLAST
col3_w = BOX_W / 2 - 0.6
col3_lx = BOX_L + 0.3
col3_rx = BOX_CX + 0.3

draw_multi_text(ax, [
    'BLAST (NCBI nt database)',
    '───',
    'Yak (n=7):',
    '  Z9853*, Z9854-56, Z9858, Z9863',
    'Cattle (n=11):',
    '  K0909, Z9850-52, Z9857, Z9859-61,',
    '  Z9864*, Z9865, Z9867',
    '  *Low coverage — ID confirmed but caution',
], col3_lx, y3t - 0.28, fontsize=8.8, line_sp=0.16)

# Right column — ML Identification Tree
draw_multi_text(ax, [
    'ML Identification Tree (IQ-TREE2)',
    '───',
    '22 seq × 15,790 bp | TIM2+F+G4 model',
    'Bovini root → Bos clade (BS=100)',
    '  Cattle clade (BS=100, 11 tips)',
    '  Yak clade   (BS=100, 7 tips)',
    'Discordant candidates identified:',
    '  Z9863 (morph=cattle) → Yak BS=100',
    '  Z9867 (morph=yak)    → Cattle BS=100',
], col3_rx, y3t - 0.28, fontsize=8.8, line_sp=0.16)

# Divider line
draw_vertical_divider(ax, BOX_CX, y3b + 0.05, y3t - 0.07, VERT_DIV)


# ═══════════════════════════════════════════════════════════
#  Box 4 — Phylogenetic Analysis (Step 2)  — HIGHLIGHT ★
# ═══════════════════════════════════════════════════════════
draw_box(ax, BOX_L, BOX_R, y4b, y4t, AMBER_BG, AMBER, 1.5)
draw_step_capsule(ax, BOX_CX, y4t, 'STEP 2  [*]', TEAL)

# ── Title ──
ax.text(BOX_CX, y4t - 0.22, 'Phylogenetic Analysis (Reference Panel & Fine-Grained Tree)',
        fontsize=11, ha='center', va='center', fontweight='bold', color=TEXT_PRI)

# ── Reference panel table (compact) ──
table_top = y4t - 0.42
table_bot = y4t - 1.35
table_cx = BOX_CX
col_xs = [BOX_L + 0.3, BOX_L + 1.3, BOX_L + 1.7, BOX_R - 0.3]
col_ws = [col_xs[1] - col_xs[0], col_xs[2] - col_xs[1], col_xs[3] - col_xs[2]]
row_h = 0.14
n_rows = 9  # header + 7 data + total

rows_data = [
    ('Source', 'n', 'Description'),
    ('REF',     '1',  'African buffalo outgroup'),
    ('CORE',    '2',  'Water buffalo + Indian I1 cattle'),
    ('LIT1',    '1',  'Ancient yak YakA (Chen 2023, Sci Adv)'),
    ('LIT2',    '9',  'Ancient YakX (Gilardet & Oppenheimer 2025)'),
    ('LIT3',    '13', 'P1a cattle ×11 + European aurochs P ×2'),
    ('LIT4',    '15', 'Ancient cattle/aurochs C (Cai & Kim 2025, Science)'),
    ('MOD',     '10', 'Modern wild yak ×5 + Qaidam T3 ×5'),
    ('TOTAL',   '51', 'Multi-source panel: 6 sources + 2 outgroups'),
]

for i, (src, n, desc) in enumerate(rows_data):
    y_row_top = table_top - i * row_h
    y_row_bot = y_row_top - row_h + 0.02
    y_mid = (y_row_top + y_row_bot) / 2

    # Header styling
    is_header = (i == 0)
    is_total  = (i == len(rows_data) - 1)
    fw = 'bold' if (is_header or is_total) else 'normal'
    fs = 7.0 if not is_header else 7.2
    fc = TEXT_PRI
    bg = None

    if is_header:
        # Draw a thin line under header instead of background
        ax.plot([col_xs[0]-0.05, col_xs[-1]+0.05], [y_row_bot, y_row_bot],
                color=NAVY, linewidth=0.5)
    if is_total:
        ax.plot([col_xs[0]-0.05, col_xs[-1]+0.05], [y_row_top, y_row_top],
                color=NAVY, linewidth=0.5, linestyle='--')
        fc = NAVY

    ax.text(col_xs[0] + 0.02, y_mid, src, fontsize=fs, ha='left', va='center',
            fontweight=fw, color=fc)
    ax.text((col_xs[1] + col_xs[2]) / 2, y_mid, n, fontsize=fs, ha='center', va='center',
            fontweight=fw, color=fc)
    ax.text(col_xs[2] + 0.02, y_mid, desc, fontsize=fs, ha='left', va='center',
            fontweight=fw, color=fc)

# Caption for table
ax.text(col_xs[0], table_bot - 0.08,
        '[*] Multi-source ancient panel -- core innovation of this study',
        fontsize=7.5, ha='left', va='top', fontweight='bold', color=AMBER)

# ── Results section between table and gray bar ──
res_top = table_bot - 0.28
res_bot = y4b + 0.35  # top of gray bar area

res_lines = [
    '69 seq (18 WJD + 51 ref)  ×  16,338 bp  |  Model: TIM2+F+I+G4',
    'Bovini root → Bos (BS=100) → Cattle clade (BS=100, 45 tips) + Yak clade (BS=92, 22 tips)',
]
# Indented sub-results
res_lines2 = [
    'Z9863 → Yak clade BS=92   |   Z9867 → Cattle clade BS=100   |   Discordant confirmed',
    'K0914 -> LIT2 YakX clade BS=99  |  Links WJD to Qinghai-Tibet Plateau ancient yak lineage',
    'Two outgroups: African buffalo (REF) + Water buffalo (CORE) — robust rooting',
]

draw_multi_text(ax, res_lines, BOX_L + 0.3, res_top,
                fontsize=8.8, line_sp=0.19)
draw_multi_text(ax, res_lines2, BOX_L + 0.3, res_top - 0.24,
                fontsize=8.5, line_sp=0.18)

# ── Gray bar at bottom: sensitivity test ──
bar_top = y4b + 0.32
bar_bot = y4b + 0.04
bar_rect = FancyBboxPatch(
    (BOX_L + 0.1, bar_bot), BOX_W - 0.2, bar_top - bar_bot,
    boxstyle="round,pad=0.02",
    facecolor=GRAY_BG, edgecolor=GRAY_TEXT, linewidth=0.8,
)
ax.add_patch(bar_rect)

ax.text(BOX_CX, (bar_top + bar_bot) / 2,
        'Sensitivity Test (remove 6 lowQ seq): Topology unchanged  |  Z9863->Yak BS=93  |  Z9867->Cattle BS=100  |  K0914->YakX BS=99 [OK]',
        fontsize=8, ha='center', va='center', color=TEXT_PRI, fontweight='bold')


# ═══════════════════════════════════════════════════════════
#  Box 5 — Haplotype Analysis (Step 3)  — dual column
# ═══════════════════════════════════════════════════════════
draw_box(ax, BOX_L, BOX_R, y5b, y5t, 'white', NAVY, 1.2)
draw_step_capsule(ax, BOX_CX, y5t, 'STEP 3', TEAL)

col5_lx = BOX_L + 0.3
col5_rx = BOX_CX + 0.3

draw_multi_text(ax, [
    'Haplotype Diversity (16 samples)',
    '───',
    'Gap/ambig removal → 10,668 bp clean alignment',
    'Total haplotypes:  13  (of 16 samples)',
    'Haplotype diversity (Hd):  0.9667',
    'Sub-site sharing:  0 shared (Tower=2, Main=11)',
], col5_lx, y5t - 0.28, fontsize=8.8, line_sp=0.16)

draw_multi_text(ax, [
    'By mtDNA Lineage (valid π)',
    '───',
    'Yak group (n=6):',
    '  6 haplotypes, Hd=1.000, π=0.21%',
    'Cattle group (n=10):',
    '  9 haplotypes, Hd=0.978, π=0.058%',
    '  (T3 lineage — low diversity, East Asian bottleneck)',
    'Network → Fig8 (MST + spring layout)',
], col5_rx, y5t - 0.28, fontsize=8.8, line_sp=0.16)

# Divider
draw_vertical_divider(ax, BOX_CX, y5b + 0.05, y5t - 0.07, VERT_DIV)


# ═══════════════════════════════════════════════════════════
#  Box 6 — Key Findings (amber fill, white text)
# ═══════════════════════════════════════════════════════════
draw_box(ax, BOX_L, BOX_R, y6b, y6t, AMBER, AMBER, 1.5)
draw_step_capsule(ax, BOX_CX, y6t, 'CONCLUSION', WHITE)

conclusions = [
    '(1)  7 Yak + 11 Cattle mtDNA lineages identified among 18 WJD samples  -  mixed pastoral system at Wenjiangduo',
    '(2)  2 discordant candidates (Z9863, Z9867) confirmed by BLAST + ML tree + RF k-mer  -  robust cross-method validation',
    '(3)  K0914 -> YakX clade (BS=99) links Wenjiangduo to ancient Qinghai-Tibet Plateau yak lineage (Gilardet & Oppenheimer 2025)',
]

for i, line in enumerate(conclusions):
    ax.text(BOX_CX, y6t - 0.28 - i * 0.22, line,
            fontsize=9.5, ha='center', va='top', color=WHITE, fontweight='normal')


# ═══════════════════════════════════════════════════════════
#  Save
# ═══════════════════════════════════════════════════════════
stem = 'fig0_workflow'
png_path = OUT_DIR / f'{stem}.png'
svg_path = OUT_DIR / f'{stem}.svg'

fig.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(svg_path, format='svg', bbox_inches='tight', facecolor='white')
plt.close(fig)

print(f'DONE: {png_path}')
print(f'DONE: {svg_path}')

# ─── Log ───
log_path = LOG_DIR / f'{stem}.log'
with open(log_path, 'w') as f:
    f.write(f'# {stem}.log\n')
    f.write(f'# Generated: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n')
    f.write(f'# Script: inline (CC)\n')
    f.write('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
    f.write('Input sources:\n')
    f.write('  - CLAUDE.md / ANALYSIS_STEPS.md\n')
    f.write('  - 02_ref/README.md (v2.0, 51-seq panel)\n')
    f.write('  - logs/step0_prepare.log\n')
    f.write('  - logs/step3_phylo.log\n')
    f.write('  - logs/step4_haplotype.log\n')
    f.write('  - 09_figures/PLOT_GUIDE.md\n')
    f.write('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
    f.write('Figure specs:\n')
    f.write('  - Size: 12×9 in (4:3)\n')
    f.write('  - DPI: 300 (PNG), vector (SVG)\n')
    f.write('  - Font: Times New Roman → Georgia → DejaVu Serif\n')
    f.write('  - Layout: 6 equal-width boxes, top-to-bottom\n')
    f.write('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
    f.write('Box contents:\n')
    f.write('  Box1 (Data Input): 18 WJD mtDNA, 16,322-16,341 bp, 2 sub-sites\n')
    f.write('  Box2 (QC): HIGH=8, MID=8, LOW=2, WARN=3\n')
    f.write('  Box3 (Species ID Step1): BLAST(7Y/11C) + ML tree(BS=100); discordant identified\n')
    f.write('  Box4 (Phylo Step2): 69 seq, TIM2+F+I+G4; 51-ref panel table; K0914→YakX BS=99\n')
    f.write('  Box5 (Haplotype Step3): 13 hap, Hd=0.967; Yak π=0.21%, Cattle π=0.058%\n')
    f.write('  Box6 (Conclusion): 3 core findings, amber fill\n')
    f.write('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
    f.write('DONE\n')

print(f'DONE: {log_path}')
