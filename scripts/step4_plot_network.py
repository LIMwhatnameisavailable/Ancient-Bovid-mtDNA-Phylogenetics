#!/usr/bin/env python3
"""
Fig8: Haplotype network for Wenjiangduo ancient bovid mitogenomes.
Output: 09_figures/fig8_haplotype_network.png / .svg
"""
import os
import math
from collections import OrderedDict
from itertools import combinations

import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE    = "/mnt/disk2/srtp2024/LIM/IVPP/analysis"
ALIGN   = os.path.join(BASE, "04_align", "final_aln.fasta")
HAP_CSV = os.path.join(BASE, "06_haplo", "haplotype_assignment.csv")
OUT_DIR = os.path.join(BASE, "09_figures")
os.makedirs(OUT_DIR, exist_ok=True)

plt.rcParams.update({
    'font.family':        'serif',
    'font.serif':         ['Times New Roman', 'Georgia', 'DejaVu Serif'],
    'mathtext.fontset':   'stix',
    'axes.unicode_minus': False,
    'font.size':          11,
    'figure.facecolor':   'white',
    'savefig.facecolor':  'white',
    'figure.dpi':         100,
    'savefig.dpi':        300,
    'axes.grid':          False,
})

LOW_COV    = {"Z9853", "Z9864"}
TOWER_IDS  = {"K0909", "K0914"}
COLOR_MAIN  = "#4A90D9"
COLOR_TOWER = "#2A9D8F"
COLOR_BOTH  = "#7B6FAB"
COLOR_LOW   = "#AAAAAA"

# ─── 1. Read sequences ───────────────────────────────────────────────────────
with open(ALIGN) as f:
    raw = f.read().strip().split(">")

seqs = OrderedDict()
for block in raw:
    if not block.strip():
        continue
    id_line, *seq_lines = block.strip().split("\n")
    sid = id_line.split()[0]
    seq = "".join(seq_lines).upper()
    if not any(p in sid for p in ["MOD_", "LIT", "CORE_", "REF_", "all_combined"]):
        seqs[sid] = seq

print(f"Read {len(seqs)} WJD sequences")

# ─── 2. Haplotype assignment ─────────────────────────────────────────────────
samp_hap = {}
with open(HAP_CSV) as f:
    next(f)
    for line in f:
        parts = line.strip().split(",")
        samp_hap[parts[0]] = parts[1]

# ─── 3. Filter clean columns (vectorised) ────────────────────────────────────
AMBIG   = set("NRYWSMKBDHV")
BAD     = set('-.N?') | AMBIG
aln_len = len(next(iter(seqs.values())))

seq_matrix = np.array([[c for c in seqs[s]] for s in seqs])
bad_mask   = np.zeros(aln_len, dtype=bool)
for ch in BAD:
    bad_mask |= (seq_matrix == ch).any(axis=0)
keep_cols = np.where(~bad_mask)[0]
print(f"Clean columns: {len(keep_cols)} / {aln_len}")

# ─── 4. Build haplotype sequences ────────────────────────────────────────────
hap_seqs     = {}
hap_site     = {}
hap_members  = {}
hap_excluded = {}
hap_label    = {}   # ← 新增：存储双行标签文本

for hid in sorted(set(samp_hap.values())):
    members = [s for s in seqs if samp_hap.get(s) == hid]
    hap_members[hid]  = members
    rep               = members[0]
    seq_arr           = np.array(list(seqs[rep]))
    hap_seqs[hid]     = "".join(seq_arr[keep_cols])
    sites = {"Tower" if m in TOWER_IDS else "Main" for m in members}
    hap_site[hid]     = "|".join(sorted(sites))
    hap_excluded[hid] = all(m in LOW_COV for m in members)
    # ── 构建标签：单倍型 ID（粗体行）+ 样本编号（细体行）────────────────
    # 多个样本按字母序排列，用 "/" 分隔
    sample_str = "/".join(sorted(members))
    hap_label[hid] = f"{hid}\n{sample_str}"

# ─── 5. Pairwise Hamming distances ───────────────────────────────────────────
hids = sorted(hap_seqs.keys())
n    = len(hids)

seq_int = np.array([[ord(c) for c in hap_seqs[h]] for h in hids], dtype=np.int32)
dist    = np.zeros((n, n), dtype=int)
for i, j in combinations(range(n), 2):
    d          = int(np.sum(seq_int[i] != seq_int[j]))
    dist[i, j] = d
    dist[j, i] = d

print(f"\nPairwise distances ({n} haplotypes):")
for i in range(n):
    for j in range(i + 1, n):
        print(f"  {hids[i]} vs {hids[j]}: {dist[i,j]} bp")

# ─── 6. Minimum spanning tree ────────────────────────────────────────────────
G = nx.Graph()
for hid in hids:
    G.add_node(hid,
               size     = len(hap_members[hid]),
               site     = hap_site[hid],
               excluded = hap_excluded[hid])

edge_list = sorted((int(dist[i, j]), hids[i], hids[j])
                   for i, j in combinations(range(n), 2))

parent = {h: h for h in hids}
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x
def union(x, y):
    px, py = find(x), find(y)
    if px != py:
        parent[px] = py
        return True
    return False

for d, h1, h2 in edge_list:
    if union(h1, h2):
        G.add_edge(h1, h2, weight=d)

# ─── 7. Layout — log-compressed weights ──────────────────────────────────────
G_layout = nx.Graph()
for hid in hids:
    G_layout.add_node(hid)
for u, v, w in G.edges(data='weight'):
    G_layout.add_edge(u, v, weight=math.log1p(w))

pos_raw = nx.kamada_kawai_layout(G_layout, weight='weight')
scale   = 3.2
pos     = {node: (xy[0] * scale, xy[1] * scale) for node, xy in pos_raw.items()}

# ─── 8. Node radii ───────────────────────────────────────────────────────────
node_radius = {hid: 0.10 + math.sqrt(G.nodes[hid]['size']) * 0.048
               for hid in hids}

# ─── 9. 8-direction label system ─────────────────────────────────────────────
_D45 = math.sqrt(2) / 2
DIRECTION_VEC = {
    'N':  ( 0.0,  1.0),
    'S':  ( 0.0, -1.0),
    'E':  ( 1.0,  0.0),
    'W':  (-1.0,  0.0),
    'NE': ( _D45,  _D45),
    'NW': (-_D45,  _D45),
    'SE': ( _D45, -_D45),
    'SW': (-_D45, -_D45),
}

LEADER_LEN = 0.38

LABEL_DIR = {
    "H01":  'S',
    "H08":  'W',
    "H07":  'N',
    "H13":  'NE',
    "H04":  'SW',
    "H06":  'W',
    "H12":  'N',
    "H14":  'S',
    "H09":  'SW',
    "H10":  'SE',
    "H02":  'NW',
    "H11":  'NW',
    "H03":  'NE',
    "H15":  'E',
    "H05":  'SE',
}

def compute_label_positions_8dir(pos, node_radius, hids, leader_len, label_dir):
    lpos = {}
    for hid in hids:
        x, y = pos[hid]
        r    = node_radius[hid]
        if hid in label_dir:
            dvec = DIRECTION_VEC[label_dir[hid]]
        else:
            rx, ry = 0.0, 0.0
            for nb in hids:
                if nb == hid:
                    continue
                dx2, dy2 = x - pos[nb][0], y - pos[nb][1]
                d2 = dx2**2 + dy2**2 + 1e-6
                rx += dx2 / d2
                ry += dy2 / d2
            norm = math.hypot(rx, ry)
            if norm > 1e-6:
                rx /= norm; ry /= norm
            best, best_dot = 'N', -999
            for dname, (dx2, dy2) in DIRECTION_VEC.items():
                dot = rx * dx2 + ry * dy2
                if dot > best_dot:
                    best_dot = dot
                    best = dname
            dvec = DIRECTION_VEC[best]
        dist_total = r + leader_len
        lpos[hid] = (x + dvec[0] * dist_total,
                     y + dvec[1] * dist_total)
    return lpos

label_pos = compute_label_positions_8dir(
    pos, node_radius, hids, LEADER_LEN, LABEL_DIR
)

# ─── 10. Shared draw function ─────────────────────────────────────────────────
def draw_network(ax):
    all_x = [p[0] for p in pos.values()] + [p[0] for p in label_pos.values()]
    all_y = [p[1] for p in pos.values()] + [p[1] for p in label_pos.values()]

    PAD_LEFT   = 0.2
    PAD_RIGHT  = 0.2
    PAD_TOP    = 0.25
    PAD_BOTTOM = 0.15

    ax.set_xlim(min(all_x) - PAD_LEFT,   max(all_x) + PAD_RIGHT)
    ax.set_ylim(min(all_y) - PAD_BOTTOM,  max(all_y) + PAD_TOP)
    ax.set_aspect('equal')
    ax.axis('off')

    # Edges
    for u, v, w in G.edges(data='weight'):
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        ax.plot([x0, x1], [y0, y1],
                color='#888888', linewidth=1.3, alpha=0.65, zorder=1,
                solid_capstyle='round')
        mx, my = (x0 + x1) / 2, (y0 + y1) / 2
        dx, dy = x1 - x0, y1 - y0
        length = max(1e-6, math.hypot(dx, dy))
        px, py = -dy / length, dx / length
        ax.text(mx + px * 0.18, my + py * 0.18,
                str(int(w)), fontsize=10, color='#333333',
                ha='center', va='center', zorder=4,
                bbox=dict(boxstyle='round,pad=0.10', facecolor='white',
                          edgecolor='none', alpha=0.80))

    # Nodes
    for hid in hids:
        x, y     = pos[hid]
        site     = G.nodes[hid]['site']
        excluded = G.nodes[hid]['excluded']
        r        = node_radius[hid]
        if excluded:
            fc, ec, lw, ls = COLOR_LOW,   '#666666', 2.0, '--'
        elif site == "Tower":
            fc, ec, lw, ls = COLOR_TOWER, '#333333', 1.5, '-'
        elif site == "Main":
            fc, ec, lw, ls = COLOR_MAIN,  '#333333', 1.5, '-'
        else:
            fc, ec, lw, ls = COLOR_BOTH,  '#333333', 1.5, '-'
        ax.add_patch(plt.Circle((x, y), r,
                                facecolor=fc, edgecolor=ec,
                                linewidth=lw, linestyle=ls,
                                alpha=0.90, zorder=5))

    # Leader lines + labels（双行：单倍型ID + 样本编号）
    label_bbox_r = 0.20   # 双行标签 bbox 更高，适当增大缩进半径
    for hid in hids:
        x, y   = pos[hid]
        lx, ly = label_pos[hid]
        r      = node_radius[hid]
        dx, dy = lx - x, ly - y
        dn     = math.hypot(dx, dy)
        # 引导线起点：节点边缘
        sx = x + dx / dn * (r + 0.03) if dn > 1e-6 else x
        sy = y + dy / dn * (r + 0.03) if dn > 1e-6 else y
        # 引导线终点：标签 bbox 边缘
        dl = math.hypot(lx - sx, ly - sy)
        if dl > label_bbox_r:
            ex = lx - (lx - sx) / dl * label_bbox_r
            ey = ly - (ly - sy) / dl * label_bbox_r
        else:
            ex, ey = lx, ly
        ax.plot([sx, ex], [sy, ey],
                color='#BBBBBB', linewidth=0.7, alpha=0.7, zorder=3)

        # ── 双行标签：用 ax.text 的 multialignment 分别渲染两行 ──────────
        label_text = hap_label[hid]   # "H01\nZ9854/Z9856/Z9858"
        ax.text(lx, ly, label_text,
                fontsize=8, color='#1a1a1a',
                ha='center', va='center', zorder=7,
                multialignment='center',
                linespacing=1.4,
                fontweight='semibold',
                bbox=dict(boxstyle='round,pad=0.30', facecolor='white',
                          edgecolor='#CCCCCC', alpha=0.85, linewidth=0.7))

    # Legend
    legend_handles = [
        plt.Line2D([0], [0], marker='o', color='w',
                   markerfacecolor=COLOR_MAIN, markeredgecolor='#333333',
                   markeredgewidth=1.2, markersize=11, label='Ancient_WJD_Main'),
        plt.Line2D([0], [0], marker='o', color='w',
                   markerfacecolor=COLOR_TOWER, markeredgecolor='#333333',
                   markeredgewidth=1.2, markersize=11, label='Ancient_WJD_Tower'),
        plt.Line2D([0], [0], marker='o', color='w',
                   markerfacecolor=COLOR_LOW, markeredgecolor='#666666',
                   markeredgewidth=1.8, markersize=11,
                   label='LOW coverage (excl. from stats)'),
    ]
    ax.legend(handles=legend_handles, loc='upper right',
              framealpha=0.95, edgecolor='#CCCCCC', fontsize=10,
              borderpad=0.8, handletextpad=0.6)
    ax.set_title('Fig8 · Haplotype network — Wenjiangduo ancient Bovini',
                 fontsize=13, fontweight='bold', pad=10, color='#222222')

# ─── 11. Save PNG ─────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(12, 10))
fig.subplots_adjust(left=0.02, right=0.98, top=0.96, bottom=0.02)
draw_network(ax)
out_png = os.path.join(OUT_DIR, "fig8_haplotype_network.png")
fig.savefig(out_png, dpi=300, bbox_inches=None, facecolor='white', format='png')
plt.close(fig)

# ─── 12. Save SVG ─────────────────────────────────────────────────────────────
fig2, ax2 = plt.subplots(figsize=(12, 10))
fig2.subplots_adjust(left=0.02, right=0.98, top=0.96, bottom=0.02)
draw_network(ax2)
out_svg = os.path.join(OUT_DIR, "fig8_haplotype_network.svg")
fig2.savefig(out_svg, format='svg', bbox_inches=None, facecolor='white')
plt.close(fig2)

# ─── 13. Verify ───────────────────────────────────────────────────────────────
for path in [out_png, out_svg]:
    size   = os.path.getsize(path)
    status = "OK" if size > 10_000 else "WARNING: file suspiciously small"
    print(f"Saved ({size/1024:.1f} KB) [{status}]: {path}")

print("DONE")
