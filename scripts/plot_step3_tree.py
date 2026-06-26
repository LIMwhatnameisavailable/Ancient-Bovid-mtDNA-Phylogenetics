#!/usr/bin/env python3
"""Generate Step 3 ML tree figure, styled like fig3_id_tree()."""
import os, copy, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from Bio import Phylo

BASE = "/mnt/disk2/srtp2024/LIM/IVPP/analysis"
TREEFILE = os.path.join(BASE, "05_phylo", "final_tree.contree")
OUT_DIR  = os.path.join(BASE, "09_figures")
os.makedirs(OUT_DIR, exist_ok=True)

plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['Times New Roman','Georgia','DejaVu Serif'],
    'mathtext.fontset': 'stix', 'font.size': 10, 'axes.titlesize': 14, 'axes.labelsize': 12,
    'xtick.labelsize': 10, 'ytick.labelsize': 10, 'legend.fontsize': 10,
})

# Colors
C = {
    'wild_yak':  '#6A4C93', 'cattle_text': '#8A6A16', 'hybrid': '#E64B35',
    'text_p': '#333333', 'text_s': '#666666', 'low_cov': '#AAAAAA',
}

DISCORDANT = {'Z9863', 'Z9867'}
QCCAUTION = {'Z9853', 'Z9859', 'Z9864', 'Z9865'}
YAK_WJD = {'K0914','Z9853','Z9854','Z9855','Z9856','Z9858','Z9863'}
CATTLE_WJD = {'K0909','Z9850','Z9851','Z9852','Z9857','Z9859','Z9860','Z9861','Z9864','Z9865','Z9867'}

print("Reading tree...")
tree = Phylo.read(TREEFILE, 'newick')
display = copy.deepcopy(tree)

print(f"Tree has {len(list(display.get_terminals()))} taxa (including outgroups)")

# Cladogram
for clade in display.find_clades():
    if clade.branch_length is not None:
        clade.branch_length = 1.0

fig, ax = plt.subplots(figsize=(14, 14))
ax.set_title('Wenjiangduo mtDNA Fine Phylogenetic Tree (Step 3)', fontsize=16, pad=12)

# Label colors
label_colors = {}
for clade in display.find_clades():
    if clade.name and clade.is_terminal():
        name = clade.name
        if name in DISCORDANT:
            label_colors[name] = C['hybrid']
        elif name in YAK_WJD or 'yak' in name.lower() or 'mutus' in name.lower() or 'grunniens' in name.lower():
            label_colors[name] = C['wild_yak']
        else:
            label_colors[name] = C['cattle_text']

def label_func(c):
    if not c.name or not c.is_terminal():
        return ''
    name = c.name
    markers = ''
    if name in DISCORDANT: markers += ' [*]'
    if name in QCCAUTION: markers += ' [+]'
    return f'{name}{markers}'

def branch_label_func(c):
    conf = getattr(c, 'confidence', None)
    if conf is None: return None
    try:
        v = float(conf)
        return f'{v:.0f}' if v >= 70 else None
    except: return None

Phylo.draw(display, axes=ax, label_func=label_func, do_show=False,
           show_confidence=False, branch_labels=branch_label_func, label_colors=label_colors)

for line in ax.get_lines():
    line.set_linewidth(1.2); line.set_color('#111111')

for txt in ax.texts:
    label = txt.get_text()
    base = label.replace(' [*]','').replace(' [+]','').strip()
    if base in label_colors:
        txt.set_color(label_colors[base])
        txt.set_fontsize(7.5)
        if base in DISCORDANT or base in QCCAUTION:
            txt.set_fontweight('bold')
    else:
        try:
            float(label)
            txt.set_fontsize(9); txt.set_color('#777777'); txt.set_alpha(0.75)
        except:
            txt.set_fontsize(6.5); txt.set_color(C['text_p'])

ax.set_xlabel(''); ax.set_ylabel(''); ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values(): s.set_visible(False)

fig.canvas.draw()
renderer = fig.canvas.get_renderer()
x_max_data = ax.get_xlim()[1]
for txt in ax.texts:
    try:
        bb = txt.get_window_extent(renderer=renderer)
        inv = ax.transData.inverted()
        x_right = inv.transform((bb.x1, bb.y0))[0]
        if x_right > x_max_data:
            x_max_data = x_right
    except Exception:
        pass
ax.set_xlim(ax.get_xlim()[0], x_max_data * 0.95) 
ylim = ax.get_ylim(); ax.set_ylim(ylim[0]+0.35, ylim[1]-0.35)

leg = [
    mpatches.Patch(color=C['wild_yak'], label='Yak-related mtDNA'),
    mpatches.Patch(color=C['cattle_text'], label='Cattle-related mtDNA'),
    mpatches.Patch(facecolor='none', edgecolor='none',
           label='[*] = morphology–mtDNA discordant'),
    mpatches.Patch(facecolor='none', edgecolor='none',
           label='[+] = QC caution'),
]
fig.legend(handles=leg, loc='lower center', ncol=4, fontsize=8,
           framealpha=0.9, bbox_to_anchor=(0.5, 0.03))
fig.text(0.5, 0.02, 'ML tree (TIM2+F+I+G4, 1000 ultrafast bootstrap). BS ≥ 70 shown. '
         'Outgroups: Syncerus caffer (REF) + Bubalus bubalis (CORE).',
         ha='center', fontsize=7.5, style='italic', color=C['text_s'])

plt.tight_layout(rect=[0.01, 0.05, 0.99, 0.96])

stem = 'fig3_step3_tree'
fig.savefig(os.path.join(OUT_DIR, f"{stem}.png"), dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(OUT_DIR, f"{stem}.svg"), format='svg', bbox_inches='tight', facecolor='white')
print(f"✓ {stem}.png + {stem}.svg")
plt.close(fig)
print("Done!")