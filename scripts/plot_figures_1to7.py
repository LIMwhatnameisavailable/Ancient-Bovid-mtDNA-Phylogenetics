#!/mnt/disk2/srtp2024/miniconda3/envs/LIM_env/bin/python3
"""
温江多遗址牛族分析 — 综合绘图脚本 (v2.4)
脚本: scripts/plot_figures_1to7.py
日志: logs/figures_revision.log

Generates all figures from Steps 0-2c data following PLOT_GUIDE.md conventions.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
import os
import copy
import textwrap
from Bio import Phylo

# ── Paths ──
BASE = "/mnt/disk2/srtp2024/LIM/IVPP/analysis"
QC_CSV      = os.path.join(BASE, "01_qc", "qc_report.csv")
META_CSV    = os.path.join(BASE, "logs", "sample_metadata.csv")
SUMMARY_CSV = os.path.join(BASE, "03_blast", "species_summary.csv")
GROUP_CSV   = os.path.join(BASE, "04_tree_identification", "sample_grouping.csv")
TREEFILE    = os.path.join(BASE, "04_tree_identification", "id_tree.treefile")
OUT_DIR     = os.path.join(BASE, "09_figures")
LOG_DIR     = os.path.join(BASE, "logs")

os.makedirs(OUT_DIR, exist_ok=True)

# ── Load data ──
qc    = pd.read_csv(QC_CSV)
meta  = pd.read_csv(META_CSV)
summ  = pd.read_csv(SUMMARY_CSV)
group = pd.read_csv(GROUP_CSV)

# ── Global style ──
def apply_global_style():
    plt.rcParams.update({
        'font.family':        'serif',
        'font.serif':         ['Times New Roman','Georgia','DejaVu Serif'],
        'mathtext.fontset':   'stix',
        'axes.unicode_minus': False,
        'font.size':          10,
        'axes.titlesize':     14,
        'axes.labelsize':     12,
        'xtick.labelsize':    10,
        'ytick.labelsize':    10,
        'legend.fontsize':    10,
        'figure.titlesize':   14,
        'axes.linewidth':     1.0,
        'axes.edgecolor':     '#444444',
        'axes.facecolor':     'white',
        'axes.labelcolor':    '#333333',
        'axes.spines.top':    False,
        'axes.spines.right':  False,
        'xtick.direction':    'out',
        'ytick.direction':    'out',
        'xtick.major.size':   4.0,
        'ytick.major.size':   4.0,
        'xtick.color':        '#333333',
        'ytick.color':        '#333333',
        'lines.linewidth':    2.0,
        'legend.frameon':     True,
        'legend.framealpha':  0.9,
        'legend.edgecolor':   '#CCCCCC',
        'legend.fancybox':    False,
        'figure.facecolor':   'white',
        'figure.dpi':         100,
        'savefig.dpi':        300,
        'savefig.bbox':       'tight',
        'savefig.facecolor':  'white',
        'axes.grid':          False,
    })

apply_global_style()

# ── Color palette ──
C = {
    "ancient_wjd_main":  "#4A90D9",
    "ancient_wjd_tower": "#2A9D8F",
    "modern_tibet_yak":  "#1A3A5C",
    "wild_yak":          "#6A4C93",
    "modern_cattle":     "#E9C46A",
    "ancient_other":     "#8B5E3C",
    "modern_other":      "#D3D3D3",
    "hybrid_candidate":  "#E64B35",
    "low_coverage":      "#AAAAAA",
    "bar_ancient":       "#4A90D9",
    "bar_modern":        "#D3D3D3",
    "threshold":         "#E9C46A",
    "bsp_median":        "#4A90D9",
    "text_primary":      "#333333",
    "text_secondary":    "#666666",
    "axis_line":         "#444444",
    "grid_line":         "#E0E0E0",
    "success":           "#2A9D8F",
    "warning":           "#E9C46A",
    "error":             "#E64B35",
    "prior_cattle":      "#E9C46A",
    "prior_yak":         "#6A4C93",
    "bg_even":           "#F5F5F5",
    "morph_bovine":      "#E8E8E8",
    "morph_uncertain":   "#FFF8DC",
    "morph_cattle":      "#E8F5E9",
    "morph_yak":         "#F3E5F5",
    "discordant_red":    "#FFE8E8",
}

def save_figure(fig, stem):
    fig.savefig(os.path.join(OUT_DIR, f"{stem}.png"), dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(os.path.join(OUT_DIR, f"{stem}.svg"), format='svg', bbox_inches='tight', facecolor='white')
    print(f"  ✓ {stem}.png + {stem}.svg")

# ── Helper constants ──
DISCORDANT_IDS = {'Z9863', 'Z9867'}
QC_CAUTION_IDS = {'Z9853', 'Z9859', 'Z9864', 'Z9865'}
ALL_FOCUS = DISCORDANT_IDS | QC_CAUTION_IDS

def get_info(name):
    rows = group[group['sample_id'] == name]
    return rows.iloc[0] if len(rows) else None

def morph_category(morph_id, morph_conf):
    """Map morphology ID to category string."""
    if morph_id == '普通牛':
        return 'Cattle'
    elif morph_id == '牦牛':
        return 'Yak'
    elif morph_id == '牛?' or morph_conf == 'UNCERTAIN':
        return 'Uncertain'
    else:
        return 'Bovine'

def parse_percent_series(s):
    """Parse percent values stored as either numeric or strings with '%'."""
    return pd.to_numeric(s.astype(str).str.replace('%', '', regex=False), errors='coerce')

def get_percent_identity_column(df):
    """Return the most likely percent identity column name."""
    for col in ['percent_identity', 'Per. ident', 'per_ident', 'pident']:
        if col in df.columns:
            return col
    raise KeyError("No percent identity column found in summary table.")

def clean_tree_label(label):
    """Remove markers added to sample labels when matching text objects."""
    if label is None:
        return ''
    return (
        str(label)
        .replace(' [*]', '')
        .replace(' [+]', '')
        .replace(' [?]', '')
        .strip()
    )

def wrap_text(text, width=42):
    return '\n'.join(textwrap.wrap(str(text), width=width, break_long_words=False))


# ════════════════════════════════════════════
# FIG 1: QC覆盖度分层
# ════════════════════════════════════════════
def fig1_qc_coverage():
    """Bar chart with coverage tiers, 20× and 100× thresholds annotated."""
    df = meta.copy().sort_values('coverage', ascending=False)
    colors = {'HIGH': C['ancient_wjd_main'], 'MID': C['ancient_wjd_tower'], 'LOW': C['low_coverage']}
    bar_colors = [colors[r['coverage_tier']] for _, r in df.iterrows()]

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(range(len(df)), df['coverage'], color=bar_colors, edgecolor='white', linewidth=0.5)

    ax.axhline(100, color=C['threshold'], linewidth=1.2, linestyle='--', alpha=0.7)
    ax.axhline(20,  color=C['low_coverage'], linewidth=1.2, linestyle=':', alpha=0.7)

    ax.set_xticks(range(len(df)))
    tick_labels = []
    for sid in df['ID']:
        if sid in DISCORDANT_IDS:
            tick_labels.append(f'{sid}*')
        else:
            tick_labels.append(sid)
    ax.set_xticklabels(tick_labels, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Coverage depth (×)')
    ax.set_xlabel('Sample')
    ax.set_title('Wenjiangduo — mtDNA Coverage Depth')

    for i, (_, row) in enumerate(df.iterrows()):
        sid = row['ID']
        if sid in DISCORDANT_IDS:
            ax.text(i, row['coverage'] + 3, f'*{sid}', ha='center', fontsize=8,
                    color=C['hybrid_candidate'], fontweight='bold')
        elif row['coverage_tier'] == 'LOW':
            ax.text(i, row['coverage'] + 2, f"{row['coverage']}×", ha='center', fontsize=8,
                    color=C['error'], fontweight='bold')

    legend_elements = [
        mpatches.Patch(color=C['ancient_wjd_main'], label='HIGH (≥100×)'),
        mpatches.Patch(color=C['ancient_wjd_tower'], label='MID (20–99×)'),
        mpatches.Patch(color=C['low_coverage'], label='LOW (<20×)'),
        Line2D([0],[0], color=C['threshold'], linestyle='--', linewidth=1.2, label='100× threshold'),
        Line2D([0],[0], color=C['low_coverage'], linestyle=':', linewidth=1.2, label='20× threshold'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', framealpha=0.9)
    ax.set_ylim(0, max(df['coverage']) * 1.18)
    save_figure(fig, 'fig1_qc_coverage')
    plt.close(fig)


# ════════════════════════════════════════════
# FIG 2: 物种鉴定组成
# ════════════════════════════════════════════
def fig2_species_composition():
    """Pie chart + horizontal bar chart with BLAST identity per sample."""
    g = group[group['sample_id'].isin(meta['ID'])].copy()
    n_yak    = (g['mtDNA_clade'] == 'yak').sum()
    n_cattle = (g['mtDNA_clade'] == 'cattle').sum()
    total = n_yak + n_cattle

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    pie_colors = [C['wild_yak'], C['modern_cattle']]
    ax1.pie(
        [n_yak, n_cattle],
        colors=pie_colors,
        startangle=90,
        textprops={'fontsize': 10, 'color': C['text_primary']}
    )
    ax1.set_title('mtDNA Clade Assignment', fontsize=13)

    pie_legend = [
        mpatches.Patch(color=C['wild_yak'], label=f'Yak-related mtDNA ({n_yak}/{total})'),
        mpatches.Patch(color=C['modern_cattle'], label=f'Cattle-related mtDNA ({n_cattle}/{total})'),
    ]
    ax1.legend(
        handles=pie_legend,
        loc='lower center',
        fontsize=9,
        framealpha=0.9,
        bbox_to_anchor=(0.5, -0.06)
    )

    df = meta.merge(summ, left_on='ID', right_on='sample_id', how='left')
    pid_col = get_percent_identity_column(df)
    df['pident_val'] = parse_percent_series(df[pid_col])
    df = df.sort_values('pident_val')
    sample_ids = df['ID'].tolist()

    bar_colors = []
    edge_colors = []
    linewidths = []

    for _, row in df.iterrows():
        sid = row['ID']
        if sid in DISCORDANT_IDS:
            bar_colors.append(C['ancient_other'])
            edge_colors.append(C['hybrid_candidate'])
            linewidths.append(2.5)
        elif sid in QC_CAUTION_IDS:
            bar_colors.append(C['low_coverage'])
            edge_colors.append(C['low_coverage'])
            linewidths.append(1.0)
        else:
            clade = group.loc[group['sample_id'] == sid, 'mtDNA_clade'].values
            if len(clade) and clade[0] == 'yak':
                bar_colors.append(C['wild_yak'])
            else:
                bar_colors.append(C['modern_cattle'])
            edge_colors.append('none')
            linewidths.append(0.5)

    ax2.barh(
        range(len(df)),
        df['pident_val'],
        color=bar_colors,
        edgecolor=edge_colors,
        linewidth=linewidths
    )

    ax2.set_yticks(range(len(df)))
    ax2.set_yticklabels(sample_ids, fontsize=8)
    ax2.set_xlabel('BLAST Percent Identity (%)')
    ax2.set_title('BLAST Top Hit Identity per Sample')
    ax2.axvline(99, color=C['threshold'], linestyle='--', linewidth=1.0, alpha=0.6)

    xmin = max(94.0, np.nanmin(df['pident_val']) - 0.4)
    ax2.set_xlim(xmin, 100.55)

    legend_elements = [
        mpatches.Patch(color=C['wild_yak'], label='Yak-related mtDNA'),
        mpatches.Patch(color=C['modern_cattle'], label='Cattle-related mtDNA'),
        mpatches.Patch(facecolor=C['ancient_other'], edgecolor=C['hybrid_candidate'],
                       linewidth=1.5, label='Discordant candidate'),
        mpatches.Patch(color=C['low_coverage'], label='QC caution'),
    ]

    ax2.legend(
        handles=legend_elements,
        fontsize=8,
        loc='lower right',
        framealpha=0.9
    )

    plt.tight_layout()
    save_figure(fig, 'fig2_species_composition')
    plt.close(fig)


# ════════════════════════════════════════════
# FIG 3: ML鉴定树
# ════════════════════════════════════════════
def fig3_id_tree():
    """Phylogenetic identification tree shown as a compact cladogram topology."""
    tree = Phylo.read(TREEFILE, 'newick')
    display_tree = copy.deepcopy(tree)

    for clade in display_tree.find_clades():
        if clade.branch_length is not None:
            clade.branch_length = 1.0

    fig, ax = plt.subplots(figsize=(12, 9))
    ax.set_title('Wenjiangduo mtDNA Identification Tree', fontsize=16, pad=12)

    cattle_text_color = '#8A6A16'

    label_colors = {}
    for clade in display_tree.find_clades():
        if clade.name and clade.is_terminal():
            info = get_info(clade.name)
            if info is not None:
                if clade.name in DISCORDANT_IDS:
                    label_colors[clade.name] = C['hybrid_candidate']
                elif info['mtDNA_clade'] == 'yak':
                    label_colors[clade.name] = C['wild_yak']
                else:
                    label_colors[clade.name] = cattle_text_color
            else:
                if 'taurus' in clade.name or 'indicus' in clade.name:
                    label_colors[clade.name] = cattle_text_color
                elif 'mutus' in clade.name or 'grunniens' in clade.name:
                    label_colors[clade.name] = C['wild_yak']
                else:
                    label_colors[clade.name] = C['text_primary']

    def label_func(clade):
        if not clade.name or not clade.is_terminal():
            return ''
        name = clade.name
        info = get_info(name)
        if info is None:
            return name
        markers = ''
        if name in DISCORDANT_IDS:
            markers = ' [*]'
        elif name in QC_CAUTION_IDS:
            markers = ' [+]'
        elif info.get('morphology_confidence') == 'UNCERTAIN':
            markers = ' [?]'
        return f'{name}{markers}'

    def branch_label_func(clade):
        conf = getattr(clade, 'confidence', None)
        if conf is None:
            return None
        try:
            conf = float(conf)
        except Exception:
            return None
        if conf >= 70:
            return f'{conf:.0f}'
        return None

    Phylo.draw(
        display_tree,
        axes=ax,
        label_func=label_func,
        do_show=False,
        show_confidence=False,
        branch_labels=branch_label_func,
        label_colors=label_colors
    )

    for line in ax.get_lines():
        line.set_linewidth(1.6)
        line.set_color('#111111')

    for txt in ax.texts:
        label = txt.get_text()
        base = clean_tree_label(label)

        if base in label_colors:
            txt.set_color(label_colors[base])
            txt.set_fontsize(11)
            if base in DISCORDANT_IDS or base in QC_CAUTION_IDS:
                txt.set_fontweight('bold')
        else:
            try:
                float(label)
                txt.set_fontsize(12)
                txt.set_color('#777777')
                txt.set_alpha(0.75)
            except ValueError:
                txt.set_fontsize(9)
                txt.set_color(C['text_primary'])

    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_yticklabels([])

    for spine in ax.spines.values():
        spine.set_visible(False)

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.set_xlim(xlim[0], xlim[1] * 1.04)

    y_pad = 0.35
    ax.set_ylim(ylim[0] + y_pad, ylim[1] - y_pad)

    legend_elements = [
        mpatches.Patch(color=C['wild_yak'],      label='Yak-related mtDNA'),
        mpatches.Patch(color=C['modern_cattle'], label='Cattle-related mtDNA'),
        mpatches.Patch(facecolor='none', edgecolor='none',
               label='[*] = morphology–mtDNA discordant'),
        mpatches.Patch(facecolor='none', edgecolor='none',
               label='[?] = morphology uncertain'),
        mpatches.Patch(facecolor='none', edgecolor='none',
               label='[+] = QC caution'),
    ]

    fig.legend(
        handles=legend_elements,
        loc='lower center',
        ncol=5,
        fontsize=9,
        framealpha=0.9,
        bbox_to_anchor=(0.5, 0.092)
    )

    fig.text(
        0.5,
        0.073,
        'Cladogram topology shown for mtDNA clade assignment; branch lengths are not proportional.',
        ha='center',
        fontsize=8.5,
        style='italic',
        color=C['text_secondary']
    )

    plt.tight_layout(rect=[0.02, 0.10, 0.98, 0.96])
    save_figure(fig, 'fig3_id_tree')
    plt.close(fig)


# ════════════════════════════════════════════
# FIG 4: 覆盖度 vs N%
# ════════════════════════════════════════════
def fig4_coverage_vs_n():
    """Scatter plot — coverage vs N%, with fixed label positions and straight leader lines."""
    df = qc.merge(meta, on='ID', suffixes=('_qc', '_meta'))
    cov_col = 'coverage_meta' if 'coverage_meta' in df.columns else 'coverage'

    fig, ax = plt.subplots(figsize=(8, 6))

    colors = []
    for _, row in df.iterrows():
        g = group[group['sample_id'] == row['ID']]
        if len(g):
            clade = g.iloc[0]['mtDNA_clade']
            colors.append(C['wild_yak'] if clade == 'yak' else C['modern_cattle'])
        else:
            colors.append(C['modern_other'])

    ax.scatter(
        df[cov_col],
        df['N%'],
        c=colors,
        s=80,
        alpha=0.85,
        edgecolors='white',
        linewidth=0.5
    )

    ax.axvline(20,  color=C['low_coverage'], linestyle=':',  linewidth=1.0, alpha=0.55)
    ax.axvline(100, color=C['threshold'],     linestyle='--', linewidth=1.0, alpha=0.55)
    ax.axhline(5,   color=C['error'],         linestyle=':',  linewidth=1.0, alpha=0.45)

    label_positions = {
        'Z9853': (2.6, 27.5),
        'Z9864': (8.9, 6.4),
        'Z9859': (30, 10),
        'Z9865': (32, 5.6),
        'Z9863': (70, 3.4),
        'Z9867': (130, 3.8),
    }

    for _, row in df.iterrows():
        sid = row['ID']
        if sid in ALL_FOCUS and sid in label_positions:
            clr = C['hybrid_candidate'] if sid in DISCORDANT_IDS else C['text_primary']
            ax.annotate(
                sid,
                xy=(row[cov_col], row['N%']),
                xytext=label_positions[sid],
                textcoords='data',
                fontsize=8,
                fontweight='bold',
                color=clr,
                ha='center',
                va='center',
                arrowprops=dict(
                    arrowstyle='->',
                    color=clr,
                    lw=1.05,
                    alpha=0.95,
                    shrinkA=3,
                    shrinkB=4,
                    connectionstyle='arc3,rad=0'
                ),
                bbox=dict(
                    boxstyle='round,pad=0.18',
                    facecolor='white',
                    edgecolor='none',
                    alpha=0.78
                )
            )

    ax.set_xlabel('Coverage depth (×)')
    ax.set_ylabel('N content (%)')
    ax.set_title('Coverage vs N Content — Sample Quality Assessment')
    ax.set_xscale('log')

    y_max = max(df['N%'].max() * 1.12, 8)
    ax.set_ylim(-0.8, y_max)

    legend_elements = [
        mpatches.Patch(color=C['wild_yak'], label='Yak clade'),
        mpatches.Patch(color=C['modern_cattle'], label='Cattle clade'),
        Line2D([0],[0], color=C['threshold'],     linestyle='--', linewidth=1.0, label='100× coverage'),
        Line2D([0],[0], color=C['low_coverage'],  linestyle=':',  linewidth=1.0, label='20× coverage'),
        Line2D([0],[0], color=C['error'],         linestyle=':',  linewidth=1.0, label='5% N threshold'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8)

    save_figure(fig, 'fig4_quality_scatter')
    plt.close(fig)


# ════════════════════════════════════════════
# FIG 5: 子遗址物种对比
# ════════════════════════════════════════════
def fig5_site_comparison():
    """Side-by-side bar charts of mtDNA clade by subsite."""
    g = group[group['sample_id'].isin(meta['ID'])].copy()
    g = g.merge(meta[['ID', 'site']], left_on='sample_id', right_on='ID', how='left')

    tower = g[g['site'] == '西南塔遗址']
    main  = g[g['site'] == '主遗址']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5.5))

    for ax, df_site, site_name in [
        (ax1, tower, 'Xinan Pagoda (n=2)\n(descriptive comparison only)'),
        (ax2, main,  'Main Site (n=16)')]:

        clade_counts = df_site['mtDNA_clade'].value_counts()
        cats = ['yak', 'cattle']
        vals = [clade_counts.get(c, 0) for c in cats]
        clrs = [C['wild_yak'], C['modern_cattle']]

        bars = ax.bar(['Yak', 'Cattle'], vals, color=clrs, edgecolor='white', linewidth=1.0)
        for bar in bars:
            h = bar.get_height()
            if h > 0:
                ax.text(bar.get_x() + bar.get_width()/2., h + 0.08,
                        f'{int(h)}', ha='center', va='bottom', fontsize=11, fontweight='bold')

        ax.set_title(site_name, fontsize=10)
        ax.set_ylabel('Number of samples')
        ax.set_ylim(0, max(max(vals), 1) * 1.4)

    fig.suptitle('mtDNA Clade Composition by Subsite', fontsize=14, y=1.01)

    legend_elements = [
        mpatches.Patch(color=C['wild_yak'], label='Yak mtDNA'),
        mpatches.Patch(color=C['modern_cattle'], label='Cattle mtDNA'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=2, fontsize=9,
               bbox_to_anchor=(0.5, -0.10))

    fig.text(0.5, -0.02, 'Discordant candidates are listed in Fig7',
             ha='center', fontsize=8, style='italic', color=C['text_secondary'])

    plt.tight_layout()
    save_figure(fig, 'fig5_site_comparison')
    plt.close(fig)


# ════════════════════════════════════════════
# FIG 6: BLAST 匹配度分布
# ════════════════════════════════════════════
def fig6_blast_identity():
    """Histogram + per-sample dot plot of BLAST identity."""
    df = summ.copy()
    pid_col = get_percent_identity_column(df)
    df['pident'] = parse_percent_series(df[pid_col])

    taxon = df['inferred_taxon'].fillna('').str.lower()
    df_yak    = df[taxon.str.contains('yak')]
    df_cattle = df[taxon.str.contains('cattle|taurus')]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    bins = np.arange(94, 100.5, 0.5)
    ax1.hist([df_yak['pident'], df_cattle['pident']],
             bins=bins, color=[C['wild_yak'], C['modern_cattle']],
             label=['Yak', 'Cattle'], alpha=0.8, edgecolor='white')
    ax1.axvline(99, color=C['threshold'], linestyle='--', linewidth=1.2, alpha=0.7, label='99% threshold')
    ax1.set_xlabel('BLAST Percent Identity (%)')
    ax1.set_ylabel('Number of samples')
    ax1.set_title('Distribution of BLAST Top Hit Identity')
    ax1.legend(fontsize=9)

    df_sorted = df.sort_values('pident')
    y_pos = range(len(df_sorted))
    dot_colors = [C['wild_yak'] if 'yak' in str(r['inferred_taxon']).lower() else C['modern_cattle']
                  for _, r in df_sorted.iterrows()]
    ax2.scatter(df_sorted['pident'], y_pos, c=dot_colors, s=60, alpha=0.8, edgecolors='white', linewidth=0.5)
    ax2.axvline(99, color=C['threshold'], linestyle='--', linewidth=1.0, alpha=0.5)
    ax2.axvline(95, color=C['error'],     linestyle=':',  linewidth=1.0, alpha=0.5)

    for idx, (_, row) in enumerate(df_sorted.iterrows()):
        sid = row['sample_id']
        if sid in ALL_FOCUS:
            marker = ' [*]' if sid in DISCORDANT_IDS else ' [+]'
            clr = C['hybrid_candidate'] if sid in DISCORDANT_IDS else C['text_secondary']
            ax2.annotate(f'{sid}{marker}', (row['pident'], idx),
                        textcoords="offset points", xytext=(5, 0),
                        fontsize=7, color=clr, fontweight='bold')

    ax2.set_yticks(list(y_pos))
    ax2.set_yticklabels(df_sorted['sample_id'], fontsize=7)
    ax2.set_xlabel('BLAST Percent Identity (%)')
    ax2.set_title('Per-Sample BLAST Identity')
    ax2.set_xlim(93, 100.5)

    save_figure(fig, 'fig6_blast_identity')
    plt.close(fig)

# ════════════════════════════════════════════
# FIG 7: 形态鉴定 × mtDNA谱系 — 矩阵图
# ════════════════════════════════════════════
def fig7_morphology_vs_mtdna():
    """
    Matrix plot:
      Morphology × mtDNA count matrix.
    """
    g = group[group['sample_id'].isin(meta['ID'])].copy()

    g['morph_cat'] = g.apply(
        lambda r: morph_category(r['morphology_id'], r['morphology_confidence']),
        axis=1
    )
    g['mtdna_cat'] = g['mtDNA_clade'].map({
        'yak': 'Yak mtDNA',
        'cattle': 'Cattle mtDNA',
        'outgroup': 'Outgroup'
    })

    morph_order = ['Bovine', 'Cattle', 'Yak', 'Uncertain']
    mtdna_order = ['Yak mtDNA', 'Cattle mtDNA']

    cell_data = {}
    for mtdna in mtdna_order:
        for morph in morph_order:
            cell_data[(mtdna, morph)] = {
                'count': 0,
                'ids': [],
                'discordant_ids': set()
            }

    for _, row in g.iterrows():
        key = (row['mtdna_cat'], row['morph_cat'])
        if key not in cell_data:
            continue
        sid = row['sample_id']
        cell_data[key]['count'] += 1
        cell_data[key]['ids'].append(sid)
        if sid in DISCORDANT_IDS:
            cell_data[key]['discordant_ids'].add(sid)

    def cell_color(mtdna, morph, data):
        if data['count'] == 0:
            return '#FAFAFA'
        if len(data['discordant_ids']) > 0:
            return C['discordant_red']
        if morph == 'Uncertain':
            return C['morph_uncertain']
        if morph == 'Bovine':
            return '#E8EEF8'
        if morph == 'Yak' and mtdna == 'Yak mtDNA':
            return C['morph_cattle']
        if morph == 'Cattle' and mtdna == 'Cattle mtDNA':
            return C['morph_cattle']
        return '#F5F5F5'

    fig, ax = plt.subplots(figsize=(7.6, 5.2))

    ax.set_title(
        'Morphology–mtDNA Concordance Matrix',
        fontsize=14,
        pad=14
    )

    n_rows = len(mtdna_order)
    n_cols = len(morph_order)

    discordant_boxes = []

    for i, mtdna in enumerate(mtdna_order):
        y = n_rows - 1 - i
        for j, morph in enumerate(morph_order):
            data = cell_data[(mtdna, morph)]
            bg = cell_color(mtdna, morph, data)

            rect = Rectangle(
                (j, y),
                1.0,
                1.0,
                facecolor=bg,
                edgecolor='#CCCCCC',
                linewidth=1.0
            )
            ax.add_patch(rect)

            if data['discordant_ids']:
                discordant_boxes.append((j, y))

            if data['count'] == 0:
                ax.text(
                    j + 0.5,
                    y + 0.5,
                    '—',
                    ha='center',
                    va='center',
                    fontsize=13,
                    color='#BBBBBB'
                )
            else:
                ax.text(
                    j + 0.5,
                    y + 0.60,
                    f'n = {data["count"]}',
                    ha='center',
                    va='center',
                    fontsize=13,
                    fontweight='bold',
                    color=C['text_primary']
                )
                if data['discordant_ids']:
                    ids = ', '.join(sorted(data['discordant_ids']))
                    ax.text(
                        j + 0.5,
                        y + 0.35,
                        ids + '*',
                        ha='center',
                        va='center',
                        fontsize=9,
                        fontweight='bold',
                        color=C['hybrid_candidate']
                    )

    for j, y in discordant_boxes:
        rect = Rectangle(
            (j + 0.012, y + 0.012),
            0.976,
            0.976,
            facecolor='none',
            edgecolor=C['hybrid_candidate'],
            linewidth=2.2,
            joinstyle='miter',
            zorder=5
        )
        ax.add_patch(rect)

    ax.set_xlim(0, n_cols)
    ax.set_ylim(0, n_rows)
    ax.set_box_aspect(0.62)

    ax.set_xticks(np.arange(n_cols) + 0.5)
    ax.set_xticklabels(morph_order, fontsize=11)

    ax.set_yticks([1.5, 0.5])
    ax.set_yticklabels(['Yak mtDNA', 'Cattle mtDNA'], fontsize=11)

    ax.set_xlabel('Morphology Identification', fontsize=12, labelpad=10)
    ax.set_ylabel('mtDNA Clade', fontsize=12, labelpad=10)

    ax.tick_params(axis='both', length=0)

    for spine in ax.spines.values():
        spine.set_visible(False)

    legend_elements = [
        mpatches.Patch(facecolor=C['discordant_red'], edgecolor=C['hybrid_candidate'],
                       linewidth=2.0, label='Discordant candidate'),
        mpatches.Patch(facecolor='#E8EEF8', edgecolor='#CCCCCC',
                       label='Generic bovine morphology'),
        mpatches.Patch(facecolor=C['morph_uncertain'], edgecolor='#CCCCCC',
                       label='Uncertain morphology'),
        mpatches.Patch(facecolor='#FAFAFA', edgecolor='#CCCCCC',
                       label='No sample'),
    ]

    ax.legend(
        handles=legend_elements,
        loc='upper center',
        bbox_to_anchor=(0.5, -0.3),
        ncol=2,
        fontsize=8.5,
        framealpha=0.9,
        columnspacing=1.6,
        handlelength=1.8,
        borderpad=0.6,
        labelspacing=0.6
    )

    plt.tight_layout(rect=[0, 0.08, 1, 0.96])
    save_figure(fig, 'fig7_morphology_vs_mtdna')
    plt.close(fig)

# ════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════
def main():
    print("=" * 60)
    print("温江多遗址牛族分析 — 最终修订绘图 v2.4")
    print(f"输出目录: {OUT_DIR}")
    print("=" * 60)

    for lbl, fn in [
        ("Fig1: QC 覆盖度分层",       fig1_qc_coverage),
        ("Fig2: 物种鉴定组成",        fig2_species_composition),
        ("Fig3: ML 鉴定树",           fig3_id_tree),
        ("Fig4: 覆盖度 vs N% 散点图", fig4_coverage_vs_n),
        ("Fig5: 子遗址物种对比",      fig5_site_comparison),
        ("Fig6: BLAST 匹配度分布",    fig6_blast_identity),
        ("Fig7: 形态×mtDNA 矩阵图",   fig7_morphology_vs_mtdna),
    ]:
        print(f"\n[{lbl}]")
        fn()

    print("\n" + "=" * 60)
    print("最终修订全部完成！")

if __name__ == '__main__':
    main()
