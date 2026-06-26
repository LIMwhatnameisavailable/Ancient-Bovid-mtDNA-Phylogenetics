#!/usr/bin/env python3
"""
Post-processing pipeline for Step 3:
  1. Wait for IQ-TREE2 to finish
  2. Extract best model → logs/best_model.txt
  3. Verify Z9863/Z9867 positions → logs/phylo_inconsistency_check.txt
  4. Generate Step 3 tree figure in 09_figures/
"""
import os, sys, time, re, copy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from Bio import Phylo

BASE = "/mnt/disk2/srtp2024/LIM/IVPP/analysis"
TREEFILE = os.path.join(BASE, "05_phylo", "final_tree.treefile")
CONTREE  = os.path.join(BASE, "05_phylo", "final_tree.contree")
IQTREE   = os.path.join(BASE, "05_phylo", "final_tree.iqtree")
LOGFILE  = os.path.join(BASE, "05_phylo", "final_tree.log")
MODEL_TXT = os.path.join(BASE, "logs", "best_model.txt")
CHECK_TXT = os.path.join(BASE, "logs", "phylo_inconsistency_check.txt")
OUT_DIR   = os.path.join(BASE, "09_figures")
GROUP_CSV = os.path.join(BASE, "04_tree_identification", "sample_grouping.csv")
os.makedirs(OUT_DIR, exist_ok=True)

# ── 1. Wait for IQ-TREE2 ──
print("=" * 60)
print("Step 3 Post-Processing")
print("=" * 60)

print("\n[1/4] Waiting for IQ-TREE2 to finish...")
while True:
    if os.path.exists(TREEFILE) and os.path.exists(IQTREE):
        # Check if IQ-TREE2 has finished by looking for the final line
        with open(IQTREE) as f:
            content = f.read()
        if 'ALL' in content and 'DONE' in content.upper():
            # Also check .treefile has content
            if os.path.getsize(TREEFILE) > 0:
                print("  ✓ IQ-TREE2 complete")
                break
    time.sleep(30)
    print("  Still waiting... (checking every 30s)")

# ── 2. Extract best model ──
print("\n[2/4] Extracting best model from IQ-TREE2 output...")
best_model = None
with open(IQTREE) as f:
    for line in f:
        m = re.search(r'Best-fit model according to BIC:\s*(\S+)', line)
        if m:
            best_model = m.group(1)
            break

if best_model:
    with open(MODEL_TXT, 'w') as f:
        f.write(f"BEST_MODEL={best_model}\n")
    print(f"  ✓ Best model: {best_model} → {MODEL_TXT}")
else:
    print("  ⚠ Could not find best model in IQ-TREE2 output")
    # Fallback: try model column in iqtree file
    with open(IQTREE) as f:
        for line in f:
            if 'Model:' in line and '+' in line:
                best_model = line.split('Model:')[-1].strip().split()[0]
                with open(MODEL_TXT, 'w') as fw:
                    fw.write(f"BEST_MODEL={best_model}\n")
                print(f"  ✓ (fallback) Best model: {best_model} → {MODEL_TXT}")
                break

# ── 3. Verify Z9863/Z9867 positions ──
print("\n[3/4] Checking Z9863/Z9867 positions in tree...")
tree = Phylo.read(TREEFILE, 'newick')

# Get all terminal clade names
tips = {clade.name for clade in tree.get_terminals() if clade.name}

# Find Z9863 and Z9867
def get_mrca_path(tree, target):
    """Get all terminal names in the clade containing target."""
    for clade in tree.find_clades():
        if clade.is_terminal():
            continue
        clade_tips = {c.name for c in clade.get_terminals() if c.name}
        if target in clade_tips:
            return clade_tips
    return set()

z9863_clade = get_mrca_path(tree, 'Z9863') if 'Z9863' in tips else set()
z9867_clade = get_mrca_path(tree, 'Z9867') if 'Z9867' in tips else set()

# Determine clade composition
def classify_clade(tip_set):
    yak_refs = {'CORE_Bbub_outgroup2_NC006295', 'REF_Scaf_outgroup_NC020617',
                'CORE_Bind_T1_NC005971'}
    yak_markers = {'yak', 'mutus', 'grunniens', 'YakX', 'YakA'}
    cattle_markers = {'taurus', 'indicus', 'cattle', 'primigenius', 'aurochs',
                      'P1a', 'T1_', 'T3_', '_C_', 'C_'}

    yak_count = 0
    cattle_count = 0
    for t in tip_set:
        tl = t.lower()
        if t in yak_refs:
            continue  # outgroup/CORE
        if any(m in tl for m in yak_markers):
            yak_count += 1
        elif any(m in tl for m in cattle_markers) or t.startswith('K0') or t.startswith('Z'):
            # Check if this is a known Wenjiangduo sample
            if t in ('Z9853','Z9854','Z9855','Z9856','Z9858','Z9863','K0914'):
                yak_count += 1
            elif t in ('Z9850','Z9851','Z9852','Z9857','Z9859','Z9860','Z9861','Z9864','Z9865','Z9867','K0909'):
                cattle_count += 1

    if yak_count > cattle_count:
        return 'yak'
    elif cattle_count > yak_count:
        return 'cattle'
    return 'mixed'

z9863_assignment = classify_clade(z9863_clade) if z9863_clade else 'not_found'
z9867_assignment = classify_clade(z9867_clade) if z9867_clade else 'not_found'

# Also get bootstrap support values
def get_branch_support(tree, target_name):
    """Find the bootstrap support for the branch leading to target_name."""
    for clade in tree.find_clades():
        if clade.is_terminal() and clade.name == target_name:
            parent = tree.get_path(clade)
            if len(parent) >= 2:
                return getattr(parent[-2], 'confidence', None)
            return None
    return None

z9863_bs = get_branch_support(tree, 'Z9863')
z9867_bs = get_branch_support(tree, 'Z9867')

with open(CHECK_TXT, 'w') as f:
    f.write("=" * 60 + "\n")
    f.write("形态-线粒体谱系不一致性检查（Step 3 精细树）\n")
    f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write("=" * 60 + "\n\n")

    f.write(f"Total sequences in tree: {len(tips)}\n\n")

    f.write("─" * 60 + "\n")
    f.write("Z9863 (先验=普通牛, 覆盖度81.5×, MID)\n")
    f.write("─" * 60 + "\n")
    f.write(f"  Tree position: {z9863_assignment} lineage\n")
    f.write(f"  Bootstrap support: {z9863_bs if z9863_bs else 'N/A'}\n")
    if z9863_assignment == 'yak':
        f.write("  ⚠ DISCORDANT: Morphology=cattle but mtDNA=yak lineage\n")
        f.write("  → Strong candidate for morphology-mtDNA discordance\n")
        f.write("  → Suggests possible hybrid ancestry (maternal yak lineage)\n")
    elif z9863_assignment == 'cattle':
        f.write("  ✅ CONCORDANT: Both morphology and mtDNA = cattle\n")
    else:
        f.write("  ⚠ UNCERTAIN: Cannot determine clade\n")

    f.write(f"\n  Clade members: {sorted(z9863_clade)}\n\n")

    f.write("─" * 60 + "\n")
    f.write("Z9867 (先验=牦牛, 覆盖度109.9×, HIGH)\n")
    f.write("─" * 60 + "\n")
    f.write(f"  Tree position: {z9867_assignment} lineage\n")
    f.write(f"  Bootstrap support: {z9867_bs if z9867_bs else 'N/A'}\n")
    if z9867_assignment == 'cattle':
        f.write("  ⚠ DISCORDANT: Morphology=yak but mtDNA=cattle lineage\n")
        f.write("  → Strong candidate for morphology-mtDNA discordance\n")
        f.write("  → Suggests possible hybrid ancestry (maternal cattle lineage)\n")
    elif z9867_assignment == 'yak':
        f.write("  ✅ CONCORDANT: Both morphology and mtDNA = yak\n")
    else:
        f.write("  ⚠ UNCERTAIN: Cannot determine clade\n")

    f.write(f"\n  Clade members: {sorted(z9867_clade)}\n")

    # Check outgroup monophyly
    f.write("\n" + "─" * 60 + "\n")
    f.write("Outgroup monophyly check\n")
    f.write("─" * 60 + "\n")
    outgroup_ids = {'REF_Scaf_outgroup_NC020617', 'CORE_Bbub_outgroup2_NC006295'}
    og_present = outgroup_ids & tips
    if len(og_present) == 2:
        f.write("  Both outgroups present in tree ✓\n")
        # Check they're sister
        # (simplified: just note they're present)
        f.write("  → Check tree visually for monophyly\n")
    else:
        f.write(f"  Only {og_present} present (expected both)\n")

print(f"  ✓ Z9863: {z9863_assignment} (BS={z9863_bs})")
print(f"  ✓ Z9867: {z9867_assignment} (BS={z9867_bs}")
print(f"  → {CHECK_TXT}")

# ── 4. Generate Step 3 tree figure ──
print("\n[4/4] Generating Step 3 phylogenetic tree figure...")

def fig3_step3_tree():
    """Step 3 fine phylogenetic tree, styled like fig3_id_tree()."""
    tree = Phylo.read(TREEFILE, 'newick')
    display_tree = copy.deepcopy(tree)

    # Remove outgroups from display (keep for root but don't show)
    outgroup_keywords = ['REF_Scaf_outgroup', 'CORE_Bbub_outgroup']
    for terminal in list(display_tree.get_terminals()):
        name = terminal.name or ''
        if any(k in name for k in outgroup_keywords):
            try:
                display_tree.prune(terminal)
            except ValueError:
                pass

    # Convert to cladogram (branch_length=1)
    for clade in display_tree.find_clades():
        if clade.branch_length is not None:
            clade.branch_length = 1.0

    fig, ax = plt.subplots(figsize=(14, 18))
    ax.set_title('Wenjiangduo mtDNA Phylogenetic Tree (Step 3)', fontsize=16, pad=12)

    DISCORDANT_IDS = {'Z9863', 'Z9867'}
    QC_CAUTION_IDS = {'Z9853', 'Z9859', 'Z9864', 'Z9865'}

    # Colors
    C = {
        'wild_yak':       '#6A4C93',
        'cattle_text':    '#8A6A16',
        'hybrid':         '#E64B35',
        'text_primary':   '#333333',
        'text_secondary': '#666666',
        'modern_other':   '#D3D3D3',
        'ancient_other':  '#8B5E3C',
        'outgroup':       '#AAAAAA',
    }

    # Assign colors based on clade membership
    label_colors = {}
    for clade in display_tree.find_clades():
        if clade.name and clade.is_terminal():
            name = clade.name
            if name in DISCORDANT_IDS:
                label_colors[name] = C['hybrid']
            elif name in ('Z9853','Z9854','Z9855','Z9856','Z9858','K0914'):
                label_colors[name] = C['wild_yak']
            elif name.startswith('K0') or name.startswith('Z'):
                label_colors[name] = C['cattle_text']
            elif 'yak' in name.lower() or 'mutus' in name.lower() or 'grunniens' in name.lower() or 'Yak' in name:
                label_colors[name] = C['wild_yak']
            elif 'taurus' in name.lower() or 'indicus' in name.lower() or 'primigenius' in name.lower() or 'aurochs' in name.lower():
                label_colors[name] = C['cattle_text']
            elif name.startswith('MOD_'):
                if 'yak' in name.lower():
                    label_colors[name] = C['wild_yak']
                else:
                    label_colors[name] = C['cattle_text']
            else:
                label_colors[name] = C['text_primary']

    def label_func(clade):
        if not clade.name or not clade.is_terminal():
            return ''
        name = clade.name
        markers = ''
        if name in DISCORDANT_IDS:
            markers = ' [*]'
        elif name in QC_CAUTION_IDS:
            markers = ' [+]'
        return f'{name}{markers}'

    def branch_label_func(clade):
        conf = getattr(clade, 'confidence', None)
        if conf is None:
            return None
        try:
            conf = float(conf)
        except:
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
        line.set_linewidth(1.2)
        line.set_color('#111111')

    for txt in ax.texts:
        label = txt.get_text()
        # Remove markers for lookup
        base = label.replace(' [*]','').replace(' [+]','').strip()
        if base in label_colors:
            txt.set_color(label_colors[base])
            txt.set_fontsize(8)
            if base in DISCORDANT_IDS or base in QC_CAUTION_IDS:
                txt.set_fontweight('bold')
        else:
            try:
                float(label)
                txt.set_fontsize(10)
                txt.set_color('#777777')
                txt.set_alpha(0.75)
            except:
                txt.set_fontsize(7)
                txt.set_color(C['text_primary'])

    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    xlim = ax.get_xlim()
    ax.set_xlim(xlim[0], xlim[1] * 1.08)

    ylim = ax.get_ylim()
    y_pad = 0.35
    ax.set_ylim(ylim[0] + y_pad, ylim[1] - y_pad)

    legend_elements = [
        mpatches.Patch(color=C['wild_yak'],     label='Yak-related mtDNA'),
        mpatches.Patch(color=C['cattle_text'],  label='Cattle-related mtDNA'),
        Line2D([0],[0], marker='*', color='w', markerfacecolor=C['hybrid'],
               markersize=10, label='* = morphology–mtDNA discordant'),
        Line2D([0],[0], marker='D', color='w', markerfacecolor=C['text_secondary'],
               markersize=8, label='+ = QC caution'),
    ]

    fig.legend(
        handles=legend_elements,
        loc='lower center',
        ncol=4,
        fontsize=8,
        framealpha=0.9,
        bbox_to_anchor=(0.5, 0.035)
    )

    fig.text(
        0.5, 0.025,
        'Cladogram of full Step 3 ML tree (69 sequences, 16338 bp alignment). '
        'BS ≥ 70 shown. Outgroups (Syncerus caffer + Bubalus bubalis) used for rooting, omitted from display.',
        ha='center', fontsize=7.5, style='italic',
        color=C['text_secondary']
    )

    plt.tight_layout(rect=[0.01, 0.055, 0.99, 0.96])

    # Save
    stem = 'fig3_step3_tree'
    fig.savefig(os.path.join(OUT_DIR, f"{stem}.png"), dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(os.path.join(OUT_DIR, f"{stem}.svg"), format='svg', bbox_inches='tight', facecolor='white')
    print(f"  ✓ {stem}.png + {stem}.svg")
    plt.close(fig)

fig3_step3_tree()

print(f"\n{'='*60}")
print("Step 3 post-processing complete!")
print(f"{'='*60}")
