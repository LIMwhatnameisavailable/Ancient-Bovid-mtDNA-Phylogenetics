#!/usr/bin/env python3
"""
Step 4: Haplotype identification for Wenjiangduo samples.
- Excludes LOW coverage samples (Z9853/Z9864) from statistics
- Removes gap/N/ambiguity columns before haplotype comparison
- Output: 06_haplo/haplotype_assignment.csv
"""
import os, sys, json
from collections import Counter, OrderedDict
from datetime import datetime

BASE = "/mnt/disk2/srtp2024/LIM/IVPP/analysis"
ALIGN_FILE = os.path.join(BASE, "04_align", "final_aln.fasta")
OUT_DIR = os.path.join(BASE, "06_haplo")
LOG_DIR = os.path.join(BASE, "logs")
os.makedirs(OUT_DIR, exist_ok=True)

LOW_COV = {"Z9853", "Z9864"}

# ── Read alignment ──
print("=" * 60)
print("Step 4: Haplotype Analysis")
print("=" * 60)
print(f"\n[1/4] Reading alignment: {os.path.basename(ALIGN_FILE)}")

seqs = OrderedDict()
with open(ALIGN_FILE) as f:
    cur_id, cur_seq = None, []
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if cur_id:
                seqs[cur_id] = "".join(cur_seq)
            cur_id = line[1:].split()[0]
            cur_seq = []
        else:
            cur_seq.append(line.upper())
    if cur_id:
        seqs[cur_id] = "".join(cur_seq)

# Retain only Wenjiangduo samples (short numeric IDs: K0909, Z9850, etc.)
wjd_ids = [sid for sid in seqs if not any(
    prefix in sid for prefix in ["MOD_", "LIT", "CORE_", "REF_", "all_combined"]
)]
# Keep original order
wjd_seqs = OrderedDict((sid, seqs[sid]) for sid in wjd_ids)

print(f"  Total WJD samples: {len(wjd_seqs)}")

# ── Identify and report excluded LOW coverage samples ──
excluded = [sid for sid in wjd_seqs if sid in LOW_COV]
included = [sid for sid in wjd_seqs if sid not in LOW_COV]
print(f"  Excluded (LOW cov — no stats): {excluded}")
print(f"  Included in statistics: {len(included)}")

# ── Remove gap/ambig columns ──
print(f"\n[2/4] Filtering sites...")
aln_len = len(list(wjd_seqs.values())[0])
print(f"  Raw alignment length: {aln_len} bp")

# First, determine which columns to keep
AMBIG_CODES = set("NRYWSMKBDHV")
keep_cols = []
for col in range(aln_len):
    bases = set(seq[col] for seq in wjd_seqs.values())
    # Exclude if any sequence has gap OR N OR ambiguity code at this position
    if '-' in bases or '?' in bases or '.' in bases:
        continue
    if 'N' in bases:
        continue
    if bases & AMBIG_CODES:
        continue
    keep_cols.append(col)

print(f"  Columns kept (no gap/N/ambig): {len(keep_cols)} / {aln_len}")

# Build cleaned sequences
cleaned_seqs = {sid: "".join(seq[col] for col in keep_cols) for sid, seq in wjd_seqs.items()}

# ── Identify haplotypes (on full wjd set, inc. LOW) ──
print(f"\n[3/4] Identifying haplotypes (no-gap cleaned sequence)...")
hap_to_ids = OrderedDict()
for sid in wjd_ids:
    hkey = cleaned_seqs[sid]
    if hkey not in hap_to_ids:
        hap_to_ids[hkey] = []
    hap_to_ids[hkey].append(sid)

# Sort haplotypes by size desc, then by min sample ID
hap_sorted = sorted(hap_to_ids.items(), key=lambda x: (-len(x[1]), min(x[1])))

haplotype_map = {}  # sid -> haplotype_id
hap_details = []   # list of dicts for CSV
for i, (hkey, members) in enumerate(hap_sorted, 1):
    hap_id = f"H{i:02d}"
    for sid in members:
        haplotype_map[sid] = hap_id
    hap_details.append({
        "haplotype": hap_id,
        "size": len(members),
        "members": "|".join(sorted(members)),
        "included_in_stats": "|".join(sorted(m for m in members if m not in LOW_COV))
    })

# ── Statistics (excluding LOW coverage) ──
print(f"\n[4/4] Computing statistics (excl. {LOW_COV})...")

# Haplotype distribution among included samples
included_hap_counts = Counter()
for sid in included:
    included_hap_counts[haplotype_map[sid]] += 1

n_included = len(included)
n_haps = len(set(haplotype_map[sid] for sid in included))
hd = (n_included / (n_included - 1)) * (1 - sum((c / n_included)**2 for c in included_hap_counts.values()))

# Site-level stats for π computation (biallelic only)
variable_sites = 0
transitions = 0
transversions = 0
pi_sum = 0
n_pairs = n_included * (n_included - 1) / 2

# Compute pairwise differences for π
# Sample up to 2000 sites if alignment is very large
n_sites = len(keep_cols)
site_range = range(n_sites)

for site_idx in site_range:
    bases_at_site = [cleaned_seqs[sid][site_idx] for sid in included]
    base_set = set(bases_at_site)
    if len(base_set) <= 1:
        continue
    variable_sites += 1

    # Count allele frequencies
    base_counts = Counter(bases_at_site)
    # For π: fraction of pairwise differences at this site
    # = 1 - sum( (n_i/n)^2 ) * n/(n-1)
    site_het = 1 - sum((c / n_included)**2 for c in base_counts.values())
    # Adjust for finite sample: site_heterozygosity * n/(n-1)
    if n_included > 1:
        site_het_adj = site_het * n_included / (n_included - 1)
    else:
        site_het_adj = 0
    pi_sum += site_het_adj

    # Transition/transversion classification (simplified: only for biallelic)
    if len(base_set) == 2:
        b1, b2 = sorted(base_set)
        ts_pairs = {('A','G'), ('G','A'), ('C','T'), ('T','C')}
        if (b1, b2) in ts_pairs:
            transitions += 1
        else:
            transversions += 1

pi = pi_sum / n_sites if (n_included > 1 and n_sites > 0) else 0.0

print(f"\n{'='*60}")
print("RESULTS SUMMARY")
print(f"{'='*60}")
print(f"  Samples included (excl LOW): {n_included}")
print(f"  Cleaned alignment length:    {n_sites} bp")
print(f"  Variable sites (biallelic):  {variable_sites}")
print(f"  Total haplotypes:            {n_haps}")
print(f"  Haplotype diversity (Hd):    {hd:.4f}")
print(f"  Nucleotide diversity (π):    {pi:.6f}")

# Count unique/shared by site
tower_ids = {"K0909", "K0914"}
main_ids = set(sid for sid in included if sid not in tower_ids)
tower_haps = set(haplotype_map[sid] for sid in tower_ids if sid in included)
main_haps = set(haplotype_map[sid] for sid in main_ids if sid in included)
shared_haps = tower_haps & main_haps
tower_only = tower_haps - main_haps
main_only = main_haps - tower_haps

print(f"\n  ── By site (excl. LOW cov) ──")
print(f"  西南塔遗址独有:    {len(tower_only)} 个")
print(f"  主遗址独有:        {len(main_only)} 个")
print(f"  共享:              {len(shared_haps)} 个")

# ── Write CSV ──
csv_path = os.path.join(OUT_DIR, "haplotype_assignment.csv")
with open(csv_path, "w") as f:
    f.write("sample_id,haplotype,site,excluded_from_stats\n")
    for sid in wjd_ids:
        site = "Ancient_WJD_Tower" if sid in tower_ids else "Ancient_WJD_Main"
        exc = "YES" if sid in LOW_COV else ""
        f.write(f"{sid},{haplotype_map[sid]},{site},{exc}\n")

print(f"\nOutput: {csv_path}")

# ── Write detailed haplotype info ──
detail_path = os.path.join(OUT_DIR, "haplotype_details.txt")
with open(detail_path, "w") as f:
    f.write(f"Step 4 Haplotype Details\n")
    f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
    f.write(f"Cleaned alignment: {n_sites} bp (removed gap/N/ambig columns)\n")
    f.write(f"{'='*60}\n\n")
    for h in hap_details:
        f.write(f"{h['haplotype']}: {h['size']} sequences\n")
        f.write(f"  Members: {h['members']}\n")
        if h['included_in_stats']:
            f.write(f"  In stats: {h['included_in_stats']}\n")
        f.write("\n")

print(f"Output: {detail_path}")
print(f"\n{'='*60}")
print("DONE")
print(f"{'='*60}")