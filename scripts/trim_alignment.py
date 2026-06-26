#!/usr/bin/env python3
"""
Trim MAFFT alignment for identification tree:
  1. Remove columns where ≥50% of sequences have gap/N
  2. Remove leading/trailing columns where any sequence has gap (ragged ends)
  3. Report stats before/after
"""
import sys
from collections import Counter

inpath  = "/mnt/disk2/srtp2024/LIM/IVPP/analysis/04_tree_identification/identification_aln.fasta"
outpath = "/mnt/disk2/srtp2024/LIM/IVPP/analysis/04_tree_identification/identification_aln_trimmed.fasta"

# Read alignment
seqs = {}
cur_id = None
cur_seq = []
with open(inpath) as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if cur_id:
                seqs[cur_id] = "".join(cur_seq)
            cur_id = line.lstrip(">").split()[0]
            cur_seq = []
        else:
            cur_seq.append(line.upper())
if cur_id:
    seqs[cur_id] = "".join(cur_seq)

names = list(seqs.keys())
nseq = len(names)
alen = len(seqs[names[0]])
print(f"Before: {nseq} sequences, {alen} bp")

# Check all same length
for n in names:
    if len(seqs[n]) != alen:
        print(f"  WARNING: {n} length={len(seqs[n])}, expected {alen}")

# ── Step 1: Compute per-column gap/ambiguity fraction ──
gap_chars = {"-", "N", "?", ".", "X", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"}
col_gap = []
for i in range(alen):
    col = [seqs[n][i] for n in names]
    gap_count = sum(1 for b in col if b in gap_chars)
    col_gap.append(gap_count / nseq)

# ── Step 2: Trim ragged ends (5' and 3' ends with any gap) ──
start = 0
while start < alen and col_gap[start] > 0:
    start += 1
end = alen
while end > start and col_gap[end-1] > 0:
    end -= 1
print(f"Ragged ends removed: 5'{'='+str(start) if start else '–'} , 3'{'='+str(alen-end) if end<alen else '–'} ")

# ── Step 3: Remove internal columns with ≥50% gaps ──
keep_cols = []
removed = 0
for i in range(start, end):
    if col_gap[i] >= 0.5:
        removed += 1
    else:
        keep_cols.append(i)

trimmed_len = len(keep_cols)
print(f"Internal columns removed (≥50% gap/N): {removed}")
print(f"After trimming: {trimmed_len} bp")

# Write trimmed alignment
with open(outpath, "w") as f:
    for n in names:
        trimmed = "".join(seqs[n][i] for i in keep_cols)
        f.write(f">{n}\n")
        # Wrap at 80 chars
        for i in range(0, len(trimmed), 80):
            f.write(trimmed[i:i+80] + "\n")

print(f"Output: {outpath}")

# ── Extra: per-sample N-content after trimming ──
print(f"\nPer-sample N/gap content after trimming:")
for n in names:
    trimmed = "".join(seqs[n][i] for i in keep_cols)
    ngap = sum(1 for b in trimmed if b in gap_chars)
    pct = ngap / len(trimmed) * 100
    if pct > 1:
        print(f"  {n}: {pct:.1f}% gaps/N ({ngap}/{len(trimmed)})")
    else:
        print(f"  {n}: {pct:.2f}% gaps/N")