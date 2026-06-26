#!/usr/bin/env python3
"""
step3_filter_alignment.py
同源区段裁剪与缺失位点过滤

用法: python3 step3_filter_alignment.py < input.fasta > output.fasta

步骤:
  1. 读取多序列比对（FASTA 格式）
  2. 裁剪两端：去除所有序列两端全为 gap 的"不完整延伸"列
  3. 列过滤：去除缺失比例 ≥ 50% 的列（gap 或 N 合计占 50% 以上）
  4. 输出过滤后的比对
"""

import sys
from collections import defaultdict

def parse_fasta(stream):
    seqs = {}
    current_id = None
    current_seq = []
    for line in stream:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id:
                seqs[current_id] = "".join(current_seq)
            current_id = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line.upper())
    if current_id:
        seqs[current_id] = "".join(current_seq)
    return seqs

def filter_alignment(seqs):
    ids = list(seqs.keys())
    seq_list = [seqs[i] for i in ids]
    n_seq = len(seq_list)
    if n_seq == 0:
        return {}, 0, 0, 0
    aln_len = len(seq_list[0])

    # Convert to columns
    columns = [[seq[i] for seq in seq_list] for i in range(aln_len)]

    # Step 1: Trim leading/trailing all-gap columns
    start = 0
    while start < aln_len and all(c == "-" for c in columns[start]):
        start += 1
    end = aln_len
    while end > start and all(c == "-" for c in columns[end - 1]):
        end -= 1

    # Step 2: Filter columns with ≥ 50% missing (gap or N)
    kept_cols = []
    dropped_positions = []
    for i in range(start, end):
        missing = sum(1 for c in columns[i] if c in ("-", "N"))
        if missing / n_seq < 0.5:
            kept_cols.append(i)
        else:
            dropped_positions.append(i)

    # Rebuild sequences
    filtered = {}
    for idx, sid in enumerate(ids):
        filtered[sid] = "".join(seq_list[idx][c] for c in kept_cols)

    trim_left = start
    trim_right = aln_len - end
    dropped = len(dropped_positions)
    return filtered, trim_left, trim_right, dropped

def main():
    seqs = parse_fasta(sys.stdin)
    if not seqs:
        print("ERROR: No sequences found", file=sys.stderr)
        sys.exit(1)

    filtered, trim_left, trim_right, dropped = filter_alignment(seqs)

    # Output
    for sid in filtered:
        print(f">{sid}")
        # Wrap at 80 chars
        s = filtered[sid]
        for i in range(0, len(s), 80):
            print(s[i:i+80])

    # Log to stderr
    orig_len = len(list(seqs.values())[0])
    new_len = len(list(filtered.values())[0])
    print(f"Filtered: {len(seqs)} seqs, {orig_len} → {new_len} bp "
          f"(trimmed {trim_left} left + {trim_right} right, "
          f"dropped {dropped} high-missing columns)",
          file=sys.stderr)

if __name__ == "__main__":
    main()
