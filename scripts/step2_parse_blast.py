#!/usr/bin/env python3
"""
Step 2a/2b · BLAST 结果解析
合并 18 个 *-HitTable.csv → results.csv
解析 18 个 *-Descriptions.csv → species_summary.csv + discordant_or_low_confidence_candidates.csv
"""
import csv, os, re, sys
from datetime import datetime

BLAST_DIR = "/mnt/disk2/srtp2024/LIM/IVPP/analysis/03_blast"
HIT_DIR   = os.path.join(BLAST_DIR, "blast-HitTable")
DESC_DIR  = os.path.join(BLAST_DIR, "blast-DescriptionsTable")
OUT_DIR   = BLAST_DIR
LOG_DIR   = "/mnt/disk2/srtp2024/LIM/IVPP/analysis/logs"

SAMPLE_IDS = ["K0909","K0914","Z9850","Z9851","Z9852","Z9853","Z9854","Z9855",
              "Z9856","Z9857","Z9858","Z9859","Z9860","Z9861","Z9863","Z9864",
              "Z9865","Z9867"]

# ──────────────────────────────────────────
# 1. Merge HitTable CSVs
# ──────────────────────────────────────────
def parse_hittable(sid):
    """Parse a single HitTable CSV; skip BLAST header/footer lines."""
    path = os.path.join(HIT_DIR, f"{sid}-HitTable.csv")
    rows = []
    with open(path, newline="") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("BLAST"):
                continue
            # The NCBI HitTable (CSV) is comma-separated, no quoted fields:
            # query,sbjct,%identity,aln-len,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore
            parts = line.split(",")
            if len(parts) < 12:
                continue
            rows.append({
                "sample_id":          parts[0].strip(),
                "subject_acc_ver":    parts[1].strip(),
                "percent_identity":   float(parts[2]),
                "alignment_length":   int(parts[3]),
                "mismatches":         int(parts[4]),
                "gap_opens":          int(parts[5]),
                "q_start":            int(parts[6]),
                "q_end":              int(parts[7]),
                "s_start":            int(parts[8]),
                "s_end":              int(parts[9]),
                "evalue":             parts[10].strip(),
                "bit_score":          float(parts[11]),
            })
    return rows

def merge_hittables():
    cols = ["sample_id","subject_acc_ver","percent_identity","alignment_length",
            "mismatches","gap_opens","q_start","q_end","s_start","s_end",
            "evalue","bit_score"]
    all_rows = []
    for sid in SAMPLE_IDS:
        all_rows.extend(parse_hittable(sid))
    out_path = os.path.join(OUT_DIR, "results.csv")
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(all_rows)
    print(f"  → {out_path}  —  {len(all_rows)} rows from {len(SAMPLE_IDS)} samples")
    return all_rows

# ──────────────────────────────────────────
# 2. Parse Descriptions CSVs
# ──────────────────────────────────────────
def parse_description_csv(sid):
    """
    Parse a single Descriptions CSV.
    The NCBI export uses quoted Description, plus a HYPERLINK formula in Accession.
    We extract: Description, Scientific Name, Common Name, Taxid,
                Max Score, Total Score, Query Cover, E value, Per. ident,
                Acc. Len, Accession.
    """
    path = os.path.join(DESC_DIR, f"{sid}-Descriptions.csv")
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        # Read raw lines — we need to handle quoted commas
        content = f.read()

    # Split into lines, skip header
    lines = content.strip().split("\n")
    if len(lines) < 2:
        return rows

    # Parse header to know column order, then parse data rows
    # The header is: Description,Scientific Name,Common Name,Taxid,Max Score,Total Score,
    #                Query Cover,E value,Per. ident,Acc. Len,Accession
    #
    # Data rows have Description in quotes, which may contain commas.
    # Use csv.reader for proper handling.

    reader = csv.reader(lines, quotechar='"')
    for i, row in enumerate(reader):
        if i == 0:
            continue  # skip header
        if len(row) < 11:
            continue

        desc      = row[0].strip('"')
        sciname   = row[1].strip()
        common    = row[2].strip()
        taxid     = row[3].strip()
        max_score = row[4].strip()
        total_score = row[5].strip()
        qcover    = row[6].strip().replace("%", "")
        evalue    = row[7].strip()
        pident    = row[8].strip().replace("%", "")
        acc_len   = row[9].strip()
        acc_raw   = row[10].strip().strip('"')

        # Extract accession from HYPERLINK formula if needed
        # Format: =HYPERLINK("...","ACCESSION")
        m = re.search(r'"([^"]+)"\s*,\s*"([^"]+)"', acc_raw)
        if m:
            accession = m.group(2)
        else:
            accession = acc_raw

        rows.append({
            "sample_id":   sid,
            "description": desc,
            "scientific_name": sciname,
            "common_name": common,
            "taxid":       taxid,
            "max_score":   float(max_score) if max_score else 0,
            "total_score":  float(total_score) if total_score else 0,
            "query_cover":  float(qcover) if qcover else 0,
            "evalue":      evalue,
            "percent_identity": float(pident) if pident else 0,
            "acc_len":     int(acc_len) if acc_len else 0,
            "accession":   accession,
        })
    return rows

# ──────────────────────────────────────────
# 3. Infer taxon & confidence
# ──────────────────────────────────────────
# YAK species (Bos mutus, Bos grunniens)
YAK_SCINAMES   = {"Bos mutus", "Bos grunniens"}
YAK_COMMON     = {"wild yak", "domestic yak", "yak"}
CATTLE_SCINAMES = {"Bos taurus", "Bos indicus"}
CATTLE_COMMON   = {"domestic cattle", "cattle", "indicine cattle", "zebu"}

def infer_taxon(sciname, common):
    sciname_l = sciname.lower()
    common_l  = common.lower()

    # Yak check
    for s in YAK_SCINAMES:
        if s.lower() in sciname_l:
            return "Bos_mutus_group (yak)"
    for c in YAK_COMMON:
        if c in common_l:
            return "Bos_mutus_group (yak)"

    # Cattle check
    for s in CATTLE_SCINAMES:
        if s.lower() in sciname_l:
            return "Bos_taurus_group (cattle)"
    for c in CATTLE_COMMON:
        if c in common_l:
            return "Bos_taurus_group (cattle)"

    return f"Other ({sciname})"

# ──────────────────────────────────────────
# 4. Main
# ──────────────────────────────────────────
def main():
    print("Step 2a/2b · BLAST 结果解析")
    print(f"时间: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print()

    # 4a. Merge HitTable
    print("─" * 50)
    print("1. 合并 HitTable → results.csv")
    merge_hittables()

    # 4b. Parse Descriptions
    print()
    print("─" * 50)
    print("2. 解析 Descriptions → species_summary.csv")

    all_desc = []
    for sid in SAMPLE_IDS:
        rows = parse_description_csv(sid)
        if not rows:
            print(f"  ⚠ {sid}: 未找到数据行")
            continue
        all_desc.append(rows[0])  # top hit (first data row)

    # Build species_summary
    summary_cols = [
        "sample_id","top_species","top_accession","top_title",
        "common_name","taxid","query_cover","percent_identity",
        "evalue","max_score","total_score","acc_len",
        "inferred_taxon","confidence","qc_note"
    ]

    summary_rows = []
    discordant_rows = []

    # Priority samples to watch
    focus_ids = {"Z9853","Z9859","Z9864","Z9865","Z9863","Z9867"}

    for r in all_desc:
        sid     = r["sample_id"]
        sciname = r["scientific_name"]
        common  = r["common_name"]
        taxon   = infer_taxon(sciname, common)
        pident  = r["percent_identity"]
        qcover  = r["query_cover"]

        # Determine confidence
        qc_note = ""
        confidence = "HIGH"

        # Rules from ANALYSIS_STEPS.md:
        # pident < 99%  → candidate
        # Also track LOW coverage and UNCERTAIN morphology

        notes = []

        if pident < 99:
            notes.append(f"pident={pident:.2f}% (<99%)")
            confidence = "LOW"

        if qcover < 80:
            notes.append(f"query_cover={qcover:.0f}% (<80%)")
            confidence = "LOW"
        elif qcover < 95:
            notes.append(f"query_cover={qcover:.0f}% (<95%)")
            if confidence != "LOW":
                confidence = "MEDIUM"

        # Check prior discordance
        if sid == "Z9863":
            notes.append("⚠ PRIOR_CATTLE but BLAST→yak; 形态-母系谱系不一致候选")
            confidence = "DISCORDANT"
        if sid == "Z9867":
            notes.append("⚠ PRIOR_YAK but BLAST→cattle; 形态-母系谱系不一致候选")
            confidence = "DISCORDANT"

        # Low coverage
        if sid == "Z9853":
            notes.append("LOW coverage (2.1×), 30% N bases")
            if confidence == "HIGH":
                confidence = "MEDIUM"
        if sid == "Z9864":
            notes.append("LOW coverage (10.9×)")
            if confidence == "HIGH":
                confidence = "MEDIUM"

        # UNCERTAIN morphology
        if sid in {"Z9850","Z9857","Z9859","Z9864"}:
            notes.append("形态鉴定='牛？'(UNCERTAIN)")

        qc_note = " | ".join(notes) if notes else ""

        row = {
            "sample_id":       sid,
            "top_species":     sciname,
            "top_accession":   r["accession"],
            "top_title":       r["description"],
            "common_name":     common,
            "taxid":           r["taxid"],
            "query_cover":     f"{qcover:.0f}%",
            "percent_identity": f"{pident:.2f}%",
            "evalue":          r["evalue"],
            "max_score":       f"{r['max_score']:.0f}",
            "total_score":     f"{r['total_score']:.0f}",
            "acc_len":         str(r["acc_len"]),
            "inferred_taxon":  taxon,
            "confidence":      confidence,
            "qc_note":         qc_note,
        }
        summary_rows.append(row)

        # Discordant / low confidence if applicable
        if confidence in ("LOW", "DISCORDANT") or sid in focus_ids:
            discordant_rows.append(row)

    # Write species_summary.csv
    sum_path = os.path.join(OUT_DIR, "species_summary.csv")
    with open(sum_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=summary_cols)
        w.writeheader()
        w.writerows(summary_rows)
    print(f"  → {sum_path}  —  {len(summary_rows)} samples")

    # Write discordant_or_low_confidence_candidates.csv
    dis_path = os.path.join(OUT_DIR, "discordant_or_low_confidence_candidates.csv")
    with open(dis_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=summary_cols)
        w.writeheader()
        w.writerows(discordant_rows)
    print(f"  → {dis_path}  —  {len(discordant_rows)} candidates")

    # ── Print summary ──
    print()
    print("─" * 50)
    print("📊 BLAST 结果汇总")
    print()

    yak_samples   = [r for r in summary_rows if "yak" in r["inferred_taxon"]]
    cattle_samples = [r for r in summary_rows if "cattle" in r["inferred_taxon"]]
    discordant    = [r for r in summary_rows if r["confidence"] == "DISCORDANT"]
    low_conf      = [r for r in summary_rows if r["confidence"] == "LOW"]

    print(f"  牦牛 (Bos mutus/grunniens) 谱系:  {len(yak_samples)} 条")
    for r in yak_samples:
        print(f"    {r['sample_id']}  —  {r['percent_identity']} id, {r['query_cover']} cov")

    print(f"\n  普通牛 (Bos taurus) 谱系:  {len(cattle_samples)} 条")
    for r in cattle_samples:
        print(f"    {r['sample_id']}  —  {r['percent_identity']} id, {r['query_cover']} cov")

    print(f"\n  ⚠ 形态-母系谱系不一致候选 (DISCORDANT):  {len(discordant)} 条")
    for r in discordant:
        print(f"    {r['sample_id']}: {r['qc_note']}")

    print(f"\n  ⚠ 低置信度 (LOW/MEDIUM):  {len(low_conf)} 条")
    for r in low_conf:
        print(f"    {r['sample_id']}: {r['qc_note']}")

    # ── Write log ──
    log_path = os.path.join(LOG_DIR, "step2_parse_blast.log")
    with open(log_path, "w", encoding="utf-8") as f:
        f.write(f"# Step 2b · BLAST 结果解析 — 操作日志\n")
        f.write(f"# 时间: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"# 解析: CC (Python script)\n\n")
        f.write("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")
        f.write("## 输入\n")
        f.write(f"  HitTable:     {HIT_DIR}/  — 18 个 CSV\n")
        f.write(f"  Descriptions: {DESC_DIR}/  — 18 个 CSV\n\n")
        f.write("## 输出\n")
        f.write(f"  1. {sum_path}\n")
        f.write(f"  2. {dis_path}\n\n")
        f.write("## 关键指标\n")
        f.write(f"  样品总数: 18\n")
        f.write(f"  牦牛谱系 (Bos mutus/grunniens): {len(yak_samples)} 条\n")
        for r in yak_samples:
            f.write(f"    {r['sample_id']} — {r['percent_identity']} id, {r['query_cover']} cov\n")
        f.write(f"  普通牛谱系 (Bos taurus): {len(cattle_samples)} 条\n")
        for r in cattle_samples:
            f.write(f"    {r['sample_id']} — {r['percent_identity']} id, {r['query_cover']} cov\n")
        f.write(f"  不一致候选 (DISCORDANT): {len(discordant)} 条\n")
        for r in discordant:
            f.write(f"    {r['sample_id']}: {r['qc_note']}\n")
        f.write(f"  低置信度: {len(low_conf)} 条\n")
        for r in low_conf:
            f.write(f"    {r['sample_id']}: {r['qc_note']}\n")
        f.write(f"\n## 先验检查\n")
        f.write(f"  Z9863 形态=普通牛 → BLAST 归属: {summary_rows[SAMPLE_IDS.index('Z9863')]['inferred_taxon']}\n")
        f.write(f"  Z9867 形态=牦牛 → BLAST 归属: {summary_rows[SAMPLE_IDS.index('Z9867')]['inferred_taxon']}\n")
        f.write(f"\n## 状态\n")
        f.write(f"  DONE\n")
    print(f"\n  → {log_path}")

if __name__ == "__main__":
    main()