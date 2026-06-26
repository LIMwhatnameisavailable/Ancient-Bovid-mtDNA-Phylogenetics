#!/usr/bin/env python3
"""
Step 6: k-mer based ML classifier for cross-validation of BLAST species assignment.
Uses k=4 (256-dim), Random Forest (200 trees) + SVM (RBF) on MOD reference panel.
"""
import os, sys, json, warnings
import numpy as np
import pandas as pd
from collections import Counter
warnings.filterwarnings('ignore')

BASE = "/mnt/disk2/srtp2024/LIM/IVPP/analysis"
OUT_DIR = os.path.join(BASE, "08_ml")
os.makedirs(OUT_DIR, exist_ok=True)

K = 4

def seq_to_kmer(seq, k=K):
    """Count k-mer frequencies, returns normalized vector."""
    seq = seq.upper().replace('-', '').replace('N', '')
    if len(seq) < k:
        return None
    counter = Counter()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if all(c in 'ACGT' for c in kmer):
            counter[kmer] += 1
    total = sum(counter.values())
    if total == 0:
        return None
    return {kmer: count/total for kmer, count in counter.items()}

def read_fasta(path):
    """Read FASTA, return {id: seq}."""
    seqs = {}
    cur_id = None
    cur_seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if cur_id:
                    seqs[cur_id] = ''.join(cur_seq)
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
    if cur_id:
        seqs[cur_id] = ''.join(cur_seq)
    return seqs

def build_feature_matrix(seq_dict, all_kmers=None):
    """Build feature matrix from sequence dictionary."""
    ids = []
    vectors = []
    valid_kmers = None
    for sid, seq in seq_dict.items():
        kmer_vec = seq_to_kmer(seq)
        if kmer_vec is None:
            continue
        ids.append(sid)
        if all_kmers is not None:
            vec = np.array([kmer_vec.get(k, 0.0) for k in all_kmers])
        else:
            if valid_kmers is None:
                valid_kmers = sorted(kmer_vec.keys())
            vec = np.array([kmer_vec.get(k, 0.0) for k in valid_kmers])
        vectors.append(vec)

    if all_kmers is None:
        all_kmers = valid_kmers

    return ids, np.array(vectors), all_kmers

print("=" * 60)
print("Step 6: k-mer ML Classifier")
print("=" * 60)

# ── 1. Load MOD reference panel (training set) ──
print("\n[1/6] Loading MOD reference panel...")
mod_dir = os.path.join(BASE, "02_ref", "renamed")
mod_seqs = {}
for fname in os.listdir(mod_dir):
    if fname.startswith("MOD_") and fname.endswith(".fa"):
        seqs = read_fasta(os.path.join(mod_dir, fname))
        mod_seqs.update(seqs)

mod_ids = list(mod_seqs.keys())
print(f"  Loaded {len(mod_ids)} MOD sequences")

# Create labels: yak vs cattle
mod_labels = []
for sid in mod_ids:
    if 'yak' in sid.lower() or 'mutus' in sid.lower():
        mod_labels.append('yak')
    else:
        mod_labels.append('cattle')

n_yak_mod = sum(1 for l in mod_labels if l == 'yak')
n_cattle_mod = sum(1 for l in mod_labels if l == 'cattle')
print(f"  Training: {n_yak_mod} yak + {n_cattle_mod} cattle")

# Build feature matrix from MOD
train_ids, X_train, all_kmers = build_feature_matrix(mod_seqs)
y_train = np.array([mod_labels[mod_ids.index(sid)] for sid in train_ids])
print(f"  Feature dim: {len(all_kmers)} k-mers (k={K})")
print(f"  Training matrix: {X_train.shape}")

# ── 2. Load Wenjiangduo samples (test set) ──
print("\n[2/6] Loading Wenjiangduo samples...")
wjd_seqs = read_fasta(os.path.join(BASE, "00_raw", "wenjiangduo_18seqs.fasta"))
test_ids, X_test, _ = build_feature_matrix(wjd_seqs, all_kmers)
print(f"  Test samples: {len(test_ids)}")

# ── 3. Train Random Forest ──
print("\n[3/6] Training Random Forest (200 trees)...")
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler

rf = RandomForestClassifier(n_estimators=200, max_depth=None,
                             random_state=42, n_jobs=-1, verbose=0)
rf.fit(X_train, y_train)
rf_pred = rf.predict(X_test)
rf_prob = rf.predict_proba(X_test)
print(f"  RF training complete")

# ── 4. Train SVM ──
print("\n[4/6] Training SVM (RBF kernel)...")
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

svm = SVC(kernel='rbf', probability=True, random_state=42, verbose=0)
svm.fit(X_train_scaled, y_train)
svm_pred = svm.predict(X_test_scaled)
svm_prob = svm.predict_proba(X_test_scaled)
print(f"  SVM training complete")

# ── 5. Load BLAST results for comparison ──
print("\n[5/6] Loading BLAST results for comparison...")
blast_path = os.path.join(BASE, "03_blast", "species_summary.csv")
blast_df = pd.read_csv(blast_path)
blast_dict = {}
for _, row in blast_df.iterrows():
    sid = row.get('sample_id', '')
    taxon = str(row.get('inferred_taxon', ''))
    if 'yak' in taxon.lower():
        blast_dict[sid] = 'yak'
    elif 'cattle' in taxon.lower() or 'taurus' in taxon.lower():
        blast_dict[sid] = 'cattle'
    else:
        blast_dict[sid] = 'unknown'

# ── 6. Compile results ──
print("\n[6/6] Compiling results...")
results = []
for i, sid in enumerate(test_ids):
    rf_label = rf_pred[i]
    rf_conf = max(rf_prob[i]) * 100
    svm_label = svm_pred[i]
    svm_conf = max(svm_prob[i]) * 100
    blast_label = blast_dict.get(sid, 'unknown')

    rf_match = 'SAME' if rf_label == blast_label else 'DIFF'
    svm_match = 'SAME' if svm_label == blast_label else 'DIFF'

    results.append({
        'sample_id': sid,
        'RF_prediction': rf_label,
        'RF_confidence': f"{rf_conf:.1f}%",
        'SVM_prediction': svm_label,
        'SVM_confidence': f"{svm_conf:.1f}%",
        'BLAST_species': blast_label,
        'RF_vs_BLAST': rf_match,
        'SVM_vs_BLAST': svm_match,
        'RF_SVM_agree': 'YES' if rf_label == svm_label else 'NO'
    })

df_out = pd.DataFrame(results)
out_path = os.path.join(OUT_DIR, "kmer_classification_result.csv")
df_out.to_csv(out_path, index=False)

# Summary stats
rf_agree = sum(1 for r in results if r['RF_vs_BLAST'] == 'SAME')
svm_agree = sum(1 for r in results if r['SVM_vs_BLAST'] == 'SAME')
rf_svm_agree = sum(1 for r in results if r['RF_SVM_agree'] == 'YES')
disagreements = [r['sample_id'] for r in results if r['RF_vs_BLAST'] == 'DIFF' or r['SVM_vs_BLAST'] == 'DIFF']

print(f"\n{'='*60}")
print(f"RESULTS SUMMARY")
print(f"{'='*60}")
print(f"RF vs BLAST agreement:  {rf_agree}/{len(results)} ({rf_agree/len(results)*100:.0f}%)")
print(f"SVM vs BLAST agreement: {svm_agree}/{len(results)} ({svm_agree/len(results)*100:.0f}%)")
print(f"RF-SVM agreement:       {rf_svm_agree}/{len(results)} ({rf_svm_agree/len(results)*100:.0f}%)")
print(f"\nDisagreements (RF or SVM ≠ BLAST): {len(disagreements)}")
for sid in disagreements:
    r = [x for x in results if x['sample_id'] == sid][0]
    print(f"  {sid}: RF={r['RF_prediction']} SVM={r['SVM_prediction']} BLAST={r['BLAST_species']}")

# Check priors
print(f"\nPrior check:")
z9863 = [r for r in results if r['sample_id'] == 'Z9863']
if z9863:
    r = z9863[0]
    print(f"  Z9863 (prior=cattle): RF={r['RF_prediction']} SVM={r['SVM_prediction']} BLAST={r['BLAST_species']}")
z9867 = [r for r in results if r['sample_id'] == 'Z9867']
if z9867:
    r = z9867[0]
    print(f"  Z9867 (prior=yak):    RF={r['RF_prediction']} SVM={r['SVM_prediction']} BLAST={r['BLAST_species']}")

# Feature importance (RF)
importances = rf.feature_importances_
top_k = 10
top_idx = np.argsort(importances)[-top_k:]
print(f"\nTop {top_k} discriminative k-mers:")
for idx in reversed(top_idx):
    print(f"  {all_kmers[idx]}: importance={importances[idx]:.4f}")

print(f"\nOutput: {out_path}")
print(f"{'='*60}")
