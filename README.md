# 🐂 Ancient-Bovid-mtDNA-Phylogenetics

![Python](https://img.shields.io/badge/Python-3.8+-blue)
![IQ-TREE2](https://img.shields.io/badge/IQ--TREE2-2.x-green)
![MAFFT](https://img.shields.io/badge/MAFFT-7.x-orange)
![License](https://img.shields.io/badge/License-MIT-lightgrey)

**Mitochondrial genome-based phylogenetic identification of ancient bovids from Wenjiangduo site (4–5th century CE, Tibet Plateau).** This project implements a complete bioinformatics pipeline integrating BLAST taxonomy, IQ-TREE2 maximum-likelihood phylogeny, haplotype network analysis, and k-mer random forest classification across 18 ancient DNA (aDNA) sequences to resolve species identity and morphology–genome discordance.

---

## 📌 Project Overview

The Wenjiangduo site (温江多遗址), located on the Tibet Plateau, yielded 18 bovid skeletal specimens dated to the 4–5th century CE. This project addresses a central question in zooarchaeology: **do morphological identifications agree with mitochondrial genomic evidence?**

Three independent computational approaches — BLAST homology search, ML phylogenetics, and k-mer machine learning — converge on a consistent species composition, while also revealing two specimens with clear **morphology–mtDNA discordance**, providing direct genomic evidence for bovid management and hybridization on the ancient Plateau.

### Key Findings

- **Species composition**: 7 yak (*Bos mutus*) + 11 taurine cattle (*Bos taurus*) confirmed across all three methods
- **Morphology–genome discordance**: Z9863 (morphologically cattle, mtDNA yak) and Z9867 (morphologically yak, mtDNA cattle)
- **Lineage resolution**: Sample K0914 assigned to the YakX lineage (bootstrap support = 99)
- **Genetic diversity**: Yak π = 0.21%; cattle π = 0.058%
- **Classifier validation**: Random forest (k-mer) reliable; SVM rejected due to overfitting (p >> n)

---

## 📂 Repository Structure

```
Ancient-Bovid-mtDNA-Phylogenetics/
├── 00_raw/                         # 18 ancient mtDNA consensus sequences (FASTA)
├── 01_qc/                          # Coverage statistics and quality-tier assignment
│   └── coverage_summary.csv        # HIGH / MID / LOW classification per sample
├── 02_ref/                         # Reference panel construction
│   ├── renamed/                    # 51 curated reference sequences (renamed FASTA)
│   │   └── all_combined.fa         # Merged reference + outgroup panel
│   └── README.md                   # Reference accession metadata
├── 03_blast/                       # BLAST-based species identification
│   ├── blast-DescriptionsTable/    # Per-sample top-hit description tables (CSV)
│   └── blast-HitTable/             # Per-sample full hit tables (CSV)
├── 04_align/                       # Multiple sequence alignment (MAFFT)
│   ├── final_aln.fasta             # Full alignment (18 ancient + 51 refs)
│   └── final_trimmed.fasta         # TrimAl-trimmed alignment for tree inference
├── 04_tree_identification/         # Species ID tree (IQ-TREE2, rapid bootstrap)
│   ├── id_tree.treefile            # ML tree (Newick)
│   └── sample_grouping.csv        # Per-sample species assignment from tree
├── 05_phylo/                       # Final phylogenetic analysis
│   ├── final_tree.treefile         # Full ML tree with all 18 samples
│   └── remove_lowQ_tree.treefile   # Sensitivity test: LOW-quality samples excluded
├── 06_haplo/                       # Haplotype network analysis
│   ├── haplotype_assignment.csv    # Per-sample haplotype ID
│   ├── diversity_by_clade.csv      # π and haplotype diversity by species
│   └── popart_input.nex            # PopART-compatible NEXUS file
├── 07_beast/                       # [Optional] BEAST2 BSP — not yet executed
├── 08_ml/                          # k-mer random forest classifier
│   └── kmer_classification_result.csv
├── 09_figures/                     # All publication-ready figures (PNG + SVG)
│   ├── fig0_workflow.svg           # Analysis workflow diagram
│   ├── fig1_qc_coverage.svg        # Coverage distribution by quality tier
│   ├── fig2_species_composition.svg
│   ├── fig3_id_tree.svg            # Species identification ML tree
│   ├── fig4_quality_scatter.svg
│   ├── fig5_site_comparison.svg
│   ├── fig6_blast_identity.svg
│   ├── fig7_morphology_vs_mtdna.svg # Discordance visualization
│   └── fig8_haplotype_network.svg
├── logs/                           # Step-by-step run logs (fully reproducible)
├── scripts/                        # 13 custom Python analysis scripts
│   ├── step2_parse_blast.py
│   ├── step3_filter_alignment.py
│   ├── step3_postprocess.py
│   ├── step4_haplotype.py
│   ├── step4_plot_network.py
│   ├── step4_to_nexus.py
│   ├── step6_kmer_classifier.py
│   ├── trim_alignment.py
│   ├── plot_fig0_workflow.py
│   ├── plot_figures_1to7.py
│   ├── plot_step3_remove_lowQ.py
│   └── plot_step3_tree.py
├── ANALYSIS_STEPS.md               # Step-by-step reproduction guide
├── SAMPLE_META.md                  # Sample metadata and morphological records
├── CLAUDE.md                       # AI-assisted development notes
└── .gitignore
```

> **Note:** Raw sequencing data (FASTQ/BAM) and BWA index files are excluded from this repository due to size constraints. Only consensus FASTA sequences and all downstream analysis files are included.

---

## 🔬 Analysis Pipeline

The pipeline consists of 9 sequential stages:

```
00_raw → 01_qc → 02_ref → 03_blast → 04_align → 05_phylo → 06_haplo → 08_ml → 09_figures
```

### Stage 1 · Quality Control (`01_qc`)
Coverage-based quality tiering of 18 ancient mtDNA sequences:
- **HIGH**: mean depth ≥ 10×, used in all analyses
- **MID**: depth 5–10×, included with caution
- **LOW**: depth < 5×, excluded in sensitivity tests

### Stage 2 · Reference Panel Construction (`02_ref`)
A curated panel of **51 reference sequences** was assembled from four literature sources:
- Modern yak (*Bos mutus*): 5 sequences (MOD_Bmut)
- Ancient yak (YakA / YakX lineages): 10 sequences (LIT1, LIT2)
- Taurine cattle (*Bos taurus*, T1/T3 haplogroups): 20 sequences (LIT3, LIT4)
- Domestic cattle × yak hybrids: 1 sequence
- Outgroups: African buffalo (*Syncerus caffer*, NC020617) + domestic water buffalo (*Bubalus bubalis*, NC006295)

### Stage 3 · BLAST Identification (`03_blast`)
Each of the 18 ancient sequences was queried against the NCBI nt database. Top hits were parsed and summarized per sample. Species assignment was based on majority-rule across top 10 hits.

### Stage 4 · Multiple Sequence Alignment (`04_align`)
All 18 ancient sequences were aligned with the 51-reference panel using **MAFFT** (FFT-NS-2 strategy). Alignment was trimmed with **TrimAl** (automated1 mode) to remove hypervariable flanking regions.

### Stage 5 · Phylogenetic Inference (`05_phylo`)
Maximum-likelihood trees were inferred with **IQ-TREE2**:
- Model selection: ModelFinder (best-fit = GTR+F+I+G4)
- Branch support: 1000 ultrafast bootstraps (UFBoot2)
- Rooting: *Syncerus caffer* (African buffalo) as outgroup
- Sensitivity test: LOW-quality samples removed and tree re-inferred

### Stage 6 · Haplotype Analysis (`06_haplo`)
- Haplotype assignment performed on the trimmed alignment
- Nucleotide diversity (π) calculated per species clade
- NEXUS file generated for PopART visualization

### Stage 7 · k-mer Machine Learning (`08_ml`)
A **Random Forest classifier** was trained on k-mer frequency profiles (k=4) of reference sequences and applied to classify the 18 ancient samples:
- Cross-validation accuracy: >97% on reference panel
- SVM classifier was tested but rejected (overfitting, p >> n)

### Stage 8 · Figure Generation (`09_figures`)
22 publication-ready figures generated in PNG (300 dpi) + SVG format using custom Python scripts (Matplotlib / ETE3).

---

## 🧬 Sample Information

| Sample ID | Site | Morphological ID | mtDNA Assignment | Discordant? |
|-----------|------|-----------------|-----------------|-------------|
| Z9850 | Wenjiangduo | Cattle | *Bos taurus* | No |
| Z9851 | Wenjiangduo | Cattle | *Bos taurus* | No |
| Z9852 | Wenjiangduo | Cattle | *Bos taurus* | No |
| Z9853 | Wenjiangduo | Cattle | *Bos taurus* | No |
| Z9854 | Wenjiangduo | Yak | *Bos mutus* | No |
| Z9855 | Wenjiangduo | Yak | *Bos mutus* | No |
| Z9856 | Wenjiangduo | Yak | *Bos mutus* | No |
| Z9857 | Wenjiangduo | Yak | *Bos mutus* | No |
| Z9858 | Wenjiangduo | Cattle | *Bos taurus* | No |
| Z9859 | Wenjiangduo | Cattle | *Bos taurus* | No |
| Z9860 | Wenjiangduo | Cattle | *Bos taurus* | No |
| Z9861 | Wenjiangduo | Cattle | *Bos taurus* | No |
| **Z9863** | Wenjiangduo | **Cattle** | ***Bos mutus*** | **⚠️ Yes** |
| Z9864 | Wenjiangduo | Cattle | *Bos taurus* | No |
| Z9865 | Wenjiangduo | Cattle | *Bos taurus* | No |
| **Z9867** | Wenjiangduo | **Yak** | ***Bos taurus*** | **⚠️ Yes** |
| K0909 | Wenjiangduo | Cattle | *Bos taurus* | No |
| K0914 | Wenjiangduo | Yak | *Bos mutus* (YakX) | No |

---

## 🛠️ Dependencies

| Tool / Package | Version | Purpose |
|----------------|---------|---------|
| Python | ≥ 3.8 | Scripting and ML |
| MAFFT | ≥ 7.490 | Multiple sequence alignment |
| TrimAl | ≥ 1.4.1 | Alignment trimming |
| IQ-TREE2 | ≥ 2.2.0 | ML phylogenetic inference |
| BLAST+ | ≥ 2.13 | Sequence homology search |
| BioPython | ≥ 1.79 | Sequence I/O |
| scikit-learn | ≥ 1.1 | k-mer random forest classifier |
| Matplotlib | ≥ 3.5 | Figure generation |
| ETE3 | ≥ 3.1 | Tree visualization |
| pandas / numpy | latest | Data processing |

---

## ▶️ Reproduction

All steps are logged in `logs/` and documented in `ANALYSIS_STEPS.md`. To reproduce the downstream analysis from the provided FASTA sequences:

```bash
# 1. Clone the repository
git clone https://github.com/LIMwhatnameisavailable/Ancient-Bovid-mtDNA-Phylogenetics.git
cd Ancient-Bovid-mtDNA-Phylogenetics

# 2. Install Python dependencies
pip install biopython scikit-learn matplotlib ete3 pandas numpy

# 3. Run alignment and trimming
mafft --auto 04_align/final.fasta > 04_align/final_aln.fasta
trimal -in 04_align/final_aln.fasta -out 04_align/final_trimmed.fasta -automated1

# 4. Infer ML tree
iqtree2 -s 04_align/final_trimmed.fasta -m TEST -B 1000 \
        --prefix 05_phylo/final_tree -T AUTO

# 5. Run haplotype analysis
python scripts/step4_haplotype.py

# 6. Run k-mer classifier
python scripts/step6_kmer_classifier.py

# 7. Generate all figures
python scripts/plot_figures_1to7.py
python scripts/plot_fig0_workflow.py
```

---

## 📊 Output Figures

| Figure | Description |
|--------|-------------|
| `fig0_workflow.svg` | Full analysis workflow diagram |
| `fig1_qc_coverage.svg` | Coverage distribution by quality tier |
| `fig2_species_composition.svg` | Pie/bar chart of species composition |
| `fig3_id_tree.svg` | Species identification ML tree |
| `fig4_quality_scatter.svg` | Quality metrics scatter plot |
| `fig5_site_comparison.svg` | Inter-site comparison |
| `fig6_blast_identity.svg` | BLAST percent identity distribution |
| `fig7_morphology_vs_mtdna.svg` | Morphology vs. mtDNA discordance matrix |
| `fig8_haplotype_network.svg` | Median-joining haplotype network |

---

## 📖 Citation & Data Sources

Reference sequences were compiled from the following publications:

- **LIT1** — Chen et al. (2023): Ancient yak YakA lineage sequences
- **LIT2** — [YakX lineage ancient sequences, AG/UP/BS series]
- **LIT3** — Taurine cattle P1a haplogroup (LC537308–LC537317; DQ124389)
- **LIT4** — Cai et al. (2025): Southwest China cattle and hybrid sequences
- **MOD** — Modern reference genomes from NCBI (NC025563, KM233417, KR106993, KY829451, MK033130, OQ513044–OQ513048)
- **Outgroups** — *Syncerus caffer* NC020617; *Bubalus bubalis* NC006295; *Capra hircus* NC020617

---

## 👤 Author

**LIMwhatnameisavailable**
Southeast University · Institute of Vertebrate Paleontology and Paleoanthropology (IVPP), CAS
ORCID: [0009-0009-8876-0154](https://orcid.org/0009-0009-8876-0154)

---

## 📄 License

This project is licensed under the MIT License.
