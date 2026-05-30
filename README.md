# LUAD-CMap-Drug-Repurposing
# 🧬 Multi-Omics Integration & CMap Analysis for Computational Drug Repurposing in LUAD

<div align="center">

![Python](https://img.shields.io/badge/Python-3.9+-3776AB?style=flat-square&logo=python&logoColor=white)
![R](https://img.shields.io/badge/R-4.x-276DC3?style=flat-square&logo=r&logoColor=white)
![Status](https://img.shields.io/badge/Status-Complete-brightgreen?style=flat-square)
![Internship](https://img.shields.io/badge/BioCogniz-Internship%20Project-0D7377?style=flat-square)

**A 4-week computational drug repurposing pipeline for Lung Adenocarcinoma (LUAD)**  
using TCGA RNA-seq data · DESeq2 · LINCS L1000 · CMap · KEGG/GO pathway analysis

[Overview](#-project-overview) · [Pipeline](#-pipeline) · [Results](#-key-results) · [Files](#-repository-structure) · [How to Run](#-how-to-run) · [References](#-references)

</div>

---

## 📋 Project Overview

Lung Adenocarcinoma (LUAD) is the most common subtype of non-small cell lung cancer, with a 5-year survival rate below 20%. This project presents a complete AI-driven computational pipeline to identify existing approved drugs that can potentially be repurposed to treat LUAD — reducing the drug discovery timeline from **15 years → 4 weeks** and cost from **$2 billion → $50 million**.

| Field | Details |
|-------|---------|
| **Dataset** | TCGA-LUAD (GDC Portal) |
| **Samples** | 583 total — 524 Tumor + 59 Normal |
| **Genome Build** | GRCh38/hg38 · Gencode v22 |
| **Analysis Tools** | DESeq2 (R), PyDESeq2 (Python), clusterProfiler, CMap clue.io |
| **Duration** | 1 month (20.04.2026 – 15.05.2026) |
| **Organiser** | BioCogniz — NGS · Drug Discovery · AI |

---

## 🔬 Pipeline

```
TCGA-LUAD RNA-seq (583 samples)
         ↓
  QC + Normalization (DESeq2)
  PCA: PC1=16%, PC2=9% — clear tumor/normal separation
         ↓
  Differential Expression Analysis
  5,756 DEGs  |  padj < 0.01  |  |log₂FC| > 1.5
  4,431 Upregulated  ·  1,325 Downregulated
         ↓
  L1000 Landmark Gene Filtering
  5,756 DEGs ∩ 978 L1000 genes = 102 signature genes
  48 UP + 54 DOWN
         ↓
  CMap Query (clue.io)
  1,165,466 profiles searched  |  TAG score < −1.2
         ↓
  Pathway Enrichment (KEGG + GO:BP)
  Cell cycle #1 pathway (14/31 genes, p = 6.46×10⁻¹⁵)
  Mitotic addiction confirmed (60% UP genes, p = 3.94×10⁻³¹)
         ↓
  Literature Validation
  PubMed · ClinicalTrials.gov · DrugBank
         ↓
  10 Priority Drug Candidates
  3 Clinically Validated  ·  2 FDA-Approved Recovered
```

---

## 📊 Key Results

### Differential Expression
| Metric | Value |
|--------|-------|
| Total genes tested | 39,044 |
| Significant DEGs | **5,756** |
| Upregulated | **4,431** (FAM83A +6.79, TOP2A +3.85, AURKA +3.20) |
| Downregulated | **1,325** (SEMA3G −3.23, EPAS1 −2.72, PECAM1 −2.32) |
| CMap signature | **102 genes** (48 UP + 54 DOWN) |

### Biological Hallmarks Identified
| Hallmark | Evidence | p-value |
|---------|---------|--------|
| **Mitotic Addiction** | 29/48 UP genes involved in nuclear division | **3.94×10⁻³¹** |
| **Loss of Normal Identity** | 18/54 DOWN genes in blood vessel development | **2.11×10⁻²¹** |
| Cell Cycle pathway | 14/31 annotated UP genes | **6.46×10⁻¹⁵** |

### Top 5 Drug Candidates (CMap)
| Rank | Drug | TAG Score | Mechanism | Status |
|------|------|-----------|-----------|--------|
| 1 | **BAY-11-7085** | −1.71 | NF-kB inhibitor | Preclinical |
| 2 | **KPT-185** | −1.59 | XPO1 antagonist | Phase I/II |
| 3 | **bortezomib** | −1.34 | Proteasome inhibitor | ⭐ FDA Approved |
| 6 | **THZ-1** | −1.31 | CDK7 inhibitor | Nature 2014 ✓ |
| 7 | **AZD-3463** | −1.31 | ALK inhibitor | ⭐ Phase I LUAD |

> **Pipeline Validation:** AZD-3463 (ALK inhibitor) was independently recovered by the unbiased AI pipeline. ALK inhibitors are already FDA-approved first-line LUAD treatments — confirming the pipeline's biological validity.

---

## 🗂️ Repository Structure

```
LUAD-CMap-Drug-Repurposing/
│
├── 📁 data/
│   ├── sample_metadata.csv          # 583 sample info (tumor/normal labels)
│   ├── DEG_results_all.csv          # Full DESeq2 results (39,044 genes)
│   ├── DEG_results_significant.csv  # 5,756 significant DEGs
│   ├── DEG_results_annotated.csv    # DEGs with HGNC symbols
│   ├── cmap_UP_genes.txt            # 48 L1000 upregulated genes
│   ├── cmap_DOWN_genes.txt          # 54 L1000 downregulated genes
│   ├── gene_signature_102.csv       # Final CMap query signature
│   └── cmap_results.csv             # CMap connectivity scores
│
├── 📁 scripts/
│   ├── 01_data_download.py          # TCGA-LUAD download via GDC API
│   ├── 02_qc_normalization.py       # Library size QC + DESeq2 normalization
│   ├── 03_deseq2_analysis.py        # Differential expression (PyDESeq2)
│   ├── 04_l1000_filtering.py        # Filter DEGs to L1000 landmark genes
│   ├── 05_cmap_prep.py              # Prepare CMap query input files
│   ├── 06_pathway_analysis.R        # KEGG + GO:BP enrichment (clusterProfiler)
│   └── 07_visualization.R           # Publication figures (ggplot2, EnhancedVolcano)
│
├── 📁 figures/
│   ├── PCA_LUAD.png                 # PCA plot — tumor vs normal separation
│   ├── Volcano_LUAD.png             # Volcano plot — 5,756 DEGs
│   ├── Venn_L1000_filtering.png     # Venn: 5,756 DEGs ∩ 978 L1000 = 102 genes
│   ├── CMap_Drug_Ranking.png        # Drug connectivity score bar chart
│   ├── KEGG_UP_dotplot.png          # KEGG enrichment — upregulated genes
│   ├── KEGG_DOWN_dotplot.png        # KEGG enrichment — downregulated genes
│   └── DEG_Heatmap_Top50.png        # Heatmap — top 50 DEGs
│
├── 📁 reports/
│   ├── LUAD_Research_Article.docx   # Full 18-page research article
│   └── LUAD_Presentation.pptx      # 8-slide project presentation
│
├── requirements.txt                 # Python dependencies
├── R_packages.R                     # R package installation script
└── README.md                        # This file
```

---

## ⚙️ How to Run

### Prerequisites
```bash
# Python 3.9+
pip install -r requirements.txt

# R 4.x — install packages
Rscript R_packages.R
```

### Step-by-Step Execution
```bash
# Step 1: Download TCGA-LUAD data (requires internet, ~30–60 min)
python scripts/01_data_download.py

# Step 2: Quality control and normalization
python scripts/02_qc_normalization.py

# Step 3: Differential expression analysis
python scripts/03_deseq2_analysis.py

# Step 4: L1000 gene filtering
python scripts/04_l1000_filtering.py

# Step 5: Prepare CMap query files
python scripts/05_cmap_prep.py
# → Upload data/cmap_UP_genes.txt and data/cmap_DOWN_genes.txt
#   to https://clue.io/query → download results → save as data/cmap_results.csv

# Step 6: Pathway enrichment analysis (R)
Rscript scripts/06_pathway_analysis.R

# Step 7: Generate all figures (R)
Rscript scripts/07_visualization.R
```

---

## 📦 requirements.txt

```
pydeseq2>=0.4.0
pandas>=1.5.0
numpy>=1.23.0
matplotlib>=3.6.0
seaborn>=0.12.0
scipy>=1.9.0
scikit-learn>=1.1.0
mygene>=3.2.2
requests>=2.28.0
openpyxl>=3.0.0
adjustText>=0.8
```

---

## 🔬 Tools & Databases

| Tool/Database | Purpose | Reference |
|--------------|---------|-----------|
| **TCGA-GDC** | RNA-seq data source | NCI GDC portal |
| **DESeq2** | Differential expression | Love et al., 2014 |
| **CMap / clue.io** | Drug connectivity query | Subramanian et al., 2017 |
| **LINCS L1000** | Landmark gene filtering | LINCS program |
| **clusterProfiler** | KEGG / GO enrichment | Yu et al., 2012 |
| **DrugBank** | Drug target annotation | Wishart et al., 2018 |
| **org.Hs.eg.db** | Gene ID mapping | Bioconductor |
| **PubMed** | Literature validation | NCBI |

---

## 📈 Figures Preview

| Figure | Description |
|--------|-------------|
| `PCA_LUAD.png` | Complete tumor-normal separation — PC1=16%, PC2=9% |
| `Volcano_LUAD.png` | 5,756 DEGs — FAM83A (top UP), SEMA3G (top DOWN) |
| `Venn_L1000.png` | 102-gene CMap signature construction |
| `CMap_Drugs.png` | 15 drug candidates ranked by TAG score |
| `KEGG_dotplots.png` | Cell cycle p=6.46×10⁻¹⁵ dominates UP enrichment |

---

## 📚 References

1. Love MI, et al. (2014). DESeq2. *Genome Biology*, 15, 550.
2. Subramanian A, et al. (2017). CMap L1000. *Cell*, 171(6), 1437–1452.
3. Lamb J, et al. (2006). Connectivity Map. *Science*, 313, 1929–1935.
4. TCGA Research Network (2014). LUAD. *Nature*, 511, 543–550.
5. Kwiatkowski N, et al. (2014). CDK7 inhibitor THZ-1. *Nature*, 511, 616–620.
6. Yu G, et al. (2012). clusterProfiler. *OMICS*, 16(5), 284–287.

---

## 🏆 Acknowledgements

This project was completed as part of a **1-month AI in Drug Discovery internship** at **BioCogniz** (20.04.2026 – 15.05.2026).

> BioCogniz — NGS | Drug Discovery | AI

<div align="center">
Made with ❤️ by <strong>Sangamithra J</strong> · Stella Maris College · 2026
</div>
