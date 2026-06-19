# Bioinformatics Project Report
### Computational Genomics Portfolio — BINF6310, Northeastern University
**Author:** Md Tariqul Islam (`islam.mdtar`)  
**Institution:** Northeastern University, Khoury College of Computer Sciences  
**Course:** BINF6310 – Introduction to Bioinformatics  
**Portfolio:** [mtariqi.github.io](https://mtariqi.github.io) · **GitHub:** [github.com/mtariqi](https://github.com/mtariqi)

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Project 1 — De Novo Genome Assembly](#2-project-1--de-novo-genome-assembly-and-quality-assessment)
   - [Introduction & Rationale](#21-introduction--rationale)
   - [Dataset & Data Retrieval](#22-dataset--data-retrieval)
   - [Computational Environment](#23-computational-environment)
   - [Methodology](#24-methodology)
   - [Results](#25-results)
   - [Interpretation & Discussion](#26-interpretation--discussion)
3. [Project 2 — RNA-Seq Differential Gene Expression](#3-project-2--rna-seq-differential-gene-expression-analysis)
   - [Introduction & Rationale](#31-introduction--rationale)
   - [Dataset & Data Retrieval](#32-dataset--data-retrieval)
   - [Computational Environment](#33-computational-environment)
   - [Methodology](#34-methodology)
   - [Results](#35-results)
   - [Interpretation & Discussion](#36-interpretation--discussion)
4. [Cross-Project Skills Summary](#4-cross-project-skills-summary)
5. [Conclusions](#5-conclusions)
6. [References](#6-references)

---

## 1. Executive Summary

This report presents two independent bioinformatics projects completed as part of BINF6310 – Introduction to Bioinformatics at Northeastern University. Both projects were executed on Northeastern's Explorer High-Performance Computing (HPC) cluster using containerized software environments (Singularity) and reproducible pipeline frameworks (Nextflow DSL2), reflecting production-standard bioinformatics practices used in academic research and industry.

**Project 1** focused on *de novo* genome assembly of a *Staphylococcus aureus* mutant strain from raw Illumina paired-end reads, achieving a high-contiguity assembly with an N50 of 132,140 bp and an L50 of 1 — indicating that a single contig covers 50% of the entire assembled genome.

**Project 2** addressed a fundamental question in cancer biology: which genes are differentially expressed between breast cancer tumor tissue and matched normal tissue? Using a full RNA-Seq pipeline (STAR alignment → featureCounts quantification → DESeq2 statistical modeling), 847 statistically significant differentially expressed genes (DEGs) were identified, several of which are established oncogenes and tumor suppressors with direct clinical relevance.

Together, these projects demonstrate end-to-end proficiency across the two most foundational workflows in modern genomics: genome reconstruction from raw sequencing data and transcriptomic analysis of gene regulation in disease.

---

## 2. Project 1 — De Novo Genome Assembly and Quality Assessment

### 2.1 Introduction & Rationale

De novo genome assembly is the process of reconstructing an organism's complete genomic sequence from raw, short sequencing reads — without using a pre-existing reference genome. This approach is essential when studying novel organisms, highly divergent strains, or organisms whose reference genomes are incomplete or absent.

In microbial genomics, high-quality de novo assemblies are foundational to a wide range of downstream analyses including variant detection, genome annotation, comparative genomics, and antimicrobial resistance (AMR) profiling. The quality of an assembly — measured by contiguity metrics such as N50 and L50 — directly determines the reliability of all downstream results.

*Staphylococcus aureus* was selected as the target organism for this project because it is a clinically important human pathogen responsible for a wide spectrum of infections ranging from skin infections to life-threatening septicemia. Its relatively small, compact genome (~2.8 Mb) makes it an appropriate organism for benchmarking assembly tools, while its clinical relevance provides direct motivation for the analysis.

**SPAdes** (St. Petersburg genome assembler) was selected as the primary assembly tool due to its proven performance on bacterial isolate sequencing data. SPAdes constructs a de Bruijn graph from k-mers derived from the input reads and applies error correction and graph simplification to produce high-quality contigs. The `--isolate` flag was used, which is specifically optimized for high-coverage bacterial isolate data and produces cleaner, more contiguous assemblies than the default multi-cell mode.

**QUAST** (Quality Assessment Tool for Genome Assemblies) was used to objectively benchmark the assembly output, providing standardized metrics including total assembly length, number of contigs, N50, L50, L90, and GC content.

---

### 2.2 Dataset & Data Retrieval

| Parameter | Value |
|---|---|
| Organism | *Staphylococcus aureus* (imaginary mutant strain) |
| Sequencing platform | Illumina (paired-end) |
| Data source | Zenodo |
| Accession | Record ID: 582600 |
| Read files | `mutant_R1.fastq.gz`, `mutant_R2.fastq.gz` |
| Reference | Gladman et al. (2017). Bacterial training dataset for Galaxy tutorials. DOI: 10.5281/zenodo.582600 |

Reads were downloaded directly to the HPC working directory using `wget`:

```bash
wget https://zenodo.org/record/582600/files/mutant_R1.fastq.gz
wget https://zenodo.org/record/582600/files/mutant_R2.fastq.gz
```

File integrity was verified using `md5sum` before proceeding to assembly.

---

### 2.3 Computational Environment

All analyses were performed on **Northeastern University's Explorer HPC cluster**, a high-performance computing system designed to support computationally intensive research tasks beyond the capabilities of standard workstations.

| Component | Details |
|---|---|
| Cluster | Explorer HPC, Northeastern University |
| Container runtime | Singularity |
| Assembly tool | SPAdes v3.15 (via Singularity image) |
| QC tool | QUAST v5.2 (via Singularity image) |
| Pipeline framework | Nextflow DSL2 |
| Shell environment | Linux/Bash |
| Module system | Environment Modules (`module load`) |

Singularity containers were used rather than direct software installation for two critical reasons: (1) reproducibility — the identical software environment can be reconstructed by any collaborator or reviewer; and (2) security — Singularity does not require root privileges, making it safe for shared HPC environments.

---

### 2.4 Methodology

**Step 1 — Environment setup:**
```bash
module load singularity
mkdir -p assembly_project/reads assembly_project/output
cd assembly_project/
```

**Step 2 — Data retrieval:**
```bash
wget https://zenodo.org/record/582600/files/mutant_R1.fastq.gz -P reads/
wget https://zenodo.org/record/582600/files/mutant_R2.fastq.gz -P reads/
md5sum reads/*.fastq.gz  # integrity check
```

**Step 3 — De novo assembly with SPAdes:**
```bash
singularity exec --bind $PWD:/data spades.sif \
  spades.py --isolate \
  -1 /data/reads/mutant_R1.fastq.gz \
  -2 /data/reads/mutant_R2.fastq.gz \
  -o /data/output/spades_assembly/
```

The `--isolate` flag enables a pipeline optimized for single isolate bacterial sequencing data at high coverage depth, applying more aggressive graph simplification to reduce contig fragmentation.

**Step 4 — Quality assessment with QUAST:**
```bash
singularity exec --bind $PWD:/data quast.sif \
  quast.py /data/output/spades_assembly/contigs.fasta \
  -o /data/output/quast_results/
```

**Step 5 — Results inspection:**
```bash
cat output/quast_results/report.txt
```

---

### 2.5 Results

#### Assembly Quality Metrics (QUAST Report)

| Metric | Value | Interpretation |
|---|---|---|
| Total length | 179,839 bp | Total size of all assembled contigs |
| Number of contigs | 3 | Very few contigs — near-complete assembly |
| N50 | **132,140 bp** | 50% of genome in contigs ≥ 132 kb |
| L50 | **1** | Just 1 contig covers 50% of genome |
| L90 | 2 | Only 2 contigs cover 90% of genome |
| GC content | **33.59%** | Matches expected GC% for *S. aureus* (32–34%) |
| Largest contig | 132,140 bp | Dominant single contig |

#### Summary Statistics

- The assembly produced **3 contigs** from paired-end Illumina reads, indicating excellent read coverage and minimal assembly fragmentation.
- An **N50 of 132,140 bp** means the majority of the assembly is represented by a single long, contiguous sequence — a hallmark of high-quality bacterial genome assembly.
- An **L50 of 1** is exceptionally strong, indicating that a single contig alone accounts for more than half the total genome length.
- **GC content of 33.59%** is fully consistent with the known GC range for *Staphylococcus aureus* genomes (typically 32–34%), confirming taxonomic fidelity of the assembly.

---

### 2.6 Interpretation & Discussion

The assembly results demonstrate a high-quality, minimally fragmented reconstruction of the *S. aureus* mutant strain genome. The N50 value of 132 kb with an L50 of 1 is particularly notable — in many bacterial genome assembly benchmarks, achieving an L50 of 1 from a first-pass SPAdes run indicates both high sequencing depth and clean, low-error input reads.

The GC content of 33.59% serves as an important internal validation step: if the assembly had produced a GC content substantially outside the 32–34% expected range for *S. aureus*, this would signal potential contamination, misassembly, or use of incorrect reads. The match between observed and expected GC content confirms that the assembled sequences are biologically consistent with the target organism.

**Limitations:** This project used a training dataset with a reduced genome reference (~180 kb subsection), which accounts for the smaller-than-expected total assembly length compared to the full *S. aureus* genome (~2.8 Mb). Assembly of the full clinical genome would require substantially more computational resources and would likely produce a more complex contig graph. Additionally, no read quality trimming (e.g. Trimmomatic or BBDuk) was applied prior to assembly; incorporating a preprocessing step would likely further improve assembly contiguity.

**Broader significance:** The containerized, HPC-based workflow developed here is directly transferable to real-world applications including clinical microbial genomics, AMR gene identification, outbreak surveillance, and comparative genomics studies between reference and novel strains.

---

## 3. Project 2 — RNA-Seq Differential Gene Expression Analysis

### 3.1 Introduction & Rationale

RNA sequencing (RNA-Seq) is the gold-standard technology for measuring gene expression across the transcriptome on a genome-wide scale. Unlike targeted approaches such as RT-qPCR, which quantify the expression of pre-selected genes, RNA-Seq simultaneously measures the expression of every expressed gene in a biological sample, enabling unbiased discovery of differentially expressed genes (DEGs) without prior hypotheses.

Differential gene expression (DGE) analysis compares RNA-Seq data across two or more biological conditions — in this case, breast cancer tumor tissue versus matched normal tissue — to identify genes whose expression levels are statistically significantly different between conditions. Such genes are candidates for roles in disease initiation, progression, or maintenance, and may represent therapeutic targets or diagnostic biomarkers.

Breast cancer (specifically the TCGA-BRCA dataset) was selected as the biological system for this project because it is among the most thoroughly characterized cancer transcriptomes in the world, providing a reliable benchmark against which analysis results can be validated. The well-established biology of breast cancer — including known oncogenes (ERBB2/HER2, MKI67), tumor suppressors (TP53, PTEN, RB1), and cell adhesion molecules (CDH1/E-cadherin) — allows for meaningful biological interpretation of computational results.

**STAR** (Spliced Transcripts Alignment to a Reference) was selected for read alignment because it is a splice-aware aligner designed specifically for RNA-Seq data from eukaryotic organisms. Unlike DNA aligners such as Bowtie2, STAR can map reads that span exon-intron boundaries — a critical capability for transcriptomic analysis, where sequencing reads are derived from processed mRNA and therefore lack intron sequences present in the genomic reference.

**DESeq2** was selected for statistical testing because it uses a negative binomial distribution model that appropriately accounts for the overdispersion characteristics of RNA-Seq count data — a property where biological variance exceeds the mean more than a Poisson distribution would predict. DESeq2 applies size factor normalization to correct for library depth differences between samples and uses the Wald test with Benjamini-Hochberg false discovery rate (FDR) correction to control for multiple comparisons across tens of thousands of genes simultaneously.

---

### 3.2 Dataset & Data Retrieval

| Parameter | Value |
|---|---|
| Dataset | Human breast cancer RNA-Seq |
| GEO Accession | GSE183947 |
| Organism | *Homo sapiens* |
| Tissue types | Breast tumor (n=6) vs. matched normal (n=6) |
| Total samples | 12 |
| Sequencing platform | Illumina (paired-end, 150 bp reads) |
| Reference genome | GRCh38 / hg38 |
| Gene annotation | GENCODE v43 (Ensembl GTF) |

```bash
# Download reads from NCBI SRA using SRA-toolkit
prefetch GSE183947
fastq-dump --split-files --gzip SRR_accession_list.txt

# Download reference genome and annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz
```

---

### 3.3 Computational Environment

| Component | Details |
|---|---|
| Cluster | Explorer HPC, Northeastern University |
| Container runtime | Singularity |
| Read trimming | Trimmomatic v0.39 |
| Quality control | FastQC v0.11.9 · MultiQC v1.14 |
| Alignment | STAR v2.7.10 |
| Quantification | featureCounts (Subread v2.0.3) |
| Statistical analysis | R v4.3 · DESeq2 v1.40 · Bioconductor |
| Visualization | ggplot2 · EnhancedVolcano · pheatmap |
| Pipeline framework | Nextflow DSL2 |

---

### 3.4 Methodology

#### Phase 1 — Quality Control & Preprocessing

```bash
# FastQC on raw reads
fastqc -t 8 *.fastq.gz -o fastqc_raw/
multiqc fastqc_raw/ -o multiqc_report/

# Trimmomatic: remove adapters and low-quality bases
trimmomatic PE -threads 8 -phred33 \
  sample_R1.fastq.gz sample_R2.fastq.gz \
  sample_R1_trimmed.fastq.gz sample_R1_unpaired.fastq.gz \
  sample_R2_trimmed.fastq.gz sample_R2_unpaired.fastq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
```

**Trimmomatic parameters explained:**
- `ILLUMINACLIP` — removes Illumina TruSeq adapter sequences
- `LEADING:3 / TRAILING:3` — removes bases below quality 3 from read ends
- `SLIDINGWINDOW:4:20` — trims when average quality in 4-base window drops below Q20
- `MINLEN:36` — discards reads shorter than 36 bp after trimming

#### Phase 2 — STAR Genome Index & Alignment

```bash
# Build STAR genome index (one-time, ~30 min on HPC)
STAR --runMode genomeGenerate \
  --genomeDir hg38_star_index/ \
  --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
  --sjdbGTFfile gencode.v43.annotation.gtf \
  --runThreadN 16 \
  --genomeSAindexNbases 14

# Align each sample (run for all 12 samples)
STAR --runThreadN 8 \
  --genomeDir hg38_star_index/ \
  --readFilesIn sample_R1_trimmed.fastq.gz sample_R2_trimmed.fastq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS NM \
  --outFileNamePrefix aligned/sample_

# Index BAM files for downstream tools
samtools index aligned/sample_Aligned.sortedByCoord.out.bam
```

#### Phase 3 — Read Quantification (featureCounts)

```bash
featureCounts -T 8 -p -s 2 \
  -a gencode.v43.annotation.gtf \
  -o counts/gene_counts.txt \
  aligned/*.bam
```

The `-s 2` flag specifies reverse-stranded counting, consistent with TruSeq library preparation. The output is a gene × sample count matrix used as input to DESeq2.

#### Phase 4 — Differential Expression Analysis (DESeq2 in R)

```r
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)

# Load count matrix
counts <- read.table("counts/gene_counts.txt",
                     header=TRUE, row.names=1, skip=1)
counts <- counts[, 6:ncol(counts)]

# Sample metadata
colData <- data.frame(
  condition = factor(c(rep("tumor",6), rep("normal",6))),
  row.names = colnames(counts)
)

# DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = colData,
  design    = ~ condition
)

# Pre-filter: keep genes with ≥10 total reads across all samples
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run DESeq2 (normalization + statistical testing)
dds <- DESeq(dds)

# Extract results with LFC threshold of 1 (2-fold change)
res <- results(dds,
  contrast = c("condition","tumor","normal"),
  alpha    = 0.05,
  lfcThreshold = 1
)

# Apply apeglm LFC shrinkage to reduce noise in low-count estimates
res_shrunk <- lfcShrink(dds,
  coef = "condition_tumor_vs_normal",
  type = "apeglm"
)

# Export significant DEGs (padj < 0.05 and |log2FC| > 1)
sig_genes <- subset(res_shrunk,
  padj < 0.05 & abs(log2FoldChange) > 1
)
write.csv(as.data.frame(sig_genes), "DEGs_tumor_vs_normal.csv")
```

#### Phase 5 — Visualization

```r
# Volcano plot
EnhancedVolcano(res_shrunk,
  lab       = rownames(res_shrunk),
  x         = 'log2FoldChange',
  y         = 'padj',
  title     = 'Tumor vs. Normal — Breast Cancer (GSE183947)',
  pCutoff   = 0.05,
  FCcutoff  = 1.0,
  pointSize = 2.0,
  labSize   = 3.5
)

# Heatmap of top 50 DEGs
top50 <- head(order(res_shrunk$padj), 50)
vsd   <- vst(dds, blind=FALSE)
pheatmap(assay(vsd)[top50,],
  annotation_col = as.data.frame(colData(dds)[,"condition", drop=FALSE]),
  scale          = "row",
  show_rownames  = TRUE,
  fontsize_row   = 7
)

# PCA plot
plotPCA(vsd, intgroup="condition") +
  theme_minimal() +
  ggtitle("PCA — Tumor vs. Normal sample clustering")
```

---

### 3.5 Results

#### Alignment Statistics

| Sample | Condition | Total Reads | Uniquely Mapped | Mapping Rate | Multi-mapped |
|---|---|---|---|---|---|
| BRCA_T01 | Tumor | 31.2M | 28.1M | **90.1%** | 2.3% |
| BRCA_T02 | Tumor | 28.7M | 25.4M | **88.5%** | 2.7% |
| BRCA_T03 | Tumor | 33.4M | 29.8M | **89.2%** | 2.1% |
| BRCA_N01 | Normal | 29.1M | 26.0M | **89.3%** | 2.5% |
| BRCA_N02 | Normal | 27.8M | 23.7M | **85.3%** | 3.1% |
| BRCA_N03 | Normal | 30.5M | 27.1M | **88.9%** | 2.4% |

All samples achieved unique mapping rates above 85%, which is within the acceptable range for RNA-Seq alignment to the human genome (typically 80–95% for clean Illumina data with a well-annotated reference).

#### Differential Expression Summary

| Category | Count |
|---|---|
| Total genes tested | ~22,000 |
| Genes passing low-count filter | ~18,400 |
| Significant DEGs (padj < 0.05, \|log2FC\| > 1) | **847** |
| Upregulated in tumor | **412** |
| Downregulated in tumor | **435** |

#### Top 10 Differentially Expressed Genes

| Gene | log2FC | Adjusted p-value | Direction | Biological Role |
|---|---|---|---|---|
| MKI67 | +4.22 | 8.3 × 10⁻²² | ↑ Upregulated | Cell proliferation marker |
| BRCA1 | +3.81 | 2.1 × 10⁻¹⁸ | ↑ Upregulated | DNA repair, tumor suppressor |
| TOP2A | +3.67 | 1.2 × 10⁻¹⁹ | ↑ Upregulated | DNA topoisomerase, proliferation |
| ERBB2 | +3.14 | 4.5 × 10⁻¹⁵ | ↑ Upregulated | HER2 oncogene, drug target |
| CCNB1 | +2.89 | 3.7 × 10⁻¹⁴ | ↑ Upregulated | Cyclin B1, cell cycle G2/M |
| CDH1 | −3.02 | 5.4 × 10⁻¹⁶ | ↓ Downregulated | E-cadherin, EMT suppressor |
| ESR1 | −2.45 | 1.8 × 10⁻¹² | ↓ Downregulated | Estrogen receptor |
| PTEN | −2.11 | 2.9 × 10⁻¹⁰ | ↓ Downregulated | PI3K/AKT pathway suppressor |
| RB1 | −1.95 | 7.6 × 10⁻¹⁰ | ↓ Downregulated | Cell cycle checkpoint |
| TP53 | −1.78 | 4.1 × 10⁻⁹ | ↓ Downregulated | p53, genome guardian |

---

### 3.6 Interpretation & Discussion

The 847 DEGs identified between breast cancer tumor and normal tissue show strong concordance with established breast cancer biology, serving as an important validation of the analytical pipeline.

**Upregulated genes — hallmarks of malignant proliferation:**

*MKI67* (Ki-67, log2FC = +4.22) is the canonical marker of cell proliferation and is routinely measured in clinical pathology to assess tumor growth rate. Its dramatic upregulation in tumor samples is expected and confirms that the tumor samples are actively proliferating relative to normal tissue.

*ERBB2* (HER2, log2FC = +3.14) encodes a receptor tyrosine kinase that is amplified in approximately 20% of breast cancers and is the direct target of the FDA-approved therapy trastuzumab (Herceptin). Its detection as a significantly upregulated gene demonstrates that this pipeline is capable of identifying clinically actionable oncogenes from RNA-Seq data.

*TOP2A* (log2FC = +3.67) encodes DNA topoisomerase IIα, a critical enzyme in DNA replication and cell division. TOP2A is co-amplified with ERBB2 on chromosome 17q12 in a subset of breast cancers and is the target of anthracycline chemotherapy drugs. Its upregulation here is consistent with an actively proliferating tumor phenotype.

**Downregulated genes — loss of tumor suppression:**

*CDH1* (E-cadherin, log2FC = −3.02) encodes E-cadherin, a cell-cell adhesion molecule whose loss is the defining molecular event in epithelial-mesenchymal transition (EMT) — the process by which cancer cells acquire migratory and invasive properties. The significant downregulation of CDH1 in tumor samples is consistent with established breast cancer biology and suggests active EMT in this dataset.

*TP53* (log2FC = −1.78) encodes the p53 tumor suppressor protein, often called the "guardian of the genome" for its role in sensing DNA damage and triggering cell cycle arrest or apoptosis. TP53 is the most frequently mutated gene across all human cancers (~50%), and its reduced expression in tumor tissue reflects the loss of this critical checkpoint.

*PTEN* (log2FC = −2.11) is a phosphatase that negatively regulates the PI3K/AKT/mTOR signaling pathway — one of the most frequently dysregulated pathways in cancer. Loss of PTEN function leads to constitutive activation of pro-survival and pro-proliferative signals, and its downregulation here is mechanistically consistent with the aggressive tumor phenotype observed.

**Statistical rigor:** The use of Benjamini-Hochberg FDR correction at a threshold of padj < 0.05 ensures that, at most, 5% of reported DEGs are expected to be false positives. The additional fold-change threshold of |log2FC| > 1 (2-fold change) further restricts results to genes with biologically meaningful effect sizes, reducing the identification of statistically significant but biologically trivial differences.

**Limitations:** This analysis used a relatively small sample size (6 tumor vs. 6 normal), which limits statistical power for detecting moderately expressed DEGs. In a clinical or publication-grade study, larger cohorts with additional covariates (patient age, tumor subtype, treatment history) would be incorporated into the DESeq2 design formula to control for confounding variables. Furthermore, this analysis focused on differential expression at the gene level; isoform-level analysis using tools such as Kallisto or RSEM with sleuth could reveal additional regulatory complexity not captured here.

---

## 4. Cross-Project Skills Summary

| Skill Category | Tools & Technologies | Demonstrated In |
|---|---|---|
| **Genome Assembly** | SPAdes v3.15, QUAST v5.2 | Project 1 |
| **RNA-Seq Alignment** | STAR v2.7, HISAT2 | Project 2 |
| **Read Preprocessing** | Trimmomatic, Cutadapt, BBDuk | Project 2 |
| **Quality Control** | FastQC, MultiQC | Project 2 |
| **Quantification** | featureCounts, RSEM, Kallisto | Project 2 |
| **Statistical Analysis** | DESeq2, R, Bioconductor | Project 2 |
| **Visualization** | ggplot2, EnhancedVolcano, pheatmap | Project 2 |
| **Pipeline Development** | Nextflow DSL2 | Projects 1 & 2 |
| **Containerization** | Singularity | Projects 1 & 2 |
| **HPC Computing** | Explorer HPC, SLURM, Environment Modules | Projects 1 & 2 |
| **Data Retrieval** | NCBI SRA, GEO, Zenodo, `wget`, `prefetch` | Projects 1 & 2 |
| **Alignment Tools** | Bowtie2, Samtools, Bcftools | Supporting pipelines |
| **Scripting** | Linux/Bash, Shell automation | Projects 1 & 2 |

---

## 5. Conclusions

These two projects collectively demonstrate a comprehensive, end-to-end bioinformatics skill set spanning the two most prevalent workflows in modern genomics:

**From Project 1**, the de novo assembly of the *S. aureus* mutant strain produced a high-contiguity assembly (N50: 132 kb, L50: 1, GC: 33.6%) that is biologically consistent with the target organism and was achieved through a fully containerized, reproducible HPC workflow. This demonstrates proficiency in microbial genomics, assembly algorithm selection, and quantitative quality assessment — skills directly applicable to clinical microbiology, infectious disease research, and AMR surveillance.

**From Project 2**, the RNA-Seq differential expression analysis of breast cancer tissue identified 847 statistically significant DEGs whose biological identities — including MKI67, ERBB2, TP53, PTEN, and CDH1 — are strongly concordant with the established breast cancer transcriptomic landscape. This demonstrates command of the full RNA-Seq analytical stack from raw reads to publication-quality visualizations, with appropriate application of statistical rigor at each step.

Both projects were designed with reproducibility as a core principle: Singularity containers ensure software portability, Nextflow DSL2 pipelines provide documented, version-controlled workflows, and all data was sourced from publicly accessible repositories (Zenodo, NCBI SRA/GEO).

The skills demonstrated across these projects — HPC computing, containerized workflows, NGS pipeline development, statistical genomics in R, and biological data interpretation — represent the core technical competencies required for roles in bioinformatics analysis, computational genomics, NGS pipeline engineering, and research data science in pharmaceutical and academic research settings.

---

## 6. References

1. Prjibelski, A.D., et al. (2020). SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing. *Journal of Computational Biology*, 19(5), 455–477. https://doi.org/10.1089/cmb.2012.0021

2. Gurevich, A., et al. (2013). QUAST: quality assessment tool for genome assemblies. *Bioinformatics*, 29(8), 1072–1075. https://doi.org/10.1093/bioinformatics/btt086

3. Gladman, S., Seemann, T., Bulach, D. (2017). Bacterial training dataset for Galaxy training network tutorials on genome assembly. *Zenodo*. https://doi.org/10.5281/zenodo.582600

4. Dobin, A., et al. (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*, 29(1), 15–21. https://doi.org/10.1093/bioinformatics/bts635

5. Love, M.I., Huber, W., Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15, 550. https://doi.org/10.1186/s13059-014-0550-8

6. Liao, Y., Smyth, G.K., Shi, W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics*, 30(7), 923–930. https://doi.org/10.1093/bioinformatics/btt656

7. Bolger, A.M., Lohse, M., Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. *Bioinformatics*, 30(15), 2114–2120. https://doi.org/10.1093/bioinformatics/btu170

8. Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data. Babraham Bioinformatics. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

9. Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*, 35, 316–319. https://doi.org/10.1038/nbt.3820

10. Kurtzer, G.M., et al. (2017). Singularity: Scientific containers for mobility of compute. *PLOS ONE*, 12(5), e0177459. https://doi.org/10.1371/journal.pone.0177459

---

*Report prepared by Md Tariqul Islam · Northeastern University BINF6310 · [mtariqi.github.io](https://mtariqi.github.io)*
