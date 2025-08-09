MCF7 RNA-Seq Analysis — Preprocessing & Analysis

This repository contains an end-to-end workflow for RNA-seq analysis of lncRNAs results in tight lncRNA–target couplings.
This is an analysis of a gene expression profiling (RNA-seq) and differentail genes expression analyses of MCF-7 cells following estrogen stimulation.

It includes two main components: 

Preprocessing (Clean.sh)
Downloads SRA runs, builds 4 sample FASTQs, performs QC, alignment, and outputs gene transcription level.

Analysis (Differential_expression.sh)
Performs differential expression analysis using edgeR, generates PCA and volcano plots, and visualizes gene expression levels.

Data

SRA BioProject: PRJNA1055006
Gene expression profiling (RNA-seq) and differentail genes expression analyses of MCF-7 cells following estrogen stimulation (after 6 hours) either alone or precultured in transwells with CAFs (for 48 hours). MCF-7 cells alone and without estrogen stimulation is the control condition.

Sample 1: GSM7984517 (SRX22974823)
cell line: MCF7
treatment: 0hr estrogen alone

Sample 2: GSM7984518 (SRX22974824)
cell line: MCF7
treatment: 0hr estrogen alone

Sample 3: GSM7984519 (SRX22974825)
cell line: MCF7
treatment: 0hr estrogen CAF 48hr

Sample 4: GSM7984520 (SRX22974826)
cell line: MCF7
treatment: 0hr estrogen CAF 48hr

Sample 5: GSM7984521 (SRX22974827)
cell line: MCF7
treatment: 6hr estrogen alone

Sample 6: GSM7984522 (SRX22974828)
cell line: MCF7
treatment: 6hr estrogen alone

Sample 7: GSM7984523 (SRX22974829)
cell line: MCF7
treatment: 6hr estrogen CAF 48hr

Sample 8: GSM7984524 (SRX22974830)
cell line: MCF7
treatment: 6hr estrogen CAF 48hr

Reference files:

Transcript annotation: gencode.v29.annotation.gtf
Alignment index: hg38
Related publication:
PMID: 38388711

Directory Output

The preprocessing script creates a structured working directory data_pre_processing/ containing:

raw/      # downloaded and processed SRA FASTQs  
fastq/    # concatenated FASTQs (one per sample)  
aligned/  # STAR-aligned BAMs  
counts/   # gene-level count matrix  
