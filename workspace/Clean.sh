#!/bin/bash

# RNA-seq mini-pipeline: QC → alignment → counting

# -------------------- Setup --------------------


# Define project structure relative to current location
PROJECT_DIR="./data_pre_processing"
mkdir -p "${PROJECT_DIR}"/{raw,fastq,aligned,counts,logs,qc,STAR_index}


cd "${PROJECT_DIR}/raw"


# Group SRA run IDs by biological sample 
ECC1_L1=(SRX24155225)   # SRX24155225
ECC1_L2=(SRX24155226)   # SRX24155226
H460_L1=(SRX24155239)    # SRX24155239	
H460_L2=(SRX24155240)    # SRX24155240
PC-3_L1=(SRX24155254)   # SRX24155254	
PC-3_L2=(SRX24155255)   # SRX24155255
PC_3_L1_siD=(SRX24155181)   # SRX24155181
PC_3_L2_siD=(SRX24155182)   # SRX24155182

# -------------------- Download & Convert --------------------

# Download .sra files
for r in "${ECC1_L1[@]}" "${ECC1_L1[@]}" "${H460_L1[@]}" "${H460_L1[@]}" "${PC-3_L1[@]}" "${PC-3_L2[@]}" "$PC_3_L2_siD[@]}" "${PC_3_L2_siD[@]}"; do
  prefetch "$r"
done

# Convert to gzipped FASTQ

for r in "${ECC1_L1[@]}" "${ECC1_L1[@]}" "${H460_L1[@]}" "${H460_L1[@]}" "${PC-3_L1[@]}" "${PC-3_L2[@]}" "$PC_3_L2_siD[@]}" "${PC_3_L2_siD[@]}"; do
  fasterq-dump -e 16 -p -O . "$r"
  gzip -f "${r}.fastq"
done

# Concatenate per-sample FASTQs
cat "${ECC1_L1[@]/%/.fastq.gz}"  > ECC1_L1.fastq.gz
cat "${ECC1_L2[@]/%/.fastq.gz}"  > ECC1_L2.fastq.gz
cat "${H460_L1[@]/%/.fastq.gz}" > H460_L1.fastq.gz
cat "${H460_L2[@]/%/.fastq.gz}" > H460_L2.fastq.gz
cat "${PC-3_L1[@]/%/.fastq.gz}"  > PC-3_L1.fastq.gz
cat "${PC-3_L2[@]/%/.fastq.gz}"  > PC-3_L2.fastq.gz
cat "${PC_3_L1_siD[@]/%/.fastq.gz}" > PC_3_L1_siD.fastq.gz
cat "${PC_3_L2_siD[@]/%/.fastq.gz}" > PC_3_L2_siD.fastq.gz

# Move to fastq/ folder
mv ECC1*.fastq.gz H460*.fastq.gz PC_3*.fastq.gz ../fastq/

# -------------------- QC --------------------

cd ../fastq
fastqc ECC1_L1.fastq.gz ECC1_L2.fastq.gz H460_L1.fastq.gz H460_L2.fastq.gz PC-3_L1.fastq.gz PC-3_L2.fastq.gz PC_3_L1_siD.fastq.gz PC_3_L2_siD.fastq.gz \
  -o ../qc --threads 16

# -------------------- Alignment (STAR) --------------------

cd ../STAR_index
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz
unzip GRCh38.primary_assembly.genome.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
unzip gencode.v29.annotation.gtf.gz

module load star
GENOMEDIR="./RNAseq/genome/"
mdkir -p $GENOMEDIR/STAR
STAR --runThreadN 23 --runMode genomeGenerate --genomeDir $GENOMEDIR/STAR --genomeFastaFiles $GENOMEDIR/GRCh38.primary_assembly.genome.fa --sjdbGTFfile $GENOMEDIR/gencode.v29.primary_assembly.annotation.gtf
cd ../trimmed
STAR --genomeDir indexes/chr10 \
      --readFilesIn ECC1_L1.fastq.gz ECC1_L2.fastq.gz H460_L1.fastq.gz H460_L2.fastq.gz PC-3_L1.fastq.gz PC-3_L2.fastq.gz PC_3_L1_siD.fastq.gz PC_3_L2_siD.fastq.gz  \
      --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode GeneCounts \
      --outFileNamePrefix alignments/

# -------------------- Quantification (featureCounts) --------------------

cd ..
curl -L -o gencode.v16.annotation.gtf.gz \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v16.annotation.gtf.gz
gunzip -f gencode.v16.annotation.gtf.gz

featureCounts -T 16 -t exon -g gene_name \
  -a gencode.v16.annotation.gtf \
  -o counts/raw_counts_gene_sym.txt aligned/*.bam \
  &> logs/featureCounts_gene_sym.log

# Format counts matrix
{ printf "GeneSymbol\t"; head -n 2 counts/raw_counts_gene_sym.txt | tail -n 1 | cut -f7-; } > counts/final_counts_symbols.tsv
tail -n +3 counts/raw_counts_gene_sym.txt | \
  awk -v OFS="\t" '{ out=$1; for(i=7;i<=NF;i++) out=out OFS $i; print out }' >> Processed_counts/final_counts_symbols.tsv

sed -i '' '1 s|aligned/||g; 1 s|\.bam||g' counts/final_counts_symbols.tsv

# Done
echo "Pipeline complete. Output saved in: ${PROJECT_DIR}/counts/final_counts_symbols.tsv"
