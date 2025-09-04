#!/bin/bash

# RNA-seq mini-pipeline: QC → alignment → counting

# -------------------- Setup --------------------


# Define project structure relative to current location
PROJECT_DIR="./data_pre_processing"
mkdir -p "${PROJECT_DIR}"/{raw,fastq,aligned,counts,logs,qc,STAR_index}


cd "${PROJECT_DIR}/raw"


# Group SRA run IDs by biological sample 
T0hrrep1=(SRX22974823)   # SRX22974823
T0hrrep2=(SRX22974824)   # SRX22974824
T0hrCAFrep1=(SRX22974825)    # SRX22974825
T0hrCAFrep2=(SRX22974826)    # SRX22974826
T6hrrep1=(SRX22974827)   # SRX22974827
T6hrrep2=(SRX22974828)   # SRX22974828
T6hrCAFrep1=(SRX22974829)   # SRX22974829
T6hrCAFrep2=(SRX22974820)   # SRX22974830

# -------------------- Download & Convert --------------------

# Download .sra files
for r in "${T0hrrep1[@]}" "${T0hrrep2[@]}" "${T0hrCAFrep1[@]}" "${T0hrCAFrep1[@]}" "${T6hrrep1[@]}" "${T6hrrep2[@]}" "${T6hrCAFrep1[@]}" "${T6hrCAFrep1[@]}"; do
  prefetch "$r"
done

# Convert to gzipped FASTQ

for r in "${T0hrrep1[@]}" "${T0hrrep2[@]}" "${T0hrCAFrep1[@]}" "${T0hrCAFrep2[@]}" "${T6hrrep1[@]}" "${T6hrrep2[@]}" "${T6hrCAFrep1[@]}" "${T6hrCAFrep1[@]}"; do
  fasterq-dump -e 16 -p -O . "$r"
  gzip -f "${r}.fastq"
done

# Concatenate per-sample FASTQs
cat "${T0hrrep1[@]/%/.fastq.gz}"  > T0hrrep1.fastq.gz
cat "${T0hrrep2[@]/%/.fastq.gz}"  > T0hrrep2.fastq.gz
cat "${T0hrCAFrep1[@]/%/.fastq.gz}" > T0hrCAFrep1.fastq.gz
cat "${T0hrCAFrep2[@]/%/.fastq.gz}" > T0hrCAFrep2.fastq.gz
cat "${T6hrrep1[@]/%/.fastq.gz}"  > T0hrrep1.fastq.gz
cat "${T6hrrep2[@]/%/.fastq.gz}"  > T0hrrep2.fastq.gz
cat "${T6hrCAFrep1[@]/%/.fastq.gz}" > T6hrCAFrep1.fastq.gz
cat "${T6hrCAFrep2[@]/%/.fastq.gz}" > T6hrCAFrep2.fastq.gz

# Move to fastq/ folder
mv T*.fastq.gz ../fastq/

# -------------------- QC --------------------

cd ../fastq
fastqc T0hrrep1.fastq.gzv T0hrrep2.fastq.gz T0hrCAFrep1.fastq.gz T0hrCAFrep2.fastq.gz T6hrrep1.fastq.gzv T6hrrep2.fastq.gz T6hrCAFrep1.fastq.gz T6hrCAFrep2.fastq.gz \
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
      --readFilesIn T0hrrep1.fastq.gzv T0hrrep2.fastq.gz T0hrCAFrep1.fastq.gz T0hrCAFrep2.fastq.gz T6hrrep1.fastq.gzv T6hrrep2.fastq.gz T6hrCAFrep1.fastq.gz T6hrCAFrep2.fastq.gz  \
      --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode GeneCounts \
      --outFileNamePrefix alignments/

# -------------------- Quantification (featureCounts) --------------------

cd ..
curl -L -o gencode.v29.annotation.gtf.gz \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_94/gencode.v29.annotation.gtf.gz
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
