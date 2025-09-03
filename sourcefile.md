# To download the data:

# Create environment
conda create -n Cancer-associated_fibroblasts -c bioconda -c conda-forge \
  sra-tools fastqc multiqc hisat2 samtools trimmomatic subread -y
conda activate Cancer-associated_fibroblasts

# Folder setup
mkdir -p ~/0_Cancer-associated_fibroblasts/{data,fastq,trimmed,aligned,counts,logs,qc}
cd ~/0_Cancer-associated_fibroblasts/data

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos8/sra-pub-zq-818/SRR028/28555/SRR28555713/SRR28555713.lite.1 # SRX22974823: RNA-Seq 0hr alone rep1

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos8/sra-pub-zq-818/SRR028/28555/SRR28555712/SRR28555712.lite.1 # SRX22974824: RNA-Seq 0hr alone rep2

  wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos8/sra-pub-zq-818/SRR028/28555/SRR28555699/SRR28555699.lite.1 # SRX22974825: RNA-Seq 0hr CAF 48hr rep1

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos8/sra-pub-zq-818/SRR028/28555/SRR28555697/SRR28555697.lite.1 # SRX22974826: RNA-Seq 0hr CAF 48hr rep2

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos8/sra-pub-zq-818/SRR028/28555/SRR28555684/SRR28555684.lite.1 # SRX22974827: RNA-Seq 6hr alone rep1

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos8/sra-pub-zq-818/SRR028/28555/SRR28555683/SRR28555683.lite.1 # SRX22974828: RNA-Seq 6hr alone rep2

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos8/sra-pub-zq-818/SRR028/28555/SRR28555757/SRR28555757.lite.1 # SRX22974829: RNA-Seq 6hr CAF 48hr rep1

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos8/sra-pub-zq-818/SRR028/28555/SRR28555756/SRR28555756.lite.1 # SRX22974830: RNA-Seq 6hr CAF 48hr rep2

for r in "${SRR[@]}"; do
  fasterq-dump -e 16 -p -O . "$r"
  gzip -f "${r}.fastq"
done
