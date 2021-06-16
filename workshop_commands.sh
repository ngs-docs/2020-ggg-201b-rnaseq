# The goal is to perform Differential Gene Expression analysis
# We will run through the workflow starting from raw files in FASTQ format all the way to generating gene list
# The input files are yeast RNAseq dataset from Schurch et al, 2016 study
# We will demonstrate the different workflow steps using a single file for concept
# Inlcuded in this binder are also commands to run all 6 files together as well as the Snakefile to run all the workflow at once

# Create conda environment from YAML file
conda env create -f rnaseq-env.yml

# Initialize conda
conda init

# Update for changes to take place
source .bashrc
## Note you could also restart the terminal for the same effect

# Activate conda environment
conda activate rnaseq-qc-map

# Create data folder
mkdir rnaseq/raw_data

# Download data file
curl -L https://osf.io/5daup/download -o rnaseq/raw_data/ERR458493.fq.gz

# Create fastqc folder for raw data
mkdir rnaseq/raw_data/fastqc

# QC on FASTQ file
fastqc rnaseq/raw_data/ERR458493.fq.gz --outdir rnaseq/raw_data/fastqc/

# Download trimmomatic adapter file
curl -L https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq2-SE.fa -o TruSeq2-SE.fa

# Create quality data folder
mkdir rnaseq/quality

# Quality trim to remove adaptor from raw FASTQ file usin trimmomatic
trimmomatic SE rnaseq/raw_data/ERR458493.fq.gz rnaseq/quality/ERR458493.qc.fq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:0:15 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2 MINLEN:25

# Create fastqc folder for trimmed files
mkdir rnaseq/quality/fastqc

# QC on trimmed file
fastqc rnaseq/quality/ERR458493.qc.fq.gz --outdir rnaseq/quality/fastqc/

# Download and index the yeast transcriptome (not run)
# curl -L ftp://ftp.ensembl.org/pub/release-99/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz -o yeast_ref/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz

# Generate salmon index
salmon index --index rnaseq/quant/sc_ensembl_index --transcripts yeast_ref/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz # --type quasi

# Quantify reads with salmon
salmon quant -i rnaseq/quant/sc_ensembl_index --libType A -r rnaseq/quality/ERR458493.qc.fq.gz -o rnaseq/quant/ERR458493_quant --seqBias --gcBias --validateMappings
