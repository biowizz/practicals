# https://rnabio.org/

# Reference fasta file
wget http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa

# Annotation file
wget http://genomedata.org/rnaseq-tutorial/annotations/GRCh38/chr22_with_ERCC92.gtf

# Adapter seq
wget http://genomedata.org/rnaseq-tutorial/illumina_multiplex.fa

# Creating index
cd $RNA_REFS_DIR
hisat2_extract_splice_sites.py $RNA_REF_GTF > $RNA_REFS_DIR/splicesites.tsv
hisat2_extract_exons.py $RNA_REF_GTF > $RNA_REFS_DIR/exons.tsv
hisat2-build -p 4 --ss $RNA_REFS_DIR/splicesites.tsv --exon $RNA_REFS_DIR/exons.tsv $RNA_REF_FASTA $RNA_REF_INDEX
ls
# On my comp
hisat2-build -p 4 --ss splicesites.tsv --exon exons.tsv chr22_with_ERCC92.fa chr22

# Raw data
wget http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar

# Trimming
export RNA_HOME=/data/results/sandeep/practical/rna
export RNA_REFS_DIR=/data/results/sandeep/practical/rna/ref
export RNA_DATA_DIR=/data/results/sandeep/practical/rna/rawdata
export RNA_DATA_TRIM_DIR=/data/results/sandeep/practical/rna/trim

cd cd $RNA_HOME
export S1=UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22
fastp -i $RNA_DATA_DIR/$S1.read1.fastq.gz -I $RNA_DATA_DIR/$S1.read2.fastq.gz -o $RNA_DATA_TRIM_DIR/$S1.read1.fastq.gz -O $RNA_DATA_TRIM_DIR/$S1.read2.fastq.gz -l 25 --adapter_fasta $RNA_REFS_DIR/illumina_multiplex.fa --trim_front1 13 --trim_front2 13 --json $RNA_DATA_TRIM_DIR/$S1.fastp.json --html $RNA_DATA_TRIM_DIR/$S1.fastp.html 2>$RNA_DATA_TRIM_DIR/$S1.fastp.log


export S2=UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22
fastp -i $RNA_DATA_DIR/$S2.read1.fastq.gz -I $RNA_DATA_DIR/$S2.read2.fastq.gz -o $RNA_DATA_TRIM_DIR/$S2.read1.fastq.gz -O $RNA_DATA_TRIM_DIR/$S2.read2.fastq.gz -l 25 --adapter_fasta $RNA_REFS_DIR/illumina_multiplex.fa --trim_front1 13 --trim_front2 13 --json $RNA_DATA_TRIM_DIR/$S2.fastp.json --html $RNA_DATA_TRIM_DIR/$S2.fastp.html 2>$RNA_DATA_TRIM_DIR/$S2.fastp.log

export S3=UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22
fastp -i $RNA_DATA_DIR/$S3.read1.fastq.gz -I $RNA_DATA_DIR/$S3.read2.fastq.gz -o $RNA_DATA_TRIM_DIR/$S3.read1.fastq.gz -O $RNA_DATA_TRIM_DIR/$S3.read2.fastq.gz -l 25 --adapter_fasta $RNA_REFS_DIR/illumina_multiplex.fa --trim_front1 13 --trim_front2 13 --json $RNA_DATA_TRIM_DIR/$S3.fastp.json --html $RNA_DATA_TRIM_DIR/$S3.fastp.html 2>$RNA_DATA_TRIM_DIR/$S3.fastp.log

export S4=HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22
fastp -i $RNA_DATA_DIR/$S4.read1.fastq.gz -I $RNA_DATA_DIR/$S4.read2.fastq.gz -o $RNA_DATA_TRIM_DIR/$S4.read1.fastq.gz -O $RNA_DATA_TRIM_DIR/$S4.read2.fastq.gz -l 25 --adapter_fasta $RNA_REFS_DIR/illumina_multiplex.fa --trim_front1 13 --trim_front2 13 --json $RNA_DATA_TRIM_DIR/$S4.fastp.json --html $RNA_DATA_TRIM_DIR/$S4.fastp.html 2>$RNA_DATA_TRIM_DIR/$S4.fastp.log

export S5=HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22
fastp -i $RNA_DATA_DIR/$S5.read1.fastq.gz -I $RNA_DATA_DIR/$S5.read2.fastq.gz -o $RNA_DATA_TRIM_DIR/$S5.read1.fastq.gz -O $RNA_DATA_TRIM_DIR/$S5.read2.fastq.gz -l 25 --adapter_fasta $RNA_REFS_DIR/illumina_multiplex.fa --trim_front1 13 --trim_front2 13 --json $RNA_DATA_TRIM_DIR/$S5.fastp.json --html $RNA_DATA_TRIM_DIR/$S5.fastp.html 2>$RNA_DATA_TRIM_DIR/$S5.fastp.log

export S6=HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22
fastp -i $RNA_DATA_DIR/$S6.read1.fastq.gz -I $RNA_DATA_DIR/$S6.read2.fastq.gz -o $RNA_DATA_TRIM_DIR/$S6.read1.fastq.gz -O $RNA_DATA_TRIM_DIR/$S6.read2.fastq.gz -l 25 --adapter_fasta $RNA_REFS_DIR/illumina_multiplex.fa --trim_front1 13 --trim_front2 13 --json $RNA_DATA_TRIM_DIR/$S6.fastp.json --html $RNA_DATA_TRIM_DIR/$S6.fastp.html 2>$RNA_DATA_TRIM_DIR/$S6.fastp.log
