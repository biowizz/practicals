# https://datacarpentry.github.io/wrangling-genomics/03-trimming.html

curl -O http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_1.fastq.gz
curl -O http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_2.fastq.gz
curl -O http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_1.fastq.gz
curl -O http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_2.fastq.gz
curl -O http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_1.fastq.gz
curl -O http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_2.fastq.gz

curl -L -o ecoli_rel606.fasta.gz http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz
gunzip ecoli_rel606.fasta.gz

## Concatenated processed files because of the long processing times
curl -L -o sub.tar.gz https://ndownloader.figshare.com/files/14418248
tar xvf sub.tar.gz

fastq=/data/results/sandeep/practical/sequencing/rawdata
flexbar=/data/results/sandeep/practical/sequencing/flexbar
adapter=/data/results/sandeep/sanjiban/illumina_multiplex.fa
hisat=/data/reference/rnaseq/gencode.v32
aligned=/data/results/sandeep/practical/sequencing/aligned
ref=/data/reference/rnaseq


##### Looping
cd "$fastq"
for infile in *_1.fastq.gz
do
base=$(basename ${infile} _1.fastq.gz)
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters "$adapter" --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 28 --zip-output GZ --reads ${infile} --reads2 ${base}_2.fastq.gz --target "$flexbar"/${base}
done


cd "$fastq"
for infile in *_1.fastq.gz
do
base=$(basename ${infile} _1.fastq.gz)
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters "$adapter" --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 28 --zip-output GZ --reads ${infile} --reads2 ${base}_2.fastq.gz --target "$flexbar"/${base}
done