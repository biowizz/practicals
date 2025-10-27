# the data subset we will use in this course
export CHRS='chr6_and_chr17'
# GATK regions string for many gatk commands that require it
# export GATK_REGIONS='-L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM'
#export GATK_REGIONS='-L chr6 -L chr17'
export GATK_REGIONS='-L chr17' 


wget http://genomedata.org/pmbio-workshop/references/genome/$CHRS/ref_genome.tar
# unpack the archive using `tar -xvf` (`x` for extract, `v` for verbose, `f` for file)
tar -xvf ref_genome.tar

# view contents
tree

# remove the archive
rm -f ref_genome.tar

# uncompress the reference genome FASTA file
gunzip ref_genome.fa.gz


# download the files
wget http://genomedata.org/pmbio-workshop/references/transcriptome/$CHRS/ref_transcriptome.gtf
wget http://genomedata.org/pmbio-workshop/references/transcriptome/$CHRS/ref_transcriptome.fa

wget -c http://genomedata.org/pmbio-workshop/fastqs/$CHRS/RNAseq_Norm.tar
wget -c http://genomedata.org/pmbio-workshop/fastqs/$CHRS/RNAseq_Tumor.tar

