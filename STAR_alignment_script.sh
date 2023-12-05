#!/bin/bash
#$ -l m_mem_free=8G
#$ -cwd
#$ -o STAR_Test.out
#$ -j y
#$ -N STAR_Testing
#$ -pe threads 4

####################
### IN CASE#  YOU WANT TO MAKE YOUR OWN INDEX FILE:
# # e.g. Human Index from https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=9606 > GRCh38.p14
# 
# # load modules
# module load EBModules-LegacyBNB
# module load Anaconda2/5.3.0
# 
# # Make a conda environment to download ncbi data
# conda create -n ncbi_datasets
# conda activate ncbi_datasets
# conda install -c conda-forge ncbi-datasets-cli
# datasets download genome accession GCF_000001405.40 --include gff3,rna,cds,protein,genome,seq-report
# unzip ncbi_dataset.zip
# rm -r ncbi_dataset.zip
# 
# ##load software modules to generate the index
# module load EBModules
# module load EBModules-LegacyBNB
# module load STAR/2.7.10a-GCC-10.3.0
# 
# STAR --runThreadN 4 \
#  --runMode genomeGenerate \
#  --genomeDir /grid/genomicscourse/home/shared/genomes/STAR_hg38_index \
#  --genomeFastaFiles /grid/genomicscourse/home/vanterve/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna \
#  --sjdbGTFfile /grid/genomicscourse/home/vanterve/ncbi_dataset/data/GCF_000001405.40/genomic.gtf \
#  --sjdbOverhang 100 #specifies the length of the  sequence around the exon-intron junctions that is used in constructing the splice junctions database. 
# # The value should be set to one less than the read length.

####################

##load software modules for FastQC / STAR
module load EBModules
module load EBModules-LegacyBNB
module load FastQC/0.11.8-Java-1.8
module load Trimmomatic/0.39-Java-11
module load STAR/2.7.10a-GCC-10.3.0
module load Subread/2.0.2-GCC-10.2.0

#create working directory
wkdir=/grid/genomicscourse/home/vanterve/STAR

#mkdir -p $wkdir

#dir with fastq files
fastq=/grid/genomicscourse/home/shared/RNA_seq_tutorial/FASTQ

# dir for STAR output
#mkdir $wkdir/STARoutput

cd $fastq

#This loop will perform the pipeline of QC, trimming, aligning and then annotating reads
for r1 in `ls *R1*.fastq.gz` 
do 
	name=${r1/_L004_R1_001.fastq.gz/} 
 	r2=${r1/R1/R2}
 	echo $r1
 	echo $name 
 	echo $r2
# 	
# 
# #QC pre trimming
# 	mkdir -p $wkdir/fastqc
# 	fastqc -o $wkdir/fastqc $fastq/$r1
# 	fastqc -o $wkdir/fastqc $fastq/$r2 
# 
	#mkdir -p $wkdir/trimmedFASTQ

#Trimming
	java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE  -threads 4 $fastq/$r1 $fastq/$r2 $wkdir/trimmedFASTQ/$name.PE.R1.fastq.gz $wkdir/trimmedFASTQ/$name.SE.R1.fastq.gz $wkdir/trimmedFASTQ/$name.PE.R2.fastq.gz $wkdir/trimmedFASTQ/$name.SE.R2.fastq.gz ILLUMINACLIP:/grid/genomicscourse/home/shared/bin/all_adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75 
	#LEADING: Removes leading low-quality bases (below quality 3).
	#TRAILING:Removes trailing low-quality bases (below quality 3).
	#SLIDINGWINDOW:4:15: Performs sliding window trimming, cutting once the average quality within the window (4 bases) falls below a threshold (15).
	#MINLEN: Drops any reads that fall below 75 bases long after trimming, lowering minlen will lower the standards
 
#QC on trimmed reads 
	fastqc -o $wkdir/fastqc $wkdir/trimmedFASTQ/$name.PE.R1.fastq.gz
	fastqc -o $wkdir/fastqc $wkdir/trimmedFASTQ/$name.PE.R2.fastq.gz
 
#Mapping the reads 
	mkdir -p $wkdir/STARoutput
 	
# running STAR
STAR --runThreadN 4 \
--genomeDir /grid/genomicscourse/home/shared/genomes/STAR_hg38_index \
--readFilesIn $wkdir/trimmedFASTQ/$name.PE.R1.fastq.gz $wkdir/trimmedFASTQ/$name.PE.R2.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix $wkdir/STARoutput/$name

# convert .sam to .bam
cd $STARoutput
samtools view -bS *.sam -o *.bam
samtools index *.bam

###Feature Counts
# install conda environment for featurecounts (only once)
#conda create -c bioconda -n featurecounts subread
featureCounts -t exon -g gene_id -a /grid/genomicscourse/home/shared/genomes/STAR_hg38_index/gencode.v30.GRCh38.ERCC.genes.collapsed_only.gtf -o feature_counts_after_STAR.txt *.bam

done