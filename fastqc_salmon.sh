#!/bin/bash
#$ -l m_mem_free=4G
#$ -cwd
#$ -o tutorial.out
#$ -j y
#$ -N salmonTutorial
#$ -pe threads 4

##load software modules that are available 
module load EBModules
module load EBModules-LegacyBNB
module load FastQC/0.11.8-Java-1.8
module load Trimmomatic/0.36-Java-1.8

#make var for local software
salmon=/grid/genomicscourse/home/shared/bin/salmon-latest_linux_x86_64/salmon

#make var for working directory to make things easier
wkdir=/grid/genomicscourse/home/vanterve
mkdir -p $wkdir

#same for dir with fastq files
fastq=/grid/genomicscourse/home/shared/RNA_seq_tutorial/FASTQ

cd $fastq

for r1 in `ls *R1*.fastq.gz`
do
name=${r1/_L004_R1_001.fastq.gz/}
r2=${r1/R1/R2}
echo $name
echo $r2
mkdir -p $wkdir/fastqc
fastqc -o $wkdir/fastqc $fastq/$r1
fastqc -o $wkdir/fastqc $fastq/$r2 

mkdir -p $wkdir/trimmedFASTQ

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar \
PE  -threads 4 \
 $fastq/$r1 \
 $fastq/$r2 \
 $wkdir/trimmedFASTQ/$
 *1.fastq.gz \
 $wkdir/trimmedFASTQ/$name.SE.R1.fastq.gz \
 $wkdir/trimmedFASTQ/$name.PE.R2.fastq.gz \
 $wkdir/trimmedFASTQ/$name.SE.R2.fastq.gz \
 ILLUMINACLIP:/grid/genomicscourse/home/shared/bin/all_adapters.fa:2:30:10 \
 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75
 
 fastqc -o $wkdir/fastqc $wkdir/trimmedFASTQ/$name.PE.R1.fastq.gz
 fastqc -o $wkdir/fastqc $wkdir/trimmedFASTQ/$name.PE.R2.fastq.gz
 
 mkdir -p $wkdir/quant
 index=/grid/genomicscourse/home/shared/genomes/salmon_hg38_index/Homo_sapiens.GRCh38.cdna.all_index
 $salmon quant -i $index -l A \
  -1 $wkdir/trimmedFASTQ/$name.PE.R1.fastq.gz \
  -2 $wkdir/trimmedFASTQ/$name.PE.R2.fastq.gz \
  -p 4
  --validateMappings -o $wkdir/quant/${name}_quant
 done
 