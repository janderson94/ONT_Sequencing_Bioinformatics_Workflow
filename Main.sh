#!/bin/#!/usr/bin/env bash


#This is the main workflow
#So we first start with the file location and need to split each fast5 into its own directory

#If anything fails lets break from main
set -e




###############################################################
#              PRE-WORKFLOW PREPARATION
###############################################################

#We have a number of directories we want to make sure exist in this folder.
#I have included in this directory so you can easily copy a main directory into a new name, drop the fast5 files in and go.
mkdir -p fast5
mkdir -p fastq
mkdir -p genotypes
mkdir -p aligned_sam

#######################################################################################
# The fast5 files from your run should be directly placed into the fast5 directory.
# e.g. there will be a fast5/SEQUENCING_RUN_0.FAST5, fast5/SEQUENCING_RUN_1.FAST5, etc.
#######################################################################################

#You can rename the run as desired. You can easily modify parts of this script to pull and store desired run info
#E.g. mean Qscore?

run='flongle_run_INFO'
touch run_summary_${run}.txt
touch commands.sh



#Other prep needed before beginning:
# barcodes.txt (you need a file that is a single column of the barcodes added (e.g. barcode01, barcode13, etc...))
# array.sh (a typical array file, an example included here. YOU WILL NEED TO CHANGE THE DIRECTORY)

# There is a template here which can be pointed to the current directory
# Feel free to make this array script fancier if desired.
tmp=`pwd`
sed -i "s|CURRENT_DIR|$tmp|g" array.sh

# Depending on your purposes, you may need to change the genome
# I will add support for this soon. Everything included uses the panubis1_preNCBI genome










###############################################################
#             BASECALLING AND DEMULTIPLEXING
###############################################################

#For each fast5 we need to generate its own directory for efficient basecalling
for f in `ls fast5/*.fast5 | sed -e 's/fast5\///g' `
do
  mkdir fast5/$f"_folder"
  mv fast5/$f fast5/$f"_folder"
done

echo 'Reads subsectioned into files'


#Great so those are sorted now.
#Now we need to parallelize basecalling. This is by far the bottleneck in this workflow.
for f in `ls fast5/ | sed -e 's/\_folder//g'`
do
  cat guppy.sh | sed -e 's/RUN/'$run'/g' | sed -e 's/FILENAME/'$f'/g' > guppy.$f.sh
  echo "sh "guppy.$f.sh >> commands.sh
done

l=`wc -l commands.sh | awk '{print $1}'`

echo 'Basecalling in parallel'

#Here we will run 50 jobs at 8GB, this should be enough if you are basecalling for 1000 reads at a time.
#You can easily increase the number of jobs based on the amount of data generated and current HARDAC usage.
#This will also demultiplex.
sbatch --mem=8G -a 1-$l%50 --wait --job-name='basecall' array.sh

echo 'Basecalling finished'

#Select the barcodes of interest -- e.g. take barcode 2 from each of the filenames above and concatenated them into one fastq file

for x in `cat barcodes.txt`
do
  for y in `ls -d fastq/*`
  do
    if [ -d "$y/$x" ]
    then
      cat $y/$x/* >> fastq/$x.fastq
    fi
  done
done

tmp3=0
for f in `cat barcodes.txt`
do
  tmp=`wc -l fastq/$f.fastq | awk '{print $1}'`
  tmp2=`expr $tmp / 4`
  echo $f $tmp2 "reads demultiplexed" >> run_summary_${run}.txt
  tmp3=`expr $tmp3 + $tmp2`
done

echo $tmp3 "reads demultiplexed in total" >> run_summary_${run}.txt






###############################################################
#               ADAPTOR TRIMMING
###############################################################

#Trim adaptor from demultiplexed reads
touch commands.sh; rm commands.sh

#Porechop is increasingly deprecated. This means that workflows in the future will need to work around this with cutadapt, TrimGalore, or trimmomatic.
for f in `cat barcodes.txt`
do
  cat porechop.sh | sed -e 's/BARCODE/'$f'/g' > porechop.$f.sh
  echo "sh "porechop.$f.sh >> commands.sh
done

l=`wc -l commands.sh | awk '{print $1}'`

echo "Trimming adaptor sequences"
#Because you can have a max of 96 samples multiplexed, this at most would be ~500GB. If you are combining multiple runs, alter ${l} accordingly.
sbatch --mem=5G -a 1-$l%$l --wait --job-name='trim' array.sh
echo "Reads trimmed successfully"







###############################################################
#               MAPPING, SORTING, AND INDEXING
###############################################################

#Map using BWA, remove unmapped reads, sort, and index
touch commands.sh
rm commands.sh

for f in `cat barcodes.txt`
do
  cat bwa.sh | sed -e 's/BARCODE/'$f'/g' > bwa.$f.sh
  echo "sh "bwa.$f.sh >> commands.sh
done

l=`wc -l commands.sh | awk '{print $1}'`

echo "Mapping reads"
sbatch --mem=10G -a 1-$l%$l --wait --job-name='map' array.sh
echo "Reads mapped"
