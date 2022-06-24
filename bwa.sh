#!/bin/bash

module load bwa
module load samtools

bwa mem -x ont2d -t 8 /data/tunglab/shared/genomes/panubis1_bwa_preNCBI/Panubis_1.0.fa \
fastq/BARCODE_trimmed.fastq | samtools view -b -F 4 | samtools sort -o aligned_sam/aligned_BARCODE.bam -T BARCODE.tmp -

samtools index aligned_sam/aligned_BARCODE.bam
