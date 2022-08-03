#!/bin/bash

module load samtools
module load bcftools

#Double check genome path
samtools mpileup -uf /data/tunglab/shared/genomes/panubis1_bwa_preNCBI/Panubis1.0.fa aligned_sam/aligned_BARCODE.bam | bcftools call -V indels -m > genotypes/BARCODE.vcf
