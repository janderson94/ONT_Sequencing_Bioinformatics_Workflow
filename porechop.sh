#!/bin/bash

module load python/3.4.1-fasrc01
module load gcc

../Porechop/porechop-runner.py -i fastq/BARCODE.fastq -o fastq/BARCODE_trimmed.fastq
