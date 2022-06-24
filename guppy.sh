#!/bin/bash

module load python

mkdir fastq/FILENAME

../guppy/ont-guppy-cpu/bin/guppy_basecaller \
--input_path fast5/FILENAME_folder \
--save_path fastq/FILENAME \
--flowcell FLO-FLG001 \
--kit SQK-LSK109 \
--cpu_threads_per_caller 4 \
--barcode_kits "EXP-PBC096"
