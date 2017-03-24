#!/bin/bash

# Get individuals
cut -d " " -f 2 ~/repo/sa_course/randomise_data/scratch/chr16.fam > inc.txt

# Get SNPs in FTO region
cat ~/data/alspac_1kg/data/snp-stats/data_chr16.snp-stats | head -n 495687 | tail -n 104026 | cut -f 2 > keepsnps 

# extract
qctool -g ${HOME}/data/alspac_1kg/data/dosage_bgen/data_chr16.bgen -s ${HOME}/data/alspac_1kg/data/data.sample -incl-rsids keepsnps -incl-samples inc.txt -og data_chr16.gen
gzip data_chr16.gen

# Make sample file
Rscript make_sample.r

# Make vcf file
qctool -g data_chr16.gen.gz -s data.sample -og data_chr16.vcf
gzip data_chr16.vcf

