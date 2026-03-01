#!/usr/bin/env bash

# SMR analysis: plasma pQTL × MS severity GWAS
# Chromosome-wise analysis (chr1–22 )

for CHR in {1..22}
do
    ./smr-1.3.1 \
        --bfile data/ldref/1000G_EUR_Phase3_plink/1000G.EUR.QC.${CHR} \
        --gwas-summary data/gwas/ms_severity.ma \
        --beqtl-summary data/pqtl/plasma/build_37/plasma_smr \
        --maf 0.01 \
        --diff-freq 0.5 \
        --out results/smr/plasma_ms_severity/chr${CHR} \
        --thread-num 30
done