#!/usr/bin/env bash

# FUSION TWAS: MS severity × whole blood expression
# Chromosome-wise analysis (chr1–22 )

for CHR in {1..22}
do
    Rscript fusion_twas/FUSION.assoc_test.R \
        --sumstats data/sumstats/MS_severity_forpwas.txt \
        --weights data/fusion/whole_blood/YFS.BLOOD.RNAARR.pos \
        --weights_dir data/fusion/whole_blood/ \
        --ref_ld_chr data/ldref/1000G.EUR. \
        --force_model enet \
        --chr ${CHR} \
        --out results/twas_fusion/whole_blood/chr${CHR}.out
done