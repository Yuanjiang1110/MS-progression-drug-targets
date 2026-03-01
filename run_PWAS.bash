#!/usr/bin/env bash

# PWAS association testing for MS severity (plasma proteome)
# Run chromosome 1–22 sequentially

for CHR in {1..22}
do
    Rscript scripts/PWAS.assoc_test.R \
        --sumstats data/sumstats/MS_severity_forpwas.txt \
        --weights data/pwas_weights/Plasma_Protein_EA_hg19.pos \
        --weights_dir data/pwas_weights/Plasma_Protein_weights_EA/ \
        --ref_ld_chr data/ldref/EUR/chr \
        --force_model enet \
        --chr ${CHR} \
        --out results/pwas/plasma/chr${CHR}.out
done
