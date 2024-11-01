#!/bin/bash

WEIGHTS_DIR="./DATA/123"

REF_LD_PREFIX="./LDREF/1000G.EUR."

BASE_OUTPUT_DIR="./result"
mkdir -p "${BASE_OUTPUT_DIR}"

task_list_file="fusion_task_list.txt"
> $task_list_file

for SUMSTATS_FILE in ./GWAS/OCD/NEW/*.sumstats; do
    SUMSTATS_BASENAME=$(basename "${SUMSTATS_FILE}" .sumstats)
    
    OUTPUT_DIR="${BASE_OUTPUT_DIR}/${SUMSTATS_BASENAME}"
    mkdir -p "${OUTPUT_DIR}"
    
    for WEIGHTS_FILE in "${WEIGHTS_DIR}"/*.pos; do
        WEIGHTS_BASENAME=$(basename "${WEIGHTS_FILE}" .pos)
        
        for CHR in {1..22}; do
            OUT_FILE="${OUTPUT_DIR}/MAMT_${WEIGHTS_BASENAME}_chr${CHR}.dat"
            
            echo "Rscript FUSION.assoc_test.R --sumstats ${SUMSTATS_FILE} \
                                            --weights ${WEIGHTS_FILE} \
                                            --weights_dir ${WEIGHTS_DIR}/ \
                                            --chr ${CHR} \
                                            --ref_ld_chr ${REF_LD_PREFIX} \
                                            --out ${OUT_FILE}" >> $task_list_file
        done
    done
done

cat $task_list_file | parallel -j 1
