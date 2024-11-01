#!/bin/bash

FUMA_DIR="./result"

OUT_DIR="./POPS_result"

mkdir -p "$OUT_DIR"

task_list_file="pops_task_list.txt"
> $task_list_file

for SUBDIR in "$FUMA_DIR"/*; do
  if [ -d "$SUBDIR" ]; then
    SUBDIR_NAME=$(basename "$SUBDIR")
    OUT_PREFIX="$OUT_DIR/$SUBDIR_NAME"
    echo "python pops.py --gene_annot_path ./gene_all.txt \
                       --feature_mat_prefix ./features_munged/pops_features \
                       --num_feature_chunks 2 \
                       --magma_prefix '$SUBDIR/$SUBDIR_NAME' \
                       --control_features_path ./features_jul17_control.txt \
                       --out_prefix '$OUT_PREFIX'" >> $task_list_file
  fi
done

cat $task_list_file | parallel -j 5
