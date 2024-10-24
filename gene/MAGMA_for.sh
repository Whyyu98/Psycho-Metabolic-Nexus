#!/bin/bash

pval_folder="./PD-MeTs"
output_folder="./result"
bfile="./g1000_eur/g1000_eur"
gene_annot="./new.genes.annot"
N=500000

task_list_file="magma_task_list.txt"
> $task_list_file

for pval_file in "$pval_folder"/*.txt; do
  file_prefix=$(basename "$pval_file" .txt)
  
  subfolder="$output_folder/$file_prefix"
  
  mkdir -p "$subfolder"
  
  output_file="$subfolder/$file_prefix"

  echo "./magma --bfile '$bfile' --pval '$pval_file' N='$N' --gene-annot '$gene_annot' --out '$output_file'" >> $task_list_file
done

cat $task_list_file | parallel -j 10
