#!/bin/bash

GWAS_DIR="./GWAS"
RESULTS_DIR="./result"
GCTA="./gcta"
BFILE="./g1000_eur/g1000_eur"
GENE_LIST="hg19.txt"
THREADS=10

mkdir -p "$RESULTS_DIR"

task_list_file="gcta_task_list.txt"
> $task_list_file

for gwas_file in "$GWAS_DIR"/*.txt; do
  file_prefix=$(basename "$gwas_file" .txt)
  
  output_file="$RESULTS_DIR/$file_prefix"

  echo "$GCTA --bfile $BFILE --mBAT-combo $gwas_file --mBAT-gene-list $GENE_LIST --out $output_file --thread-num $THREADS" >> $task_list_file
done

cat $task_list_file | parallel -j 10
