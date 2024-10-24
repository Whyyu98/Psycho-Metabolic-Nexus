#!/bin/bash

input_folder="./GWAS/"
input_files=(${input_folder}*.txt)

eqtl_folder="./Tissue/"
eqtl_prefixes=($(ls ${eqtl_folder}*.besd | sed 's/\.besd$//' | xargs -n 1 basename))

task_list_file="task_list.txt"
> $task_list_file

for input_file in "${input_files[@]}"; do
    filename=$(basename -- "${input_file}")
    filename_no_ext="${filename%.*}"

    for eqtl_prefix in "${eqtl_prefixes[@]}"; do
        output_file="./result/${filename_no_ext}_${eqtl_prefix}"

        echo "./smr-1.3.1 --bfile ./g1000_eur/g1000_eur --gwas-summary '${input_file}' --beqtl-summary '${eqtl_folder}${eqtl_prefix}' --diff-freq-prop 0.90 --out '${output_file}' --thread-num 10" >> $task_list_file
    done
done

cat $task_list_file | parallel -j 20
