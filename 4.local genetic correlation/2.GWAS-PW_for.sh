#!/bin/bash

input_dir="DATA"
output_dir="PD-MeTs"

mkdir -p "$output_dir"

for input_file in "$input_dir"/*.txt; do
    filename=$(basename "$input_file" .txt)
    
    output_file="$output_dir/$filename"
    
    ./gwas-pw -i "$input_file" -bed example_data/fourier_ls-all_modified.bed -phenos PD MeTs -o "$output_file"
done
