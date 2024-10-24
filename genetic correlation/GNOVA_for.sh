#!/bin/bash

AD_folder="./AD"
PD_folder="./PD"
output_folder="./ALL-LDSC"

mkdir -p $output_folder

AD_files=("$AD_folder"/*.gz)

PD_files=("$PD_folder"/*.gz)

total_files=$(( ${#AD_files[@]} * ${#PD_files[@]} ))

count=0

for ad_file in "${AD_files[@]}"; do
  for pd_file in "${PD_files[@]}"; do
    ad_base=$(basename "${ad_file%.*}")
    pd_base=$(basename "${pd_file%.*}")

    output_file="$output_folder/${ad_base}-${pd_base}.results"

    ((count++))
    progress=$((count * 100 / total_files))
    echo -ne "Progress: $progress%    \r"

    python gnova.py ./"$ad_file" "$pd_file" --bfile ./EUR1000g/EUR.QC --out "$output_file"
  done
done

echo -e "\nAnalysis completed!"

