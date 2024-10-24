#!/bin/bash

folder_path="/mnt/e/GXN/MTAG_CPASSOC/MTAG_CPASSOC" 

for txt_file in "$folder_path"/*.txt; do
  file_prefix=$(basename "$txt_file" .txt)
  
  ./plink --bfile ./g1000_eur/g1000_eur --clump-p1 1 --clump-r2 0.2 --clump-kb 500 --clump "$txt_file" --clump-snp-field SNP --clump-field P.MTAG --out ./SNV/"$file_prefix"
done
