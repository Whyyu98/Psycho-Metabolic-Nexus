#!/bin/bash

AD_folder="./AD"
output_folder="./s-LDSC/result/AD-Franke"

mkdir -p $output_folder

AD_files=("$AD_folder"/*.gz)

for ad_file in "${AD_files[@]}"; do
 
    ad_base=$(basename "${ad_file%.*}")

    output_file="$output_folder/${ad_base}.results"

    ./ldsc.py --h2-cts "$ad_file" --ref-ld-chr ./s-LDSC/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD. --w-ld-chr ./s-LDSC/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. --ref-ld-chr-cts ./Franke1.ldcts --out "$output_file"
  done
done

