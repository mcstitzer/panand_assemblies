#!/bin/bash

# Define the taxon names
taxonnames=("zTIL11" "zmB735" "zTIL01" "zTIL25" "zTIL18" "zmhuet" "zluxur" "znicar" "zdmomo" "zdgigi" 
            "tzopol" "tdacs1" "tdacs2" "tdacn2" "tdacn1" "tdactm" "udigit" "vcuspi" "rrottb" "rtuber" 
            "hcompr" "etrips" "sscopa" "smicro" "avirgi" "achine" "agerar" "crefra" "ccitra" "hconto" 
            "ttrian" "blagur" "ppanic" "sbicol" "irugos" "snutan" "atenui" "telega" "cserru" "pvagin")

# Specify the directory containing the input files
input_dir="/data1/users/shengkaihsu/p_panAndOGASR/output/anchorAln_omega_v2/"

# Loop over each file in the directory
for file in "$input_dir"/*
do
  echo "Processing file: $file"
  
  # Read the file once and use a loop to grep for each taxon
  while IFS= read -r line; do
    for taxon in "${taxonnames[@]}"; do
      # Check if the line contains the taxon twice
      if [[ $line =~ $taxon.*$taxon ]]; then
        echo "$line" >> "${taxon}_matches.txt"
      fi
    done
  done < "$file"
done

echo "Processing of all files in $input_dir completed."
