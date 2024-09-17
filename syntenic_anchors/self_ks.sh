#!/bin/bash

# Define the taxon names
taxonnames=("zTIL11" "zmB735" "zTIL01" "zTIL25" "zTIL18" "zmhuet" "zluxur" "znicar" "zdmomo" "zdgigi" 
            "tzopol" "tdacs1" "tdacs2" "tdacn2" "tdacn1" "tdactm" "udigit" "vcuspi" "rrottb" "rtuber" 
            "hcompr" "etrips" "sscopa" "smicro" "avirgi" "achine" "agerar" "crefra" "ccitra" "hconto" 
            "ttrian" "blagur" "ppanic" "sbicol" "irugos" "snutan" "atenui" "telega" "cserru" "pvagin")

# Specify the input file
input_file="/data1/users/shengkaihsu/p_panAndOGASR/output/anchorAln_omega_v2/*"

# Loop over each taxon name and use grep to find lines with the pattern "taxon.*taxon"
for taxon in "${taxonnames[@]}"
do
  # Use grep to search for lines where the taxon appears twice with any characters between
  grep -E "${taxon}.*${taxon}" $input_file > "${taxon}_matches.txt"
  
  # Optionally, print a message
  echo "Searching for ${taxon} and saving results to ${taxon}_matches.txt"
done
