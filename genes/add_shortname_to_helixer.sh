#!/bin/bash

## oops ran helixer without name field
## and then when i put it in for orthofinder i added an extra underscore :(
## so these will be avirgi__ctgX_000001

# Define input file with mappings
mapping_file="../panand_sp_ploidy.txt"
mkdir -p for_maizegdb
# Loop through all GFF files
for gff_file in *_helixer.gff; do
  # Extract the prefix from the GFF filename
  file_prefix=$(basename "$gff_file" "_helixer.gff")

  # Find the corresponding species code from the mapping file
  species_code=$(awk -v prefix="$file_prefix" '$1 ~ prefix {print $2}' "$mapping_file")

  # Skip if no matching species code is found
  if [ -z "$species_code" ]; then
    echo "No matching species code for $gff_file"
    continue
  fi

  # Update the ID and Parent fields using sed and save to a new file
  output_file="for_maizegdb/${gff_file}"
  sed -E "s/ID=/ID=${species_code}_/g; s/Parent=/Parent=${species_code}_/g" "$gff_file" > "$output_file"

  echo "Processed and saved: $output_file"
done