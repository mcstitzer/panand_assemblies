
# 
# bedtools flank -i genes.gff -g genome.sizes -l 20000 -r 20000 > flanked_genes.bed
# bedtools makewindows -b flanked_genes.bed -w 100 > windows.bed
# 
# bedtools coverage -a windows.bed -b tes.gff > coverage.txt

## chat gpt made me really angry on this one... it couldn't keep track of strand.

#!/bin/bash

# Read plant information file
while read -r line; do
  # Extract genome ID
  GENOME_ID=$(echo "$line" | awk '{print $1}')
  
  # Define file paths
  GENE_FILE="../../helixer_annotations/${GENOME_ID}_helixer.gff"
  TE_FILE="trash/${GENOME_ID}_EDTAandTandemRepeat.gff3"
  GENOME_SIZE_FILE="../genomes/${GENOME_ID}.fasta.fai"
  OUTPUT_DIR="./coverage_results"

  # Create output directory if it doesn't exist
  mkdir -p "$OUTPUT_DIR"

  # Check if files exist
  if [[ -f "$GENE_FILE" && -f "$TE_FILE" && -f "$GENOME_SIZE_FILE" ]]; then

    # Filter gene entries and keep gene name
    FILTERED_GENE_FILE="$OUTPUT_DIR/${GENOME_ID}_filtered_genes.gff"
    awk '$3 == "gene" {print $1, $4, $5, $9, $6, $7}' OFS="\t" "$GENE_FILE" > "$FILTERED_GENE_FILE"

    # Generate output filenames
    FLANKED_GENES="$OUTPUT_DIR/${GENOME_ID}_flanked_genes.bed"
    WINDOWS="$OUTPUT_DIR/${GENOME_ID}_windows.bed"
    COVERAGE="$OUTPUT_DIR/${GENOME_ID}_coverage.txt"

    echo "Processing $GENOME_ID..."

    # Run bedtools commands
    bedtools flank -i "$FILTERED_GENE_FILE" -g "$GENOME_SIZE_FILE" -l 20000 -r 20000 > "$FLANKED_GENES"
#     awk 'BEGIN {OFS="\t"} {gene_id=$4; print $1, $2, $3, gene_id}' "$FLANKED_GENES" > "$FLANKED_GENES.tmp"
#     mv "$FLANKED_GENES.tmp" "$FLANKED_GENES"

#     bedtools makewindows -b "$FLANKED_GENES" -w 100 -i src | \
#       awk 'BEGIN {OFS="\t"} {
#         if ($4 != prev_gene) {
#           num_windows_upstream = 200;
#           num_windows_downstream = 200;
#           window_count = 0;
#           prev_gene = $4;
#         }
#         window_index = (++window_count <= num_windows_upstream) ? window_count - num_windows_upstream - 1 : window_count - num_windows_upstream;
#         print $0, window_index;
#       }' > "$WINDOWS"

	bedtools makewindows -b "$FLANKED_GENES" -w 100 -i srcwinnum  > "$WINDOWS"

    bedtools coverage -a "$WINDOWS" -b "$TE_FILE" > "$COVERAGE"
  else
    echo "Missing files for $GENOME_ID, skipping."
  fi

done < "../../panand_sp_ploidy.txt"

echo "Coverage calculations complete!"

## then I'll just make sure cach gene has 400 windows, and ddrop if not




## redo for syntenic
  
  while read -r line; do
  # Extract genome ID
  GENOME_ID=$(echo "$line" | awk '{print $1}')
  echo $GENOME_ID
  # Define file paths
  SYNT_FILE="../syntenic/${GENOME_ID}.syntenicAnchors.gff3"
  TE_FILE="trash/${GENOME_ID}_EDTAandTandemRepeat.gff3"
  GENOME_SIZE_FILE="../genomes/${GENOME_ID}.fasta.fai"
  OUTPUT_DIR="./coverage_results"
  echo $SYNT_FILE
  echo $TE_FILE
  echo $GENOME_SIZE_FILE
  echo $OUTPUT_DIR
  # Create output directory if it doesn't exist
  mkdir -p "$OUTPUT_DIR"

  # Check if files exist
  if [[ -f "$SYNT_FILE" && -f "$TE_FILE" && -f "$GENOME_SIZE_FILE" ]]; then

    # Generate output filenames
    FLANKED_SYNT="$OUTPUT_DIR/${GENOME_ID}_flanked_synt.bed"
    WINDOWS="$OUTPUT_DIR/${GENOME_ID}_syntwindows.bed"
    COVERAGE="$OUTPUT_DIR/${GENOME_ID}_syntcoverage.txt"

    echo "Processing $GENOME_ID..."

    # Run bedtools commands
    bedtools flank -i "$SYNT_FILE" -g "$GENOME_SIZE_FILE" -l 20000 -r 20000 > "$FLANKED_SYNT"

	bedtools makewindows -b "$FLANKED_SYNT" -w 100 -i srcwinnum  > "$WINDOWS"

    bedtools coverage -a "$WINDOWS" -b "$TE_FILE" > "$COVERAGE"
  else
    echo "Missing files for $GENOME_ID, skipping."
  fi

done < "../../panand_sp_ploidy.txt"
  
  
  