#!/bin/bash

## if conda's being annoying, run this first
source /home/$USER/miniconda3/bin/activate

# Usage: ./workflow.sh Cr-AUB069-DRAFT-PanAnd-1.0 crefraCHR
# This script assumes all software and dependencies are correctly installed and in the PATH.

# Input arguments
GENOME_BASE=$1  # e.g., Cr-AUB069-DRAFT-PanAnd-1.0
PREFIX=$2       # e.g., crefraCHR (output prefix)
AGP_FILE=$3     # e.g., crefra.agp
PLOIDY=$4       # e.g., 2 for tetraploid, as these get multiplied by 2

# Paths to input/output files (update these paths as needed)
GENOME_FILE="../genomes/${GENOME_BASE}.fasta"
#AGP_FILE="crefra.agp"  # Ensure the AGP file is supplied and resides in this directory
TRASH_GFF="../repeats/trash/${GENOME_BASE}_EDTAandTandemRepeat.gff3"
HELIXER_GFF="../../helixer_annotations/${GENOME_BASE}_helixer.gff"
CDS_FILE="udigit/Pv.CDS.fa"
OUTPUT_SCAFFOLD="${GENOME_BASE}.chrSuperScaf.fa"
OUTPUT_TRASH_GFF="${GENOME_BASE}.chrSuperScaf.EDTAandTandemRepeat.gff3"
OUTPUT_HELIXER_GFF="${GENOME_BASE}.chrSuperScaf.helixer.gff"
REF_GENOME="../genomes/Pvaginatum_672_v3.0.fa"
REF_GFF="udigit/Pvaginatum_672_v3.1.gene.gff3"
#PLOIDY=2
PLOIDY=$PLOIDY*2

# Exit script if any command fails
set -e

source activate scaffolding
# Step 1: Use agptools to assemble the genome
echo "Step 1: Assembling genome with agptools..."
agptools assemble "$GENOME_FILE" "$AGP_FILE" > "$OUTPUT_SCAFFOLD"

# Step 2: Uplift the repeat annotation GFF using the AGP file
echo "Step 2: Uplifting repeat annotations..."
python uplift_gff.py "$AGP_FILE" "$TRASH_GFF" "$OUTPUT_TRASH_GFF"

# Step 3: Uplift the helixer GFF using the AGP file
echo "Step 3: Uplifting helixer annotations..."
python uplift_gff.py "$AGP_FILE" "$HELIXER_GFF" "$OUTPUT_HELIXER_GFF"

# Step 4: Align sequences and perform AnchorWave analysis
echo "Step 4: Running alignment and AnchorWave..."

# Activate your conda environment (CBSU-specific setup)
source /home/$USER/miniconda3/bin/activate
source activate anchorwave_new

# Run minimap2 to align the assembled genome against CDS
echo "Running minimap2..."
minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 "$OUTPUT_SCAFFOLD" "$CDS_FILE" > "${PREFIX}-Pv.CDS.sam"

## remove previous anchors
touch "${PREFIX}-Pv-${PLOIDY}" 
rm "${PREFIX}-Pv-${PLOIDY}" 

# Run AnchorWave with the alignment
echo "Running AnchorWave..."
anchorwave proali -t 10 \
  -i "$REF_GFF" \
  -as "$CDS_FILE" \
  -r "$REF_GENOME" \
  -a "${PREFIX}-Pv.CDS.sam" \
  -ar udigit/Pv.CDS.sam \
  -s "$OUTPUT_SCAFFOLD" \
  -n "${PREFIX}-Pv-${PLOIDY}" \
  -R "${PLOIDY}" \
  -Q 1 \
  -ns

# Step 5: Generate TE and gene visualization using R
echo "Step 5: Plotting TE and gene distributions..."
Rscript plot_te_genes_on_new_chr.R "$PREFIX" "$OUTPUT_TRASH_GFF" "$OUTPUT_HELIXER_GFF" "$AGP_FILE"
Rscript plot_te_genes_newChrCombined.R "$PREFIX" "$OUTPUT_TRASH_GFF" "$OUTPUT_HELIXER_GFF" "$AGP_FILE" "$PLOIDY"

#Rscript plot_te_genes_newChrCombined.R cserruCHR Cs-KelloggPI219580-DRAFT-PanAnd-1.0.chrSuperScaf.EDTAandTandemRepeat.gff3 Cs-KelloggPI219580-DRAFT-PanAnd-1.0.chrSuperScaf.helixer.gff cserru.agp 2

## eventually should also plot the dotplot!!!

# Step 6: Plot Dotplot
echo "Step 6: Plotting dotplot ..."
Rscript plot_anchor_dotplot.R "$PREFIX" "$AGP_FILE" "$PLOIDY"
#Rscript plot_anchor_dotplot.R cserruCHR cserru.agp 2

echo "Workflow completed successfully!"




## 