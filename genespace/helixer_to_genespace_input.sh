




### gffread to get aa fastas
input_file="panand_sp_ploidy.paper.txt"
while IFS=$'\t' read -r genome_name rest; do
gff_file="/data4/users/zrm22/HelixerRuns/annotations/${genome_name}_helixer.gff";     fasta_file="genomes/${genome_name}.fasta"    
output_protein_fasta="${genome_name}_proteins.fasta"
if [[ -f "$gff_file" && -f "$fasta_file" ]]; then         echo "Processing $genome_name..."
gffread -g "$fasta_file" -y "$output_protein_fasta" "$gff_file"         echo "Amino acid FASTA saved to $output_protein_fasta";     else         echo "Missing GFF or FASTA file for $genome_name. Skipping.";     fi; done < "$input_file"


## whoops go back and get cds so i can look for convergence!!!!
### gffread to get cds fastas
input_file="panand_sp_ploidy.paper.txt"
while IFS=$'\t' read -r genome_name rest; do
gff_file="/data4/users/zrm22/HelixerRuns/annotations/${genome_name}_helixer.gff";     
fasta_file="genomes/${genome_name}.fasta"   ; 
output_cds_fasta="${genome_name}_cds.fasta";
if [[ -f "$gff_file" && -f "$fasta_file" ]]; 
then         
echo "Processing $genome_name..."
gffread -g "$fasta_file" -x "$output_cds_fasta" "$gff_file"         
echo "CDS FASTA saved to $output_protein_fasta";     
else         
echo "Missing GFF or FASTA file for $genome_name. Skipping.";     
fi; 
done < "$input_file"



# sixcode="achine"
# awk -F'\t' -v sixcode="$sixcode" '
#     $3 == "gene" {
#         split($9, a, ";");
#         split(a[1], b, "=");
#         gene_id = sixcode b[2];
#         print $1, $4, $5, gene_id
#     }
# ' OFS='\t' Ac-Pasquet1232-DRAFT-PanAnd-1.0_helixer.gff >../../bed/achine.bed
# 
# sed -E "s/^>(.*)\.1$/>${sixcode}\1/" Ac-Pasquet1232-DRAFT-PanAnd-1.0_proteins.fasta > ../../peptide/achine.fasta
# 

### pvagin will work fine with phytozome preset within genespace

# Read panand_sp_ploidy.txt and loop through each line
while read -r gff_name sixcode ploidy; do
    # Define the input and output files
    fasta_file="${sixcode}/${gff_name}_proteins.fasta"
    gff_file="${sixcode}/${gff_name}_helixer.gff"
    bed_file=../bed/"${sixcode}.bed"
    updated_fasta_file=../peptide/"${sixcode}.fa"

    echo "Processing $gff_file and $fasta_file with sixcode $sixcode..."

    # Generate BED file using awk
    awk -F'\t' -v sixcode="$sixcode" '
        $3 == "gene" {
            split($9, a, ";");
            split(a[1], b, "=");
            gene_id = sixcode "_" b[2];
            print $1, $4, $5, gene_id
        }
    ' OFS='\t' "$gff_file" > "$bed_file"

    # Update FASTA file using sed
    sed -E "s/^>(.*)\.1$/>${sixcode}_\1/" "$fasta_file" > "$updated_fasta_file"

    echo "Finished processing for $gff_name. Output: $bed_file and $updated_fasta_file"

done < ../../panand_sp_ploidy.txt





#### get cdses so i can test for selection!!!!
##also fix cds names to match!!!!
while read -r gff_name sixcode ploidy; do
    # Define the input and output files
    fasta_file="helixer_protein_cds/${gff_name}_cds.fasta"
    updated_fasta_file="helixer_protein_cds/${sixcode}_cds.fa"


    # Update FASTA file using sed
    sed -E "s/^>(.*)\.1$/>${sixcode}_\1/" "$fasta_file" > "$updated_fasta_file"

    echo "Finished processing for $fasta_file. Output:  $updated_fasta_file"
done < panand_sp_ploidy.paper.txt

## fix paspalum
sed -E 's/(\.1).*//g' Pvaginatum_672_v3.1.cds_primaryTranscriptOnly.fa > pvagin.fa

## combine and index so i can extract as needed
cat *.fa > all_panand_helixer.cds.fa
samtools faidx all_panand_helixer.cds.fa


### make cds for each one!
grep '^>' ../orthofinder/Results_Nov14/Orthogroup_Sequences/OG0000016.fa | sed 's/>//' | xargs samtools faidx CDS_sequences.fa > OG0000016_CDS.fa
# Loop through all OG*.fa files
for og_file in ../orthofinder/Results_Nov14/Orthogroup_Sequences/OG*.fa; do
  # Extract the base filename (e.g., OG0000016)
  og_base=$(basename "$og_file" .fa)

  # Extract the headers, remove '>', and use samtools faidx to get CDS
  grep '^>' "$og_file" | sed 's/>//' | xargs samtools faidx ../helixer_protein_cds/all_panand_helixer.cds.fa > "${og_base}_CDS.fa"
  
  echo "Processed $og_file -> ${og_base}_CDS.fa"
done





