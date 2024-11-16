

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
