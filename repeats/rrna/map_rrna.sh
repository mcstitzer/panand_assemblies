

### simple minimap2
## cbsuxm01:/workdir/mcs368/panand_assemblies/rrna/
while read genome shortname _; do
    output_file="${shortname}_rrna.fa"
    
    # Check if the output file exists, and skip if it does
    if [[ -f "$output_file" ]]; then
        echo "Skipping ${genome}.fasta -> ${output_file} (already exists)"
        continue
    fi

    echo "Processing ${genome}.fasta -> ${output_file}"
    
    minimap2 -x asm20 ../genomes/${genome}.fasta one_its.fa | \
    awk '{print $6 "\t" $8 "\t" $9}' | \
    bedtools getfasta -fi ../genomes/${genome}.fasta -bed - | \
    awk -v sn="${shortname}" '/^>/ {print ">" sn "_" substr($0,2)} !/^>/ {print}' > "$output_file"

done < 


while read genome shortname _; do
    output_file="${shortname}_rrna.fa"
    
    # Check if the output file exists, and skip if it does
    if [[ -f "$output_file" ]]; then
        echo "Skipping ${genome}.fasta -> ${output_file} (already exists)"
        continue
    fi

    echo "Processing ${genome}.fasta -> ${output_file}"
    
    minimap2 -x asm20 -f 0.0 ../genomes/${genome}.fasta one_its.fa | \
    awk '{print $6 "\t" $8 "\t" $9}' | \
    bedtools getfasta -fi ../genomes/${genome}.fasta -bed - | \
    awk -v sn="${shortname}" '/^>/ {print ">" sn "_" substr($0,2)} !/^>/ {print}' > "$output_file"

done < ../../panand_sp_ploidy.txt