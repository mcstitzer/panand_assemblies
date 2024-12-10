


for genome in $(cut -f2 ../panand_sp_ploidy.txt); do     output_file="${genome}_blastresults.xml" if [ -f "$output_file" ]; then echo "Skipping $genome: $output_file already exists." else /programs/diamond/diamond blastp --db /home2/shared/genome_db/uniref90 --query /workdir/mcs368/genespace/peptide/${genome}.fa --outfmt 5 --max-target-seqs 100 --max-hsps 1 --evalue 1e-10 -t ./ --block-size 10 --index-chunks 1 -o "$output_file" -p 128 echo "Completed DIAMOND run for $genome." fi done

## on cbsumm10
for genome in $(cut -f2 ../panand_sp_ploidy.txt); do
    output_file="${genome}_blastresults.xml"
    annot_file="${genome}_result.annot"
    log_file="${genome}_annotatelogfile"
    result_prefix="${genome}_result"
    
    # Check if blast results exist
    if [ -f "$output_file" ]; then
        echo "Running Blast2GO for $genome..."
        /usr/local/blast2go/blast2go_cli.run \
            -properties annotation.prop \
            -useobo go.obo \
            -loadblast "$output_file" \
            -mapping -annotation -statistics all \
            -saveannot "$result_prefix" \
            -saveseqtable "$result_prefix" \
            -savereport "$result_prefix" \
            -tempfolder ./ >& "$log_file"
        echo "Completed Blast2GO for $genome. Results saved with prefix $result_prefix."
    else
        echo "Skipping Blast2GO for $genome: $output_file does not exist."
        continue
    fi

    # Check if annotation file was created
    if [ -f "$annot_file" ]; then
        echo "Processing annotation for $genome..."
        /shared_data/annotation2019_2/toTopGO.py "$annot_file"
        echo "Completed conversion for $genome. Files generated:"
        echo "  - ${result_prefix}.annot.topgoAnnot"
        echo "  - ${result_prefix}.annot.geneDesc"
    else
        echo "Skipping annotation conversion for $genome: $annot_file does not exist."
    fi
done

