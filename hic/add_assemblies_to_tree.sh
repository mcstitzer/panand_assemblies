

#!/bin/bash


OUTFA="Pavag03G039800_new.fa"
TMPBED="tmp_coords.bed"
> "$OUTFA"
> "$TMPBED"

### first column is anchorwave anchor output
##  second column is relative path to fasta...
# Loop over two-column file
while read shortname fastapath; do
    # Annotation file matching the shortname
    annofile="${shortname}"   # adjust if needed

    if [[ ! -f "$annofile" ]]; then
        echo "Annotation file for $shortname not found: $annofile" >&2
        continue
    fi
    if [[ ! -f "$fastapath" ]]; then
        echo "FASTA file for $shortname not found: $fastapath" >&2
        continue
    fi
samtools faidx $fastapath
    # Extract the Pavag03G039800 lines and convert to BED
grep 'Pavag03G039800' "$annofile" |
awk -v samp="$shortname" '{
    contig=$4; start=$5; end=$6;
    # bedtools uses 0-based start -> subtract 1
    header=samp"_"contig"_"start"-"end;
    print contig, start-1, end, header;
}' OFS='\t' > "$TMPBED"

bedtools getfasta -fi "$fastapath" -bed "$TMPBED" -nameOnly |
sed 's/^\(>[^:]\+\)::.*$/\1/' >> "$OUTFA"

done < samples_to_tree.txt

rm "$TMPBED"
echo "Combined sequences written to $OUTFA"


# Step 1: Align new sequences to existing alignment
/programs/mafft-7.525-with-extensions/bin/mafft --adjustdirection --add $OUTFA ../../Pavag03G039800.aln.fa > Pavag03G039800.updated.aln.fa

sed -E 's/^(>[^()]+)\([^)]*\)$/\1/' Pavag03G039800.updated.aln.fa > Pavag03G039800.cleaned.fa

# Step 2: Build ML tree
/programs/raxml-ng_v1.2.0/raxml-ng --all \
  --msa Pavag03G039800.cleaned.fa \
  --model GTR+G \
  --threads 4 \
  --seed 42 \
  --bs-trees 10 --redo
  
