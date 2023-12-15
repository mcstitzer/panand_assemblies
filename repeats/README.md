
Generate Tandem Repeat classifications and filter EDTA with them:

1. identify tandem repeats (`run_trash.sh`)
2. generate tandem repeat consensuses (`combine_trash.R`)
3. repeatmask genomes with tandem repeat consensuses (`repeatmask_with_tandems.sh`)
4. make repeatmasker gff look nicer (`make_repeatmasker_gff_friendlier.R`)
5. intersect EDTA repeats with TRASH tandem repeats (`intersect_EDTA_with_TRASH.sh`)

Decide on telomere repeat cutoffs (Aimee actually ran tidk) `filter_telomeres.R`
