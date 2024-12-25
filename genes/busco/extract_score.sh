#!/bin/bash

# Create CSV header
echo "Sample,Complete,Single_copy,Multi_copy,Fragmented,Missing,n_markers" > busco_summary.csv

# Find all short_summary.json files recursively within run_poales_odb10 folders
find ./helixer_protein_aa/ -path "*/run_poales_odb10/short_summary.json" | while read file
do
    # Extract sample name (3 levels up from the JSON file)
    sample=$(basename $(dirname $(dirname "$file")))

    # Parse values directly under "results"
    complete=$(jq '.results.Complete' "$file")
    single=$(jq '.results."Single copy"' "$file")
    multi=$(jq '.results."Multi copy"' "$file")
    fragmented=$(jq '.results.Fragmented' "$file")
    missing=$(jq '.results.Missing' "$file")
    markers=$(jq '.results.n_markers' "$file")

    # Handle null or missing values by replacing with 0
    complete=${complete:-0}
    single=${single:-0}
    multi=${multi:-0}
    fragmented=${fragmented:-0}
    missing=${missing:-0}
    markers=${markers:-0}

    # Append data to CSV
    echo "$sample,$complete,$single,$multi,$fragmented,$missing,$markers" >> busco_summary.csv
done

echo "Summary saved in busco_summary.csv"
