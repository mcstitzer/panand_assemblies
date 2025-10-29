import sys
import csv

def parse_agp(agp_file):
    """
    Parses the AGP file and constructs a mapping from scaffolds to chromosomes.
    """
    mapping = {}
    with open(agp_file) as agp:
        reader = csv.reader(agp, delimiter='\t')
        for row in reader:
            if row[4] == "W":  # Only deal with 'W' (contig rows)
                scaffold = row[5]
                mapping[scaffold] = {
                    'chromosome': row[0],
                    'chrom_start': int(row[1]),
                    'chrom_end': int(row[2]),
                    'scaf_start': int(row[6]),
                    'scaf_end': int(row[7]),
                    'strand': row[8]
                }
    return mapping

def uplift_gff(mapping, gff_file, output_file):
    """
    Uplifts GFF coordinates from scaffold space into chromosome space based on AGP mapping.
    """
    with open(gff_file) as gff, open(output_file, "w") as output:
        for line in gff:
            if line.startswith("#"):
                # Pass through comments and headers
                output.write(line)
                continue

            fields = line.strip().split("\t")
            scaffold = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]

            if scaffold in mapping:
                agp = mapping[scaffold]
                chrom = agp['chromosome']
                chrom_start = agp['chrom_start']
                chrom_end = agp['chrom_end']
                scaf_start = agp['scaf_start']
                scaf_end = agp['scaf_end']
                agp_strand = agp['strand']

                # Calculate new coordinates
                if agp_strand == "+":
                    uplifted_start = chrom_start + (start - scaf_start)
                    uplifted_end = chrom_start + (end - scaf_start)
                else:  # Reverse orientation
                    uplifted_start = chrom_end - (start - scaf_start)
                    uplifted_end = chrom_end - (end - scaf_start)
                    # Flip strand
                    strand = "-" if strand == "+" else "+"

                # Ensure coordinates are ordered correctly
                uplifted_start, uplifted_end = sorted([uplifted_start, uplifted_end])

                # Write updated GFF record
                fields[0] = chrom
                fields[3] = str(uplifted_start)
                fields[4] = str(uplifted_end)
                fields[6] = strand
                output.write("\t".join(fields) + "\n")
            else:
                # Skip scaffolds not found in AGP
                sys.stderr.write(f"Warning: Scaffold {scaffold} not found in AGP file\n")

def main():
    if len(sys.argv) != 4:
        print("Usage: python agp_to_gff.py <AGP_FILE> <GFF_FILE> <OUTPUT_FILE>")
        sys.exit(1)

    agp_file = sys.argv[1]
    gff_file = sys.argv[2]
    output_file = sys.argv[3]

    # Parse AGP file and uplift GFF file
    mapping = parse_agp(agp_file)
    uplift_gff(mapping, gff_file, output_file)

    print(f"Uplifted GFF file saved to {output_file}")

if __name__ == "__main__":
    main()