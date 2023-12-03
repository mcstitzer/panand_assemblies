### Generate syntenic anchors relative to Paspalum

Downloaded Pvaginatum_672 v3.0 from Phytozome, v3.1 gene annotation from Phytozome

1. `anchorwave_get_panand_anchors.sh` to use AnchorWave to get longest syntenic path, guided by max ploidy of sample
2. `make_anchor_table.R` to count how many times a Paspalum gene is in syntenic blocks in each taxon
3. `extract_regions_for_genetrees_from_anchors.R` to get fastas of each region, in parallel (96cpu)
   - outputs to `gene_tree_beds/` for each gene
