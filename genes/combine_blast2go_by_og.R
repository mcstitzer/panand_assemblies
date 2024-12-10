# Load necessary libraries
library(dplyr)
library(tidyr)
library(topGO)

# Step 1: Generate `final_result` (Aggregating GO Terms by Orthogroup)

# Example input data
orthogroups_data <- read.delim("Orthogroups/Orthogroups.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(orthogroups_data) <- c("orthogroup", "achine_genes", "agerar_genes", "atenui_genes", "avirgi_genes", 
                                "blagur_genes", "ccitra_genes", "crefra_genes", "cserru_genes", "etrips_genes",
                                "hcompr_genes", "hconto_genes", "irugos_genes", "ppanic_genes", "Pavag_genes",
                                "rrottb_genes", "rtuber_genes", "sbicol_genes", "smicro_genes", "snutan_genes",
                                "sscopa_genes", "tdacn1_genes", "tdacn2_genes", "tdacs1_genes", "tdacs2_genes",
                                "telega_genes", "ttrian_genes", "udigit_genes", "vcuspi_genes", "zTIL01_genes",
                                "zTIL11_genes", "zTIL18_genes", "zTIL25_genes", "zdgigi_genes", "zdmomo_genes",
                                "zluxur_genes", "zmB735_genes", "zmhuet_genes", "znicar_genes")

# Reshape orthogroups data into long format
orthogroups_long <- orthogroups_data %>%
  pivot_longer(cols = starts_with("achine_genes"):starts_with("znicar_genes"), 
               names_to = "species", values_to = "gene_list") %>%
  mutate(genes = strsplit(gene_list, ",")) %>%
  unnest(genes) %>%
  select(orthogroup, genes)

# Example GO term data
go_data <- read.delim("vcuspi_result.annot.topgoAnnot", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(go_data) <- c("gene_id", "go_terms")

# Merge orthogroups with GO terms
combined_data <- orthogroups_long %>%
  left_join(go_data, by = c("genes" = "gene_id"))

# Aggregate GO terms by orthogroup
final_result <- combined_data %>%
  group_by(orthogroup) %>%
  summarize(combined_go = paste(unique(unlist(strsplit(go_terms, ","))), collapse = ","))

# Save final_result as a tab-delimited file
write.table(final_result, file = "combined_go_terms_by_orthogroup.tsv", sep = "\t", row.names = FALSE, quote = FALSE)




# Step 2: Prepare Gene-to-GO Mapping for `topGO`
geneID2GO <- readMappings(file = "combined_go_terms_by_orthogroup.tsv")

# Step 3: Define Gene List for Enrichment Analysis
# All orthogroups in the mapping
all_orthogroups <- names(geneID2GO)

# Define orthogroups of interest
orthogroups_of_interest <- c("OG0000001", "OG0000003")  # Replace with your specific IDs

# Create a binary vector: 1 for orthogroups of interest, 0 for others
gene_list <- factor(as.integer(all_orthogroups %in% orthogroups_of_interest), levels = c(0, 1))
names(gene_list) <- all_orthogroups

# Step 4: Run `topGO`
# Create the topGOdata object
GOdata <- new(
  "topGOdata",
  ontology = "BP",  # Choose "BP", "MF", or "CC"
  allGenes = gene_list,
  annot = annFUN.gene2GO,
  gene2GO = geneID2GO
)

# Perform Fisher's exact test for enrichment
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

# Generate Results Table
allRes <- GenTable(
  GOdata,
  classicFisher = resultFisher,
  orderBy = "classicFisher",
  ranksOf = "classicFisher",
  topNodes = 10  # Adjust for more results
)

# Step 5: View and Save Results
print(allRes)  # Display the top results
write.table(allRes, file = "topGO_results.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
