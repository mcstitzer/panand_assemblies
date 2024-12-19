# Load necessary libraries
library(dplyr)
library(tidyr)
library(topGO)

# Step 1: Generate `final_result` (Aggregating GO Terms by Orthogroup)

# Example input data
orthogroups_data <- read.delim("Orthogroups.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(orthogroups_data) <- c("orthogroup", "achine_genes", "agerar_genes", "atenui_genes", "avirgi_genes", 
                                "blagur_genes", "ccitra_genes", "crefra_genes", "cserru_genes", "etrips_genes",
                                "hcompr_genes", "hconto_genes", "irugos_genes", "ppanic_genes", "Pavag_genes",
                                "rrottb_genes", "rtuber_genes", "sbicol_genes", "smicro_genes", "snutan_genes",
                                "sscopa_genes", "tdacn1_genes", "tdacn2_genes", "tdacs1_genes", "tdacs2_genes",
                                "telega_genes", "ttrian_genes", "udigit_genes", "vcuspi_genes", "zTIL01_genes",
                                "zTIL11_genes", "zTIL18_genes", "zTIL25_genes", "zdgigi_genes", "zdmomo_genes",
                                "zluxur_genes", "zmB735_genes", "zmhuet_genes", "znicar_genes")

# Reshape orthogroups data into long format
# Reshape orthogroups data into long format
orthogroups_long <- orthogroups_data %>%
  pivot_longer(cols = starts_with("achine_genes"):starts_with("znicar_genes"), 
               names_to = "species", values_to = "gene_list") %>%
  mutate(genes = strsplit(gene_list, ",")) %>%
  unnest(genes) %>%
  dplyr::select(orthogroup, genes)  # Explicitly use dplyr::select


# Example GO term data

go_files <- list.files(pattern = "_result.annot.topgoAnnot$", full.names = TRUE)
# Function to read a single file
read_go_file <- function(file_path) {
  read.delim(file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE) %>%
    setNames(c("gene_id", "go_terms"))  # Ensure consistent column names
}

# Combine all files into a single data frame
#go_data_combined <- do.call(rbind, lapply(go_files, read_go_file))

#go_data <- read.delim("vcuspi_result.annot.topgoAnnot", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
#colnames(go_data) <- c("gene_id", "go_terms")

go_data_combined <- do.call(rbind, lapply(go_files, read_go_file))

## also combine gene descirptions
# List all files matching the pattern
desc_files <- list.files(pattern = "_result.annot.geneDesc$", full.names = TRUE)

# Function to read a single file
read_desc_file <- function(file_path) {
  read.delim(file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE) %>%
    setNames(c("gene_id", "gene_description"))  # Ensure consistent column names
}

# Combine all files into a single data frame
desc_data_combined <- do.call(rbind, lapply(desc_files, read_desc_file))
# Map genes to their descriptions and orthogroups
orthogroup_to_desc <- orthogroups_long %>%
  left_join(desc_data_combined, by = c("genes" = "gene_id")) %>%  # Link genes to descriptions
  filter(!is.na(gene_description))  # Remove rows with missing descriptions
final_desc_result <- orthogroup_to_desc %>%
  group_by(orthogroup) %>%
  summarize(gene_descriptions = paste(unique(gene_description), collapse = "; "))  # Combine descriptions
write.table(final_desc_result, file = "final_gene_descriptions_by_orthogroup.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

#### There is a warning on one row, I'm not super concerned but will check it out
#### Warning message:
#### In left_join(., desc_data_combined, by = c(genes = "gene_id")) :
####   Detected an unexpected many-to-many relationship between `x` and `y`.
#### ℹ Row 262518 of `x` matches multiple rows in `y`.
#### ℹ Row 1786293 of `y` matches multiple rows in `x`.
#### ℹ If a many-to-many relationship is expected, set `relationship = "many-to-many"` to silence this warning.

## this is the one bothriochloa gene i messed up...
## > orthogroups_long[262518,]
## # A tibble: 1 × 2
##   orthogroup genes              
##   <chr>      <chr>              
## 1 OG0000550  _scaf_1355_000004.1
## > desc_data_combined[1786293,]
##                     gene_id             gene_description
## 1786293 _scaf_1355_000004.1 DNA-directed RNA polymerase 

geneID2GO <- go_data_combined %>%
  filter(!is.na(go_terms)) %>%  # Remove rows with missing GO terms
  mutate(go_terms = strsplit(go_terms, ",")) %>%  # Split GO terms into lists
  unnest(go_terms) %>%  # Expand each GO term into a separate row
  group_by(gene_id) %>%  # Group by gene_id
  summarize(go_terms = list(unique(go_terms))) %>%  # Aggregate unique GO terms
  { setNames(.$go_terms, .$gene_id) }  # Convert to named list


# Merge orthogroups with GO terms
combined_data <- orthogroups_long %>%
  left_join(go_data_combined, by = c("genes" = "gene_id"))

# Aggregate GO terms by orthogroup
final_result <- combined_data %>%
  group_by(orthogroup) %>%
  summarize(
    combined_go = paste(unique(unlist(strsplit(go_terms[!is.na(go_terms)], ","))), collapse = ",")
  )
# Save final_result as a tab-delimited file
write.table(final_result, file = "combined_go_terms_by_orthogroup.tsv", sep = "\t", row.names = FALSE, quote = FALSE)




# Step 2: Prepare Gene-to-GO Mapping for `topGO`
geneID2GO <- readMappings(file = "combined_go_terms_by_orthogroup.tsv")

# Step 3: Define Gene List for Enrichment Analysis
# All orthogroups in the mapping
all_orthogroups <- names(geneID2GO)

# Define orthogroups of interest
orthogroups_of_interest <- c("OG0000001", "OG0000003")  # Replace with your specific IDs


tethex=c("OG0001629", "OG0002099", "OG0001589", "OG0001754", "OG0001911", 
"OG0002143", "OG0002187", "OG0001506", "OG0001702", "OG0002008", 
"OG0001403", "OG0002095", "OG0001695", "OG0001268", "OG0001900", 
"OG0002387", "OG0002086", "OG0001653", "OG0001537", "OG0002142", 
"OG0002607", "OG0001582", "OG0002021", "OG0001848", "OG0001945", 
"OG0002715", "OG0001672", "OG0001563", "OG0001675", "OG0001868", 
"OG0001823", "OG0002278", "OG0001482", "OG0001712", "OG0002582", 
"OG0002261", "OG0002335", "OG0001931", "OG0002050", "OG0002491", 
"OG0002116", "OG0002507", "OG0001701", "OG0002318", "OG0001723", 
"OG0002210", "OG0001692", "OG0001904", "OG0001347", "OG0002910", 
"OG0001719", "OG0001725", "OG0002020", "OG0001417", "OG0002024", 
"OG0001586", "OG0001583", "OG0001921", "OG0002647", "OG0002308", 
"OG0001419", "OG0002298", "OG0001797", "OG0002188", "OG0002392", 
"OG0002471", "OG0002406", "OG0001996", "OG0002664", "OG0002157", 
"OG0002534", "OG0002610", "OG0002179", "OG0002617", "OG0001323", 
"OG0007213", "OG0002098", "OG0001761", "OG0001517", "OG0001726", 
"OG0001951", "OG0002710", "OG0001495", "OG0001682", "OG0001858", 
"OG0001630", "OG0002845", "OG0002037", "OG0002207", "OG0002135", 
"OG0001895", "OG0003052", "OG0001873", "OG0001451", "OG0001990", 
"OG0002726", "OG0002390", "OG0001609", "OG0001728", "OG0002942"
)
orthogroups_of_interest=tethex

annual=c("OG0014705", "OG0017420", "OG0017854", "OG0018138", "OG0015291", 
"OG0013545", "OG0013576", "OG0013898", "OG0014005", "OG0014228", 
"OG0014669", "OG0015001", "OG0015032", "OG0015040", "OG0015045", 
"OG0015179", "OG0015410", "OG0015932", "OG0015933", "OG0016453", 
"OG0016536", "OG0016829", "OG0016832", "OG0016873", "OG0017021", 
"OG0017037", "OG0017102", "OG0017248", "OG0017266", "OG0018363", 
"OG0018393", "OG0018763", "OG0018771", "OG0018938", "OG0018941", 
"OG0018968", "OG0018976", "OG0019133", "OG0019243", "OG0019286", 
"OG0019510", "OG0020222", "OG0014364", "OG0017885", "OG0011857", 
"OG0013882", "OG0013943", "OG0014793", "OG0015600", "OG0015682", 
"OG0015891", "OG0016301", "OG0016425", "OG0016877", "OG0016993", 
"OG0017411", "OG0017831", "OG0017848", "OG0018775", "OG0019325", 
"OG0019544", "OG0019793", "OG0019944", "OG0012556", "OG0012830", 
"OG0013396", "OG0013448", "OG0013768", "OG0013826", "OG0013886", 
"OG0013896", "OG0013987", "OG0014236", "OG0014529", "OG0014663", 
"OG0014707", "OG0014754", "OG0014757", "OG0014765", "OG0014766", 
"OG0014896", "OG0014995", "OG0015033", "OG0015034", "OG0015048", 
"OG0015231", "OG0015358", "OG0015413", "OG0015592", "OG0015602", 
"OG0015607", "OG0015622", "OG0015887", "OG0015890", "OG0015912", 
"OG0015921", "OG0015923", "OG0015929", "OG0015938", "OG0015947"
)
orthogroups_of_interest=annual


annualish=c("OG0007238", "OG0008173", "OG0008985", "OG0012253", "OG0013534", 
"OG0010263", "OG0006539", "OG0008412", "OG0004671", "OG0002812", 
"OG0007458", "OG0009859", "OG0015589", "OG0010799", "OG0009935", 
"OG0012206", "OG0006952", "OG0007670", "OG0007706", "OG0010112", 
"OG0012204", "OG0011945", "OG0011778", "OG0009289", "OG0011257", 
"OG0009260", "OG0010496", "OG0011793", "OG0008506", "OG0009463", 
"OG0009029", "OG0011366", "OG0012813", "OG0010061", "OG0011780", 
"OG0008482", "OG0003928", "OG0010965", "OG0010835", "OG0002866", 
"OG0015990", "OG0008411", "OG0008835", "OG0013381", "OG0015302", 
"OG0012260", "OG0008397", "OG0009931", "OG0009237", "OG0009679", 
"OG0009379", "OG0006729", "OG0002257", "OG0009058", "OG0009718", 
"OG0006772", "OG0007234", "OG0007757", "OG0007770", "OG0007883", 
"OG0007509", "OG0006053", "OG0010402", "OG0010076", "OG0010464", 
"OG0011214", "OG0010619", "OG0010827", "OG0007836", "OG0002924", 
"OG0011803", "OG0009715", "OG0006618", "OG0010206", "OG0009126", 
"OG0003167", "OG0007153", "OG0010173", "OG0010403", "OG0010627", 
"OG0007822", "OG0010792", "OG0016237", "OG0008580", "OG0008722", 
"OG0009515", "OG0006008", "OG0010298", "OG0011125", "OG0007675", 
"OG0007938", "OG0006597", "OG0009416", "OG0005223", "OG0013598", 
"OG0009915", "OG0004634", "OG0006435", "OG0009008", "OG0006983"
)
orthogroups_of_interest=annualish

polyploid=c("OG0001629", "OG0001506", "OG0001702", "OG0001695", "OG0001848", 
"OG0001563", "OG0001403", "OG0001653", "OG0001589", "OG0001931", 
"OG0001582", "OG0001672", "OG0001725", "OG0001268", "OG0001719", 
"OG0001754", "OG0001701", "OG0001537", "OG0001921", "OG0001675", 
"OG0001583", "OG0001823", "OG0001900", "OG0001419", "OG0001482", 
"OG0001728", "OG0001347", "OG0001712", "OG0001517", "OG0002135", 
"OG0001586", "OG0001951", "OG0001495", "OG0001723", "OG0001976", 
"OG0001451", "OG0002207", "OG0001990", "OG0001323", "OG0001973", 
"OG0001641", "OG0001930", "OG0001692", "OG0002240", "OG0001417", 
"OG0001609", "OG0001505", "OG0001651", "OG0001726", "OG0001953", 
"OG0001420", "OG0001415", "OG0001993", "OG0001873", "OG0001805", 
"OG0001708", "OG0001846", "OG0001870", "OG0001408", "OG0001758", 
"OG0001730", "OG0001945", "OG0001928", "OG0001927", "OG0001992", 
"OG0001798", "OG0002187", "OG0002015", "OG0001374", "OG0001568", 
"OG0001851", "OG0002131", "OG0001829", "OG0001660", "OG0001895", 
"OG0002091", "OG0001412", "OG0001911", "OG0001349", "OG0002172", 
"OG0002045", "OG0001535", "OG0001509", "OG0001529", "OG0001682", 
"OG0001827", "OG0001868", "OG0002132", "OG0001700", "OG0001706", 
"OG0004184", "OG0001910", "OG0003551", "OG0001861", "OG0002067", 
"OG0001473", "OG0001987", "OG0001460", "OG0001626", "OG0002152"
)
orthogroups_of_interest=polyploid

tannins=c("OG0006062", "OG0004244", "OG0008122", "OG0009316", "OG0012258", 
"OG0010087", "OG0010684", "OG0007611", "OG0012528", "OG0015438", 
"OG0009894", "OG0011572", "OG0012074", "OG0007709", "OG0008247", 
"OG0007890", "OG0012863", "OG0011309", "OG0007189", "OG0013521", 
"OG0014284", "OG0012349", "OG0013537", "OG0006940", "OG0010842", 
"OG0006154", "OG0010354", "OG0008124", "OG0002715", "OG0009992", 
"OG0009696", "OG0009703", "OG0008724", "OG0003998", "OG0012309", 
"OG0011672", "OG0008222", "OG0011694", "OG0011308", "OG0009340", 
"OG0010568", "OG0012862", "OG0012871", "OG0013434", "OG0006092", 
"OG0013520", "OG0010840", "OG0006291", "OG0008658", "OG0008738", 
"OG0008789", "OG0009376", "OG0009962", "OG0015060", "OG0012782", 
"OG0005903", "OG0012377", "OG0012866", "OG0005140", "OG0008115", 
"OG0010206", "OG0009126", "OG0011436", "OG0009503", "OG0006396", 
"OG0003072", "OG0010937", "OG0010659", "OG0011065", "OG0011242", 
"OG0009705", "OG0012174", "OG0005730", "OG0008745", "OG0009320", 
"OG0008989", "OG0010406", "OG0010600", "OG0010843", "OG0011310", 
"OG0009900", "OG0007458", "OG0008021", "OG0008249", "OG0010785", 
"OG0010802", "OG0012311", "OG0007123", "OG0009502", "OG0011770", 
"OG0010143", "OG0011314", "OG0009700", "OG0009243", "OG0007305", 
"OG0011777", "OG0013530", "OG0006339", "OG0003588", "OG0007968"
)
orthogroups_of_interest=tannins

selfing=c("OG0021131", "OG0012356", "OG0010192", "OG0013886", "OG0020275", 
"OG0020583", "OG0021455", "OG0012810", "OG0018065", "OG0019955", 
"OG0020408", "OG0020410", "OG0021866", "OG0009172", "OG0012098", 
"OG0021868", "OG0012436", "OG0016225", "OG0010645", "OG0020743", 
"OG0009016", "OG0004313", "OG0022122", "OG0006735", "OG0009212", 
"OG0009624", "OG0013946", "OG0013947", "OG0014412", "OG0014531", 
"OG0016361", "OG0017398", "OG0018157", "OG0019351", "OG0019463", 
"OG0013691", "OG0015143", "OG0015909", "OG0016063", "OG0019396", 
"OG0019715", "OG0019727", "OG0019728", "OG0019929", "OG0019931", 
"OG0019952", "OG0020256", "OG0020280", "OG0020289", "OG0020417", 
"OG0020585", "OG0020873", "OG0020893", "OG0010011", "OG0016817", 
"OG0003176", "OG0012018", "OG0008533", "OG0007096", "OG0011453", 
"OG0014718", "OG0007945", "OG0022198", "OG0015536", "OG0014618", 
"OG0015257", "OG0016263", "OG0016513", "OG0019586", "OG0014019", 
"OG0015978", "OG0017591", "OG0017812", "OG0020877", "OG0021788", 
"OG0022235", "OG0014485", "OG0015401", "OG0016199", "OG0016365", 
"OG0017335", "OG0018084", "OG0019194", "OG0019391", "OG0004501", 
"OG0003577", "OG0012593", "OG0020100", "OG0021003", "OG0011087", 
"OG0010159", "OG0022057", "OG0009569", "OG0011573", "OG0012478", 
"OG0021777", "OG0002018", "OG0010185", "OG0006797", "OG0011913"
)
orthogroups_of_interest=selfing



# Create gene list (binary vector) from orthogroups_of_interest
all_orthogroups <- names(geneID2GO)  # geneID2GO is your gene-to-GO mapping
gene_list <- factor(as.integer(all_orthogroups %in% orthogroups_of_interest), levels = c(0, 1))
names(gene_list) <- all_orthogroups  # Ensure names match orthogroups

# Step 2: Define a Function to Run Enrichment Analysis for One Ontology
run_topGO <- function(ontology, allGenes, gene2GO, nodeSize = 5) {
  # Create the topGOdata object
  tgd <- new(
    "topGOdata",
    ontology = ontology,
    allGenes = allGenes,
    nodeSize = nodeSize,
    annot = annFUN.gene2GO,
    gene2GO = gene2GO
  )
  
  # Run enrichment tests
  res_classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher")
  res_weight01 <- runTest(tgd, algorithm = "weight01", statistic = "Fisher")
  
  # Generate results table
  GO_res_table <- GenTable(
    tgd,
    Fisher.classic = res_classic,
    Fisher.weight01 = res_weight01,
    orderBy = "Fisher.weight01",
    ranksOf = "Fisher.classic",
    topNodes = length(res_classic@score),
    numChar = 100  # Adjust description length if needed
  )
  
  # Add ontology for reference
  GO_res_table$ontology <- ontology
  
  return(GO_res_table)
}

# Step 3: Run Enrichment for BP, MF, and CC
res_BP <- run_topGO("BP", gene_list, geneID2GO)
res_MF <- run_topGO("MF", gene_list, geneID2GO)
res_CC <- run_topGO("CC", gene_list, geneID2GO)

# Step 4: Combine and Sort Results
combined_res <- rbind(res_BP, res_MF, res_CC)  # Combine results
combined_res <- combined_res[order(as.numeric(combined_res$Fisher.weight01)), ]  # Sort by p-value
combined_res$Fisher.weight01=as.numeric(combined_res$Fisher.weight01)  # make numeric


## see what genes are in the enriched category...


  tgdbp <- new(
    "topGOdata",
    ontology = 'BP',
    allGenes = gene_list,
    nodeSize = 5,
    annot = annFUN.gene2GO,
    gene2GO = geneID2GO
  )

  tgdmf <- new(
    "topGOdata",
    ontology = 'MF',
    allGenes = gene_list,
    nodeSize = 5,
    annot = annFUN.gene2GO,
    gene2GO = geneID2GO
  )
  tgdcc <- new(
    "topGOdata",
    ontology = 'CC',
    allGenes = gene_list,
    nodeSize = 5,
    annot = annFUN.gene2GO,
    gene2GO = geneID2GO
  )


# Define significant GO terms for each ontology
significant_terms_list <- list(
  BP = c("GO:0048766", 'GO:0046777'),  # Replace with actual significant GO terms for BP
  CC = c("GO:0042788", 'GO:0022625'),                # Replace with actual significant GO terms for CC
  MF = c("GO:0009931", 'GO:0004683')                 # Replace with actual significant GO terms for MF
)

## do automatically
filtered_combined_res <- combined_res[combined_res$Fisher.weight01 < 0.05 , ]

significant_terms_list <- split(filtered_combined_res$GO.ID, filtered_combined_res$ontology)




# Create a list of topGOdata objects
topGOdata_list <- list(
  BP = tgdbp,
  CC = tgdcc,
  MF = tgdmf
)

# Extract enriched orthogroups for each ontology
combined_enriched_orthogroups <- do.call(rbind, lapply(names(topGOdata_list), function(ontology) {
  GOdata <- topGOdata_list[[ontology]]
  significant_terms <- significant_terms_list[[ontology]]
  
  # Extract orthogroups for each significant GO term
  enriched_orthogroups <- lapply(significant_terms, function(go_term) {
    genes <- genesInTerm(GOdata, go_term)
    significant_orthogroups <- intersect(genes[[1]], names(gene_list[gene_list == 1]))
    return(data.frame(Ontology = ontology, GO = go_term, Orthogroup = significant_orthogroups))
  })
  
  # Combine results for this ontology
  do.call(rbind, enriched_orthogroups)
}))

# View the combined results
print(combined_enriched_orthogroups)

## done above
# final_desc_result <- orthogroup_to_desc %>%
#   group_by(orthogroup) %>%
#   summarize(gene_descriptions = paste(unique(gene_description), collapse = "; "))

# Add gene descriptions to combined_enriched_orthogroups
combined_enriched_with_descriptions <- combined_enriched_orthogroups %>%
  left_join(final_desc_result, by = c("Orthogroup" = "orthogroup"))


genegroup='annual'
write.table(combined_res, paste0('enrichedGO_', genegroup, '.txt'), row.names=F, col.names=T, sep='\t', quote=F)
write.table(combined_enriched_with_descriptions, paste0('enrichedGO_eachOGdescription_', genegroup, '.txt'), row.names=F, col.names=T, sep='\t', quote=F)

genegroup='polyploid'
write.table(combined_res, paste0('enrichedGO_', genegroup, '.txt'), row.names=F, col.names=T, sep='\t', quote=F)
write.table(combined_enriched_with_descriptions, paste0('enrichedGO_eachOGdescription_', genegroup, '.txt'), row.names=F, col.names=T, sep='\t', quote=F)
