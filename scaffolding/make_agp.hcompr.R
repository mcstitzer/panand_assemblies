# Load libraries (if necessary)
library(dplyr)
library(readr)
library(tibble)

## on my laptop

all=read.table('../panand_sp_ploidy.txt')
all=all[!all$V2 %in% c('pprate', 'tdactm', 'tzopol', 'osativ', 'bdista', 'agerjg', 'svirid', 'eophiu'),]

species='hcompr'
ploidy=3 ## haploid

# Step 1: Load the fai file
fai <- read_delim(paste0("fais/", all$V1[all$V2==species], '.fasta.fai'), delim = "\t", col_names = FALSE)
colnames(fai) <- c("contig", "length", "other1", "other2", "other3") # Adjust based on the `.fai` file format

fai <- as.data.frame(fai)

contig_lengths <- fai %>%
  dplyr::select(contig, length) %>%
  deframe() # Converts to a named vector with `contig` as keys and `length` as values

# Step 2: Load the data frame with directions, chromosomes, and ordering
# Simulated input dataset
filepath=paste0('../syntenic_anchors/anchors/', species, '-Pv-', all$V3[all$V2==species]*2)
#filepath='~/Downloads/crefra-Av-2'

color_palette=muted_colors
minBlock=10
title=''
refChrs=c(paste0('Chr0', 1:9), 'Chr10')
#refChrs=c(paste0('chr', 1:10))
queryChrs=''
queryChrtoFlip=''
ylabelspecies=''
ploidy=all$V3[all$V2==species]*2


#process_anchors_to_dotplot_SUPP <- function(filepath, color_palette=muted_colors, minBlock=10, title='', refChrs=c(paste0('Chr0', 1:9), 'Chr10'), queryChrs='', 
#                                              queryChrtoFlip='', ylabelspecies='',ploidy='') {
  # Load data
  data <- read.table(filepath, header = TRUE)
  data <- data[data$gene != 'interanchor', ]
  crscaf=unique(data$queryChr)
  ## get queryChrs if they aren't supplied
  if(queryChrs[1]==''){
    queryChrs=unique(data$queryChr)
  }
  
  
  # Reduce to blocks and calculate stats
  data <- data %>%
    filter(refChr%in%refChrs)%>%
    group_by(blockIndex) %>%
    mutate(blockLength = dplyr::n()) %>%
    group_by(queryChr) %>%
    mutate(freqStrand = names(which.max(table(strand))),
           maxChr = max(queryStart),
           freqRef = names(which.max(table(refChr)))) 
           
## manually flip strand! alt-scaf_27
 data$freqStrand[data$queryChr=='alt-scaf_27']=ifelse(data$freqStrand[data$queryChr=='alt-scaf_27']=='+', '-', '+')
  
  # Filter data based on block length
  data <- data[data$blockLength > minBlock, ]
  data=data[data$refChr==data$freqRef,]

#  data$refChr <- factor(data$refChr, levels = c(paste0('Chr0', 1:9), 'Chr10'))
  data$refChr <- factor(data$refChr, levels =refChrs)
  
  # Reverse strand calculations
  data <- data %>%
    arrange(freqRef, referenceStart, queryStart)
  data$queryChr <- factor(data$queryChr, levels = rev(data$queryChr[!duplicated(data$queryChr)]))
  data$revQueryStart <- data$queryStart
  data$revQueryStart[data$freqStrand == '-'] <- abs(data$queryStart - data$maxChr)[data$freqStrand == '-']
  
# Reorder data based on your query-chr ordering logic (specific to your case)
  ordered_data <- data %>%
    arrange(freqRef, referenceStart, queryStart)
  ordered_data$Chr=ordered_data$freqRef
  ordered_data$scaf=ordered_data$queryChr



## because this is a polyploid, we have to deal with each copy of the chromosome...


filtered_data <- ordered_data %>%
  mutate(geneNumber= as.numeric(substr(gene, 9, 14)))%>%
  # Group by 'Chr'
  group_by(Chr, scaf, freqStrand) %>%
  # Count occurrences of each 'Chr' and 'scaf' combination
  summarise(count = n(), start_gene = min(geneNumber, na.rm = TRUE),  # Adjust based on extracted numeric gene value
      end_gene = max(geneNumber, na.rm = TRUE),   # Same here
     .groups = "drop_last") %>%

  # Within each 'Chr', get the top ploidy combinations based on count
  slice_max(order_by = count, n = ploidy) %>%
  # Optionally arrange by Chr for readability
  arrange(Chr, desc(count))





# ordered_data <- data %>%
#   arrange(freq, Chr, scaf) # Adjust sorting based on your actual requirements

# filtered_data <- ordered_data %>%
#   distinct(Chr, scaf, .keep_all = TRUE) 

### now let's fix as needed
## fixes here are based on translocations/fusions in contigs
filtered_data$Chr=as.character(filtered_data$Chr)

### manually call finished contigs??

hcontochr=c(  #chr1
            'scaf_3', 'scaf_4', 'scaf_2',#chr2
#            ,#chr3
            'scaf_6',#chr4
            'scaf_1', 'scaf_5',#chr5
            'scaf_11', 'scaf_13',#chr6
            'scaf_14',#chr7
 #           ,#chr8
            #chr9
 #           ,# chr10
 			'scaf_12', ## chr5-chr9 fusion
            'scaf_19', 'scaf_15', 'scaf_18', 'scaf_17', 'scaf_23', 'scaf_28' ## chr8-chr10 fusions
            )
## may need to be careful to name homeologous chromosomes...

## need to fix a misrepair in crefra - how to do programaticaly??? arising because of flip??
#filtered_data[c(8,9),]=filtered_data[c(9,8),]
# Function to find matching contigs based on gene proximity

## determine these from gene trees that are arm-specific....
# gene_tree_matches=data.frame(
#                   contig1=c('alt-scaf_57', 'scaf_93','scaf_33','scaf_32',
#                   'scaf_38', 'alt-scaf_36'), ## chr2
#                   contig2=c('scaf_39', 'scaf_31','scaf_8', 'scaf_30',
#                   'scaf_26', 'alt-scaf_25') ## chr2
#                   )


forced_contig_matches=list()
for(i in 1:nrow(gene_tree_matches)){
forced_contig_matches[[length(forced_contig_matches) + 1]] <- list(
        contig1 = gene_tree_matches$contig1[i],
        contig2 = gene_tree_matches$contig2[i],
        gap_size = NA
      )
}

## maybe this will help??
gene_tree_matches=data.frame(contig1='blah',contig2='duh')

## eventually, add in gene tree matching as first combo of scaffolds (do gene trees first!)


find_matching_contigs <- function(data, max_gene_distance = 100000, forced_contig_matches=forced_contig_matches) {## 50000 is 50 genes

  # Calculate start and end positions based on geneNumber
#   contig_breaks <- data %>%
#     group_by(queryChr) %>%
#     summarise(
#       start_gene = min(geneNumber, na.rm = TRUE),  # Adjust based on extracted numeric gene value
#       end_gene = max(geneNumber, na.rm = TRUE),   # Same here
#       freqRef = first(freqRef),                  # Reference chromosome
#       strand = first(freqStrand)                 # Strand orientation
#     )

## sort by start gene!
   contig_breaks=data%>%arrange(start_gene)
   contig_breaks$visited=F
   contig_breaks$forced=F
   contig_breaks$forced[contig_breaks$scaf%in%c(gene_tree_matches$contig1, gene_tree_matches$contig2)]=T
  # Identify matching pairs based on proximity
  matches <- list()
  for (i in seq_len(nrow(contig_breaks))) {
    current <- contig_breaks[i, ]
    if(current$forced){
    matches[[length(matches)+1]]<- list(
     contig1=as.character(current$scaf),
     contig2=gene_tree_matches$contig2[gene_tree_matches$contig1==current$scaf],
     gap_size=NA
    )
    contig_breaks$visited[contig_breaks$scaf%in%c(as.character(current$scaf), gene_tree_matches$contig2[gene_tree_matches$contig1==current$scaf])]=T ## get rid of contigs from original!!!

    }else if(current$visited==F){
    # Look for matching scaffolds based on proximity
    potential_matches <- contig_breaks %>%
      filter(
        Chr == current$Chr,             # Same chromosome
        scaf != as.character(current$scaf),           # Avoid self-match
        abs(start_gene - current$end_gene) <= max_gene_distance, # Check proximity
        visited==F,
        forced==F
      ) %>%
      mutate(gap_size = abs(start_gene - current$end_gene))

    if (nrow(potential_matches) > 0) {
      # Sort potential matches by smallest gap size
      potential_matches <- potential_matches %>% arrange(gap_size)

      # Save the best match (smallest gap)
      matches[[length(matches) + 1]] <- list(
        contig1 = as.character(current$scaf),
        contig2 = potential_matches$scaf[1],
        gap_size = potential_matches$gap_size[1]
      )
      contig_breaks$visited[contig_breaks$scaf%in%c(current$scaf, potential_matches$scaf[1])]=T ## get rid of contigs from original!!!
    }
    }
  }

  return(matches)
}

# Preprocess strand-aware gene adjustments

# Separate full-length and partial scaffolds
## cahnged this to filtered_daata ?!?!?!
full_length_scaffolds <- filtered_data %>% filter(scaf %in% hcontochr)
partial_scaffolds <- filtered_data %>% filter(!scaf %in% hcontochr)

# Initialize output rows and processed scaffolds tracker
agp_rows <- list()
processed_contigs <- c()  # Track scaffolds already processed

# Process each reference chromosome
for (chr in unique(ordered_data$freqRef)) {
  # Process full-length scaffolds
  chr_full <- filter(full_length_scaffolds, Chr == chr)
  current_start <- 1
  part_number <- 1

  # Add full-length scaffolds with underscore naming
  for (i in seq_len(nrow(chr_full))) {
    contig_name <- as.character(chr_full$scaf[i])
    contig_length <- contig_lengths[contig_name]
    orientation <- chr_full$freqStrand[i]

    # Ensure scaffold hasn't already been processed
    if (!(contig_name %in% processed_contigs)) {
      full_row <- data.frame(
        object = paste0(chr, "_", part_number),  # Add underscore to full scaffold name
        object_beg = current_start,
        object_end = current_start + contig_length - 1,
        part_number = part_number,
        component_type = "W",
      component_id = as.character(contig_name),    # Contig name as a character (W row)
      component_beg = as.character(1),                           # Start of contig in FASTA (1-based indexing)
      component_end = as.character(contig_length),               # End of contig in FASTA
      orientation = orientation,                   # Orientation of the contig
        stringsAsFactors = FALSE
      )

      # Append row and update position
      agp_rows <- append(agp_rows, list(full_row))
      processed_contigs <- c(processed_contigs, contig_name)  # Mark as processed
#      current_start <- current_start + contig_length ### DO NOT DO THIS! THIS CHROMOSOME IS FINISHED?!!!
      part_number <- part_number + 1
    }
  }

  # Process partial scaffolds: Find matching pairs
  chr_partial <- filter(partial_scaffolds, Chr == chr & count>50)%>%arrange(start_gene)
  matches <- find_matching_contigs(chr_partial, forced_contig_matches=forced_contig_matches)

## check here if length(matches)>ploidy
## then could do the joining of more than 2 contigs per chromosomes!!!
# if(length(matches)>ploidy){
# }

  # Process matched pairs first
  for (match in matches) {
    current_start=1 ## restart each pair!!
    scaf_part_number=1 ## restart this will be for this Chr_X which parts
    contig1 <- as.character(match$contig1)
    contig2 <- as.character(match$contig2)
    gap_size <- match$gap_size
    
    
    # Only process pairs if neither scaffold has been processed yet
    if (!(contig1 %in% processed_contigs || contig2 %in% processed_contigs)) {
      # Process contig1
      contig_length1 <- contig_lengths[contig1]
      orientation1 <- chr_partial %>% filter(scaf == contig1) %>% pull(freqStrand)

      row1 <- data.frame(
        object = paste0(chr, "_", part_number),
        object_beg = current_start,
        object_end = current_start + contig_length1 - 1,
        part_number = scaf_part_number,
        component_type = "W",
        component_id = as.character(contig1),
        component_beg = 1,
        component_end = as.character(contig_length1),
        orientation = orientation1,
        stringsAsFactors = FALSE
      )

      agp_rows <- append(agp_rows, list(row1))
      processed_contigs <- c(processed_contigs, contig1)  # Mark as processed
#      chr_partial=chr_partial[chr_partial$scaf!=contig1,] ## be more aggressive and remove it here??
      current_start <- current_start + contig_length1
      scaf_part_number <- scaf_part_number+1 

      # Add gap row between contig1 and contig2
      gap_row <- data.frame(
        object = paste0(chr, "_", part_number),
        object_beg = current_start,
        object_end = current_start + 100 - 1,
        part_number = scaf_part_number,
        component_type = "U",
        component_id = as.character(100), ## was putting in gene distance here, but don't do that!!!
        component_beg = "contig",
        component_end = "no",
        orientation = "na",
        stringsAsFactors = FALSE
      )

      agp_rows <- append(agp_rows, list(gap_row))
      current_start <- current_start + 100
      scaf_part_number <- scaf_part_number +1

      # Process contig2
      contig_length2 <- contig_lengths[contig2]
      orientation2 <- chr_partial %>% filter(scaf == contig2) %>% pull(freqStrand)

      row2 <- data.frame(
        object = paste0(chr, "_", part_number),
        object_beg = current_start,
        object_end = current_start + contig_length2 - 1,
        part_number = scaf_part_number,
        component_type = "W",
        component_id = contig2,
        component_beg = 1,
        component_end = as.character(contig_length2),
        orientation = orientation2,
        stringsAsFactors = FALSE
      )

      agp_rows <- append(agp_rows, list(row2))
      processed_contigs <- c(processed_contigs, contig2)  # Mark as processed
#     chr_partial=chr_partial[chr_partial$scaf!=contig2,] ## be more aggressive and remove it here??

      current_start <- current_start + contig_length2
      part_number <- part_number + 1 ## this is fine, now we're done!
    }
  }
## for this don't have unmatched?? iun reality, do this iteratively in a loop - join the big chunks, then add on the next...



#   # Process unmatched partial scaffolds
#   unmatched_contigs <- setdiff(chr_partial$scaf, processed_contigs)
# 
#   for (contig in unmatched_contigs) {
#     contig_length <- contig_lengths[contig]
#     orientation <- chr_partial %>% filter(scaf == contig) %>% pull(freqStrand)
# 
#     unmatched_row <- data.frame(
#       object = paste0(chr, "_", part_number),
#       object_beg = current_start,
#       object_end = current_start + contig_length - 1,
#       part_number = part_number,
#       component_type = "W",
#       component_id = contig,
#       component_beg = 1,
#       component_end = contig_length,
#       orientation = orientation,
#       stringsAsFactors = FALSE
#     )
# 
#     agp_rows <- append(agp_rows, list(unmatched_row))
#     processed_contigs <- c(processed_contigs, contig)  # Mark as processed
#     current_start <- current_start + contig_length
#     part_number <- part_number + 1
#   }
}

# Combine rows into final AGP DataFrame
agp_df <- do.call(rbind, agp_rows)





## now rename "chromosomes" by length

# Preprocessing: Summing component lengths for all "object"
agp_df_summed <- agp_df %>%
  group_by(object) %>%
  summarise(sum_length = sum(as.numeric(component_end[component_end != 'no']), na.rm = TRUE))

# Extract the base group (e.g., "Chr01" from "Chr01_1") and find the MAX length for each group
grouped_lengths <- agp_df_summed %>%
  mutate(base_group = case_when(
    object == "chr2-4fusion_1" ~ "Chr02", # Special case for chr2-4fusion_1
    object == "chr4-2fusion_1" ~ "Chr04", # Special case for chr4-2fusion_1
    TRUE ~ str_extract(object, "^[^_]+")  # Extract base group for other objects
  )) %>% # Handle special cases and extract group base (e.g., Chr01)
  group_by(base_group) %>%
  summarise(max_length = max(sum_length)) %>% # Find the longest object in each group
  arrange(desc(max_length)) %>%
  mutate(new_chr = paste0("chr", row_number())) # Assign new `chrX` names based on ranking of longest object
  
  
  
# Join back to assign new `chr` name to each base group
agp_df_with_new_chr <- agp_df_summed %>%
  mutate(base_group = case_when(
    object == "chr2-4fusion_1" ~ "Chr02", # Special case for chr2-4fusion_1
    object == "chr4-2fusion_1" ~ "Chr04", # Special case for chr4-2fusion_1
    TRUE ~ str_extract(object, "^[^_]+")  # Extract base group for other objects
  )) %>% # Handle special cases and extract group base (e.g., Chr01)
  left_join(grouped_lengths, by = "base_group") %>%
  group_by(base_group) %>%
  arrange(desc(sum_length)) %>% # Rank individual objects within each group
  mutate(object_rank = row_number(), # Rank objects (e.g., _1, _2, ...)
         new_object = paste0(new_chr, "_", object_rank)) %>% # Combine new chr name with rank
  ungroup()

# Update the original `agp_df` with the new object names
final_agp_df <- agp_df %>%
  left_join(agp_df_with_new_chr %>% dplyr::select(object, new_object), by = "object") %>%
  mutate(object = new_object) %>%
  dplyr::select(-new_object)



# Step 1: Identify unplaced scaffolds
unplaced_scaffolds <- setdiff(names(contig_lengths), final_agp_df$component_id)

# Step 2: Add unplaced scaffolds as individual records
# Initialize a list for unplaced scaffold rows
unplaced_rows <- list()

for (i in seq_along(unplaced_scaffolds)) {
  scaffold_name <- unplaced_scaffolds[i]
  scaffold_length <- contig_lengths[scaffold_name]
  
  # Create a single AGP entry for the unplaced scaffold
  unplaced_row <- data.frame(
    object = as.character(scaffold_name),               # Give each unplaced scaffold a unique name
    object_beg = 1,                                # Start always at 1
    object_end = scaffold_length,                 # End is the full length of the scaffold
    part_number = 1,                               # Only one part for unplaced scaffolds
    component_type = "W",                          # W indicates contig placement
    component_id = as.character(scaffold_name),    # Contig name
    component_beg = "1",                           # Start of the contig
    component_end = as.character(scaffold_length), # End of the contig
    orientation = "+"                              # Default orientation
  )
  
  # Append to the unplaced scaffold list
  unplaced_rows <- append(unplaced_rows, list(unplaced_row))
}

# Step 3: Combine the unplaced rows with the existing AGP data
unplaced_df <- do.call(rbind, unplaced_rows)
final_agp_df <- bind_rows(final_agp_df, unplaced_df)


# Save AGP file
write_delim(final_agp_df, paste0(species, ".agp"), delim = "\t", col_names = FALSE)

# AGP file is now generated and saved as "output.agp"


## source /programs/miniconda3/bin/activate ragtag
## ragtag.py agp2fa crefra.agp ../genomes/Cr-AUB069-DRAFT-PanAnd-1.0.fasta > Cr-AUB069-DRAFT-PanAnd-1.0.chrSuperScaf.fa
## will need to add back in scaffolds that aren't in the agp....

# 
# source /home/$USER/miniconda3/bin/activate ## on cbsu
#  source activate anchorwave_new
#  reffa=../genomes/Pvaginatum_672_v3.0.fa
#  ls
#  refgff=udigit/Pvaginatum_672_v3.1.gene.gff3 
#  ref=Pv
#  ##
#  ploidy=2
#  minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 Cr-AUB069-DRAFT-PanAnd-1.0.chrSuperScaf.fa udigit/Pv.CDS.fa > crefraCHR-Pv.CDS.sam
#  anchorwave proali -t 10 -i $refgff -as  udigit/Pv.CDS.fa -r $reffa -a crefraCHR-Pv.CDS.sam -ar udigit/Pv.CDS.sam -s Cr-AUB069-DRAFT-PanAnd-1.0.chrSuperScaf.fa -n crefraCHR-Pv-2 -R 2 -Q 1 -ns





