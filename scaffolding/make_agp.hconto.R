# Load libraries (if necessary)
library(dplyr)
library(readr)
library(tibble)



all=read.table('../panand_sp_ploidy.txt')
all=all[!all$V2 %in% c('pprate', 'tdactm', 'tzopol', 'osativ', 'bdista', 'agerjg', 'svirid', 'eophiu'),]

species='hconto'

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
  
  # Filter data based on block length
  data <- data[data$blockLength > minBlock, ]
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
  # Group by 'Chr'
  group_by(Chr, scaf, freqStrand) %>%
  # Count occurrences of each 'Chr' and 'scaf' combination
  summarise(count = n(), .groups = "drop_last") %>%
  # Within each 'Chr', get the top 4 combinations based on count
  slice_max(order_by = count, n = 4) %>%
  # Optionally arrange by Chr for readability
  arrange(Chr, desc(count))


# ordered_data <- data %>%
#   arrange(freq, Chr, scaf) # Adjust sorting based on your actual requirements

# filtered_data <- ordered_data %>%
#   distinct(Chr, scaf, .keep_all = TRUE) 

### now let's fix as needed
## fixes here are based on translocations/fusions in contigs
filtered_data$Chr=as.character(filtered_data$Chr)
filtered_data$Chr[filtered_data$scaf=='scaf_9']='chr2-4fusion' ## this is t2t so leave it be
filtered_data$Chr[filtered_data$scaf=='scaf_40']='chr4-2fusion' ## not T2T, has to merge with scaf_41
filtered_data$Chr[filtered_data$scaf=='scaf_41']='chr4-2fusion' ## not T2T, has to merge with scaf_41

### manually call finished contigs??

hcontochr=c('scaf_11', 'scaf_7', 'scaf_6', 'scaf_5',#chr1
            'scaf_4', 'scaf_10', 'scaf_8',#chr2
            'scaf_13', 'scaf_15', 'scaf_14', 'scaf_12',#chr3
            'scaf_20', 'scaf_21', 'scaf_19',#chr4
            'scaf_1', 'scaf_2', 'scaf_3',#chr5
            'scaf_27', 'scaf_29', 'scaf_35',#chr6
            'scaf_30', 'scaf_26', 'scaf_18', 'scaf_33',#chr7
            'scaf_16', 'scaf_22', 'scaf_17',#chr8
            'scaf_25', 'scaf_32', 'scaf_34',#chr9
            'scaf_31', 'scaf_23', 'scaf_28',# chr10
            'scaf_9' # fusion 2-4
            )
## may need to be careful to name homeologous chromosomes...

## need to fix a misrepair in crefra - how to do programaticaly??? arising because of flip??
#filtered_data[c(8,9),]=filtered_data[c(9,8),]

agp_rows <- list()
# Loop over unique Chromosomes


# Process each chromosome individually
for (chr in unique(filtered_data$Chr)) {
  # Filter data for the current chromosome
  chr_data <- filter(filtered_data, Chr == chr)
  chrname=
  current_start <- 1  # Reset position at the start of every chromosome
  part_number <- 1    # Reset part number at the start of every chromosome
  part_number_nocontig <- 1


## this loop is working just fine...
if(nrow(chr_data)==ploidy){ ## then we just make lines for each (add in new name for object??)
  for (i in seq_len(nrow(chr_data))) {
    # Extract details for this contig
    contig_name <- as.character(chr_data$scaf[i])
    contig_length <- contig_lengths[contig_name]  # Get contig length from `.fai`
    orientation <- chr_data$freqStrand[i]          # Get contig orientation
    
    # Generate a 'W' row for the contig placement
    agp_row <- data.frame(
      object = paste0(chr, '_', i),                # Chromosome name (needs to be different for each "subgenome")
      object_beg = current_start,                  # Start of the contig
      object_end = current_start + contig_length - 1,  # End of the contig
      part_number = part_number,                   # Increment part number
      component_type = "W",                        # W indicates contig placement
      component_id = as.character(contig_name),    # Contig name as a character (W row)
      component_beg = as.character(1),                           # Start of contig in FASTA (1-based indexing)
      component_end = as.character(contig_length),               # End of contig in FASTA
      orientation = orientation,                   # Orientation of the contig
      stringsAsFactors = FALSE
    )
#    print(agp_row)
    # Append this contig row to the list
    agp_rows <- append(agp_rows, list(agp_row))
}
}else{ ## then we need to combine things maybe... but how to know which ones are done?? sort by length?

  # Generate rows for the current chromosome
  for (i in seq_len(nrow(chr_data))) {
    # Extract details for this contig
    contig_name <- as.character(chr_data$scaf[i])
    contig_length <- contig_lengths[contig_name]  # Get contig length from `.fai`
    orientation <- chr_data$freqStrand[i]          # Get contig orientation
    count = chr_data$count[i]
    
    if(chr_data$scaf[i]%in%hcontochr){
     
    # Generate a 'W' row for the contig placement
    agp_row <- data.frame(
      object = paste0(chr, '_', i),                                 # Chromosome name
      object_beg = current_start,                  # Start of the contig
      object_end = current_start + contig_length - 1,  # End of the contig
      part_number = part_number,                   # Increment part number
      component_type = "W",                        # W indicates contig placement
      component_id = as.character(contig_name),    # Contig name as a character (W row)
      component_beg = as.character(1),                           # Start of contig in FASTA (1-based indexing)
      component_end = as.character(contig_length),               # End of contig in FASTA
      orientation = orientation,                   # Orientation of the contig
      stringsAsFactors = FALSE
    )
    
    # Append this contig row to the list
    agp_rows <- append(agp_rows, list(agp_row))
    
    part_number_nocontig=i+1
    
    }else{
    
        # Generate a 'W' row for the contig placement
    agp_row <- data.frame(
      object = paste0(chr, '_', part_number_nocontig),                                 # Chromosome name
      object_beg = current_start,                  # Start of the contig
      object_end = current_start + contig_length - 1,  # End of the contig
      part_number = part_number,                   # Increment part number
      component_type = "W",                        # W indicates contig placement
      component_id = as.character(contig_name),    # Contig name as a character (W row)
      component_beg = as.character(1),                           # Start of contig in FASTA (1-based indexing)
      component_end = as.character(contig_length),               # End of contig in FASTA
      orientation = orientation,                   # Orientation of the contig
      stringsAsFactors = FALSE
    )
    
    # Append this contig row to the list
    agp_rows <- append(agp_rows, list(agp_row))

    ## make a gap row to continue it...
    
    current_start <- current_start + contig_length
    part_number <- part_number + 1
    # Add a 'N' scaffold gap row (unless it's the last row of the chromosome)
 #   if (i < nrow(chr_data)) {
 
    if(i < nrow(chr_data[chr_data$Chr==chr,])){
      gap_size <- 100  # Set the gap size (adjust as needed)

      # Explicitly assign gap length as a character
      gap_row <- data.frame(
        object = paste0(chr, '_', part_number_nocontig),                              # Chromosome name
        object_beg = current_start,               # Start of the gap
        object_end = current_start + gap_size - 1,  # End of the gap
        part_number = part_number,                # Increment the part number
        component_type = "U",                     # N indicates a gap
        component_id = as.character(gap_size),    # Gap length as a character
        component_beg = "contig",               # Scaffold gap type
        component_end = "no",                    # Linkage evidence
        orientation = 'na',                         # Not applicable for gaps
        stringsAsFactors = FALSE
      )
      
      # Append the gap row to the list
      agp_rows <- append(agp_rows, list(gap_row))
      
      # Update start position for the next contig after the gap
      current_start <- current_start + gap_size
      part_number <- part_number + 1
#    part_number_nocontig <- part_number_nocontig +1
    
    }
    

        
    # Update start position for the next row (after the contig ends)

### NOW, NEED TO DECIDE WHETHER IT'S A HOMOLOG THAT COMES NEXT...
   

    
    


    }
  }
}
}









# Combine all rows into a single data frame
agp_df <- do.call(rbind, agp_rows)  # Merges rows for all chromosomes


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





