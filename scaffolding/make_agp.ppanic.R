# Load libraries (if necessary)
library(dplyr)
library(readr)
library(tibble)

# Step 1: Load the fai file
fai <- read_delim("fais/Pi-Clark-DRAFT-PanAnd-1.0.fasta.fai", delim = "\t", col_names = FALSE)
colnames(fai) <- c("contig", "length", "other1", "other2", "other3") # Adjust based on the `.fai` file format

fai <- as.data.frame(fai)

contig_lengths <- fai %>%
  dplyr::select(contig, length) %>%
  deframe() # Converts to a named vector with `contig` as keys and `length` as values

# Step 2: Load the data frame with directions, chromosomes, and ordering
# Simulated input dataset
filepath='../syntenic_anchors/anchors/ppanic-Pv-2'
#filepath='~/Downloads/crefra-Av-2'

color_palette=muted_colors
minBlock=10
title=''
refChrs=c(paste0('Chr0', 1:9), 'Chr10')
#refChrs=c(paste0('chr', 1:10))
queryChrs=''
queryChrtoFlip=''
ylabelspecies=''
ploidy=''


# process_anchors_to_dotplot_SUPP(filepath = filepath, minBlock=3, queryChrtoFlip = 'chr9', ylabelspecies = 'Pogonatherum paniceum',  ploidy='Diploid', color_palette = muted_colors)+
#  new_scale_color()+geom_hline(data=tr[tr$queryChr%in%unique(data$queryChr),], aes(yintercept=V4/1e6, color=rp), lty='dotted')+scale_color_brewer(palette='Set1')+theme(legend.position='bottom')
# #process_anchors_to_dotplot_SUPP <- function(filepath, color_palette=muted_colors, minBlock=10, title='', refChrs=c(paste0('Chr0', 1:9), 'Chr10'), queryChrs='', 
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



# ordered_data <- data %>%
#   arrange(freq, Chr, scaf) # Adjust sorting based on your actual requirements

filtered_data <- ordered_data %>%
  distinct(Chr, scaf, .keep_all = TRUE) 

## need to fix a misrepair in crefra - how to do programaticaly??? arising because of flip??
#filtered_data[c(8,9),]=filtered_data[c(9,8),]


## need to make afusion

filtered_data$Chr=as.character(filtered_data$Chr)
filtered_data$Chr[filtered_data$scaf=='ctg_16986']='chr1-7fusion' ## this is t2t so leave it be
filtered_data$Chr[filtered_data$scaf%in%c('ctg_42969', 'ctg_64743')]='chr7-1fusion' ## this is t2t so leave it be

## and make a forced order
filtered_data$Chr[filtered_data$scaf %in%c('ctg_150187', 'ctg_171461')]='chr2-4fusion'

## manually call finished contigs
finishedchr=c('ctg_1',#chr10
            'ctg_49407'#chr9
            )
filtered_data$Chr[filtered_data$Chr%in%filtered_data$Chr[filtered_data$scaf%in%finishedchr]]=NA
filtered_data$Chr[filtered_data$scaf%in%finishedchr]=finishedchr ## this might be out of order but shouldn't matter after renaming


agp_rows <- list()
# Loop over unique Chromosomes


# Process each chromosome individually
for (chr in unique(filtered_data$Chr)) {
  # Filter data for the current chromosome
  chr_data <- filter(filtered_data, Chr == chr)
  chrname=
  current_start <- 1  # Reset position at the start of every chromosome
  part_number <- 1    # Reset part number at the start of every chromosome

  # Generate rows for the current chromosome
  for (i in seq_len(nrow(chr_data))) {
    # Extract details for this contig
    contig_name <- as.character(chr_data$scaf[i])
    contig_length <- contig_lengths[contig_name]  # Get contig length from `.fai`
    orientation <- chr_data$freqStrand[i]          # Get contig orientation
    
    # Generate a 'W' row for the contig placement
    agp_row <- data.frame(
      object = chr,                                 # Chromosome name
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
    
    # Update start position for the next row (after the contig ends)
    current_start <- current_start + contig_length
    part_number <- part_number + 1
    
    # Add a 'N' scaffold gap row (unless it's the last row of the chromosome)
    if (i < nrow(chr_data) & !contig_name%in%finishedchr) {
      gap_size <- 100  # Set the gap size (adjust as needed)

      # Explicitly assign gap length as a character
      gap_row <- data.frame(
        object = chr,                              # Chromosome name
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
    }
  }
}

# Combine all rows into a single data frame
agp_df <- do.call(rbind, agp_rows)  # Merges rows for all chromosomes

agp_df <- agp_df %>%
  group_by(object) %>%
  summarise(sum_component_end = sum(as.numeric(component_end[component_end!='no']), na.rm=T)) %>%
  arrange(desc(sum_component_end)) %>%
  mutate(new_object = paste0("chr", row_number())) %>%
  dplyr::select(object, new_object) %>%
  right_join(agp_df, by = "object") %>%
  mutate(object = new_object) %>%
  dplyr::select(-new_object)


# Step 1: Identify unplaced scaffolds
unplaced_scaffolds <- setdiff(names(contig_lengths), agp_df$component_id)

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
final_agp_df <- bind_rows(agp_df, unplaced_df)


# Save AGP file
write_delim(final_agp_df, "ppanic.agp", delim = "\t", col_names = FALSE)

# AGP file is now generated and saved as "output.agp"
## on cbsu, run
##./agp_to_validationplots.sh Pi-Clark-DRAFT-PanAnd-1.0 ppanicCHR ppanic.agp 1

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
# 
# 



