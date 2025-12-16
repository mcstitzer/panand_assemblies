library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())



process_anchors_to_karyotype <- function(filepath, alleles=1, color_palette=muted_colors, minBlock=10, title='', refChrs=c(paste0('Chr0', 1:9), 'Chr10'), queryChrs='') {
  # Load data
  data <- read.table(filepath, header = TRUE)
  data <- data[data$gene != 'interanchor', ]
  
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
  data$refChr <- factor(data$refChr, levels = c(paste0('Chr0', 1:9), 'Chr10'))
  
  # Reverse strand calculations
  data <- data %>%
    arrange(freqRef, referenceStart, queryStart)
  data$queryChr <- factor(data$queryChr, levels = rev(data$queryChr[!duplicated(data$queryChr)]))
  data$revQueryStart <- data$queryStart
  data$revQueryStart[data$freqStrand == '-'] <- abs(data$queryStart - data$maxChr)[data$freqStrand == '-']
  
  
  data=data[ data$refChr %in% names(color_palette) & data$refChr%in%refChrs & data$queryChr%in%queryChrs, ]
  
  data %>%
  # 1. Summarize to get the segment length (maxChr) and count (n)
  group_by(refChr, queryChr, maxChr) %>%
  summarize(n = n(), .groups = "drop") %>%
  # 2. Assign each queryChr uniquely to the refChr with the longest segment
  group_by(queryChr) %>%
  slice_max(maxChr, n = 1, with_ties = FALSE) %>%
  # 3. From the unique assignments, keep only the six longest assignments for each refChr
  group_by(refChr) %>%
  slice_max(maxChr, n = alleles, with_ties = FALSE) %>%
  # 4. Final arrangement for clear presentation
  arrange(refChr, desc(maxChr))
  
  
  
  

# --- Define the Top 6 Unique Assignments (The Y-Axis Skeleton) ---
# This part is the same: it defines which queryChr belongs to which refChr (the facet)
# and its Y-position (rank 1-6).
assignment_skeleton <- data %>%
  group_by(refChr, queryChr, maxChr) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(queryChr) %>%
  slice_max(maxChr, n = 1, with_ties = FALSE) %>%
  group_by(refChr) %>%
  slice_max(maxChr, n = alleles, with_ties = FALSE) %>%
  arrange(refChr, desc(maxChr)) %>%
  mutate(rank = row_number()) %>% # Local rank (1-6)
  ungroup() %>%
  # FIX: Ensure maxChr is included here so it's available for the next step
  select(refChr, queryChr, maxChr, rank)

# --- Determine the Correct Facet Order ---
refChr_order <- assignment_skeleton %>%
  filter(rank == 1) %>% # Keep only the longest (Rank 1) segment for each refChr group
  # The maxChr column is now available here to perform the ordering:
  arrange(desc(maxChr)) %>% # Sort these Rank 1 segments by length, descending
  pull(refChr) # Extract the refChr names in the desired order

# --- Prepare Plot Data for Points (using the new order) ---
plot_data_points <- data %>%
  # Join with the skeleton (which now contains maxChr and rank)
  inner_join(assignment_skeleton, by = "queryChr", suffix = c(".x", ".y")) %>%
  rename(refChr_original = refChr.x, refChr_assigned = refChr.y) %>%
  # Convert refChr_assigned to a factor with the derived order
  mutate(refChr_assigned = factor(refChr_assigned, levels = refChr_order)) %>%
  select(refChr_assigned, refChr_original, queryChr, rank, queryStart)
  
  
# --- Define Label Data ---
label_data <- assignment_skeleton

#pdf('~/transfer/blagur_karyo_try.pdf',14,4)
ggplot() +
  # 1. Plot the points
  geom_point(
    data = plot_data_points,
    aes(
      x = queryStart, # X-position is the start of the alignment on the query contig
      y = rank,       # Y-position is the rank (1-6)
      color = refChr_original # Color is the original refChr (to show translocations)
    ),
    size = 2.5,
    alpha = 0.8 # Use transparency if points overlap a lot
  ) +
  # 3. Use facet_grid for the assigned refChr
  facet_grid(. ~ refChr_assigned, scales = "free_x") +
  # 4. Apply custom colors
  scale_color_manual(values = muted_colors) +
  # 5. Customize the Y-Axis
  scale_y_continuous(
    name = "Allele",
    breaks = 1:alleles,
    trans = "reverse", # Rank 1 (longest) is at the top
    expand = expansion(mult = c(0.1, 0.1))
  ) +
  # 6. Customize X-axis
  scale_x_continuous(
    name = "Position (bp)",
    labels = scales::unit_format(unit = "M", scale = 1e-6)
  ) +
  # 7. Labels and theme
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linetype = "dotted", color = "grey90"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom"
  )
  
#  dev.off()
  
  
 
}




plot_contigs_by_length <- function(filepath, color_palette=muted_colors, minBlock=10, minContigLengthMb=1) {
  
  # --- 1. Initial Data Loading and Processing ---
  data <- read.table(filepath, header = TRUE)
  data <- data[data$gene != 'interanchor', ]
  
  # Calculate stats (blockLength, freqStrand, maxChr, freqRef)
  data <- data %>%
    group_by(blockIndex) %>%
    mutate(blockLength = dplyr::n()) %>%
    group_by(queryChr) %>%
    mutate(freqStrand = names(which.max(table(strand))),
           # maxChr represents the total length of the queryChr segment involved in the alignment
           maxChr = max(queryStart), 
           freqRef = names(which.max(table(refChr)))) %>%
    ungroup()
  
  # Filter data based on minimum block size
  data <- data[data$blockLength > minBlock, ]
  
  # --- 2. Determine Contig Order and Filter by Length ---
  
  # Find the total length of each queryChr and rank them
  contig_order_data <- data %>%
    group_by(queryChr) %>%
    summarise(total_length = max(maxChr, na.rm = TRUE)) %>%
    ungroup() %>%
    # Filter for contigs greater than 1 Mb (user-defined threshold)
    filter(total_length >= minContigLengthMb * 1e6) %>%
    arrange(desc(total_length)) %>%
    # Create the X-axis factor levels
    mutate(queryChr_f = factor(queryChr, levels = queryChr))
  
  # --- 3. Prepare Final Plot Data ---
  
  plot_data_final <- data %>%
    # Keep only the query contigs that passed the length filter
    inner_join(contig_order_data %>% select(queryChr, queryChr_f), by = "queryChr") %>%
    # Assign a unique numeric Y-position to each queryChr for plotting
    mutate(y_pos = as.numeric(queryChr_f))
  
  # --- 4. Prepare Contig Line Data (Vertical lines) ---
  
  # We need one vertical line per queryChr, positioned at its numeric rank
  contig_line_data <- plot_data_final %>%
    distinct(queryChr_f, y_pos)
  
  # --- 5. Generate the Plot ---
  
  ggplot() +
    # A. Vertical line for each contig (X-axis position, Y-axis span)
    geom_segment(
      data = contig_line_data,
      aes(
        x = y_pos, 
        xend = y_pos, 
        y = 0, # Start line at Y=0 (minimum position)
        yend = max(plot_data_final$queryStart, na.rm = TRUE) * 1.05 # End line slightly beyond max position
      ),
      linewidth = 0.5,
      color = "grey80"
    ) +
    # B. Alignment Points
    geom_point(
      data = plot_data_final,
      aes(
        x = y_pos, # X-position is the unique contig rank
        y = queryStart, # Y-position is the alignment position
        color = refChr # Color by the original refChr match
      ),
      size = 2,
      alpha = 0.8
    ) +
    # C. Aesthetics and Labels
    scale_color_manual(values = color_palette) +
    scale_x_continuous(
      name = "Scaffolds longer than 1Mb",
      breaks = 1:n_distinct(contig_order_data$queryChr),
      labels = contig_order_data$queryChr # Label the vertical lines with contig names
    ) +
    scale_y_continuous(
      name = "Position (bp)",
      labels = scales::unit_format(unit = "M", scale = 1e-6),
      expand = expansion(mult = c(0, 0.05)) # Expand y-axis a bit
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8), # Rotate X labels
      panel.grid.major.x = element_blank(), # Remove major X grid lines
      legend.position = "bottom"
    )
}


## get data
process_contigs_data_local <- function(filepath, minBlock=10, minContigLengthMb=1, haploidFilter=NA) {
  
  # --- 1. Initial Data Loading and Processing (Same as before) ---
  data <- read.table(filepath, header = TRUE)
  data <- data[data$gene != 'interanchor', ]
  
  if(!is.na(haploidFilter)){
  
    data <- data[data$queryChr %in% haploidFilter, ]
 
  }
  
  
  data <- data %>%
    group_by(blockIndex) %>%
    mutate(blockLength = dplyr::n()) %>%
    group_by(queryChr) %>%
    mutate(freqStrand = names(which.max(table(strand))),
           maxChr = max(queryStart),
           freqRef = names(which.max(table(refChr)))) %>%
    ungroup()
  
  data <- data[data$blockLength > minBlock, ]
  
  # --- 2. Determine Contig Order and Filter by Length (LOCAL) ---
  
  contig_order_data_local <- data %>%
    group_by(queryChr) %>%
    summarise(total_length = max(maxChr, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(total_length >= minContigLengthMb * 1e6) %>%
    arrange(desc(total_length)) %>%
    # Create the X-axis factor levels (LOCAL order for this file)
    mutate(x_id_local = row_number(),
           queryChr_f = factor(queryChr, levels = queryChr))
  
  # --- 3. Prepare Final Plot Data ---
  
  plot_data_final <- data %>%
    inner_join(contig_order_data_local %>% select(queryChr, queryChr_f, x_id_local), by = "queryChr") %>%
    # Assign local X-position (x_id_local)
    mutate(x_pos = x_id_local,
           sample_file = basename(filepath))
    
  # --- 4. Return the Processed Data ---
  return(plot_data_final)
}



######
muted_colors <- c("#b34064", "#459abf", "#68b488", "#b3ac40", "#8d4cba", 
                  "#bf9140", "#ae459a", "#99aabf", "#409f90", "#405973")

names(muted_colors) <- c(paste0('Chr0', 1:9), 'Chr10')

ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')

taxonnames=c("Zea mays ssp. parviglumis TIL11", "Zea mays ssp. parviglumis TIL01", "Zea mays ssp. mays B73v5", "Zea mays ssp. mexicana TIL25", "Zea mays ssp. mexicana TIL18", "Zea mays ssp. huehuetenangensis", 
"Zea luxurians", "Zea nicaraguensis", "Zea diploperennis Momo", "Zea diploperennis Gigi", "Tripsacum zoloptense", "Tripsacum dactyloides FL", "Tripsacum dactyloides Southern Hap2", 
"Tripsacum dactyloides Northern Hap2", "Tripsacum dactyloides KS", "Tripsacum dactyloides tetraploid", "Urelytrum digitatum", "Vossia cuspidata", "Rhytachne rottboellioides", "Rottboellia tuberculosa", 
"Hemarthria compressa", "Elionurus tripsacoides", "Schizachyrium scoparium", "Schizachyrium microstachyum", "Anatherum virginicum", "Andropogon chinensis", "Andropogon gerardi", 
"Cymbopogon refractus", "Cymbopogon citratus", "Heteropogon contortus", "Themeda triandra", "Bothriochloa laguroides", "Pogonatherum paniceum", "Sorghum bicolor", 
"Ischaemum rugosum", "Sorghastrum nutans", '"Andropogon" burmanicus', "Thelepogon elegans", "Chrysopogon serrulatus", "Paspalum vaginatum")
names(taxonnames)=c("zTIL11", "zTIL01", "zmB735", "zTIL25", "zTIL18", "zmhuet", 
"zluxur", "znicar", "zdmomo", "zdgigi", "tzopol", "tdacs1", "tdacs2", 
"tdacn2", "tdacn1", "tdactm", "udigit", "vcuspi", "rrottb", "rtuber", 
"hcompr", "etrips", "sscopa", "smicro", "avirgi", "achine", "agerar", 
"crefra", "ccitra", "hconto", "ttrian", "blagur", "ppanic", "sbicol", 
"irugos", "snutan", "atenui", "telega", "cserru", "pvagin")



library(purrr) # Or use lapply/do.call(rbind, ...)
library(ggplot2)

# Define your file list
file_paths <- c(    'etripsHap1aggressive-Pv-6',     
    'smicroHap1aggressive-Pv-2',   'rtuberHap1aggressive-Pv-2',   'ttrianDIPLOIDqv_aggressive-Pv-6', 
'tdacs1-Pv-2','tdacn1-Pv-2',  'avirgi-Pv-2',  'udigitHap1aggressive-Pv-6', 'telegaHap1aggressive-Pv-6', 'sbicol-Pv-2',
  '../scaffolding/crefraCHR-Pv-2', '../scaffolding/ppanicCHR-Pv-2',
'zdgigi-Pv-4','zdmomo-Pv-4','zluxur-Pv-4','zmB735-Pv-4','zmhuet-Pv-4','znicar-Pv-4','zTIL01-Pv-4','zTIL11-Pv-4','zTIL18-Pv-4','zTIL25-Pv-4',
'rrottbHap1chr18aggressive-Pv-6'
)

haploidifyfiles=c('../scaffolding/cserruCHR-Pv-2','../scaffolding/hcontoCHR-Pv-4',
'blaguraggressive-Pv-6','sscopaBothHapsaggressive-Pv-6',"achineBothHapsaggressive-Pv-6",
'irugosBothHaps_aggressivecorrection-Pv-4','ttrianBothHapsUnfilteredaggressive-Pv-6','snutanaggressive-Pv-4')

cserrucontigs=paste0('chr',1:10,'_1')
hcontocontigs=c(paste0('chr',1:10,'_2'), paste0('chr',1:10,'_4'))
blagurcontigs=c(paste0('scaffold_', c(16,4,18,1,20,23,22,24,5,31,29,25,21,9,36,38,41,44,7,42,47,50,10,32,53,33,55,56,60,58)))
sscopacontigs=c(paste0('scaffold_',c(1,3,6,7,10,11,13,15,17,19,21,23,25,27,29,33,28,35,37,40)))
achinecontigs=c(paste0('scaffold_',c(1,3,5,7,9,17,15,11,21,19,27,25,14,37,30,33,31,36,24,39)))
irugoscontigs=c(paste0('scaffold_',c(16,4,5,7,9,11,15,1,18,21)))
ttriancontigs=c(paste0('scaffold_',c(9,19,3,11,1,17,14,15,8,5)))
snutancontigs=c(paste0('scaffold_',c(8,2,25,6,28,31,4,12,13,14,19,33,20,24,18,11,22,38,41,40)))

haploid_contigs=list(cserrucontigs,hcontocontigs,blagurcontigs,sscopacontigs,achinecontigs,irugoscontigs,ttriancontigs, snutancontigs)

# 1. Process all files into a list of data frames
list_of_data <- map(file_paths, ~process_contigs_data_local(filepath = .x, minContigLengthMb = 15))
combined_data_local <- bind_rows(list_of_data)

haploid_list=lapply(1:8, function(x){
                        process_contigs_data_local(filepath=haploidifyfiles[x], minContigLengthMb=15, haploidFilter=haploid_contigs[[x]])
})
combined_hap=bind_rows(haploid_list)

combined_data_local=rbind(combined_data_local, combined_hap)


# 3. Create the unique X-axis label data for each facet
#    This table ensures each facet gets the correct, locally ordered labels (queryChr names)
local_label_data <- combined_data_local %>%
  distinct(sample_file, x_id_local, queryChr) %>%
  arrange(sample_file, x_id_local)
  
 combined_data_local$six=substr(combined_data_local$sample_file,1,6) 
combined_data_local$six=factor(combined_data_local$six, levels=names(taxonnames))
  
  
pdf('~/transfer/all_combined_karyotype.pdf', 10,20)
# 4. Generate the Plot with Shared Scales
ggplot() +
  # Use geom_vline to show where each contig begins (X-axis spacing is equal)
  geom_vline(
    data = combined_data_local %>% distinct(sample_file, x_id_local),
    aes(xintercept = x_id_local),
    color = "grey85",
    linewidth = 0.5
  ) +
    geom_vline(xintercept=c(10,20,30), color='slategray', linetype='dashed')+

  # Alignment Points
  geom_point(
    data = combined_data_local,
    aes(
      x = x_id_local,            # X-position is the unique contig ID (equal spacing)
      y = queryStart,      # Y-position is the alignment position (shared scale)
      color = refChr      # Color by the original refChr match
    ),
    size = 6,
    alpha = 0.7, shape='-'
  ) +
  # Facet by sample to show the comparison
#  facet_grid(sample_file ~ ., scales = "free_x") + # Free X scale can be tricky for comparison
  facet_grid(six ~ ., scales = "free_x") + # Free X scale can be tricky for comparison

  # OPTION A: Shared X-scale (Better for comparison)
  # facet_grid(sample_file ~ ., scales = "fixed") + 
  
  # Aesthetics and Labels
  scale_color_manual(values = muted_colors) + # Ensure muted_colors is defined
  scale_y_continuous(
    name = "Position (bp)",
    labels = scales::unit_format(unit = "M", scale = 1e-6)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linetype = "dotted", color = "grey90"), # Make Y-grid lighter/more transparent
    legend.position = "bottom",
    strip.text.y = element_text(angle = 0) # Set angle to 0 degrees (horizontal)
  )
dev.off()



# Example of usage
#process_anchors_to_karyotype('../syntenic_anchors/anchors/agerar-Pv-6', alleles=6, minBlock=10, refChrs='Chr01')

genotype='blagur'
ploidy=6
alleles=6
pdf(paste0('~/transfer/', genotype, '_karyotype_try.pdf'), 14,4)
process_anchors_to_karyotype(paste0(genotype, 'aggressive-Pv-', as.character(as.numeric(ploidy))), alleles=alleles, minBlock=10)
plot_contigs_by_length(paste0(genotype, 'aggressive-Pv-', as.character(as.numeric(ploidy))), minBlock=10)
dev.off()



genotype='achineBothHaps'
ploidy=6
alleles=4
pdf(paste0('~/transfer/', genotype, '_karyotype_try.pdf'), 14,4)
process_anchors_to_karyotype(paste0(genotype, 'aggressive-Pv-', as.character(as.numeric(ploidy))), alleles=alleles, minBlock=10)
plot_contigs_by_length(paste0(genotype, 'aggressive-Pv-', as.character(as.numeric(ploidy))), minBlock=10)
dev.off()

genotype='irugosBothHaps'
ploidy=4
alleles=4
pdf(paste0('~/transfer/', genotype, '_karyotype_try.pdf'), 14,4)
process_anchors_to_karyotype(paste0(genotype, '_aggressivecorrection-Pv-', as.character(as.numeric(ploidy))), alleles=alleles, minBlock=10)
plot_contigs_by_length(paste0(genotype, '_aggressivecorrection-Pv-', as.character(as.numeric(ploidy))), minBlock=10)

dev.off()


genotype='snutan'
ploidy=4
alleles=4
pdf(paste0('~/transfer/', genotype, '_karyotype_try.pdf'), 14,4)
process_anchors_to_karyotype(paste0(genotype, 'aggressive-Pv-', as.character(as.numeric(ploidy))), alleles=alleles, minBlock=10)
plot_contigs_by_length(paste0(genotype, 'aggressive-Pv-', as.character(as.numeric(ploidy))), minBlock=10)

dev.off()


genotype='smicroHap1'
ploidy=2
alleles=1
pdf(paste0('~/transfer/', genotype, '_karyotype_try.pdf'), 14,4)
process_anchors_to_karyotype(paste0(genotype, 'aggressive-Pv-', as.character(as.numeric(ploidy))), alleles=alleles, minBlock=10)
plot_contigs_by_length(paste0(genotype, 'aggressive-Pv-', as.character(as.numeric(ploidy))), minBlock=10)

dev.off()


genotype='etripsHap1'
ploidy=6
alleles=2
pdf(paste0('~/transfer/', genotype, '_karyotype_try.pdf'), 14,4)
process_anchors_to_karyotype(paste0(genotype, 'aggressive-Pv-', as.character(as.numeric(ploidy))), alleles=alleles, minBlock=10)
plot_contigs_by_length(paste0(genotype, 'aggressive-Pv-', as.character(as.numeric(ploidy))), minBlock=10)

dev.off()

genotype='rtuberHap1'
ploidy=2
alleles=2
pdf(paste0('~/transfer/', genotype, '_karyotype_try.pdf'), 14,4)
process_anchors_to_karyotype(paste0(genotype, 'aggressive-Pv-', as.character(as.numeric(ploidy))), alleles=alleles, minBlock=10)
plot_contigs_by_length(paste0(genotype, 'aggressive-Pv-', as.character(as.numeric(ploidy))), minBlock=10)

dev.off()


genotype='ttrianDIPLOIDqv_'
ploidy=6
alleles=2
pdf(paste0('~/transfer/', genotype, '_karyotype_try.pdf'), 14,4)
process_anchors_to_karyotype(paste0(genotype, 'aggressive-Pv-', as.character(as.numeric(ploidy))), alleles=alleles, minBlock=10)
plot_contigs_by_length(paste0(genotype, 'aggressive-Pv-', as.character(as.numeric(ploidy))), minBlock=10)

dev.off()

genotype='tdacs1'
ploidy=2
alleles=2
pdf(paste0('~/transfer/', genotype, '_karyotype_try.pdf'), 14,4)
process_anchors_to_karyotype(paste0(genotype, '-Pv-', as.character(as.numeric(ploidy))), alleles=alleles, minBlock=10)
plot_contigs_by_length(paste0(genotype, '-Pv-', as.character(as.numeric(ploidy))), minBlock=10)

dev.off()

genotype='zTIL01'
ploidy=4
alleles=1
pdf(paste0('~/transfer/', genotype, '_karyotype_try.pdf'), 14,4)
process_anchors_to_karyotype(paste0(genotype, '-Pv-', as.character(as.numeric(ploidy))), alleles=alleles, minBlock=10)
plot_contigs_by_length(paste0(genotype, '-Pv-', as.character(as.numeric(ploidy))), minBlock=10)

dev.off()

genotype='sscopaBothHaps'
ploidy=6
alleles=4
pdf(paste0('~/transfer/', genotype, '_karyotype_try.pdf'), 14,4)
process_anchors_to_karyotype(paste0(genotype, 'aggressive-Pv-', as.character(as.numeric(ploidy))), alleles=alleles, minBlock=10)
plot_contigs_by_length(paste0(genotype, 'aggressive-Pv-', as.character(as.numeric(ploidy))), minBlock=10)

dev.off()


genotype='telegaHap1'
ploidy=6
alleles=2
pdf(paste0('~/transfer/', genotype, '_karyotype_try.pdf'), 14,4)
process_anchors_to_karyotype(paste0(genotype, 'aggressive-Pv-', as.character(as.numeric(ploidy))), alleles=alleles, minBlock=10)
plot_contigs_by_length(paste0(genotype, 'aggressive-Pv-', as.character(as.numeric(ploidy))), minBlock=10)

dev.off()


genotype='rrottbHap1chr18'
ploidy=6
alleles=2
pdf(paste0('~/transfer/', genotype, '_karyotype_try.pdf'), 14,4)
process_anchors_to_karyotype(paste0(genotype, 'aggressive-Pv-', as.character(as.numeric(ploidy))), alleles=alleles, minBlock=10)
plot_contigs_by_length(paste0(genotype, 'aggressive-Pv-', as.character(as.numeric(ploidy))), minBlock=10)

dev.off()


#../scaffolding/cserruCHR-Pv-2
genotype='cserru'
ploidy=2
alleles=2
pdf(paste0('~/transfer/', genotype, '_karyotype_try.pdf'), 14,4)
process_anchors_to_karyotype(paste0('../scaffolding/cserruCHR', '-Pv-', as.character(as.numeric(ploidy))), alleles=alleles, minBlock=10)
plot_contigs_by_length(paste0('../scaffolding/cserruCHR', '-Pv-', as.character(as.numeric(ploidy))), minBlock=10)

dev.off()


# '../scaffolding/hcontoCHR-Pv-4', '../scaffolding/crefraCHR-Pv-2', '../scaffolding/ppanicCHR-Pv-2',
genotype='hconto'
ploidy=4
alleles=4
pdf(paste0('~/transfer/', genotype, '_karyotype_try.pdf'), 14,4)
process_anchors_to_karyotype(paste0('../scaffolding/hcontoCHR', '-Pv-', as.character(as.numeric(ploidy))), alleles=alleles, minBlock=10)
plot_contigs_by_length(paste0('../scaffolding/hcontoCHR', '-Pv-', as.character(as.numeric(ploidy))), minBlock=10)

dev.off()
genotype='crefra'
ploidy=2
alleles=2
pdf(paste0('~/transfer/', genotype, '_karyotype_try.pdf'), 14,4)
process_anchors_to_karyotype(paste0('../scaffolding/crefraCHR', '-Pv-', as.character(as.numeric(ploidy))), alleles=alleles, minBlock=10)
plot_contigs_by_length(paste0('../scaffolding/crefraCHR', '-Pv-', as.character(as.numeric(ploidy))), minBlock=10)

dev.off()
genotype='ppanic'
ploidy=2
alleles=2
pdf(paste0('~/transfer/', genotype, '_karyotype_try.pdf'), 14,4)
process_anchors_to_karyotype(paste0('../scaffolding/ppanicCHR', '-Pv-', as.character(as.numeric(ploidy))), alleles=alleles, minBlock=10)
plot_contigs_by_length(paste0('../scaffolding/ppanicCHR', '-Pv-', as.character(as.numeric(ploidy))), minBlock=10)

dev.off()
