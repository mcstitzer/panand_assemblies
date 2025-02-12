library(rtracklayer)

knobbin=read.table('../repeats/knobs/mcclintock_knob_bins.txt', header=F, colClasses = c('character', 'character'))
bins=import.gff3('../repeats/knobs/knob_bin_positions_B73v3.gff3')

binsK=bins[bins$Name%in%knobbin$V2,]
binsK$knob=knobbin$V1[match(binsK$Name, knobbin$V2)]


bins$knob=knobbin$V1[match(bins$Name, knobbin$V2)]



df_bins <- data.frame(bins) %>%
  mutate(knob_status = if_else(is.na(knob), "No Knob", "Knob"))

# 3. Basic ggplot: plot each chromosome as a segment from start to end
#    Color the segment by whether knob is present (knob_status)
ggplot(df_bins, 
       aes(x = start,
           xend = end,
           y = factor(seqnames, levels=rev(paste0('Chr',1:10))),
           yend = factor(seqnames, levels=rev(paste0('Chr',1:10))),
           color = knob_status)) +
  geom_segment(size = 4) +
  scale_color_manual(values = c("No Knob" = "gray70", "Knob" = "red")) +
  theme(
    axis.title.y = element_blank(),
    legend.title = element_blank()
  ) +
  labs(x = "Genomic Position")
  
  
  
  
  
  
  

#### also now count knobs
countKnobsCent <- function(filepath, bins, genome='', distance=10000,color_palette=muted_colors, minBlock=20) {
  # Load data
  anchors <- read.table(filepath, header = TRUE)
  anchors <- anchors[anchors$gene != 'interanchor', ]
  
  # Reduce to blocks and calculate stats
  anchors <- anchors %>%
    group_by(blockIndex) %>%
    mutate(blockLength = n()) %>%
    group_by(queryChr) %>%
    mutate(freqStrand = names(which.max(table(strand))),
           maxChr = max(queryStart),
           freqRef = names(which.max(table(refChr))))
  
  # Filter data based on block length
  anchors <- anchors[anchors$blockLength > minBlock, ]
 # Identify fusions and calculate query positions
  fusion_data <- anchors %>%
    filter(blockLength > minBlock) %>%
    group_by(queryChr) %>%
    arrange(queryChr, queryStart) %>% # Ensure sorted order by queryStart
    mutate(
      ref_transition = refChr != lag(refChr, default = refChr[1]), # Detect transitions in refChr
      fusion_point = ifelse(ref_transition, queryStart, NA) # Record query position of the transition
    ) %>%
    filter(ref_transition) %>% # Retain only rows with transitions
    select(queryChr, fusion_point, refChr)%>% tidyr::drop_na(fusion_point)
#     summarize(
#       queryChr = dplyr::first(queryChr),
#       fusion_positions = paste(na.omit(fusion_point), collapse = ", "),
#       involved_refChrs = paste(unique(refChr), collapse = ", ")
#     ) %>% tidyr::unnest(fusion_positions)
#   
  # Print fusion information
  print(fusion_data)

breakpoints_gr=GRanges(seqnames=gsub('c', 'C', fusion_data$queryChr), ranges=IRanges(start=fusion_data$fusion_point, end=fusion_data$fusion_point))

overlaps=findOverlaps(breakpoints_gr,bins, ignore.strand=T)


overlaps=subsetByOverlaps(bins, breakpoints_gr, ignore.strand=T)
print(overlaps)
mcols(bins)[,genome]=bins$Name%in%overlaps$Name
return(bins)
}

spbins=countKnobsCent('../syntenic_anchors/anchors/zdgigi-Pv-4', bins, genome='zdgigi')
spbins=countKnobsCent('../syntenic_anchors/anchors/zdmomo-Pv-4', spbins, genome='zdmomo')
spbins=countKnobsCent('../syntenic_anchors/anchors/zluxur-Pv-4', spbins, genome='zluxur')
spbins=countKnobsCent('../syntenic_anchors/anchors/znicar-Pv-4', spbins, genome='znicar')
spbins=countKnobsCent('../syntenic_anchors/anchors/zmhuet-Pv-4', spbins, genome='zmhuet')
spbins=countKnobsCent('../syntenic_anchors/anchors/zmB735-Pv-4', spbins, genome='zmB735')
spbins=countKnobsCent('../syntenic_anchors/anchors/zTIL01-Pv-4', spbins, genome='zTIL01')
spbins=countKnobsCent('../syntenic_anchors/anchors/zTIL11-Pv-4', spbins, genome='zTIL11')
spbins=countKnobsCent('../syntenic_anchors/anchors/zTIL18-Pv-4', spbins, genome='zTIL18')
spbins=countKnobsCent('../syntenic_anchors/anchors/zTIL25-Pv-4', spbins, genome='zTIL25')

mcols(spbins)$br=rowSums(data.frame(mcols(spbins)[,8:17]))

# 3. Basic ggplot: plot each chromosome as a segment from start to end
#    Color the segment by whether knob is present (knob_status)
ggplot(df_bins, 
       aes(x = start,
           xend = end,
           y = factor(seqnames, levels=rev(paste0('Chr',1:10))),
           yend = factor(seqnames, levels=rev(paste0('Chr',1:10))),
           color = knob_status)) +
  geom_segment(size = 4) +
  scale_color_manual(values = c("No Knob" = "gray70", "Knob" = "red")) +
  theme(
    axis.title.y = element_blank(),
    legend.title = element_blank()
  ) +
  labs(x = "Genomic Position")+
  geom_segment()
  













