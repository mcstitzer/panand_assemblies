library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(ggimage)
library(pdftools)
library(gridExtra)
library(magick)
library(GenomicRanges)
library(grid)

## generalize!
# Set the theme and color palette
theme_set(theme_cowplot())
# Define color palette
#muted_colors <- c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf')
muted_colors <- c("#b34064", "#459abf", "#68b488", "#b3ac40", "#8d4cba", 
                  "#bf9140", "#ae459a", "#99aabf", "#409f90", "#405973")

names(muted_colors) <- c(paste0('Chr0', 1:9), 'Chr10')

ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')



process_anchors_to_dotplot <- function(filepath, color_palette=muted_colors, minBlock=10, title='', refChrs=c(paste0('Chr0', 1:9), 'Chr10'), queryChrs='') {
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
  
  
  # Create the plot
  ggplot(data[ data$refChr %in% names(color_palette) & data$refChr%in%refChrs & data$queryChr%in%queryChrs, ],
         aes(x = referenceStart / 1e6, y = revQueryStart / 1e6, color = refChr)) +
    geom_point() +
    facet_grid(queryChr ~ refChr, scales = 'free', space = 'free') +
    scale_color_manual(values = color_palette) +
    theme(legend.position = 'none') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=9), 
          strip.placement.y = "outside" , 
          strip.text = element_text(size = 8, color = "darkblue", face = "bold"),
          strip.text.y=element_text(angle=0),
          strip.background = element_rect(fill = "lightblue", color = "darkblue", linewidth = 1),
          #      strip.text.y = element_blank(),
  
          axis.text.y=element_text(size=5),
          panel.spacing = unit(0.1, 'lines')
    ) +
    geom_hline(aes(yintercept=maxChr/1e6), lty='dashed', color='gray') +
    ggtitle(title)
}




# Example of usage
process_anchors_to_dotplot('../syntenic_anchors/anchors/agerar-Pv-6', minBlock=20, refChrs='Chr01')

## count rearrangments
diploids=c('cserru', 'irugos', 'sbicol', 'ppanic', 'ttrian', 'crefra', 'avirgi', 'smicro', 'rtuber')
tetraploids=c('snutan', 'hconto', 'ccitra', 'achine', 'sscopa', 'etrips',  'vcuspi')
hexaploids=c('udigit', 'agerar', 'hcompr', 'blagur')

countRearrangements <- function(filepath, color_palette=muted_colors, minBlock=10) {
  # Load data
  data <- read.table(filepath, header = TRUE)
  data <- data[data$gene != 'interanchor', ]
  
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
  
  temp=data %>% filter(blockLength>minBlock) %>% group_by(queryChr)%>% summarize(n=length(unique(refChr)))
  sum(temp$n>1)
}

for(i in diploids){
  print(countRearrangements(Sys.glob(paste0('anchors/', i, '-Pv-*'))))
}

for(i in tetraploids){
  print(countRearrangements(Sys.glob(paste0('anchors/', i, '-Pv-*'))))
}
for(i in hexaploids){
  print(countRearrangements(Sys.glob(paste0('anchors/', i, '-Pv-*')), minBlock=50))
}                              
for(i in c('zmB735')){
  print(countRearrangements(Sys.glob(paste0('anchors/', i, '-Pv-*')), minBlock=50))
}  





process_anchors_to_dotplot('../syntenic_anchors/anchors/avirgi-Pv-2', minBlock=50)


process_anchors_to_dotplot_Tripsacinae <- function(filepath, color_palette=muted_colors, minBlock=10, title='', refChrs=c(paste0('Chr0', 1:9), 'Chr10'), queryChrs='', 
                                              queryChrtoFlip='', ylabelspecies='', pathtokaryotype='', pathtoalluvial='', ploidy='Paleotetraploid', origChrOrientation=F) {
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
  
  ## flip querychr that look funny
  if(queryChrtoFlip!=''){
    data$revQueryStart[data$queryChr%in%queryChrtoFlip]=abs(data$queryStart - data$maxChr)[data$queryChr%in%queryChrtoFlip]
  }
  
  
  # Create the plot
  if(origChrOrientation){
    p=ggplot(data[ data$refChr %in% names(color_palette) & data$refChr%in%refChrs & data$queryChr%in%queryChrs, ],
             aes(x = referenceStart / 1e6, y = queryStart / 1e6, color = refChr))
  }else{
   p= ggplot(data[ data$refChr %in% names(color_palette) & data$refChr%in%refChrs & data$queryChr%in%queryChrs, ],
           aes(x = referenceStart / 1e6, y = revQueryStart / 1e6, color = refChr))
  }
  p= p+
    geom_point(size=0.5) +
    facet_grid(queryChr ~ refChr, scales = 'free', space='free') +
    scale_color_manual(values = color_palette) +
    theme(legend.position = 'none') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    geom_hline(aes(yintercept=maxChr/1e6), lty='dashed', color='gray') +
    geom_vline(xintercept=0, color='gray95')+
    ggtitle(paste0(ylabelspecies,', ', ploidy)) + xlab('P. vaginatum position (Mb)') +
    ylab(paste0(ylabelspecies, ' position (Mb)')) + 
 #   ylab('Position (Mb)') + 
        theme(#strip.text.y = element_text(angle = 0, hjust = 0), 
      strip.placement.y = "outside" , 
      strip.text = element_text(size = 8, color = "darkblue", face = "bold"),
      strip.text.y=element_text(angle=0),
      strip.background = element_rect(fill = "lightblue", color = "darkblue", linewidth = 1),
#      strip.text.y = element_blank(),
      axis.text.x=element_text(size=9),
      axis.text.y=element_text(size=5),
      panel.spacing = unit(0.1, 'lines'),
      plot.title=element_text(color=ploidycolors[ploidy], size=10))+
    scale_x_continuous(breaks = scales::breaks_pretty(n = 2), expand=c(0,0)) +  # Automatically choose 3 breaks for x-axis
    scale_y_continuous(breaks = scales::breaks_pretty(n = 2), expand=c(0,0))    # Automatically choose 3 breaks for y-axis
  
#  strip.background = element_rect(fill = "lightblue", color = "darkblue", size = 1),
#  strip.text = element_text(color = "white", face = "bold"),
#  panel.spacing = unit(1, "lines")
  
  #  p + annotation_custom(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
  p
}




process_anchors_to_dotplot_ZeaTrip <- function(filepath, 
                                               paspalum_tripsacum,
                                               paspalum_maize,
                                               color_palette=muted_colors, minBlock=10, title='', refChrs=c(paste0('Chr0', 1:9), 'Chr10'), queryChrs='', 
                                                   queryChrtoFlip='', ylabelspecies='', pathtokaryotype='', pathtoalluvial='', ploidy='Paleotetraploid', origChrOrientation=F) {
  # Load data
  data <- read.table(filepath, header = TRUE)
  data <- data[data$gene != 'interanchor', ]
  
   pasptrip=read.table(paspalum_tripsacum, header=T)
   pasptrip=pasptrip[pasptrip$gene!='interanchor',]
   paspzea=read.table(paspalum_maize, header=T)
   paspzea=paspzea[paspzea$gene!='interanchor',]
   
   
   dataGRtrip=GRanges(seqnames=data$queryChr, IRanges(start=data$queryStart, end=data$queryEnd))
   dataGRzea=GRanges(seqnames=data$refChr, IRanges(start=data$referenceStart, end=data$referenceEnd))
   
   pasptripGR=GRanges(seqnames=pasptrip$queryChr, IRanges(start=pasptrip$queryStart, end=pasptrip$queryEnd))
   pasptripGR$paspChr=pasptrip$refChr
   paspzeaGR=GRanges(seqnames=paspzea$queryChr, IRanges(start=paspzea$queryStart, end=paspzea$queryEnd))
   paspzeaGR$paspChr=paspzea$refChr
   
   overlaps <- findOverlaps(dataGRtrip, pasptripGR)
   
   # Step 2: Create a vector for paspChr with NA values as default
   paspChr <- rep(NA, length(dataGRtrip))
   
   # Step 3: Assign paspChr values based on the overlaps
   paspChr[queryHits(overlaps)] <- pasptripGR$paspChr[subjectHits(overlaps)]
   
   # Step 4: Add the paspChr column to dataGRtrip
   dataGRtrip$paspChr <- paspChr
   
   data$tripPaspChr=paspChr
   ##zea
   overlaps <- findOverlaps(dataGRzea, paspzeaGR)
   # Step 2: Create a vector for paspChr with NA values as default
   paspChr <- rep(NA, length(dataGRzea))
   # Step 3: Assign paspChr values based on the overlaps
   paspChr[queryHits(overlaps)] <- paspzeaGR$paspChr[subjectHits(overlaps)]
   # Step 4: Add the paspChr column to dataGRtrip
   dataGRzea$paspChr <- paspChr
   data$zeaPaspChr=paspChr
   
   data=data %>% group_by(blockIndex) %>%  mutate(commonTripPasp = as.character(names(which.max(table(tripPaspChr, useNA="always")))[1]),
                                                  commonZeaPasp = as.character(names(which.max(table(zeaPaspChr, useNA="always")))[1]))
   
   data=data %>%
     mutate(commonPaspChr = case_when(
       is.na(commonTripPasp) ~ commonZeaPasp,             # If col1 is NA, use col2
       is.na(commonZeaPasp) ~ commonTripPasp,             # If col2 is NA, use col1
       commonTripPasp == commonZeaPasp ~ commonTripPasp,            # If both are the same, use that value
       commonTripPasp != commonZeaPasp  ~ NA_character_                  # If they are different, use NA
     ))
   
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
  data$refChr <- factor(data$refChr, levels = paste0('chr', 1:10))
  
  # Reverse strand calculations
  data <- data %>%
    arrange(freqRef, referenceStart, queryStart)
  data$queryChr <- factor(data$queryChr, levels = rev(data$queryChr[!duplicated(data$queryChr)]))
  data$revQueryStart <- data$queryStart
  data$revQueryStart[data$freqStrand == '-'] <- abs(data$queryStart - data$maxChr)[data$freqStrand == '-']
  
  ## flip querychr that look funny
  if(queryChrtoFlip!=''){
    data$revQueryStart[data$queryChr%in%queryChrtoFlip]=abs(data$queryStart - data$maxChr)[data$queryChr%in%queryChrtoFlip]
  }
  
  
  
  ##### NEVERMIND! The repeats i care about aren't in b73, since the chromosome collisions occurred sooo sooo long ago
#   ## get repeats to plot knobs and centromeres
#   teb73=import.gff('~/Downloads/Zm-B73-REFERENCE-NAM-5.0.TE.gff3.gz')
#   teb73=teb73[teb73$Classification%in% c('Cent/CentC', 'knob/knob180', 'knob/TR-1', 'rDNA/spacer'),]
#   teb73_split <- split(teb73, teb73$Classification)
#   merged_teb73 <- GRanges()
#   
#   # Loop over each split GRanges object, reduce and concatenate manually
#   for(i in 1:4){
#     reduced_gr <- GenomicRanges::reduce(teb73_split[[i]], min.gapwidth = 10000, ignore.strand=T)
#     reduced_gr$Classification=teb73_split[[i]]$Classification[1]
#     merged_teb73 <- c(merged_teb73, reduced_gr)
#   }
#   
# merged_teb73$refChr=seqnames(merged_teb73)
# merged_teb73$referenceStart=start(merged_teb73)
# merged_teb73$referenceEnd=end(merged_teb73)
# merged_teb73=merged_teb73[seqnames(merged_teb73)%in%paste0('chr',1:10),]
#   
  # Create the plot
  if(origChrOrientation){
    p=ggplot(data[ data$refChr%in%refChrs & data$queryChr%in%queryChrs,],# & data$pasptrip==data$paspzea & !is.na(data$pasptrip), ],
             aes(x = referenceStart / 1e6, y = queryStart / 1e6, color=commonPaspChr))
  }else{
    p= ggplot(data[ data$refChr%in%refChrs & data$queryChr%in%queryChrs,],# & data$pasptrip==data$paspzea & !is.na(data$pasptrip), ],
              aes(x = referenceStart / 1e6, y = revQueryStart / 1e6, color=commonPaspChr))
  }
  p= p+
    geom_point(size=0.5) +
    facet_grid(queryChr ~ refChr, scales = 'free', space='free') +
    scale_color_manual(values = color_palette) +
    theme(legend.position = 'none') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    geom_hline(aes(yintercept=maxChr/1e6), lty='dashed', color='gray') +
    geom_vline(xintercept=0, color='gray95')+
#    geom_rect(data=data.frame(merged_teb73), aes(xmin = start/1e6, xmax = end/1e6, y=-Inf, ymin = -Inf, ymax = Inf, fill = Classification), alpha = 0.9, color=NA)+
 #   scale_fill_manual(values=c("Cent/CentC" = "#ff5733", "knob/knob180" = "#33ff57", "knob/TR-1" = "#3357ff", "rDNA/spacer" = "#ff33a8")) + 
 #   ggtitle(paste0(ylabelspecies,', ', ploidy)) + 
    xlab('Z. m. ssp. mays B73 position (Mb)') +
       ylab(paste0(ylabelspecies, ' position (Mb)')) + 
 #   ylab('Position (Mb)') + 
    theme(strip.text.y = element_text(angle = 0, hjust = 0), 
      strip.placement.y = "outside" , 
      strip.text = element_text(size = 8, color = "darkblue", face = "bold"),
      strip.background = element_rect(fill = "lightblue", color = "darkblue", linewidth = 1),
      axis.text.x=element_text(size=9),
      axis.text.y=element_text(size=5),
      panel.spacing = unit(0.1, 'lines'),
      plot.title=element_text(color=ploidycolors[ploidy], size=10))+
    scale_x_continuous(breaks = scales::breaks_pretty(n = 2), expand=c(0,0)) +  # Automatically choose 3 breaks for x-axis
    scale_y_continuous(breaks = scales::breaks_pretty(n = 2), expand=c(0,0))    # Automatically choose 3 breaks for y-axis
  
  
  
  #  p + annotation_custom(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
  p
}



process_anchors_to_dotplot_FIGURE <- function(filepath, color_palette=muted_colors, minBlock=10, title='', refChrs=c(paste0('Chr0', 1:9), 'Chr10'), queryChrs='', 
                                              queryChrtoFlip='', ylabelspecies='', pathtokaryotype='', pathtoalluvial='', ploidy='') {
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
  
  ## flip querychr that look funny
  if(queryChrtoFlip!=''){
    data$revQueryStart[data$queryChr%in%queryChrtoFlip]=abs(data$queryStart - data$maxChr)[data$queryChr%in%queryChrtoFlip]
  }
  
  ## clean up future facet labels - this is very specific to chr level s
#  data$refLabel=gsub('Chr0', '', data$refChr)
#  data$refLabel=gsub('Chr', '', data$refLabel)
#  data$queryLabel=gsub('chr', '', data$queryLabel)
  
  ## set up the image
#  im=data.frame(refChr='Chr01', queryChr=NA, referenceStart=1, revQueryStart=1, path=pathtokaryotype)
#  im$queryChr=names(which.max(table(data$queryChr[data$refChr=='Chr01'])))
  
  pdf_file <- pathtokaryotype
  image_raster <- image_read_pdf(pdf_file)
  
  # Convert the image to a rasterGrob
  img <- rasterGrob(as.raster(image_raster), 
 #                   x = 0.9, y = 0.12, width = 0.25, height = 0.25, just = c("right", "bottom"))
  x = 0.90, y = 0.15, width = 0.25, height = 0.25, just = c("right", "bottom"))

  # Create a black border (rectangular grob)
#  border <- rectGrob(gp = gpar(col = "black", fill = NA, lwd = 2))  # lwd sets the thickness of the border
  
  # Combine the image and the border into a single grob
#  img <- grobTree(border, img)
  
  
  ## do same for alluvial plot
  pdf_fileA <- pathtoalluvial
  image_rasterA <- image_read_pdf(pdf_fileA)
  
  # Convert the image to a rasterGrob
  imgA <- rasterGrob(as.raster(image_rasterA), 
                    x = 0.1, y = 0.88, width = 0.35, height = 0.15, just = c("left", "top"))
  
  # Combine the image and the border into a single grob
#  imgA <- grobTree(border, imgA)
  
  ystriptextsize=ifelse(length(unique(queryChrs))>50, 6,9)
  
  # Create the plot
  p=ggplot(data[ data$refChr %in% names(color_palette) & data$refChr%in%refChrs & data$queryChr%in%queryChrs, ],
         aes(x = referenceStart / 1e6, y = revQueryStart / 1e6, color = refChr)) +
    geom_point(size=0.5) +
    facet_grid(queryChr ~ refChr, scales = 'free', space='free') +
    scale_color_manual(values = color_palette) +
    theme(legend.position = 'none') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    geom_hline(aes(yintercept=maxChr/1e6), lty='dashed', color='gray') +
    geom_vline(xintercept=0, color='gray95')+
    ggtitle(paste0(ylabelspecies,', ', ploidy)) + xlab('P. vaginatum position (Mb)') +
    ylab(paste0(ylabelspecies, ' position (Mb)')) + 
  #  ylab('Position (Mb)') +
    theme(#strip.text.y = element_text(angle = 0, hjust = 0), 
          strip.placement.y = "outside" , 
          strip.text = element_text(size = 8, color = "darkblue", face = "bold"),
          strip.background = element_rect(fill = "lightblue", color = "darkblue", linewidth = 1),
 #         strip.text.y = element_text(angle=0, size=ystriptextsize),
          strip.text.y=element_blank(),
          strip.background.y=element_blank(),
          axis.text.x=element_text(size=9),
          axis.text.y=element_text(size=5),
          panel.spacing = unit(0.1, 'lines'),
          plot.title=element_text(color=ploidycolors[ploidy], size=10))+
    scale_x_continuous(breaks = scales::breaks_pretty(n = 2), expand=c(0,0)) +  # Automatically choose 3 breaks for x-axis
    scale_y_continuous(breaks = scales::breaks_pretty(n = 2), expand=c(0,0))    # Automatically choose 3 breaks for y-axis
  
 
#  p + annotation_custom(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
  ggdraw(p) + 
    draw_grob(img) + draw_grob(imgA)
  }


av=process_anchors_to_dotplot_FIGURE(filepath = '../syntenic_anchors/anchors/avirgi-Pv-2', minBlock=20, queryChrtoFlip = 'chr9', ylabelspecies = 'A. virginicum', pathtokaryotype = '~/Downloads/avirginicum_karyotype.pdf', pathtoalluvial = '~/Downloads/evoday_alluvial.pdf', ploidy='Diploid')
av
process_anchors_to_dotplot_FIGURE(filepath = '../syntenic_anchors/anchors/snutan-Pv-4', minBlock=50, ylabelspecies = 'S. nutans', pathtokaryotype = '~/Downloads/snutans_karyotype.pdf', pathtoalluvial = '~/Downloads/blagur_alluvial.pdf', ploidy='Tetraploid')
process_anchors_to_dotplot_FIGURE(filepath = '../syntenic_anchors/anchors/hconto-Pv-4', minBlock=50, ylabelspecies = 'H. contortus', pathtokaryotype = '~/Downloads/hcontortus_karyotype.pdf', pathtoalluvial = '~/Downloads/blagur_alluvial.pdf', ploidy='Tetraploid')

bl=process_anchors_to_dotplot_FIGURE(filepath = '../syntenic_anchors/anchors/blagur-Pv-6', minBlock=20, ylabelspecies = 'B. laguroides', pathtokaryotype = '~/Downloads/blaguroides_karyotype.pdf', pathtoalluvial = '~/Downloads/blagur_alluvial.pdf', ploidy='Hexaploid')
bl

td=process_anchors_to_dotplot_Tripsacinae(filepath = '../syntenic_anchors/anchors/tdacs1-Pv-2', minBlock=20, ylabelspecies = 'T. dactyloides FL')
td
zm=process_anchors_to_dotplot_Tripsacinae(filepath = '../syntenic_anchors/anchors/zmB735-Pv-4', minBlock=20, ylabelspecies='Z. mays subsp. mays B73')
zm

tz=plot_grid(td, zm, ncol=2, align='hv', labels=c('C', 'D'))


### new trick - color zea/tripsacum by paspalum chromosomes???? Or 

tzdp=process_anchors_to_dotplot_ZeaTrip('../syntenic_anchors/anchors/tdacs1-zB73v5-2', paspalum_tripsacum = '../syntenic_anchors/anchors/tdacs1-Pv-2', paspalum_maize = '../syntenic_anchors/anchors/zmB735-Pv-4', color_palette=muted_colors, minBlock=10, refChrs = paste0('chr',1:10), queryChrs=paste0('chr',1:18), ylabelspecies='T. dactyloides KS ')
tzdp+ ylab('T. dactyloides KS position (Mb)')


### look up not1 and tga1
process_anchors_to_dotplot_ZeaTrip('anchors/tdacs1-zB73v5-2', paspalum_tripsacum = 'anchors/tdacs1-Pv-2', paspalum_maize = 'anchors/zmB735-Pv-4', color_palette=muted_colors, minBlock=3, refChrs = paste0('chr',c(4)), queryChrs=paste0('chr',1:18), ylabelspecies='T. dactyloides KS ', origChrOrientation = T)+ geom_vline(data=data.frame(refChr='chr4', refStart=c(46654396/1e6,46929803/1e6)), aes(xintercept=refStart))+
  geom_hline(data=data.frame(queryChr='chr17', queryStart=c(75078278,74824806)), aes(yintercept=queryStart/1e6)) + ylab('T. dactyloides KS position (Mb)')
process_anchors_to_dotplot_Tripsacinae(filepath = '../syntenic_anchors/anchors/tdacs1-Pv-2', minBlock=3, queryChrs=c('chr17', 'chr11'),ylabelspecies = 'T. dactyloides FL', origChrOrientation = T) + geom_hline(data=data.frame(queryChr='chr17', queryStart=c(75078278,74824806)), aes(yintercept=queryStart/1e6))+
  geom_vline(data=data.frame(refChr='Chr07', referenceStart=40239652), aes(xintercept = referenceStart/1e6))
process_anchors_to_dotplot_Tripsacinae(filepath = '../syntenic_anchors/anchors/zmB735-Pv-4', minBlock=3, queryChrs=c('chr4', 'chr1'),refChrs = 'Chr07', ylabelspecies = 'Z. mays subsp. mays B73', origChrOrientation = T) + geom_hline(data=data.frame(queryChr='chr4', queryStart=c(46654396,46929803)), aes(yintercept=queryStart/1e6))+
  geom_vline(data=data.frame(refChr='Chr07', referenceStart=40239652), aes(xintercept = referenceStart/1e6))




### fig_chrrearr from count_rearrangements.R
pdf('~/Downloads/rule1.pdf',14,10)
#plot_grid(
plot_grid(plot_grid(av, tz, 
                    plot_grid(tzdp + ylab('T. dactyloides KS position (Mb)'), fig_chrrearr + theme(legend.position='bottom'), labels=c('E','F'), align='v'),  labels=c('A', '',''), align='v', ncol=1, rel_heights = c(1,1,1)),
          bl, ncol=2, labels=c('', 'B'), align='h')#,
#tz, ncol=1, rel_heights=c(1,0.5))
dev.off()


plot_grid(plot_grid(av, td , ncol=1, align='hv', rel_heights = c(1,2)), bl, ncol=2, align='h')
#                    plot_grid(tzdp + ylab('T. dactyloides KS position (Mb)'), fig_chrrearr + theme(legend.position='bottom'), labels=c('E','F'), align='v'),  labels=c('A', '',''), align='v', ncol=1, rel_heights = c(1,1,1)),
#          bl, ncol=2, labels=c('', 'B'), align='h')#,


plot_grid(plot_grid(av, fig_chrrearr + theme(legend.position='bottom') , rel_heights=c(2,1),ncol=1, align='h', labels=c('A', 'C')), bl, ncol=2, align='h')


## dynamically determine plot dimenstions


facet_sums <- av %>%
  group_by(PANEL) %>%  # PANEL corresponds to facets
  summarise(total_y = sum(y), total_x=sum(x))

# Total y sum across all facets
total_y_sum <- sum(facet_sums$total_y)
total_x_sum=sum(facet_sums$total_x)
# Set dynamic height based on total sum
plot_height <- total_y_sum / 100  # Adjust scale as needed




### each species dotplot supp

process_anchors_to_dotplot_SUPP <- function(filepath, color_palette=muted_colors, minBlock=10, title='', refChrs=c(paste0('Chr0', 1:9), 'Chr10'), queryChrs='', 
                                              queryChrtoFlip='', ylabelspecies='',ploidy='') {
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
  
  ## flip querychr that look funny
  if(queryChrtoFlip!=''){
    data$revQueryStart[data$queryChr%in%queryChrtoFlip]=abs(data$queryStart - data$maxChr)[data$queryChr%in%queryChrtoFlip]
  }
  
  ## clean up future facet labels - this is very specific to chr level s
  #  data$refLabel=gsub('Chr0', '', data$refChr)
  #  data$refLabel=gsub('Chr', '', data$refLabel)
  #  data$queryLabel=gsub('chr', '', data$queryLabel)
  
  ## set up the image
  #  im=data.frame(refChr='Chr01', queryChr=NA, referenceStart=1, revQueryStart=1, path=pathtokaryotype)
  #  im$queryChr=names(which.max(table(data$queryChr[data$refChr=='Chr01'])))
  
  # Combine the image and the border into a single grob
  #  imgA <- grobTree(border, imgA)
  
  ystriptextsize=ifelse(length(unique(queryChrs))>50, 6,9)
  
  # Create the plot
  p=ggplot(data[ data$refChr %in% names(color_palette) & data$refChr%in%refChrs & data$queryChr%in%queryChrs, ],
           aes(x = referenceStart / 1e6, y = revQueryStart / 1e6, color = refChr)) +
    geom_point(size=0.5) +
    facet_grid(queryChr ~ refChr, scales = 'free', space='free') +
    scale_color_manual(values = color_palette) +
    theme(legend.position = 'none') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    geom_hline(aes(yintercept=maxChr/1e6), lty='dashed', color='gray') +
    geom_vline(xintercept=0, color='gray95')+
    ggtitle(paste0(ylabelspecies,', ', ploidy)) + xlab('P. vaginatum position (Mb)') +
    ylab(paste0(ylabelspecies, ' position (Mb)')) + 
    #  ylab('Position (Mb)') +
    theme(#strip.text.y = element_text(angle = 0, hjust = 0), 
      strip.placement.y = "outside" , 
      strip.text = element_text(size = 8, color = "darkblue", face = "bold"),
      strip.background = element_rect(fill = "lightblue", color = "darkblue", linewidth = 1),
      strip.text.y = element_text(angle=0, size=ystriptextsize),
      axis.text.x=element_text(size=9),
      axis.text.y=element_text(size=5),
      panel.spacing = unit(0.1, 'lines'),
      plot.title=element_text(color=ploidycolors[ploidy], size=10))+
    scale_x_continuous(breaks = scales::breaks_pretty(n = 2), expand=c(0,0)) +  # Automatically choose 3 breaks for x-axis
    scale_y_continuous(breaks = scales::breaks_pretty(n = 2), expand=c(0,0))    # Automatically choose 3 breaks for y-axis
  
p  

}

## diploids
for(i in asize$V2[asize$ploidy=='Diploid']){
  pdf(paste0('../figures/dotplots_by_sp/', i, '.pdf'),14,10)
  if(!i%in%c('rtuber', 'telega')){
  dp=process_anchors_to_dotplot_SUPP(filepath = paste0('../syntenic_anchors/anchors/', i, '-Pv-2'), minBlock=20, queryChrtoFlip = 'chr9', ylabelspecies = asize$speciesLabel[asize$V2==i],  ploidy='Diploid')
}else(
  dp=process_anchors_to_dotplot_SUPP(filepath = paste0('../syntenic_anchors/anchors/', i, '-Pv-4'), minBlock=20, queryChrtoFlip = 'chr9', ylabelspecies = asize$speciesLabel[asize$V2==i],  ploidy='Diploid')
)
  if(i=='telega'){
    dp=process_anchors_to_dotplot_SUPP(filepath = paste0('../syntenic_anchors/anchors/', i, '-Pv-4'), minBlock=40, queryChrtoFlip = 'chr9', ylabelspecies = asize$speciesLabel[asize$V2==i],  ploidy='Diploid')
  }
  
  print(dp)
  dev.off()  
}

## tetraploids
for(i in asize$V2[asize$ploidy=='Tetraploid']){
pdf(paste0('../figures/dotplots_by_sp/', i, '.pdf'),14,10)
dp=process_anchors_to_dotplot_SUPP(filepath = paste0('../syntenic_anchors/anchors/', i, '-Pv-4'), minBlock=20, queryChrtoFlip = 'chr9', ylabelspecies = asize$speciesLabel[asize$V2==i],  ploidy='Tetraploid')
print(dp)
dev.off()  
}

# A. tenuifolius has way too many contigs - plot it separately :(
i='atenui'
pdf(paste0('../figures/dotplots_by_sp/', i, '.pdf'),14,10)
dp=process_anchors_to_dotplot_SUPP(filepath = paste0('../syntenic_anchors/anchors/', i, '-Pv-4'), minBlock=60, queryChrtoFlip = 'chr9', ylabelspecies = asize$speciesLabel[asize$V2==i],  ploidy='Tetraploid')
print(dp)
dev.off()  


## hexaploids
for(i in asize$V2[asize$ploidy=='Hexaploid']){
  pdf(paste0('../figures/dotplots_by_sp/', i, '.pdf'),14,10)
  dp=process_anchors_to_dotplot_SUPP(filepath = paste0('../syntenic_anchors/anchors/', i, '-Pv-6'), minBlock=20, queryChrtoFlip = 'chr9', ylabelspecies = asize$speciesLabel[asize$V2==i],  ploidy='Hexaploid')
  print(dp)
  dev.off()  
}

## paleotetraploid
for(i in asize$V2[asize$ploidy=='Paleotetraploid']){
  pdf(paste0('../figures/dotplots_by_sp/', i, '.pdf'),14,10)
  if(!i%in%c('tdacn1', 'tdacs1')){
  dp=process_anchors_to_dotplot_SUPP(filepath = paste0('../syntenic_anchors/anchors/', i, '-Pv-4'), minBlock=20, queryChrtoFlip = 'chr9', ylabelspecies = asize$speciesLabel[asize$V2==i],  ploidy='Paleotetraploid')
  }else{
    dp=process_anchors_to_dotplot_SUPP(filepath = paste0('../syntenic_anchors/anchors/', i, '-Pv-2'), minBlock=20, queryChrtoFlip = 'chr9', ylabelspecies = asize$speciesLabel[asize$V2==i],  ploidy='Paleotetraploid')
  }
    print(dp)
  dev.off()  
}


## todo
## XXXX - axis title with species name, colored by ploidy
###. now italicize appropriately??
## XXXXX remove species from y axis then,
## figure out spacing - force square??
## make genespace ewww - three/four plots to bottom... OR genespace on top????
## names for sp can be T. dact. KS?? Z. m. hue. ? some shortened abbreviation
## make tripsacum genespace maybe??

### facet labels at top as chr numbers with paspalum colors - for zea/trip comparison, back to black!!
#### how to do tripsacum chromosome labels that aren't so wide?????
##### color trip/zea segments by paspalum color to make chromsome fusions more clear...
### decide whether to add knob positions to this figure or jsut to supplement