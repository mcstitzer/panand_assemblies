## maize knobs vs Tripsacum

library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(plyranges)

## generalize!
# Set the theme and color palette
theme_set(theme_cowplot())
# Define color palette
muted_colors <- c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf')
names(muted_colors) <- c(paste0('Chr0', 1:9), 'Chr10')

  data <- read.table('../syntenic_anchors/anchors/zmB735-Pv-4', header = TRUE)
  data <- read.table('../syntenic_anchors/anchors/tdacs1-zB73v5-2', header = TRUE)
  
    data <- data[data$gene != 'interanchor', ]
  
  ## get queryChrs if they aren't supplied
#  if(queryChrs[1]==''){
    queryChrs=unique(data$queryChr)
#  }
  
#pasp  refChrs=c(paste0('Chr0', 1:9), 'Chr10')
  refChrs=c(paste0('chr',1:10))
  # Reduce to blocks and calculate stats
  data <- data %>%
    group_by(blockIndex) %>%
    dplyr::mutate(blockLength = dplyr::n()) %>%
    group_by(queryChr) %>%
    mutate(freqStrand = names(which.max(table(strand))),
           maxChr = max(queryStart),
           freqRef = names(which.max(table(refChr))))
  
  # Filter data based on block length
  data <- data[data$blockLength > minBlock, ]
 # data$refChr <- factor(data$refChr, levels = c(paste0('Chr0', 1:9), 'Chr10'))
  data$refChr <- factor(data$refChr, levels = paste0('chr', 1:10))
  
  # Reverse strand calculations
  data <- data %>%
    arrange(freqRef, referenceStart, queryStart)
  data$queryChr <- factor(data$queryChr, levels = rev(data$queryChr[!duplicated(data$queryChr)]))
  data$revQueryStart <- data$queryStart
  data$revQueryStart[data$freqStrand == '-'] <- abs(data$queryStart - data$maxChr)[data$freqStrand == '-']
  
  
  ## get tandem repeat positions
  te=import.gff('~/Downloads/Zm-B73-REFERENCE-NAM-5.0.TE.gff3.gz')
  tandems=te[te$type %in% c('rDNA_intergenic_spacer_element', 'centromeric_repeat', 'subtelomere', 'knob'),]
  tandems=data.frame((tandems+500) %>% group_by(Classification) %>% reduce_ranges())
  tandems$refChr=tandems$seqnames
  tandems=tandems[tandems$seqnames%in%refChrs,]
  # Create the plot
  ggplot(data[  data$refChr%in%refChrs & data$queryChr%in%queryChrs, ],
         aes(x = referenceStart / 1e6, y = revQueryStart / 1e6) ) +##/, color = refChr)) +
    geom_point() +
    facet_grid(queryChr ~ refChr, scales = 'free', space = 'free') +
 #   scale_color_manual(values = muted_colors) +
 #   theme(legend.position = 'none') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    geom_hline(aes(yintercept=maxChr/1e6), lty='dashed', color='gray') + geom_vline(data=tandems[tandems$width>5e3,], aes(xintercept = start/1e6, color=Classification))


## just single repeat
  ggplot(data[ data$refChr %in% names(muted_colors) & data$refChr%in%refChrs & data$queryChr%in%queryChrs, ],
         aes(x = referenceStart / 1e6, y = revQueryStart / 1e6, color = refChr)) +
    geom_point() +
    facet_grid(queryChr ~ refChr, scales = 'free', space = 'free') +
    scale_color_manual(values = muted_colors) +
    theme(legend.position = 'none') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    geom_hline(aes(yintercept=maxChr/1e6), lty='dashed', color='gray') + geom_hline(data=tandems[tandems$width>5e3 & tandems$Classification=='rDNA/spacer',], aes(yintercept = start/1e6))
  
  
  ## knbo
  ggplot(data[ data$refChr %in% names(muted_colors) & data$refChr%in%refChrs & data$queryChr%in%queryChrs, ],
         aes(x = referenceStart / 1e6, y = revQueryStart / 1e6, color = refChr)) +
    geom_point() +
    facet_grid(queryChr ~ refChr, scales = 'free', space = 'free') +
    scale_color_manual(values = muted_colors) +
    theme(legend.position = 'none') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    geom_hline(aes(yintercept=maxChr/1e6), lty='dashed', color='gray') + geom_hline(data=tandems[tandems$width>1e4 & tandems$Classification=='knob/TR-1',], aes(yintercept = start/1e6))
  
  
  ## cent
  ggplot(data[ data$refChr %in% names(muted_colors) & data$refChr%in%refChrs & data$queryChr%in%queryChrs, ],
         aes(x = referenceStart / 1e6, y = revQueryStart / 1e6, color = refChr)) +
    geom_point() +
    facet_grid(queryChr ~ refChr, scales = 'free', space = 'free') +
    scale_color_manual(values = muted_colors) +
    theme(legend.position = 'none') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    geom_hline(aes(yintercept=maxChr/1e6), lty='dashed', color='gray') + geom_hline(data=tandems[tandems$width>1e4 & tandems$Classification=='Cent/CentC',], aes(yintercept = start/1e6))
  
  
  