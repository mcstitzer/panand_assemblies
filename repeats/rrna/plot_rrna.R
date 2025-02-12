library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
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


process_anchors_to_dotplot_nor <- function(filepath, color_palette=muted_colors, minBlock=10, title='', refChrs=c(paste0('Chr0', 1:9), 'Chr10'), queryChrs='', 
                                              queryChrtoFlip='', ylabelspecies='',ploidy='',
                                              nor=data.frame(queryChr='', queryStart='', queryEnd='')) {
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
    scale_y_continuous(breaks = scales::breaks_pretty(n = 2), expand=c(0,0)) +   # Automatically choose 3 breaks for y-axis
      geom_hline(data=nor, aes(yintercept=queryStart/1e6))
p  

}

a=fread('grep -h ">" *_rrna.fa', header=F)

a <- a %>% 
  mutate(
    queryChr = sub("^.{8}([^:]+):.*", "\\1", V1),          # Remove first 8 characters and extract until the colon
    queryStart = as.integer(sub("^.*:(\\d+)-.*", "\\1", V1)),   # Extract start position
    end = as.integer(sub("^.*-(\\d+)$", "\\1", V1))        # Extract end position
  )
  
  
  asize=fread('../panand_assembly_sizes.txt', header=T, quote="", fill=T)
asize$ploidy=factor(asize$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))

  diploids=c('cserru', 'irugos', 'sbicol', 'ppanic', 'ttrian', 'crefra', 'avirgi', 'smicro', 'rtuber')
tetraploids=c('snutan', 'hconto', 'ccitra', 'achine', 'sscopa', 'etrips',  'vcuspi')
hexaploids=c('udigit', 'agerar', 'hcompr', 'blagur')

  
  process_anchors_to_dotplot_nor('~/transfer/anchors/zmB735-Pv-4', nor=a[a$genome=='zmB735',])

  for(i in asize$V2){
print( process_anchors_to_dotplot_nor(Sys.glob(paste0('~/transfer/anchors/',i, '-Pv-*'))[1], nor=a[a$genome==i,], ylabelspecies=i)
)
  }
  

  
  dev.off()
  
  
    pdf('~/transfer/try_nor_plot.trip.pdf',15,15)

   process_anchors_to_dotplot_nor('~/transfer/anchors/tdacs1-zB73v5-2', nor=a[a$genome=='tdacs1',], refChrs=c(paste0('chr', 1:10)))

  dev.off()
  
  
  