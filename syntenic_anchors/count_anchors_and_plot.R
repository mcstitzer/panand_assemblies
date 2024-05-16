library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

## generalize!
# Set the theme and color palette
theme_set(theme_cowplot())
# Define color palette
muted_colors <- c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf')
names(muted_colors) <- c(paste0('Chr0', 1:9), 'Chr10')

process_anchors_to_dotplot <- function(filepath, color_palette=muted_colors, minBlock=10, title='', refChrs=c(paste0('Chr0', 1:9), 'Chr10'), queryChrs='', subgenomepath='') {
  # Load data
  data <- read.table(filepath, header = TRUE)
  data <- data[data$gene != 'interanchor', ]
  
  ## get queryChrs if they aren't supplied
  if(queryChrs[1]==''){
    queryChrs=unique(data$queryChr)
  }
  
  ## assign each queryChr to a subgenome, give shapes to them for plotting
  if(subgenomepath!=''){
    sg=read.table(subgenomepath, header=F)
    data$subgenome=sg$V2[match(data$queryChr, sg$V1)]
    data$subgenome[is.na(data$subgenome)]='none'
  }else{
    data$subgenome='none'
  }

  
  # Reduce to blocks and calculate stats
  data <- data %>%
    group_by(blockIndex) %>%
    mutate(blockLength = n()) %>%
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
         aes(x = referenceStart / 1e6, y = revQueryStart / 1e6, color = refChr, 
             shape=subgenome)) +
    geom_point(alpha=0.3, size=3) +
    facet_grid(queryChr ~ refChr, scales = 'free', space = 'free') +
    scale_color_manual(values = color_palette) +
    theme(legend.position = 'none') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    geom_hline(aes(yintercept=maxChr/1e6), lty='dashed', color='gray') +
    theme(strip.text.y = element_text(angle = 0)) + 
    ggtitle(title) + scale_shape_manual(values=c(SG1='1', SG2='2', SG3='3', none='o'))
}


# Example of usage
process_anchors_to_dotplot('../syntenic_anchors/anchors/agerar-Pv-6', minBlock=20, refChrs='Chr01')

process_anchors_to_dotplot('../syntenic_anchors/anchors/udigit-Pv-6', minBlock=5, subgenomepath='../subgenomes/udigitk17_q50_f2.chrom-subgenome.tsv', refChrs = c('Chr01', 'Chr02'))

pdf('~/Downloads/udigit_subgenome_dotplot.pdf', 7,9)
process_anchors_to_dotplot('../syntenic_anchors/anchors/udigit-Pv-6', minBlock=10, subgenomepath='../subgenomes/udigitk17_q50_f2.chrom-subgenome.tsv', refChrs=c('Chr01', 'Chr02', 'Chr04', 'Chr07'))
#process_anchors_to_dotplot('../syntenic_anchors/anchors/udigit-Pv-6', minBlock=10, subgenomepath='../subgenomes/udigitk17_q50_f2.chrom-subgenome.tsv')

dev.off()


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
    mutate(blockLength = n()) %>%
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
  print(countRearrangements(Sys.glob(paste0('~/Downloads/', i, '-Pv-*'))))
}

for(i in tetraploids){
  print(countRearrangements(Sys.glob(paste0('~/Downloads/', i, '-Pv-*'))))
}
for(i in hexaploids){
  print(countRearrangements(Sys.glob(paste0('~/Downloads/', i, '-Pv-*')), minBlock=50))
}                              
for(i in c('zmB735')){
  print(countRearrangements(Sys.glob(paste0('~/Downloads/', i, '-Pv-*')), minBlock=50))
}  