library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())




## actually, incorporated into script, so need to use like this
## Rscript plot_te_genes_on_new_chr.R "$PREFIX" "$OUTPUT_TRASH_GFF" "$OUTPUT_HELIXER_GFF" "$AGP_FILE"
args <- commandArgs(trailingOnly = TRUE)
genotype <- args[1]       # The PREFIX value
agpfile <- args[2]        # The agp file
ploidy <- args[3]

# Print the Inputs for Debugging
cat("Genotype Prefix: ", genotype, "\n")
cat("AGP File: ", agpfile, "\n")
cat("Ploidy: ", ploidy, "\n")



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

muted_colors <- c("#b34064", "#459abf", "#68b488", "#b3ac40", "#8d4cba", 
                  "#bf9140", "#ae459a", "#99aabf", "#409f90", "#405973")

names(muted_colors) <- c(paste0('Chr0', 1:9), 'Chr10')

ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')



# Example of usage
#process_anchors_to_dotplot('../syntenic_anchors/anchors/agerar-Pv-6', minBlock=20, refChrs='Chr01')

pdf(paste0('~/transfer/', genotype, '_dotplot.pdf'), 14,14)
process_anchors_to_dotplot(paste0(genotype, '-Pv-', as.character(as.numeric(ploidy))), minBlock=3, title=genotype)
dev.off()
