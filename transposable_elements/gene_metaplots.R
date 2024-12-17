library(ggplot2)
library(cowplot)
library(rtracklayer)
library(reshape2)
library(tidyverse)
library(data.table)
library(plyranges)
library(ggrepel)
theme_set(theme_cowplot())

### redone by chatgpt..... 
## work on figuring out sbicolor, since I don't have it in here now !!~!

# Input data
all <- read.table('../panand_sp_ploidy.txt')
all <- all[!all$V2 %in% c('pprate', 'tdactm', 'tzopol', 'osativ', 'bdista', 'agerjg', 'svirid', 'eophiu', 'tdacs2', 'tdacn2'),]

gff_types_to_keep <- c("Gypsy_LTR_retrotransposon", "LTR_retrotransposon", "Copia_LTR_retrotransposon",
                       "CACTA_TIR_transposon", "helitron", "hAT_TIR_transposon", "PIF_Harbinger_TIR_transposon",
                       "Tc1_Mariner_TIR_transposon", "Mutator_TIR_transposon", 'tandem_repeat') 

# Function to process a genome
process_genome <- function(genome, all_data, gff_types) {
  gff_file <- paste0('trash/', all_data$V1[all_data$V2 == genome], '_EDTAandTandemRepeat.gff3')
  a <- import.gff3(gff_file)
  a <- a[a$type %in% gff_types]
  a$genome <- genome
  return(a)
}

# Import and filter GFF files
genome_list <- lapply(all$V2, process_genome, all_data = all, gff_types = gff_types_to_keep)
names(genome_list) <- all$V2

# Prepare syntenic regions
syn <- fread('~/transfer/anchorPositions_AnchorWave_Paspalum.txt.gz')
syngr <- GRanges(seqnames = syn$queryChr, IRanges(start = syn$queryStart, end = syn$queryEnd), strand = syn$strand)
mcols(syngr) <- syn[, .(genome, quickgene)]
geneflanks <- promoters(syngr, upstream = 20000, downstream = 1000)

# Function to compute TE overlap stats
compute_TE_window_stats <- function(te_data, flank_ranges, genome) {
  tewindow <- join_overlap_intersect(unstrand(flank_ranges), te_data)
  posplot <- data.frame(tewindow[tewindow$ogstrand == '+',])[, c('window', 'width', 'partition')]
  posplot <- posplot %>% complete(partition, window, fill = list(width = 0))
  summary_stats <- posplot %>% group_by(window) %>% summarize(medianTE = median(width))
  return(summary_stats)
}

# Combine TE window stats for all genomes
metaplot <- data.frame(window = unique(geneflanks$window))
for (genome in names(genome_list)) {
  tes <- genome_list[[genome]]
  stats <- compute_TE_window_stats(tes, geneflanks, genome)
  metaplot[, genome] <- stats$medianTE
}

# Melt data for plotting
metaplot_melt <- melt(metaplot, id.vars = 'window')

# Load assembly sizes and annotate TE proportions
gs <- read.table('~/transfer/panand_assembly_sizes.txt', header = TRUE, sep = '\t')
gs$teprop <- gs$haploidRepeatSize / (gs$haploidAssemblySize - gs$haploidNCount)
metaplot_melt$teprop <- gs$teprop[match(metaplot_melt$variable, gs$V2)]

# Generalized Plot Function
plot_TE_distribution <- function(data, y_label, title = NULL) {
  ggplot(data, aes(group = variable, x = window, y = value, color = teprop)) +
    geom_vline(xintercept = seq(1, max(data$window), by = 10), color = 'whitesmoke') +
    geom_line() +
    scale_color_viridis_c(option = 'inferno') +
    xlab('Base pairs away from TSS/TTS') +
    ylab(y_label) +
    geom_vline(xintercept = c(200, 210, 220), lty = c('dashed', 'solid', 'dashed')) +
    scale_x_continuous(breaks = c(1, 200, 210, 220, max(data$window)),
                       labels = c('-20000', 'TranslationSS', 'X', 'TranslTermS', '+20000')) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
    ggtitle(title)
}

# Generate plots
pdf('~/transfer/te_metaplot_syntenicGenes.pdf', 12, 6)
print(plot_TE_distribution(metaplot_melt, 'Median TE base pairs in 100bp window'))
dev.off()
