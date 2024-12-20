library(GenomicRanges)
library(rtracklayer)

### chatgpt suggests and i think it's right, that i should give up and switch to bedtools!

all=read.table('../../panand_sp_ploidy.txt')
all=all[!all$V2 %in% c('pprate', 'tdactm', 'tzopol', 'osativ', 'bdista', 'agerjg', 'svirid', 'eophiu', 'tdacs2', 'tdacn2'),]

## f i lost sorghum
all=all[all$V2!='sbicol',]

genotype='achine'

#genes <- import("genes.gff", format = "gff")
syn=fread('~/transfer/anchorPositions_AnchorWave_Paspalum.txt.gz')
syn=syn[!syn$genome %in% c('pprate', 'svirid', 'tdacn2', 'tdacs2', 'tdactm', 'tzopol', 'bdista', 'agerjg', 'eophiu', 'osativ'),]

syngr=GRanges(seqnames=syn$queryChr, IRanges(start=syn$queryStart, end=syn$queryEnd), strand=syn$strand)
mcols(syngr)$genome=syn$genome
mcols(syngr)$quickgene=syn$quickgene

genes=syngr[syngr$genome==genotype,]

tes=import.gff3(paste0('trash/', all$V1[all$V2==genotype], '_EDTAandTandemRepeat.gff3'))
tes=tes[tes$type %in% gfftypestokeep,]
tes=unstrand(tes)

#tes <- import("tes.gff", format = "gff")

# Function to create sliding windows
create_windows <- function(region, range_start, range_end, step = 100, width = 100) {
  starts <- seq(range_start, range_end - width, by = step)
  ends <- starts + width - 1
  windows <- GRanges(seqnames = seqnames(region),
                     ranges = IRanges(start = start(region) + starts,
                                      end = start(region) + ends),
                     strand = strand(region))
  return(windows)
}

# Generate upstream and downstream regions
upstream <- flank(genes, width = 20000, start = TRUE)
downstream <- flank(genes, width = 20000, start = FALSE)

# Create 100 bp tiles for upstream and downstream
upstream_windows <- unlist(tile(upstream, width = 100))
downstream_windows <- unlist(tile(downstream, width = 100))
gene_windows <- unlist(tile(genes, width = 100))



# Function to calculate coverage
get_window_coverage <- function(windows, tes) {
  overlaps <- findOverlaps(windows, tes)
  covered_bp <- tapply(width(intersect(windows[queryHits(overlaps)], tes[subjectHits(overlaps)])),
                       queryHits(overlaps), sum, default = 0)
  
  total_bp <- width(windows)
  proportion <- covered_bp / total_bp
  return(proportion)
}

# Calculate proportions for each region
upstream_coverage <- get_window_coverage(upstream_windows, tes)
downstream_coverage <- get_window_coverage(downstream_windows, tes)
gene_coverage <- get_window_coverage(unlist(gene_windows), tes)


# Combine into a single data frame
results <- data.frame(
  WindowType = c(rep("Upstream", length(upstream_coverage)),
                 rep("GeneBody", length(gene_coverage)),
                 rep("Downstream", length(downstream_coverage))),
  WindowID = c(seq_along(upstream_coverage),
               seq_along(gene_coverage),
               seq_along(downstream_coverage)),
  ProportionCovered = c(upstream_coverage, gene_coverage, downstream_coverage)
)


write.csv(results, "te_coverage_windows.csv", row.names = FALSE)
