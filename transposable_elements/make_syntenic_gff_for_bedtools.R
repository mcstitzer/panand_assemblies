
library(GenomicRanges)
library(rtracklayer)
library(dplyr)

### chatgpt suggests and i think it's right, that i should give up and switch to bedtools!

all=read.table('../../panand_sp_ploidy.txt')
all=all[!all$V2 %in% c('pprate', 'tdactm', 'tzopol', 'osativ', 'bdista', 'agerjg', 'svirid', 'eophiu', 'tdacs2', 'tdacn2'),]


genotype='achine'

#genes <- import("genes.gff", format = "gff")
syn=fread('~/transfer/anchorPositions_AnchorWave_Paspalum.txt.gz')
syn=syn[!syn$genome %in% c('pprate', 'svirid', 'tdacn2', 'tdacs2', 'tdactm', 'tzopol', 'bdista', 'agerjg', 'eophiu', 'osativ'),]

syngr=GRanges(seqnames=syn$queryChr, IRanges(start=syn$queryStart, end=syn$queryEnd), strand=syn$strand)
mcols(syngr)$genome=syn$genome
mcols(syngr)$quickgene=syn$quickgene


# Iterate through each row of the 'all' data frame
for (i in seq_len(nrow(all))) {
  # Get the genome and output file name
  genome_id <- all$V2[i]
  output_file <- paste0('../syntenic/', all$V1[i], ".syntenicAnchors.gff3")

  # Filter 'syngr' based on matching genome
  filtered_syngr <- syngr[syngr$genome == genome_id]

  # Prepare data for GFF3 format
  gff_data <- data.frame(
    seqnames = as.character(seqnames(filtered_syngr)),
    source = ".",
    type = "gene",
    start = start(filtered_syngr),
    end = end(filtered_syngr),
    score = ".",
    strand = as.character(strand(filtered_syngr)),
    phase = ".",
    attributes = paste0("ID=", filtered_syngr$quickgene)
  )

  # Write to GFF3 file
  write.table(
    gff_data,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  cat("Written:", output_file, "\n")
}




## wait bed might be better to keep track of gene name
# Iterate through each row of the 'all' data frame
for (i in seq_len(nrow(all))) {
  # Get the genome and output file name
  genome_id <- all$V2[i]
  output_file <- paste0('../syntenic/', all$V1[i], ".syntenicAnchors.bed")

  # Filter 'syngr' based on matching genome
  filtered_syngr <- syngr[syngr$genome == genome_id]

  # Prepare data for BED format
  bed_data <- data.frame(
    seqnames = as.character(seqnames(filtered_syngr)),
    start = start(filtered_syngr) - 1,  # BED format uses 0-based start
    end = end(filtered_syngr),
    name = filtered_syngr$quickgene,
    score = ".",
    strand = as.character(strand(filtered_syngr))
  )

  # Write to BED file
  write.table(
    bed_data,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  cat("Written:", output_file, "\n")
}