library(dplyr)
library(stringr)

a=read.table('panand_v2_transcripts_overlap_with_anchor.txt', header=T)
a$genenotranscript=str_split_fixed(a$transcripts, '_', 2)[,1]

b=read.table('~/transfer/anchorPositions_AnchorWave_Paspalum.txt.gz', header=T)

## keep just gene models
d=merge(a, b, by.x=c('species', 'gene', 'blockIdx', 'refChr'), by.y=c('genome', 'gene', 'blockIndex', 'refChr'), all.x=T)

## keep anchors and gene models
e=merge(a, b, by.x=c('species', 'gene', 'blockIdx', 'refChr'), by.y=c('genome', 'gene', 'blockIndex', 'refChr'), all=T)

write.table(e, '~/transfer/panand_anchorPositions_v2transcripts.txt', row.names=F, col.names=T, sep='\t', quote=F)

## species = six digit code
## gene = paspalum gene id (can use as grouping variable for anchor ID)
## genenotranscript = other species gene ID without transcript (you may want to filter for unique)
## queryChr = other species chromosome
## queryStart = other species start (0 based I think)
## queryEnd = other species end
## percentage = proportion of anchor in the gene model
