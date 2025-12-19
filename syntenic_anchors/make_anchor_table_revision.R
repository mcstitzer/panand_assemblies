library(ggplot2) ## on cbsu
library(cowplot)
theme_set(theme_cowplot())
library(rtracklayer)
library(stringr)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(reshape2)


all=read.table('../panand_sp_ploidy_revision.txt', header=F, sep='\t', quote='', fill=T)

all=all[!all$V2%in%c('tdacn2','tdacs2'),]

anchors=lapply(all$V4[all$V4!=''], function(x) {
print(x)
 a=read.table(Sys.glob(paste0('../hic/',x, '*','-Pv-*')), header=T)
 a$genome=all$V2[all$V4==x][1]
 return(a)
 })


ab=Reduce(function(...) merge(..., all=T), anchors)
## for maize, remobe the transcript
ab$gene=str_split_fixed(ab$gene, '_', 2)[,1]
          
b=ab[ab$gene!='interanchor',] ## keep only genes

## have to summarize how many in each genotype, then pivot
geneanchors=b[,-c(1:7,9:10)] %>% group_by(gene, genome) %>% dplyr::summarize(count=n()) %>% pivot_wider(names_from=genome, values_from=count, values_fill=0)

## see which have a copy in all
head(geneanchors[rowSums(geneanchors[,-1]==0)==0,]) %>% data.frame()

## assign back which genes have given copy presence across genotypes
b$genomeCount=rowSums(geneanchors[,-1]!=0)[match(b$gene, geneanchors$gene)]

          
          
write.table(geneanchors, 'gene_anchorTable_AnchorWave_Paspalum_revision.txt', row.names=F, col.names=T, sep='\t', quote=F)
write.table(b, 'gene_anchors_AnchorWave_Paspalum_revision.txt', row.names=F, col.names=T, sep='\t', quote=F)

## zea chromosomes only because of alleles
geneanchorszc=rbind(b[substr(b$genome,1,1)=='z'&grepl('chr', b$queryChr),-c(1:7,9:10)], b[substr(b$genome,1,1)!='z',-c(1:7,9:10)]) %>% group_by(gene, genome) %>% dplyr::summarize(count=n()) %>% pivot_wider(names_from=genome, values_from=count, values_fill=0)
write.table(geneanchorszc, 'gene_anchorTable_AnchorWave_Paspalum_zeachronly_revision.txt', row.names=F, col.names=T, sep='\t', quote=F)


png('~/transfer/geneanchors_paspalum_revision.png',20,10, unit='in', res=300)
ggplot(melt(geneanchors), aes(y=variable, x=gene, fill=factor(value))) + geom_tile() + facet_wrap(~substr(gene,6,7), nrow=1, scale='free_x') + scale_fill_manual(values=c('white',brewer.pal(n=10, name='Paired')[c(1,2,3,4,9,10,7,8)]))#values=c('white', 'pink', 'purple', 'darkblue'))
dev.off()


          
geneanchors$refChr=b$refChr[match(geneanchors$gene, b$gene)]
