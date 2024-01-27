
library(ggtree) ## /programs/R-4.2.1-r9/bin/R
library(ggplot2) 
library(cowplot)
theme_set(theme_cowplot())
library(rtracklayer)
library(reshape2)
library(RColorBrewer)
library(tidypaleo) ## facet species names in italics!!!!!
library(data.table)
library(dplyr)


all=read.table('panand_sp_ploidy.txt')
all=all[!all$V2 %in% c('pprate', 'tdactm', 'tzopol', 'osativ', 'bdista', 'agerjg', 'svirid', 'eophiu'),]
all=all[!all$V2 %in% c('zluxur', 'sbicol'),]


for(x in 1:nrow(all)){
  a=import.gff(Sys.glob(paste0('genes/',all$V1[x], '*.gff3')))
b=findOverlaps(unstrand(a[a$type=='gene',]))
  all$anyoverlap[x]=length(unique(c(queryHits(b[queryHits(b)!=subjectHits(b),]), subjectHits(b[queryHits(b)!=subjectHits(b),]))))
  all$overlapcount[x]=length(unique(c(subjectHits(b[queryHits(b)!=subjectHits(b),]))))

  all$genecount[x]=length(unique(queryHits(b)))
  all$nonoverlapcount[x]=all$genecount[x]-all$overlapcount[x]
  
  }



gs=read.table('~/transfer/panand_assembly_sizes.txt', header=T, sep='\t')
## this drops sorghum, since it's not our annotation
g=read.csv('supplement_table_annotations.csv')
g$six=c('atenui', 'achine', 'agerar', 'avirgi', 'blagur', 'ccitra', 'crefra', 'cserru', 'etrips', 'hcompr', 'hconto', 'irugos', 'ppanic', 'rrottb', 'rtuber', 'smicro', 'snutan', 'sscopa', 'tdacs1','tdacs2', 'tdacn1', 'tdacn2', 'telega', 'ttrian', 'udigit', 'vcuspi', 'zdgigi', 'zdmomo', 'zmhuet', 'znicar', 'zTIL01', 'zTIL11', 'zTIL18', 'zTIL25', 'zluxur')
gg=merge(gs, g, by.x='V2', by.y='six')
gg$genecount=as.numeric(gsub(',','',gg$Genes))
gg$meangenelength=as.numeric(gsub(',','',gg$Avg.Gene.length..bp.))
gg$totalcds=as.numeric(gsub(',','',gg$Total.CDS.region..bp.))
gg$meancdslength=as.numeric(gsub(',','',gg$Avg.CDS.length..bp.))
gg$meanintronlength=gg$meangenelength-gg$meancdslength


ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')


gg$ploidy=factor(gg$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))
gg$nonoverlapgenecount=all$nonoverlapcount[match(gg$V2, all$V2)]
gg$doublednonoverlapgenecount=gg$nonoverlapgenecount
gg$doublednonoverlapgenecount[gg$haploid]=gg$doublednonoverlapgenecount[gg$haploid]*2
gg$doubledgenecount=gg$genecount
gg$doubledgenecount[gg$haploid]=gg$doubledgenecount[gg$haploid]*2
gg$doubledtotalcds=gg$totalcds
gg$doubledtotalcds[gg$haploid]=gg$doubledtotalcds[gg$haploid]*2





## gettin real sloppy here - this is from aum, which is from counting genes from gene table for fractionated/resistant genes
#> dput(colSums(a[,-1]))
syntcount=c(achine = 30037, agerar = 62404, avirgi = 14549, blagur = 79949, 
ccitra = 52037, crefra = 14712, cserru = 23638, etrips = 26922, 
hcompr = 65110, hconto = 54530, irugos = 16094, ppanic = 14777, 
rrottb = 27766, sbicol = 15053, smicro = 14200, sscopa = 42060, 
tdacn1 = 18980, tdacn2 = 18784, tdacs1 = 19151, tdacs2 = 19130, 
telega = 28664, ttrian = 23530, udigit = 40017, vcuspi = 41068, 
zTIL01 = 17394, zTIL11 = 17057, zTIL18 = 17311, zTIL25 = 17338, 
zdgigi = 22451, zdmomo = 18246, zluxur = 17052, zmB735 = 17032, 
zmhuet = 17500, znicar = 18980, atenui = 49008, rtuber = 14986, 
snutan = 35298) ## don't worry about b73 name, since we don't want to compare to the b73 gene annotation anyways
gg$synteniccount=syntcount[match(gg$V2, names(syntcount))]
gg$doubledsyntenic=gg$synteniccount
gg$doubledsyntenic[gg$haploid]=gg$doubledsyntenic[gg$haploid]*2

ks=fread('~/transfer/ks_to_look_for_mixtures.txt', sep='\t', quote='')
ksp=ks


mks=ks %>% group_by(genome) %>% dplyr::summarize(median(ks))

gg$mks=mks$`median(ks)`[match(gg$V2, mks$genome)]
cor(gg$doublednonoverlapgenecount, gg$mks, use='complete')

gg %>% group_by(ploidy) %>% summarize(mean(doublednonoverlapgenecount, na.rm=T), mean(mks, na.rm=T), cor(doublednonoverlapgenecount, mks, use='complete'),
                                      cor(doubledgenecount, mks, use='complete'))







pdf(paste0('~/transfer/gene_ploidy.', Sys.Date(), '.pdf'), 6,6)
ggplot(gg, aes(x=ploidy, y=doubledgenecount/2, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=120000) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('Annotated Genes') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(gg, aes(x=ploidy, y=meangenelength, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=5000) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('Average Gene Length (bp)') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(gg, aes(x=ploidy, y=doubledtotalcds/1e6/2, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=210000000/1e6) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('Total CDS Length (Mbp)') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(gg, aes(x=ploidy, y=meancdslength, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=1300) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('Average CDS Length (bp)') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(gg, aes(x=ploidy, y=meangenelength-meancdslength, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=3700) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('Average Intron Length (bp)') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(gg, aes(x=meangenelength, y=meangenelength-meancdslength, color=ploidy))  + geom_point()     +      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Average Gene length (bp)') + ylab('Average Intron Length (bp)') 
ggplot(gg, aes(x=ploidy, y=doublednonoverlapgenecount/2, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=120000) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('Annotated Genes') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggplot(gg, aes(x=mks, y=doublednonoverlapgenecount/2, color=ploidy))  + geom_point()     +      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Median Ks') + ylab('Nonoverlapping Gene Count (haploid)') 
ggplot(gg, aes(x=mks, y=doubledsyntenic/2, color=ploidy))  + geom_point()     +      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Median Ks') + ylab('Syntenic Gene Count (haploid)') 
ggplot(gg, aes(x=mks, y=doublednonoverlapgenecount/2-doubledsyntenic/2, color=ploidy))  + geom_point()     +      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Median Ks') + ylab('Nonsyntenic Gene Count (haploid)') 

dev.off()
