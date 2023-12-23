library(ggtree) ## /programs/R-4.2.1-r9/bin/R
library(ggplot2) 
library(cowplot)
theme_set(theme_cowplot())
library(rtracklayer)
library(reshape2)
library(RColorBrewer)
library(tidypaleo) ## facet species names in italics!!!!!



## this drops sorghum, since it's not our annotation
g=read.csv('../supplement_table_annotations.csv')
g$six=c('atenui', 'achine', 'agerar', 'avirgi', 'blagur', 'ccitra', 'crefra', 'cserru', 'etrips', 'hcompr', 'hconto', 'irugos', 'ppanic', 'rrottb', 'rtuber', 'smicro', 'snutan', 'sscopa', 'tdacs1','tdacs2', 'tdacn1', 'tdacn2', 'telega', 'ttrian', 'udigit', 'vcuspi', 'zdgigi', 'zdmomo', 'zmhuet', 'znicar', 'zTIL01', 'zTIL11', 'zTIL18', 'zTIL25', 'zluxur')
gg=merge(gs, g, by.x='V2', by.y='six')
gg$genecount=as.numeric(gsub(',','',gg$Genes))
gg$meangenelength=as.numeric(gsub(',','',gg$Avg.Gene.length..bp.))
gg$totalcds=as.numeric(gsub(',','',gg$Total.CDS.region..bp.))
gg$meancdslength=as.numeric(gsub(',','',gg$Avg.CDS.length..bp.))


ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')


gg$ploidy=factor(gg$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))

## come on you have to adjust for haploid assemblies!@!
gg$doubledgenecount=gg$genecount
gg$doubledgenecount[gg$haploid]=gg$doubledgenecount[gg$haploid]*2
gg$doubledtotalcds=gg$totalcds
gg$doubledtotalcds[gg$haploid]=gg$doubledtotalcds[gg$haploid]*2
#### then reduce them back down in plot, because people are used to seeing haploid values

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

dev.off()






