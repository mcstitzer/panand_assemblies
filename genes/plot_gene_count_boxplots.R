library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(rtracklayer)

## count genes as boxplot, show relaxed selection


ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')



b=read.table('../syntenic_anchors/gene_anchorTable_AnchorWave_Paspalum.txt', header=T)
b=b[,-which(colnames(b)%in%c('svirid',"pprate", "osativ", "eophiu", "bdista", "agerjg", 'tdactm', 'tzopol'))]
b=b[,-which(colnames(b)%in%c('tdacs2', 'tdacn2'))]

## themeda might loose abn arm like contortus?
## serrulatus is nanopore
#b=b[,-which(colnames(b)%in%c('ttrian', 'cserru'))]


table(rowSums(b[,-1]>0))
bb=b[rowSums(b[,-1]>0)>=32,]

## output these fun genes for hybpiper in the short reads
write.table(gsub('.1.v3.1', '', bb$gene), '~/Downloads/sharedSyntenicAnchors.txt', quote=F, row.names=F, col.names=F)
table(substr(gsub('.1.v3.1', '', bb$gene),1,7))


asize=fread('../general_summaries/panand_assembly_sizes.txt', header=T, quote="", fill=T)

asize$syntAnchors=colSums(bb[,-1])[match(asize$V2, names(colSums(bb[,-1])))]
asize$syntAnchorsCount=colSums(bb[,-1]>0)[match(asize$V2, names(colSums(bb[,-1]>0)))]
asize$syntAnchors=ifelse(asize$haploid, asize$syntAnchors, asize$syntAnchors/2)


asize$ploidy=factor(asize$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))
asize$diploidEquivalent=ifelse(asize$ploidy=='Diploid', asize$syntAnchors, ifelse(asize$ploidy%in%c('Tetraploid', 'Paleotetraploid'), asize$syntAnchors/2, ifelse(asize$ploidy=='Hexaploid', asize$syntAnchors/3, NA)))


## boxplot of syntenic gene count
ggplot(asize, aes(x=ploidy, y=syntAnchors, group=ploidy, color=ploidy)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(seed = 1), size=3)+ 
  scale_color_manual(values=ploidycolors) + ylim(0,30000) + 
  ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=30000) + 
  geom_hline(yintercept=c(median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
  annotate("text", x = 4.45, y = median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  annotate("text", x = 4.45, y = median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  annotate("text", x = 4.45, y = median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  xlab('Ploidy') + ylab('Syntenic Anchor Genes') + theme(legend.position='NULL')


## boxplot of syntenic gene count --diploid equivalent
ggplot(asize, aes(x=ploidy, y=diploidEquivalent, group=ploidy, color=ploidy)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(seed = 1), size=3)+ 
  scale_color_manual(values=ploidycolors) + ylim(0,12000) + 
  ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=12000) + 
  geom_hline(yintercept=c(median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
#  annotate("text", x = 4.45, y = median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  xlab('Ploidy') + ylab('Diploid Equivalent\nSyntenic Anchor Genes') + theme(legend.position='NULL')


### count gene annotation

## function to add taxon name
getGenes <- function(filepath, id='') {
  # Load data
  data <- import.gff3(filepath)
  data=data[data$type=='gene',]
  data$genome=id
  return(data)
}

## file from google sheets that has gff3 name
genepath=read.table('../genes/panand_gene_annotations.txt', header=F, sep='\t')

genelist=lapply(asize$V2[asize$V2!='zluxur'], function(x) getGenes(paste0('~/Downloads/Annotations_v4.1/', genepath$V4[genepath$V2==x]), x))
genes=do.call(c, genelist)

asize$geneCount=sapply(asize$V2, function(x) sum(genes$genome==x))
asize$haploidGeneCount=ifelse(asize$haploid, asize$geneCount, asize$geneCount/2)
asize$diploidEquivalentGeneCount=ifelse(asize$ploidy=='Diploid', asize$haploidGeneCount, ifelse(asize$ploidy%in%c('Tetraploid', 'Paleotetraploid'), asize$haploidGeneCount/2, ifelse(asize$ploidy=='Hexaploid', asize$haploidGeneCount/3, NA)))

## boxplot of annotated gene count
ggplot(asize, aes(x=ploidy, y=haploidGeneCount, group=ploidy, color=ploidy)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(seed = 1), size=3)+ 
  scale_color_manual(values=ploidycolors) + ylim(0,110000) + 
  ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=110000) + 
  geom_hline(yintercept=c(median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
  annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  xlab('Ploidy') + ylab('Genes') + theme(legend.position='NULL')


ggplot(asize, aes(x=syntAnchors, y=haploidGeneCount, group=ploidy, color=ploidy)) + geom_point() + scale_color_manual(values=ploidycolors)
  

## diploid equivalent gene count
## boxplot of annotated gene count
ggplot(asize, aes(x=ploidy, y=diploidEquivalentGeneCount, group=ploidy, color=ploidy)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(seed = 1), size=3)+ 
  scale_color_manual(values=ploidycolors) + ylim(0,50000) + 
  ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=50000) + 
  geom_hline(yintercept=c(median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
#  annotate("text", x = 4.45, y = median(asize$diploidEquivalentGeneCount[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
#  annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
#  annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  xlab('Ploidy') + ylab('Diploid Equivalent Genes') + theme(legend.position='NULL')



ks=fread('../general_summaries/ks_to_look_for_mixtures.txt', header=T, quote='')
ks=data.frame(ks)[ks$ks>0.001,]
mks=ks[ks$ploidy %in%c('Tetraploid', 'Paleotetraploid', 'Hexaploid'),] %>% group_by(genome, ploidy, species, haploid) %>% dplyr::summarize(median=median(ks, na.rm=T), nonallelic=median(ks[ks>0.005], na.rm=T))
mks$medianNonAllelicCorr=ifelse(mks$nonallelic-mks$median>0.01, mks$nonallelic, mks$median)


mu=6.5e-9
mks$mya=mks$medianNonAllelicCorr/2/mu/1e6
ks$mya=ks$ks/2/mu/1e6
ggplot(data.frame(ks)[ks$ploidy %in% c('Tetraploid', 'Paleotetraploid', 'Hexaploid') & !ks$genome%in%c('zdmomo', 'zdgigi', 'zdnica'),], aes(x=mya, color=ploidy, group=genome)) + geom_density() + scale_color_manual(values=ploidycolors) + facet_wrap(~ploidy, ncol=1, scales='free_y') + geom_vline(data=mks[!mks$genome %in% c('zdmomo', 'zdgigi', 'zdnica'),], aes(xintercept=mya))


asize$medianKs=mks$medianNonAllelicCorr[match(asize$V2, mks$genome)]
asize$mya=mks$mya[match(asize$V2, mks$genome)]
asize$mya[asize$ploidy=='Diploid']=0
asize$medianKs[asize$ploidy=='Diploid']=0


##
ggplot(asize, aes(x=medianKs, y=diploidEquivalent, group=ploidy, color=ploidy)) + 
#  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(seed = 1), size=3)+ 
  scale_color_manual(values=ploidycolors) + #ylim(0,110000) + 
 # ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=110000) + 
  geom_hline(yintercept=c(median(asize$diploidEquivalent[asize$ploidy=='Diploid'], na.rm=T)*c(1)), lty='dotted', color='darkgray')+
#  annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
 # annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
 # annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  xlab('Ks') + ylab('anchors') + theme(legend.position='NULL') + geom_text(aes(label=V2))

polyploidy_data <- data.frame(
  Polyploidy_event = c("Arabidopsis thaliana", "Ancestral grass", "Ancestral eudicot", 
                       "Ancestral grass", "Arabidopsis thaliana", "Musa acuminata", 
                       "Solanum lycopersicum", "Glycine max", "Gossypium raimondii", 
                       "Brassica rapa", "Populus trichocarpa", "Malus domestica", 
                       "Zea mays", "Glycine max"),
  Event = c("b", "s", "g", "r", "a", "a/b", "", "B", "", "", "", "", "", "A"),
  Age_MY = c("40-70", "~130", "~117", "70-90", "~40", "~65", "52-91", "~59", 
             "13-20", "5-9", "60-65", "30-65", "5-12", "~13"),
  Ks = c(2.27, 1.72, 1.29, 0.94, 0.78, 0.60, 0.60, 0.59, 0.50, 0.25, 0.25, 0.20, 0.18, 0.13),
  Genome_Average = c("14%", "17%", "18%", "18%", "14%", "35%", "19%", "26%", 
                     "26%", "43%", "35%", "47%", "14%", "51%"),
  Meiotic_All = c("2%", "9%", "12%", "11%", "12%", "11%", "11%", "15%", "22%", "35%", 
                  "32%", "53%", "8%", "58%"),
  Meiotic_Recombination = c("0%", "0%", "3%", "8%", "0%", "3%", "3%", "8%", 
                            "15%", "23%", "23%", "39%", "5%", "45%")
)  


bbd=bb
bbd[,asize$V2[asize$haploid]]=bbd[,asize$V2[asize$haploid]]*2

## meiotic genes
mg=fread('../genes/phytozome_meiotic_to_paspalum.txt')
mg$gene=paste0(mg$`Ortholog Gene Name`, '.1.v3.1')


bbd[bbd$gene %in% mg$gene,]


asize$ploidyNumber=ifelse(asize$ploidy=='Diploid', 2, ifelse(asize$ploidy%in%c('Tetraploid', 'Paleotetraploid'), 4, ifelse(asize$ploidy=='Hexaploid', 6, NA)))
asize$meiosisRetention=sapply(asize$V2, function(x) sum(bbd[bbd$gene %in% mg$gene,x]> asize$ploidyNumber[asize$V2==x]-2)/nrow(bbd[bbd$gene %in% mg$gene,]))
asize$geneRetention=sapply(asize$V2, function(x) sum(bbd[,x]> asize$ploidyNumber[asize$V2==x]-2)/nrow(bbd))


ggplot(asize, aes(x=geneRetention, y=meiosisRetention, color=ploidy)) + geom_point() + geom_text(aes(label=V2))+ scale_color_manual(values=ploidycolors) + geom_abline(slope=1)

ggplot(asize, aes(x=medianKs, y=geneRetention, group=ploidy, color=ploidy)) + 
  #  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(seed = 1), size=3)+ 
  geom_point(data=polyploidy_data, inherit.aes=F, aes(x=Ks, y=as.numeric(gsub('%', '', Genome_Average))/100)) + 
  scale_color_manual(values=ploidycolors) + #ylim(0,110000) + 
  # ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=110000) + 
#  geom_hline(yintercept=c(median(asize$diploidEquivalent[asize$ploidy=='Diploid'], na.rm=T)*c(1)), lty='dotted', color='darkgray')+
  #  annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  # annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  # annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  xlab('Ks') + ylab('Proportion Gene Retention\n(syntenic anchor in genome, one allele allowed to be lost)') + theme(legend.position='NULL') + geom_text(aes(label=V2))


ggplot(asize, aes(x=medianKs, y=meiosisRetention, group=ploidy, color=ploidy)) + 
  #  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(seed = 1), size=3)+ 
  geom_point(data=polyploidy_data, inherit.aes=F, aes(x=Ks, y=as.numeric(gsub('%', '', Meiotic_All))/100)) + 
  scale_color_manual(values=ploidycolors) + #ylim(0,110000) + 
  # ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=110000) + 
  #  geom_hline(yintercept=c(median(asize$diploidEquivalent[asize$ploidy=='Diploid'], na.rm=T)*c(1)), lty='dotted', color='darkgray')+
  #  annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  # annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  # annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  xlab('Ks') + ylab('Proportion Meiotic Gene Retention\n(syntenic anchor in genome, one allele allowed to be lost)') + theme(legend.position='NULL') + geom_text(aes(label=V2))


asize$meiosisEnrichment=asize$meiosisRetention-asize$geneRetention
polyploidy_data$meiosisEnrichment=as.numeric(gsub('%', '', polyploidy_data$Meiotic_All))/100 -as.numeric(gsub('%', '', polyploidy_data$Genome_Average))/100 
## offset between genome-wide and meiotic!
ggplot(asize, aes(x=medianKs, y=meiosisEnrichment, group=ploidy, color=ploidy)) + 
  #  geom_boxplot(outlier.shape = NA) + 
  geom_hline(yintercept=0, color='gray', lty='dashed') +
  geom_point(position = position_jitter(seed = 1), size=3)+ 
  geom_point(data=polyploidy_data, inherit.aes=F, aes(x=Ks, y=meiosisEnrichment)) + 
  scale_color_manual(values=ploidycolors) + #ylim(0,110000) + 
  # ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=110000) + 
  #  geom_hline(yintercept=c(median(asize$diploidEquivalent[asize$ploidy=='Diploid'], na.rm=T)*c(1)), lty='dotted', color='darkgray')+
  #  annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  # annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  # annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  xlab('Ks') + ylab('Meiotic Gene Retention - Overall Gene Retention') + theme(legend.position='NULL') + geom_text(aes(label=V2))

  




ggplot(asize, aes(x=medianKs, y=meiosisEnrichment, group=ploidy, color=ploidy)) + 
  #  geom_boxplot(outlier.shape = NA) + 
  geom_hline(yintercept=0, color='gray', lty='dashed') +
  geom_point(position = position_jitter(seed = 1), size=3)+ 
  scale_color_manual(values=ploidycolors) + #ylim(0,110000) + 
  # ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=110000) + 
  #  geom_hline(yintercept=c(median(asize$diploidEquivalent[asize$ploidy=='Diploid'], na.rm=T)*c(1)), lty='dotted', color='darkgray')+
  #  annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  # annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  # annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  xlab('Ks') + ylab('Meiotic Gene Retention - Overall Gene Retention') + theme(legend.position='NULL') + geom_text(aes(label=V2))

## what could meiotic enrichment be related to??
ggplot(asize, aes(x=(haploidAssemblySize-haploidNCount)/1e6, y=meiosisEnrichment, group=ploidy, color=ploidy)) + 
  #  geom_boxplot(outlier.shape = NA) + 
  geom_hline(yintercept=0, color='gray', lty='dashed') +
  geom_point(position = position_jitter(seed = 1), size=3)+ 
  scale_color_manual(values=ploidycolors) + #ylim(0,110000) + 
  # ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=110000) + 
  #  geom_hline(yintercept=c(median(asize$diploidEquivalent[asize$ploidy=='Diploid'], na.rm=T)*c(1)), lty='dotted', color='darkgray')+
  #  annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  # annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  # annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  xlab('Haploid Assembly Size (Mb)') + ylab('Meiotic Gene Retention - Overall Gene Retention') + theme(legend.position='NULL') + geom_text(aes(label=V2))
######NOPE

ggplot(asize, aes(x=haploidRepeatSize/1e6, y=meiosisEnrichment, group=ploidy, color=ploidy)) + 
  #  geom_boxplot(outlier.shape = NA) + 
  geom_hline(yintercept=0, color='gray', lty='dashed') +
  geom_point(position = position_jitter(seed = 1), size=3, aes(shape=haploid))+ 
  scale_color_manual(values=ploidycolors) + #ylim(0,110000) + 
  # ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=110000) + 
  #  geom_hline(yintercept=c(median(asize$diploidEquivalent[asize$ploidy=='Diploid'], na.rm=T)*c(1)), lty='dotted', color='darkgray')+
  #  annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  # annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  # annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  xlab('Haploid Repeat Size (Mb)') + ylab('Meiotic Gene Retention - Overall Gene Retention') + theme(legend.position='NULL') + geom_text(aes(label=V2))
### NOPE





## look up interesting genes
meiotic_gene_list <- data.frame(
  name = c("AtML1/AML4", "AtML2", "AtML3/AtML5", "SWI1=DYAD", "MEI1/MCD1"),
  species = c("Ath (AtML1)", "Ath", "Ath (AtML3)", "Ath", "Ath"),
  locus_identifier = c("At5G61960", "At2G42890", "At4G18120", "At5g51330", "At1g77320"),
  reference = c("", "Kaur et al., 2006", "", "Mercier et al., 2001", "Grelon et al., 2003"),
  Protein_feature = c("", "", "5 related genes with overlapping functions", "meiosis set up", "meiosis-specific DNA repair events independent of SPO11-induced DSB"),
  Function = c("", "", "chromatin organization", "", ""),
  orthologues = c("", "", "", "Zm AM1, OsAM1", ""),
  Meiosis_Specific_Somatic_Function = c("somatic", "somatic", "somatic", "Meiosis Spe", "Meiosis Spe"),
  Meiotic_Recombination = c("", "", "", "", "")
)




## at msh7 which is wheat ph2!
## Pavag03G039800

bbd[bbd$gene=='Pavag03G039800.1.v3.1',]
