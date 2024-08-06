

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
genepath=read.table('panand_gene_annotations.txt', header=F, sep='\t')

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





##
ggplot(asize, aes(x=mya, y=syntAnchors, group=ploidy, color=ploidy)) + 
#  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(seed = 1), size=3)+ 
  scale_color_manual(values=ploidycolors) + #ylim(0,110000) + 
 # ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=110000) + 
  geom_hline(yintercept=c(median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
#  annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
 # annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
 # annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  xlab('mya') + ylab('anchors') + theme(legend.position='NULL') + geom_text(aes(label=V2))


