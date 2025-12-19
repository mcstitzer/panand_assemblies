library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(rtracklayer)


ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Pentaploid', 'Paleotetraploid')



a=read.table('panand_assembly_stats.txt', header=T, sep='\t')
a$simpleploidy=gsub('Reduced','',a$ploidy)
a$haploidAssemblySize=ifelse(a$haploid, a$assemblysize-a$assemblyn, (a$assemblysize-a$assemblyn)/2)
a$haploidChrCount=ifelse(a$haploid, a$chrcount, a$chrcount/2)
a$haploidHelixer=ifelse(a$haploid, a$helixer, a$helixer/2)

ploidyorder=c('Diploid', 'DiploidReduced', 'Tetraploid','TetraploidReduced','Paleotetraploid', 'Pentaploid','Hexaploid','HexaploidReduced')

## revised
  
chrcountsp=ggplot(a, aes(x=simpleploidy, y=haploidChrCount, group=simpleploidy, color=simpleploidy)) + 
#  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values=ploidycolors) + ylim(0,32) + 
  ggpubr::stat_compare_means(aes(group=simpleploidy, x=simpleploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=31.5) + 
  geom_hline(yintercept=c(median(a$haploidChrCount[a$simpleploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
  annotate("text", x = 8.45, y = median(a$haploidChrCount[a$simpleploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  annotate("text", x = 8.45, y = median(a$haploidChrCount[a$simpleploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  annotate("text", x = 8.45, y = median(a$haploidChrCount[a$simpleploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  geom_point(position = position_jitter(w=0.3, h=0,seed = 1), size=3)+ 
  xlab('Ploidy') + ylab('Chromosome Number') + theme(legend.position='NULL')
chrcountep=ggplot(a[a$ploidy%in%ploidyorder,], aes(x=factor(ploidy, levels=ploidyorder), y=haploidChrCount, group=ploidy, color=simpleploidy)) + 
#  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values=ploidycolors) + ylim(0,32) + 
  ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=31.5) + 
  geom_hline(yintercept=c(median(a$haploidChrCount[a$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
  annotate("text", x = 8.45, y = median(a$haploidChrCount[a$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  annotate("text", x = 8.45, y = median(a$haploidChrCount[a$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  annotate("text", x = 8.45, y = median(a$haploidChrCount[a$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  geom_point(position = position_jitter(w=0.3, h=0,seed = 1), size=3)+ 
  xlab('Ploidy') + ylab('Chromosome Number') + theme(legend.position='NULL', axis.text.x = element_text(angle = 45, size=8,hjust = 1, vjust = 1))+
  geom_text(aes(label=six)) + 
  scale_x_discrete(labels=c('Diploid','Diploid\n(Chr Reduced)', 'Tetraploid', 'Tetraploid\n(Chr Reduced)', 'Paleotetraploid','Pentaploid','Hexaploid','Hexaploid\n(Chr Reduced)'))

### genome size

gssp=ggplot(a, aes(x=simpleploidy, y=haploidAssemblySize/1e6, group=simpleploidy, color=simpleploidy)) + 
#  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values=ploidycolors) + ylim(0,5000) + 
  ggpubr::stat_compare_means(aes(group=simpleploidy, x=simpleploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=4800) + 
  geom_hline(yintercept=c(median(a$haploidAssemblySize[a$simpleploidy=='Diploid']/1e6, na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
  annotate("text", x = 8.45, y = median(a$haploidAssemblySize[a$simpleploidy=='Diploid']/1e6, na.rm=T), label = "\u00D71", vjust = -0.5) + 
  annotate("text", x = 8.45, y = median(a$haploidAssemblySize[a$simpleploidy=='Diploid']/1e6, na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  annotate("text", x = 8.45, y = median(a$haploidAssemblySize[a$simpleploidy=='Diploid']/1e6, na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  geom_point(position = position_jitter(w=0.3, h=0,seed = 1), size=3)+ 

  xlab('Ploidy') + ylab('Assembly Size') + theme(legend.position='NULL')

gsep=ggplot(a[a$ploidy%in%ploidyorder,], aes(x=factor(ploidy, levels=ploidyorder), y=haploidAssemblySize/1e6, group=ploidy, color=simpleploidy)) + 
#  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values=ploidycolors) + ylim(0,5000) + 
  ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=4800) + 
  geom_hline(yintercept=c(median(a$haploidAssemblySize[a$simpleploidy=='Diploid']/1e6, na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
  annotate("text", x = 8.45, y = median(a$haploidAssemblySize[a$simpleploidy=='Diploid']/1e6, na.rm=T), label = "\u00D71", vjust = -0.5) + 
  annotate("text", x = 8.45, y = median(a$haploidAssemblySize[a$simpleploidy=='Diploid']/1e6, na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  annotate("text", x = 8.45, y = median(a$haploidAssemblySize[a$simpleploidy=='Diploid']/1e6, na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  geom_point(position = position_jitter(w=0.3, h=0,seed = 1), size=3)+ 
  xlab('Ploidy') + ylab('Assembly Size') + theme(legend.position='NULL', axis.text.x = element_text(angle = 45, size=8, hjust = 1, vjust = 1))+
  geom_text(aes(label=six)) + 
  scale_x_discrete(labels=c('Diploid','Diploid\n(Chr Reduced)', 'Tetraploid', 'Tetraploid\n(Chr Reduced)', 'Paleotetraploid','Pentaploid','Hexaploid','Hexaploid\n(Chr Reduced)'))

## helixer genes
hgsp=ggplot(a, aes(x=simpleploidy, y=haploidHelixer, group=simpleploidy, color=simpleploidy)) + 
#  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values=ploidycolors) + ylim(0,150000) + 
  ggpubr::stat_compare_means(aes(group=simpleploidy, x=simpleploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=140000) + 
  geom_hline(yintercept=c(median(a$haploidHelixer[a$simpleploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
  annotate("text", x = 8.45, y = median(a$haploidHelixer[a$simpleploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  annotate("text", x = 8.45, y = median(a$haploidHelixer[a$simpleploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  annotate("text", x = 8.45, y = median(a$haploidHelixer[a$simpleploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  geom_point(position = position_jitter(w=0.3, h=0,seed = 1), size=3)+ 
  xlab('Ploidy') + ylab('Helixer Genes') + theme(legend.position='NULL')

hgep=ggplot(a[a$ploidy%in%ploidyorder,], aes(x=factor(ploidy, levels=ploidyorder), y=haploidHelixer, group=ploidy, color=simpleploidy)) + 
#  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values=ploidycolors) + ylim(0,150000) + 
  ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=140000) + 
  geom_hline(yintercept=c(median(a$haploidHelixer[a$simpleploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
  annotate("text", x = 8.45, y = median(a$haploidHelixer[a$simpleploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  annotate("text", x = 8.45, y = median(a$haploidHelixer[a$simpleploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  annotate("text", x = 8.45, y = median(a$haploidHelixer[a$simpleploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  geom_point(position = position_jitter(w=0.3, h=0,seed = 1), size=3)+ 
  xlab('Ploidy') + ylab('Helixer Genes') + theme(legend.position='NULL', axis.text.x = element_text(angle = 45,size=8, hjust = 1, vjust = 1))+
  geom_text(aes(label=six)) + 
  scale_x_discrete(labels=c('Diploid','Diploid\n(Chr Reduced)', 'Tetraploid', 'Tetraploid\n(Chr Reduced)', 'Paleotetraploid','Pentaploid','Hexaploid','Hexaploid\n(Chr Reduced)'))



## syntenic genes
b=read.table('../syntenic_anchors/gene_anchorTable_AnchorWave_Paspalum_zeachronly_revision.txt', header=T)
bb=b[rowSums(b[,-1]>0)>=29,] ## 90%... change as asseblies finish
a$syntAnchors=colSums(bb[,-1])[match(a$six, names(colSums(bb[,-1])))]
a$syntAnchorsCount=colSums(bb[,-1]>0)[match(a$six, names(colSums(bb[,-1]>0)))]
a$syntAnchors=ifelse(a$haploid, a$syntAnchors, a$syntAnchors/2)

sysp=ggplot(a, aes(x=simpleploidy, y=syntAnchors, group=simpleploidy, color=simpleploidy)) + 
#  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values=ploidycolors) + ylim(0,25000) + 
  ggpubr::stat_compare_means(aes(group=simpleploidy, x=simpleploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=23000) + 
  geom_hline(yintercept=c(median(a$syntAnchors[a$simpleploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
  annotate("text", x = 8.45, y = median(a$syntAnchors[a$simpleploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  annotate("text", x = 8.45, y = median(a$syntAnchors[a$simpleploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  annotate("text", x = 8.45, y = median(a$syntAnchors[a$simpleploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  geom_point(position = position_jitter(w=0.3, h=0,seed = 1), size=3)+ 
  xlab('Ploidy') + ylab('Syntenic Genes') + theme(legend.position='NULL')

syep=ggplot(a[a$ploidy%in%ploidyorder,], aes(x=factor(ploidy, levels=ploidyorder), y=syntAnchors, group=ploidy, color=simpleploidy)) + 
#  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values=ploidycolors) + ylim(0,25000) + 
  ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=23000) + 
  geom_hline(yintercept=c(median(a$syntAnchors[a$simpleploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
  annotate("text", x = 8.45, y = median(a$syntAnchors[a$simpleploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  annotate("text", x = 8.45, y = median(a$syntAnchors[a$simpleploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  annotate("text", x = 8.45, y = median(a$syntAnchors[a$simpleploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  geom_point(position = position_jitter(w=0.3, h=0,seed = 1), size=3)+ 
  xlab('Ploidy') + ylab('Syntenic Genes') + theme(legend.position='NULL', axis.text.x = element_text(angle = 45,size=8, hjust = 1, vjust = 1))+
  geom_text(aes(label=six)) + 
  scale_x_discrete(labels=c('Diploid','Diploid\n(Chr Reduced)', 'Tetraploid', 'Tetraploid\n(Chr Reduced)', 'Paleotetraploid','Pentaploid','Hexaploid','Hexaploid\n(Chr Reduced)'))

sypsp=ggplot(a, aes(x=simpleploidy, y=haploidHelixer, group=simpleploidy, color=simpleploidy)) + 
#  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values=ploidycolors) + ylim(0,150000) + 
  ggpubr::stat_compare_means(aes(group=simpleploidy, x=simpleploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=140000) + 
  geom_hline(yintercept=c(median(a$haploidHelixer[a$simpleploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
  annotate("text", x = 8.45, y = median(a$haploidHelixer[a$simpleploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  annotate("text", x = 8.45, y = median(a$haploidHelixer[a$simpleploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  annotate("text", x = 8.45, y = median(a$haploidHelixer[a$simpleploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  geom_point(position = position_jitter(w=0.3, h=0,seed = 1), size=3)+ 
  xlab('Ploidy') + ylab('Helixer Genes') + theme(legend.position='NULL')

sypep=ggplot(a[a$ploidy%in%ploidyorder,], aes(x=factor(ploidy, levels=ploidyorder), y=haploidHelixer, group=ploidy, color=simpleploidy)) + 
#  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values=ploidycolors) + ylim(0,150000) + 
  ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=140000) + 
  geom_hline(yintercept=c(median(a$haploidHelixer[a$simpleploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
  annotate("text", x = 8.45, y = median(a$haploidHelixer[a$simpleploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  annotate("text", x = 8.45, y = median(a$haploidHelixer[a$simpleploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  annotate("text", x = 8.45, y = median(a$haploidHelixer[a$simpleploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  geom_point(position = position_jitter(w=0.3, h=0,seed = 1), size=3)+ 
  xlab('Ploidy') + ylab('Helixer Genes') + theme(legend.position='NULL', axis.text.x = element_text(angle = 45,size=8, hjust = 1, vjust = 1))+
  geom_text(aes(label=six)) + 
  scale_x_discrete(labels=c('Diploid','Diploid\n(Chr Reduced)', 'Tetraploid', 'Tetraploid\n(Chr Reduced)', 'Paleotetraploid','Pentaploid','Hexaploid','Hexaploid\n(Chr Reduced)'))





plot_grid(chrcountep, hgep, syep, gsep, ncol=4, align='hv', axis='l')


ggplot(a[a$ploidy%in%ploidyorder,], aes(x=similarsubgenomechr, y=haploidHelixer, color=simpleploidy)) + 
geom_point()+scale_color_manual(values=ploidycolors)
ggplot(a[a$ploidy%in%ploidyorder,], aes(x=similarsubgenomechr, y=haploidAssemblySize, color=simpleploidy)) + 
geom_point()+scale_color_manual(values=ploidycolors)

ggplot(a[a$ploidy%in%ploidyorder,], aes(x=similarsubgenomechr, y=haploidChrCount, color=simpleploidy)) + 
geom_point()+scale_color_manual(values=ploidycolors)


## use old parental divergence values
old=read.table('../general_summaries/panand_assembly_sizes.txt', header=T, sep='\t')

a$mya=old$mya[match(a$six, old$V2)]

ggplot(a[a$ploidy%in%ploidyorder,], aes(x=similarsubgenomechr, y=mya, color=simpleploidy)) + 
geom_point()+scale_color_manual(values=ploidycolors)



ggplot(a[a$ploidy%in%ploidyorder,], aes(y=haploidChrCount, shape=grepl('Reduced', ploidy), x=mya, color=simpleploidy)) + 
geom_point()+scale_color_manual(values=ploidycolors)

ggplot(a[a$ploidy%in%ploidyorder,], aes(y=haploidHelixer, shape=grepl('Reduced', ploidy), x=mya, color=simpleploidy)) + 
geom_point()+scale_color_manual(values=ploidycolors)

