library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(rtracklayer)



asize=fread('../general_summaries/panand_assembly_sizes.txt', header=T, quote="", fill=T)
## factor so correct order
asize$ploidy=factor(asize$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))

ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')

## A is ks distribution for polyploids
## generated in syntenic_anchors/polyploidy_plotting.R

ksplot

## B is chromosome number by ploidy boxplots

chrs=read.table('../general_summaries/panand_chr_counts.txt', header=T)

asize$chrcount=chrs$haploidchr[match(asize$V2, chrs$sixlettercode)]


cpl=ggplot(asize, aes(x=ploidy, y=chrcount, group=ploidy, color=ploidy)) + 
#  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values=ploidycolors) + ylim(0,32) + 
  ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=31.5) + 
  geom_hline(yintercept=c(median(asize$chrcount[asize$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
  annotate("text", x = 4.45, y = median(asize$chrcount[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  annotate("text", x = 4.45, y = median(asize$chrcount[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  annotate("text", x = 4.45, y = median(asize$chrcount[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  geom_point(position = position_jitter(w=0.3, h=0,seed = 1), size=3)+ 

  xlab('Ploidy') + ylab('Chromosome Number') + theme(legend.position='NULL')


## C is gene or anchor count by ploidy boxplots

b=read.table('../syntenic_anchors/gene_anchorTable_AnchorWave_Paspalum.txt', header=T)
b=b[,-which(colnames(b)%in%c('svirid',"pprate", "osativ", "eophiu", "bdista", "agerjg", 'tdactm', 'tzopol'))]
b=b[,-which(colnames(b)%in%c('tdacs2', 'tdacn2'))]
bb=b[rowSums(b[,-1]>0)>=32,]
asize$syntAnchors=colSums(bb[,-1])[match(asize$V2, names(colSums(bb[,-1])))]
asize$syntAnchorsCount=colSums(bb[,-1]>0)[match(asize$V2, names(colSums(bb[,-1]>0)))]
asize$syntAnchors=ifelse(asize$haploid, asize$syntAnchors, asize$syntAnchors/2)

spl=ggplot(asize, aes(x=ploidy, y=syntAnchors, group=ploidy, color=ploidy)) + 
#  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(seed = 1), size=3)+ 
  scale_color_manual(values=ploidycolors) + ylim(0,30000) + 
  ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=30000) + 
  geom_hline(yintercept=c(median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
  annotate("text", x = 4.45, y = median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  annotate("text", x = 4.45, y = median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  annotate("text", x = 4.45, y = median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  xlab('Ploidy') + ylab('Syntenic Genes') + theme(legend.position='NULL')

### helixer genes
ogs=read.table('~/Downloads/genespace_dotplots/Orthogroups.GeneCount.tsv', header=T)
asize$helixerCount=colSums(ogs[,-1])[asize$V2]
asize$helixerCount=ifelse(asize$haploid, asize$helixerCount, asize$helixerCount/2)

hpl=ggplot(asize, aes(x=ploidy, y=helixerCount, group=ploidy, color=ploidy)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(seed = 1), size=3)+ 
  scale_color_manual(values=ploidycolors) + ylim(0,135000) + 
  ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=130000) + 
  geom_hline(yintercept=c(median(asize$helixerCount[asize$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
  annotate("text", x = 4.45, y = median(asize$helixerCount[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  annotate("text", x = 4.45, y = median(asize$helixerCount[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  annotate("text", x = 4.45, y = median(asize$helixerCount[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  xlab('Ploidy') + ylab('Helixer Genes') + theme(legend.position='NULL')


## D is TE/repeat bp by ploidy boxplots

tpl=ggplot(asize, aes(x=ploidy, y=haploidRepeatSize/1e6, group=ploidy, color=ploidy)) + 
#  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(seed = 1), size=3)+ 
  scale_color_manual(values=ploidycolors) + ylim(0,3.8e3) + 
  ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=3.75e3) + 
  geom_hline(yintercept=c(median(c(asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3,4,5,6,7)), lty='dotted', color='darkgray')+
  annotate("text", x = 4.45, y = median(c(asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  annotate("text", x = 4.45, y = median(c(asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  annotate("text", x = 4.45, y = median(c(asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  annotate("text", x = 4.45, y = median(c(asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T)*4, label = "\u00D74", vjust = -0.5) + 
  annotate("text", x = 4.45, y = median(c(asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T)*5, label = "\u00D75", vjust = -0.5) + 
  annotate("text", x = 4.45, y = median(c(asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T)*6, label = "\u00D76", vjust = -0.5) + 
  annotate("text", x = 4.45, y = median(c(asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T)*7, label = "\u00D77", vjust = -0.5) + 
  xlab('Ploidy') + ylab('Repeat Megabases (Mb)') + theme(legend.position='NULL')




## bcd
threerules=plot_grid(cpl+ theme(axis.text.x=element_blank(), axis.title.x=element_blank()), spl+ theme(axis.text.x=element_blank(), axis.title.x=element_blank()), tpl+ theme(axis.text.x=element_text(size = 9)), ncol=1, labels=c('B', 'C', 'D'), align='v')

figure2=plot_grid(ksplot, threerules, ncol=2, labels=c('A',''), rel_widths=c(1,0.8), align='v', axis='b')


pdf('../figures/figure2_so-many-polyploids.pdf',18,8)
figure2
dev.off()


