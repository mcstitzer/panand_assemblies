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
  scale_color_manual(values=ploidycolors) + ylim(0,30000) + 
  ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=30000) + 
  geom_hline(yintercept=c(median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
  annotate("text", x = 4.45, y = median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  annotate("text", x = 4.45, y = median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  annotate("text", x = 4.45, y = median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  geom_point(position = position_jitter(seed = 1), size=3)+ 
  xlab('Ploidy') + ylab('Syntenic Genes') + theme(legend.position='NULL')

### helixer genes
ogs=read.table('~/Downloads/genespace_dotplots/Orthogroups.GeneCount.tsv', header=T)
asize$helixerCount=colSums(ogs[,-1])[asize$V2]
asize$helixerCount=ifelse(asize$haploid, asize$helixerCount, asize$helixerCount/2)

hpl=ggplot(asize, aes(x=ploidy, y=helixerCount, group=ploidy, color=ploidy)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values=ploidycolors) + ylim(0,135000) + 
  ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=130000) + 
  geom_hline(yintercept=c(median(asize$helixerCount[asize$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
  annotate("text", x = 4.45, y = median(asize$helixerCount[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  annotate("text", x = 4.45, y = median(asize$helixerCount[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  annotate("text", x = 4.45, y = median(asize$helixerCount[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  geom_point(position = position_jitter(seed = 1), size=3)+ 

  xlab('Ploidy') + ylab('Helixer Genes') + theme(legend.position='NULL')


## D is TE/repeat bp by ploidy boxplots

tpl=ggplot(asize, aes(x=ploidy, y=haploidRepeatSize/1e6, group=ploidy, color=ploidy)) + 
#  geom_boxplot(outlier.shape = NA) + 
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
  geom_point(position = position_jitter(seed = 1), size=3)+ 

  xlab('Ploidy') + ylab('Repeat Megabases (Mb)') + theme(legend.position='NULL')



## add in efg for ks vs bcd cats
## het from ??
asize$ks=het$ks[match(asize$V2, het$genome)]


cplks=ggplot(asize, aes(x=ks, y=chrcount, color=ploidy)) + 
#  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values=ploidycolors) + ylim(0,32) + 
  geom_hline(yintercept=c(median(asize$chrcount[asize$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
  annotate("text", x = 0.2, y = median(asize$chrcount[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  annotate("text", x = 0.2, y = median(asize$chrcount[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  annotate("text", x = 0.2, y = median(asize$chrcount[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  geom_point( size=3)+ 

  xlab('Median Ks of Syntenic Homologs') + ylab('Chromosome Number') + theme(legend.position='NULL')


splks=ggplot(asize, aes(x=ks, y=syntAnchors, color=ploidy)) + 
#  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values=ploidycolors) + ylim(0,30000) + 
  geom_hline(yintercept=c(median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
  annotate("text", x = 0.2, y = median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  annotate("text", x = 0.2, y = median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  annotate("text", x = 0.2, y = median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  geom_point(size=3)+ 
  xlab('Median Ks of Syntenic Homologs') + ylab('Syntenic Genes') + theme(legend.position='NULL')

tplks=ggplot(asize, aes(x=ks, y=haploidRepeatSize/1e6, color=ploidy)) + 
#  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values=ploidycolors) + ylim(0,3.8e3) + 
  geom_hline(yintercept=c(median(c(asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3,4,5,6,7)), lty='dotted', color='darkgray')+
  annotate("text", x = 0.2, y = median(c(asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  annotate("text", x = 0.2, y = median(c(asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  annotate("text", x = 0.2, y = median(c(asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  annotate("text", x = 0.2, y = median(c(asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T)*4, label = "\u00D74", vjust = -0.5) + 
  annotate("text", x = 0.2, y = median(c(asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T)*5, label = "\u00D75", vjust = -0.5) + 
  annotate("text", x = 0.2, y = median(c(asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T)*6, label = "\u00D76", vjust = -0.5) + 
  annotate("text", x = 0.2, y = median(c(asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T)*7, label = "\u00D77", vjust = -0.5) + 
  geom_point( size=3)+ 
  xlab('Ploidy') + ylab('Repeat Megabases (Mb)') + theme(legend.position='NULL')


## aapbar from ../gene_trees/count_topologies_genetrees.R
aapbar


## bcd
threerules=plot_grid(cpl+ theme(axis.text.x=element_blank(), axis.title.x=element_blank()), spl+ theme(axis.text.x=element_blank(), axis.title.x=element_blank()), tpl+ theme(axis.text.x=element_text(size = 9)), ncol=1, labels=c('B', 'C', 'D'), align='v')
threerulesabc=plot_grid(cpl+ theme(axis.text.x=element_blank(), axis.title.x=element_blank()), spl+ theme(axis.text.x=element_blank(), axis.title.x=element_blank()), tpl+ theme(axis.text.x=element_text(size = 9)), ncol=1, labels=c('A', 'B', 'C'), align='v')

threeks=plot_grid(cplks+ theme(axis.text.x=element_blank(), axis.title.x=element_blank()), splks+ theme(axis.text.x=element_blank(), axis.title.x=element_blank()), tplks+ theme(axis.text.x=element_text(size = 9)), ncol=1, labels=c('E', 'F', 'G'), align='v')

figure2=plot_grid(ksplot, threerules, ncol=2, labels=c('A',''), rel_widths=c(1,0.8), align='v', axis='b')

### i don't think i like this - leave this relationship for the later rules!
figure2ks=plot_grid( threerulesabc, ksplot,threeks, ncol=3, labels=c('','D', ''), rel_widths=c(0.7,1,0.7), align='v', axis='b')
figure2ks=plot_grid( ksplot,threerules, threeks, ncol=3, labels=c('A','', ''), rel_widths=c(1,0.7,0.7), align='v', axis='b')

figure2right=plot_grid(threerulesabc, 
             plot_grid(ksplotright, aapbar+theme(axis.text.y=element_blank()), ncol=2, rel_widths=c(1,0.4),labels=c('D', 'E'), align='h', axis='tb'),
             ncol=2, rel_widths=c(0.4,1), align='v', axis='tb')




pdf('../figures/figure2_so-many-polyploids.pdf',14,9)
figure2right
dev.off()



#################################################
## generate subfigures highlighting each taxon

for(i in asize$V2){
#ksplot
ksdat=final_merged_data %>% filter(ploidy!='Diploid', !is.na(shortspeciesLabel)) %>% group_by(queryChr, refChr, paspChr, genome,ploidy, speciesLabel, shortspeciesLabel) %>%summarize(ks=mean(V17), n=length(unique(gene))) %>% filter(n>10) 
pdf(paste0('../figures/fig2_by_sp/', i, '_fig2.pdf'), 18,8)


ksplotsp= ggplot(data=ksdat)+ 
  geom_vline(xintercept=seq(from=0,to=0.25,by=0.01), color='gray', lty='dashed', alpha=0.2) + 
  geom_histogram(aes(x=ks,weight=n), fill='gray', color='gray',binwidth=0.001) + 
  geom_histogram(data=ksdat[ksdat$genome==i,],aes(x=ks, color=ploidy, fill=ploidy, weight=n), binwidth=0.001) + 
  facet_wrap(~shortspeciesLabel, ncol=1, scales='free_y', strip.position = 'left', drop = TRUE, labeller=purrr::partial(label_species, dont_italicize=c('subsp.', 'ssp.', 'TIL11', 'TIL01', 'TIL25', 'TIL18', 'Momo', 'Gigi', 'Southern Hap1', 'Northern Hap1', 'FL', 'KS',  '\\*', '\\"', 'B73v5'))) + scale_color_manual(values=ploidycolors)+ scale_fill_manual(values=ploidycolors)+
  theme(strip.background = element_blank(), 
        strip.text.y.left = element_text(angle = 0, hjust=1), 
        panel.spacing = unit(3, "pt"), 
        axis.text.y = element_blank(), 
        axis.text.x = element_text(size = 9),
        legend.position = 'NULL',
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank()) +
  ylab(label = '')+
#  xlim(0,0.25)
  ## from gerardi plot
  scale_x_continuous( "Mean Ks of Syntenic Homologs per Block", limits=c(0,0.25),  sec.axis = sec_axis(~ . /6.5e-9/2/1e6, name = "Syntenic Homolog Divergence\n(million years)"))+
  theme(axis.text.x.top=element_text(color='blue'), axis.title.x.top=element_text(color='blue'))



cplsp=ggplot(asize, aes(x=ploidy, y=chrcount, group=ploidy)) + 
#  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values=c(ploidycolors, 'gray'='gray')) +  scale_fill_manual(values=c(ploidycolors, 'gray'='gray'))+ ylim(0,32) + 
  ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=31.5) + 
  geom_hline(yintercept=c(median(asize$chrcount[asize$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
  annotate("text", x = 4.45, y = median(asize$chrcount[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  annotate("text", x = 4.45, y = median(asize$chrcount[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  annotate("text", x = 4.45, y = median(asize$chrcount[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
#  geom_point(position = position_jitter(w=0.3, h=0,seed = 1), size=3, color='gray')+ 
  geom_point(data=asize %>% mutate(newcolor=ifelse(V2==i,as.character(ploidy),'gray'), size=ifelse(V2==i, 4,3)), position = position_jitter(w=0.3, h=0,seed = 1), aes(size=size, color=newcolor))+ 
  xlab('Ploidy') + ylab('Chromosome Number') + theme(legend.position='NULL')

splsp=ggplot(asize, aes(x=ploidy, y=syntAnchors, group=ploidy)) + 
#  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values=c(ploidycolors, 'gray'='gray')) +  scale_fill_manual(values=c(ploidycolors, 'gray'='gray'))+ ylim(0,30000) + 
  ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=30000) + 
  geom_hline(yintercept=c(median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
  annotate("text", x = 4.45, y = median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  annotate("text", x = 4.45, y = median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  annotate("text", x = 4.45, y = median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
#  geom_point(position = position_jitter(seed = 1), size=3, color='gray')+ 
  geom_point(data=asize %>% mutate(newcolor=ifelse(V2==i,as.character(ploidy),'gray'), size=ifelse(V2==i, 4,3)), position = position_jitter(w=0.3, h=0,seed = 1), aes(size=size, color=newcolor))+ 
  xlab('Ploidy') + ylab('Syntenic Genes') + theme(legend.position='NULL')

tplsp=ggplot(asize, aes(x=ploidy, y=haploidRepeatSize/1e6, group=ploidy, color=ploidy)) + 
#  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values=c(ploidycolors, 'gray'='gray')) +  scale_fill_manual(values=c(ploidycolors, 'gray'='gray'))+ ylim(0,3.8e3) + 
  ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=3.75e3) + 
  geom_hline(yintercept=c(median(c(asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3,4,5,6,7)), lty='dotted', color='darkgray')+
  annotate("text", x = 4.45, y = median(c(asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  annotate("text", x = 4.45, y = median(c(asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  annotate("text", x = 4.45, y = median(c(asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  annotate("text", x = 4.45, y = median(c(asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T)*4, label = "\u00D74", vjust = -0.5) + 
  annotate("text", x = 4.45, y = median(c(asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T)*5, label = "\u00D75", vjust = -0.5) + 
  annotate("text", x = 4.45, y = median(c(asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T)*6, label = "\u00D76", vjust = -0.5) + 
  annotate("text", x = 4.45, y = median(c(asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T)*7, label = "\u00D77", vjust = -0.5) + 
#  geom_point(position = position_jitter(seed = 1), size=3, color='gray')+ 
  geom_point(data=asize %>% mutate(newcolor=ifelse(V2==i,as.character(ploidy),'gray'), size=ifelse(V2==i, 4,3)), position = position_jitter(w=0.3, h=0,seed = 1), aes(size=size, color=newcolor))+ 
  xlab('Ploidy') + ylab('Repeat Megabases (Mb)') + theme(legend.position='NULL')

threerulessp=plot_grid(cplsp+ theme(axis.text.x=element_blank(), axis.title.x=element_blank()), splsp+ theme(axis.text.x=element_blank(), axis.title.x=element_blank()), tplsp+ theme(axis.text.x=element_text(size = 9)), ncol=1, labels=c('B', 'C', 'D'), align='v')
figure2sp=plot_grid(ksplotsp, threerulessp, ncol=2, labels=c('A',''), rel_widths=c(1,0.8), align='v', axis='b')
print(figure2sp)
dev.off()

}

