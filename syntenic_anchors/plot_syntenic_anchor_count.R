library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())



asize=fread('../general_summaries/panand_assembly_sizes.txt', header=T, quote="", fill=T)
## factor so correct order
asize$ploidy=factor(asize$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))

ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')




### similar to figure 2

b=read.table('../syntenic_anchors/gene_anchorTable_AnchorWave_Paspalum.txt', header=T)
b=b[,-which(colnames(b)%in%c('svirid',"pprate", "osativ", "eophiu", "bdista", "agerjg", 'tdactm', 'tzopol'))]
b=b[,-which(colnames(b)%in%c('tdacs2', 'tdacn2'))]
bb=b[rowSums(b[,-1]>0)>=32,]
asize$syntAnchors=colSums(bb[,-1])[match(asize$V2, names(colSums(bb[,-1])))]
asize$syntAnchorsCount=colSums(bb[,-1]>0)[match(asize$V2, names(colSums(bb[,-1]>0)))]
asize$syntAnchors=ifelse(asize$haploid, asize$syntAnchors, asize$syntAnchors/2)

asize$diploidEquivalentsyntAnchors=ifelse(asize$ploidy=='Diploid', asize$syntAnchors, ifelse(asize$ploidy%in%c('Tetraploid', 'Paleotetraploid'), asize$syntAnchors/2, ifelse(asize$ploidy=='Hexaploid', asize$syntAnchors/3, NA)))



asizezt=asize
asizezt$V2zt=ifelse(asize$V2 %in% c('tdacn1', 'tdacs1'), 'trips', asize$V2)
asizezt$V2zt=ifelse(asize$V2 %in% c("zdgigi", "zdmomo", "zluxur", "zmhuet", "zTIL18", "zTIL25", "zTIL01", "zTIL11", "znicar", "zmB735"), 'zea', asize$V2)
asizezt=asizezt%>% group_by(V2zt, ploidy) %>% summarize(syntAnchorsCount=median(syntAnchorsCount), syntAnchors=median(syntAnchors), diploidEquivalentsyntAnchors=median(diploidEquivalentsyntAnchors), mya=median(mya))

cor.test(asizezt$diploidEquivalentsyntAnchors, asizezt$mya)
summary(lm(asizezt$diploidEquivalentsyntAnchors~asizezt$mya))
## negative -0.5081658 
cor.test(asizezt$diploidEquivalentsyntAnchors[!asizezt$ploidy%in%c('Diploid', 'Paleotetraploid')], asizezt$mya[!asizezt$ploidy%in%c('Diploid', 'Paleotetraploid')])
## positive 0.5994824 
summary(lm(asizezt$diploidEquivalentsyntAnchors[!asizezt$ploidy%in%c('Diploid', 'Paleotetraploid')]~asizezt$mya[!asizezt$ploidy%in%c('Diploid', 'Paleotetraploid')])))


fig_fracage = ggplot(asize, aes(x=mya, y=diploidEquivalentsyntAnchors, color=ploidy)) + geom_vline(xintercept=c(2,4,6,8,10,12,14), color='gray90', lty='dotted', alpha=0.3) + geom_vline(xintercept=c(1,3,5,7,9,11,13), color='gray80', lty='dashed', alpha=0.3)+
  geom_point(size=3) + scale_color_manual(values=ploidycolors)+ ylim(0,10000)+ 
  theme(legend.position='none') + 
  xlab('Divergence between\nParental Subgenomes (Mya)') + 
  ylab('Diploid Equivalent\nSyntenic Copies') +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=9))+
  stat_smooth(data=asizezt, method='lm', aes(group=NA), se=F, color='gray80')+
  stat_smooth(data=asizezt[!asizezt$ploidy%in%c('Diploid', 'Paleotetraploid'),], method='lm', aes(group=NA), se=F, color='gray80', lty='longdash')


