





### gene-gene distance
## from ../genes/gene-gene_distance.R
rpgg ## repeat proportion
gsgg ## genome size

###v  te prop vs age colored by polyploidy
## no burst!!

asize$repeatProp=asize$haploidRepeatSize/(asize$haploidAssemblySize-asize$haploidNCount)

## upside down triangles for the points :)
teploidy=ggplot(asize, aes(x=ploidy, y=repeatProp, group=ploidy, color=ploidy, fill=ploidy)) + 
#  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values=ploidycolors) + scale_fill_manual(values=ploidycolors) + ylim(0,1) + 
  ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=0.98) + 
  geom_hline(yintercept=c(median(asize$repeatProp[asize$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
#  annotate("text", x = 4.35, y = median(asize$repeatProp[asize$ploidy=='Diploid']-0.07, na.rm=T), label = "Diploid\nMedian", vjust = -0.5) + 
#  annotate("text", x = 4.45, y = median(asize$repeatProp[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
#  annotate("text", x = 4.45, y = median(asize$repeatProp[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  geom_point(position = position_jitter(w=0.3, h=0,seed = 1), size=3, pch=25)+ 
  xlab('Ploidy') + ylab('Repeat Proportion') + theme(legend.position='NULL')

### haploid bp vs time since polyploidy

asizezt=asize
asizezt$V2zt=ifelse(asize$V2 %in% c('tdacn1', 'tdacs1'), 'trips', asize$V2)
asizezt$V2zt=ifelse(asize$V2 %in% c("zdgigi", "zdmomo", "zluxur", "zmhuet", "zTIL18", "zTIL25", "zTIL01", "zTIL11", "znicar", "zmB735"), 'zea', asize$V2)
asizezt=asizezt%>% group_by(V2zt, ploidy) %>% summarize(haploidRepeatSize=median(haploidRepeatSize), medFrac=median(medFrac), meanFrac=median(meanFrac), syntAnchorsCount=median(syntAnchorsCount), syntAnchors=median(syntAnchors), diploidEquivalentsyntAnchors=median(diploidEquivalentsyntAnchors), mya=median(mya))

summary(lm(asizezt$haploidRepeatSize~asizezt$mya))
## so 110189996 bp per million years!!

teploidyage=ggplot(asize, aes(x=mya, y=haploidRepeatSize, color=ploidy)) + geom_vline(xintercept=c(2,4,6,8,10,12,14), color='gray90', lty='dotted', alpha=0.3) + geom_vline(xintercept=c(1,3,5,7,9,11,13), color='gray80', lty='dashed', alpha=0.3)+ 
  geom_point(size=3) + 
  scale_color_manual(values=ploidycolors)+ xlab('Divergence between\nParental Subgenomes (Mya)') + 
#  stat_smooth(data=asizezt, method='lm', aes(group=NA), alpha=0.1, se=F, color='gray90')+
  stat_smooth(data=asize[asize$ploidy!='Paleotetraploid' & !asize$V2%in%lowQualAssemblies,], method='lm', lty='longdash', color='gray', se=F)+ 
  stat_smooth(data=asizezt, method='lm',color='gray', se=F) + 
  ylab('Repeat Base Pairs')  + theme(legend.position='NULL')
teploidyage



#### types of TEs no burst??

### 

plot_grid(rpgg, teploidy, teploidyage, labels='AUTO', align='hv', ncol=1)