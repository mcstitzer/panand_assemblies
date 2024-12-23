





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
  xlab('Ploidy') + ylab('Repeat Proportion') + theme(legend.position='NULL') +
  scale_x_discrete(labels = c("Diploid" = "Dip.", "Tetraploid" = "Tet.", "Hexaploid" = "Hex.", "Paleotetraploid" = "Paleo."))

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


### distance to genes. 


combined_data_helixer=read.table('../transposable_elements/combined_tedist_helixer_genes.txt', header=T, sep='\t')
combined_data_helixer$ploidy=factor(combined_data_helixer$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))

### upstream
uph=ggplot(combined_data_helixer[combined_data_helixer$windowadj%in%1:200,], aes(group=genome_id, x=(windowadj-200)*100, y=mean, color=ploidy))+ geom_vline(xintercept=c(0:20*10*-100), color='whitesmoke', lty='dotted')+ geom_hline(yintercept=c(0.2,0.4,0.6,0.8), color='seashell2', lty='longdash') + geom_line() + scale_color_manual(values=ploidycolors)  + xlab('Base pairs away from TSS') + ylab('Mean TE base pairs in 100bp window') + theme(legend.position='NULL')
### downstream
downh=ggplot(combined_data_helixer[combined_data_helixer$windowadj%in%201:400,], aes(group=genome_id, x=(windowadj-200)*100, y=mean, color=ploidy))+ geom_vline(xintercept=c(0:20*10*100), color='whitesmoke', lty='dotted')+ geom_hline(yintercept=c(0.2,0.4,0.6,0.8), color='seashell2', lty='longdash') + geom_line() + scale_color_manual(values=ploidycolors)  + xlab('Base pairs away from TTS') + ylab('Mean TE base pairs in 100bp window') + theme(legend.position='NULL')

genedists=plot_grid(uph, downh, align='hv', axis='tb', ncol=2, labels=c('D', 'E'))



#### types of TEs no burst??

### 

rowte=plot_grid(rpgg, teploidy, teploidyage, labels='AUTO', align='h', ncol=3)

plot_grid(rowte, genedists, ncol=1, rel_heights=c(0.7,1), align='hv')

pdf('../figures/figure5_tes-are-unpredictable.pdf', 12,8)

plot_grid(rowte, genedists, ncol=1, rel_heights=c(0.7,1), align='hv')

dev.off()
