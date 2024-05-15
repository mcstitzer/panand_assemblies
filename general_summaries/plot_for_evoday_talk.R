


ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')


b=read.table('../syntenic_anchors/gene_anchorTable_AnchorWave_Paspalum.txt', header=T)
b=b[,-which(colnames(b)%in%c('zdmomo',"pprate", "osativ", "eophiu", "bdista", "agerjg"))]

## themeda might loose abn arm like contortus?
## serrulatus is nanopore
b=b[,-which(colnames(b)%in%c('ttrian', 'cserru'))]


table(rowSums(b[,-1]>0))
bb=b[rowSums(b[,-1]>0)>=35,]

asize=fread('../general_summaries/panand_assembly_sizes.txt', header=T, quote="", fill=T)

asize$syntAnchors=colSums(bb[,-1])[match(asize$V2, names(colSums(bb[,-1])))]
asize$syntAnchorsCount=colSums(bb[,-1]>0)[match(asize$V2, names(colSums(bb[,-1]>0)))]
asize$syntAnchors=ifelse(asize$haploid, asize$syntAnchors, asize$syntAnchors/2)
## dibn't change this one! this is already corrected for haploidy implicitly
#asize$syntAnchorsCount=ifelse(asize$haploid, asize$syntAnchorsCount, asize$syntAnchorsCount/2)

ggplot(asize, aes(x=ploidy, y=syntAnchors, group=ploidy, color=ploidy)) + geom_boxplot(outlier.shape = NA) + geom_jitter() + scale_color_manual(values=ploidycolors) 
ggplot(asize, aes(x=ploidy, y=syntAnchorsCount, group=ploidy, color=ploidy)) + geom_boxplot(outlier.shape = NA) + geom_jitter() + scale_color_manual(values=ploidycolors) 


ggplot(asize, aes(x=factor(ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid')), y=syntAnchors, group=ploidy, color=ploidy)) + geom_boxplot(outlier.shape = NA) + geom_jitter() + scale_color_manual(values=ploidycolors) + ylim(0,24000) + geom_hline(yintercept=c(median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')

asize$ploidy=factor(asize$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))
asize$diploidEquivalent=ifelse(asize$ploidy=='Diploid', asize$syntAnchors, ifelse(asize$ploidy%in%c('Tetraploid', 'Paleotetraploid'), asize$syntAnchors/2, ifelse(asize$ploidy=='Hexaploid', asize$syntAnchors/3, NA)))

set.seed(1)
pdf('~/Downloads/syntenic_anchors_evoday.pdf',8,4)
ggplot(asize, aes(x=ploidy, y=syntAnchors, group=ploidy, color=ploidy)) + geom_boxplot(outlier.shape = NA) +geom_point(position = position_jitter(seed = 1), size=3) + scale_color_manual(values=ploidycolors) + ylim(0,25000) + ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=25000) #+ geom_hline(yintercept=c(median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')
ggplot(asize, aes(x=ploidy, y=syntAnchors, group=ploidy, color=ploidy)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(seed = 1), size=3)+ scale_color_manual(values=ploidycolors) + ylim(0,25000) + ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=25000) + geom_hline(yintercept=c(median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')
ggplot(asize[asize$V2!='zdgigi',], aes(x=ploidy, y=diploidEquivalent, group=ploidy, color=ploidy)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(seed = 1), size=3)+ scale_color_manual(values=ploidycolors) + ylim(0,10000) + ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=10000) + geom_hline(yintercept=c(median(asize$syntAnchors[asize$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')

dev.off()



## now add ks to these values
ks=fread('../general_summaries/ks_to_look_for_mixtures.txt', header=T, quote='')

ks=data.frame(ks)[ks$ks>0.001,]

mks=ks[ks$ploidy %in%c('Tetraploid', 'Paleotetraploid', 'Hexaploid'),] %>% group_by(genome, ploidy, species, haploid) %>% dplyr::summarize(median=median(ks, na.rm=T), nonallelic=median(ks[ks>0.005], na.rm=T))

mks$medianNonAllelicCorr=ifelse(mks$nonallelic-mks$median>0.01, mks$nonallelic, mks$median)

## wei-yun says mode could work, but median is more real and trustworthy - I'll go with that!

ggplot(data.frame(ks)[ks$ploidy %in% c('Tetraploid', 'Paleotetraploid', 'Hexaploid')& !ks$genome%in%c('zdmomo', 'zdgigi', 'zdnica'),], aes(x=ks, color=ploidy, group=genome)) + geom_density() + scale_color_manual(values=ploidycolors)


ggplot(data.frame(ks)[ks$ploidy %in% c('Tetraploid', 'Paleotetraploid', 'Hexaploid') & !ks$genome%in%c('zdmomo', 'zdgigi', 'zdnica'),], aes(x=ks, color=ploidy, group=genome)) + geom_density() + scale_color_manual(values=ploidycolors) + facet_wrap(~ploidy, ncol=1, scales='free_y') + geom_vline(data=mks[!mks$genome %in% c('zdmomo', 'zdgigi', 'zdnica'),], aes(xintercept=medianNonAllelicCorr))
ggplot(data.frame(ks)[ks$ploidy %in% c('Tetraploid', 'Paleotetraploid', 'Hexaploid') & !ks$genome%in%c('zdmomo', 'zdgigi', 'zdnica'),], aes(x=ks, color=ploidy, group=genome)) + geom_density() + scale_color_manual(values=ploidycolors) + facet_wrap(~ploidy, ncol=1, scales='free_y') + geom_vline(data=mks[!mks$genome %in% c('zdmomo', 'zdgigi', 'zdnica'),], aes(xintercept=median))

ggplot(data.frame(ks)[ks$ploidy %in% c('Tetraploid', 'Paleotetraploid', 'Hexaploid') & !ks$genome%in%c('zdmomo', 'zdgigi', 'zdnica'),], aes(x=ks, color=ploidy, group=genome)) + geom_density() + scale_color_manual(values=ploidycolors) + facet_wrap(~ploidy, ncol=1, scales='free_y') + geom_vline(data=mks[!mks$genome %in% c('zdmomo', 'zdgigi', 'zdnica'),], aes(xintercept=median), color='darkgray', alpha=0.5)

## examples= 'zmB735', 'udigit', 'achine'
ggplot(data.frame(ks)[ks$ploidy %in% c('Tetraploid', 'Paleotetraploid', 'Hexaploid') & !ks$genome%in%c('zdmomo', 'zdgigi', 'zdnica') & ks$genome%in%c('zmB735', 'udigit', 'vcuspi'),], aes(x=ks, color=ploidy, group=genome)) + geom_density(lwd=2) + scale_color_manual(values=ploidycolors) + facet_wrap(~ploidy, ncol=1, scales='free_y') + geom_vline(data=mks[!mks$genome %in% c('zdmomo', 'zdgigi', 'zdnica') & mks$genome%in%c('zmB735', 'udigit', 'vcuspi'),], aes(xintercept=median), color='darkgray', alpha=0.5) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.background = element_blank(), legend.position='none', strip.text.x = element_blank()) + ylab('')

ggplot(data.frame(ks)[ks$ploidy %in% c('Tetraploid', 'Paleotetraploid', 'Hexaploid') & !ks$genome%in%c('zdmomo', 'zdgigi', 'zdnica') ,], aes(x=ks, color=ploidy, group=genome)) + geom_density(aes(lwd=ifelse(genome%in%c('zmB735', 'udigit', 'vcuspi'), 2,1))) + scale_color_manual(values=ploidycolors) + facet_wrap(~ploidy, ncol=1, scales='free_y') + geom_vline(data=mks[!mks$genome %in% c('zdmomo', 'zdgigi', 'zdnica'),], aes(xintercept=median), color='darkgray', alpha=0.5) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.background = element_blank(), legend.position='none', strip.text.x = element_blank()) + ylab('')
ggplot(data.frame(ks)[ks$ploidy %in% c('Tetraploid', 'Paleotetraploid', 'Hexaploid') & !ks$genome%in%c('zdmomo', 'zdgigi', 'zdnica') ,], aes(x=ks, color=ploidy, group=genome)) + geom_density() + scale_color_manual(values=ploidycolors) + facet_wrap(~ploidy, ncol=1, scales='free_y') + geom_vline(data=mks[!mks$genome %in% c('zdmomo', 'zdgigi', 'zdnica'),], aes(xintercept=median), color='darkgray', alpha=0.5) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.background = element_blank(), legend.position='none', strip.text.x = element_blank()) + ylab('')


pdf('~/Downloads/ks_distributions.pdf', 8,4)
ggplot(data.frame(ks)[ks$ploidy %in% c('Tetraploid', 'Paleotetraploid', 'Hexaploid') & !ks$genome%in%c('zdmomo', 'zdgigi', 'zdnica') & ks$genome%in%c('zmB735', 'udigit', 'vcuspi'),], aes(x=ks, color=ploidy, group=genome)) + geom_density(lwd=2) + scale_color_manual(values=ploidycolors) + facet_wrap(~ploidy, ncol=1, scales='free_y') + geom_vline(data=mks[!mks$genome %in% c('zdmomo', 'zdgigi', 'zdnica') & mks$genome%in%c('zmB735', 'udigit', 'vcuspi'),], aes(xintercept=median), color='darkgray', alpha=0.5) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.background = element_blank(), legend.position='none', strip.text.x = element_blank()) + ylab('')
ggplot(data.frame(ks)[ks$ploidy %in% c('Tetraploid', 'Paleotetraploid', 'Hexaploid') & !ks$genome%in%c('zdmomo', 'zdgigi', 'zdnica') & ks$genome%in%c('zmB735', 'udigit', 'vcuspi'),], aes(x=mya, color=ploidy, group=genome)) + geom_density(lwd=2) + scale_color_manual(values=ploidycolors) + facet_wrap(~ploidy, ncol=1, scales='free_y') + geom_vline(data=mks[!mks$genome %in% c('zdmomo', 'zdgigi', 'zdnica') & mks$genome%in%c('zmB735', 'udigit', 'vcuspi'),], aes(xintercept=mya), color='darkgray', alpha=0.5) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.background = element_blank(), legend.position='none', strip.text.x = element_blank()) + ylab('')
ggplot(data.frame(ks)[ks$ploidy %in% c('Tetraploid', 'Paleotetraploid', 'Hexaploid') & !ks$genome%in%c('zdmomo', 'zdgigi', 'zdnica') ,], aes(x=ks, color=ploidy, group=genome)) + geom_density() + scale_color_manual(values=ploidycolors) + facet_wrap(~ploidy, ncol=1, scales='free_y') + geom_vline(data=mks[!mks$genome %in% c('zdmomo', 'zdgigi', 'zdnica'),], aes(xintercept=median), color='darkgray', alpha=0.5) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.background = element_blank(), legend.position='none', strip.text.x = element_blank()) + ylab('')
ggplot(data.frame(ks)[ks$ploidy %in% c('Tetraploid', 'Paleotetraploid', 'Hexaploid') & !ks$genome%in%c('zdmomo', 'zdgigi', 'zdnica') ,], aes(x=mya, color=ploidy, group=genome)) + geom_density() + scale_color_manual(values=ploidycolors) + facet_wrap(~ploidy, ncol=1, scales='free_y') + geom_vline(data=mks[!mks$genome %in% c('zdmomo', 'zdgigi', 'zdnica'),], aes(xintercept=mya), color='darkgray', alpha=0.5) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), strip.background = element_blank(), legend.position='none', strip.text.x = element_blank()) + ylab('')
dev.off()


mu=6.5e-9
mks$mya=mks$medianNonAllelicCorr/2/mu/1e6
ks$mya=ks$ks/2/mu/1e6
ggplot(data.frame(ks)[ks$ploidy %in% c('Tetraploid', 'Paleotetraploid', 'Hexaploid') & !ks$genome%in%c('zdmomo', 'zdgigi', 'zdnica'),], aes(x=mya, color=ploidy, group=genome)) + geom_density() + scale_color_manual(values=ploidycolors) + facet_wrap(~ploidy, ncol=1, scales='free_y') + geom_vline(data=mks[!mks$genome %in% c('zdmomo', 'zdgigi', 'zdnica'),], aes(xintercept=mya))


asize$medianKs=mks$medianNonAllelicCorr[match(asize$V2, mks$genome)]
asize$mya=mks$mya[match(asize$V2, mks$genome)]
asize$mya[asize$ploidy=='Diploid']=0
asize$medianKs[asize$ploidy=='Diploid']=0


ggplot(asize[asize$V2!='zdgigi',], aes(x=mya, y=diploidEquivalent, color=ploidy)) + geom_point(size=3) + scale_color_manual(values=ploidycolors)+ ylim(0,10000)+ geom_text(aes(label=V2))

summary(lm(diploidEquivalent~mya,data=asize[asize$V2!='zdgigi' & asize$ploidy%in%c('Diploid', 'Tetraploid', 'Hexaploid'),]))
summary(lm(diploidEquivalent~mya,data=asize[asize$V2!='zdgigi' ,]))



pdf('~/Downloads/gene_vs_age.pdf',8,4)
ggplot(asize[asize$V2!='zdgigi',], aes(x=mya, y=diploidEquivalent, color=ploidy)) + geom_point(size=3) + scale_color_manual(values=ploidycolors)+ ylim(0,10000) + theme(legend.position='none') #+ stat_smooth(aes(group=1), color='gray10', method='lm', se=F) 
ggplot(asize[asize$V2!='zdgigi',], aes(x=mya, y=diploidEquivalent, color=ploidy)) + geom_point(size=3) + scale_color_manual(values=ploidycolors)+ ylim(0,10000)+ theme(legend.position='none') + stat_smooth(aes(group=1), color='gray30', method='lm', se=F) 
dev.off()
## see who
ggplot(asize[asize$V2!='zdgigi',], aes(x=mya, y=diploidEquivalent, color=ploidy)) + geom_point(size=3) + geom_text(aes(label=V2))+ scale_color_manual(values=ploidycolors)+ ylim(0,10000)+ theme(legend.position='none') + stat_smooth(aes(group=1), color='gray30', method='lm', se=F) 


### tes

asize$diploidEquivalentTE=ifelse(asize$ploidy=='Diploid', asize$haploidRepeatSize, ifelse(asize$ploidy%in%c('Tetraploid', 'Paleotetraploid'), asize$haploidRepeatSize/2, ifelse(asize$ploidy=='Hexaploid', asize$haploidRepeatSize/3, NA)))


pdf('~/Downloads/te_boxplots_evoday.pdf',8,4)
ggplot(asize, aes(x=ploidy, y=haploidRepeatSize/1e6, group=ploidy, color=ploidy)) + geom_boxplot(outlier.shape=NA) +geom_point(position = position_jitter(seed = 1), shape=25, size=2)+ scale_color_manual(values=ploidycolors) + theme(legend.position='none')+ ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=4000)
ggplot(asize, aes(x=ploidy, y=haploidRepeatSize/1e6, group=ploidy, color=ploidy)) + geom_boxplot(outlier.shape=NA) +geom_point(position = position_jitter(seed = 1), shape=25, size=2)+ scale_color_manual(values=ploidycolors) + theme(legend.position='none')+ ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=4000)+ geom_hline(yintercept=c(median((asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')
ggplot(asize, aes(x=ploidy, y=diploidEquivalentTE/1e6, group=ploidy, color=ploidy)) + geom_boxplot(outlier.shape=NA) +geom_point(position = position_jitter(seed = 1), shape=25, size=2)+ scale_color_manual(values=ploidycolors) + theme(legend.position='none')+ ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=2000)+ geom_hline(yintercept=c(median((asize$haploidRepeatSize/1e6)[asize$ploidy=='Diploid'], na.rm=T)*c(1)), lty='dotted', color='darkgray')
dev.off()




ggplot(asize, aes(x=mya, y=diploidEquivalentTE/1e6, color=ploidy)) + geom_point() + scale_color_manual(values=ploidycolors)

# Creating a data frame in R with the precise values
ploidy <- c("Diploid", "Diploid", "Diploid", "Diploid", "Diploid", "Diploid", "Diploid",
            "Tetraploid", "Tetraploid", "Tetraploid", "Tetraploid", "Tetraploid", "Tetraploid", "Tetraploid",
            "Paleotetraploid", "Paleotetraploid", "Paleotetraploid", "Paleotetraploid", "Paleotetraploid",
            "Hexaploid", "Hexaploid", "Hexaploid", "Hexaploid", "Hexaploid")

point_values <- c(0.840, 0.850, 0.856, 0.862, 0.866, 0.869, 0.880,
                  0.843, 0.846, 0.852, 0.856, 0.858, 0.861, 0.864,
                  0.845, 0.848, 0.852, 0.854, 0.858,
                  0.840, 0.845, 0.848, 0.853, 0.858)

# Create the data frame
data_frame <- data.frame(ploidy = ploidy, point = point_values)

ggplot(data_frame, aes(x=ploidy, y=point, group=ploidy, color=ploidy)) + geom_boxplot(outlier.shape=NA) +geom_point(position = position_jitter(seed = 1),size=2)+ scale_color_manual(values=ploidycolors) + theme(legend.position='none')+ ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=0.9)
