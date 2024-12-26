
library(ggtext)




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



#### barplot!!
## generated in ../transposable_elements/plot_te_summaries_along_chr.R
#### WATCH oUYT
### changing this here so parv is shortened!!!!! for TILs
shorttaxonnames=c("Z. mays ssp. parv. TIL11",  "Z. mays ssp. parv. TIL01", "Z. mays ssp. mays B73v5", "Z. mays ssp. mex. TIL25", "Z. mays ssp. mex. TIL18", "Z. mays ssp. huehue.", 
                  "Z. luxurians", "Z. nicaraguensis", "Z. diploperennis Momo", "Z. diploperennis Gigi", "T. zoloptense", "T. dactyloides FL", "T. dactyloides Southern Hap2", 
                  "T. dactyloides Northern Hap2", "T. dactyloides KS", "T. dactyloides tetraploid", "U. digitatum", "V. cuspidata", "R. rottboellioides", "R. tuberculosa", 
                  "H. compressa", "E. tripsacoides", "S. scoparium", "S. microstachyum", "A. virginicum", "A. chinensis", "A. gerardi", 
                  "C. refractus", "C. citratus", "H. contortus", "T. triandra", "B. laguroides", "P. paniceum", "S. bicolor", 
                  "I. rugosum", "S. nutans", '"A." burmanicus', "T. elegans", "C. serrulatus", "P. vaginatum")
names(shorttaxonnames)=c("zTIL11",  "zTIL01", "zmB735", "zTIL25", "zTIL18", "zmhuet", 
                         "zluxur", "znicar", "zdmomo", "zdgigi", "tzopol", "tdacs1", "tdacs2", 
                         "tdacn2", "tdacn1", "tdactm", "udigit", "vcuspi", "rrottb", "rtuber", 
                         "hcompr", "etrips", "sscopa", "smicro", "avirgi", "achine", "agerar", 
                         "crefra", "ccitra", "hconto", "ttrian", "blagur", "ppanic", "sbicol", 
                         "irugos", "snutan", "atenui", "telega", "cserru", "pvagin")

# List of substrings to exclude from italics
exceptions <- c('subsp.', 'ssp.', 'TIL11', 'TIL01', 
                'TIL25', 'TIL18', 'Momo', 'Gigi', 
                'Southern Hap1', 'Northern Hap1', 
                'FL', 'KS', '\\*', '\\"', 'B73v5')

# Custom function to apply italics with exceptions
format_labels <- function(labels, exceptions) {
  sapply(labels, function(label) {
    # Apply italics to everything initially
    parts <- unlist(strsplit(label, " "))  # Split label into words
    formatted_parts <- sapply(parts, function(word) {
      if (word %in% exceptions) {
        return(word)  # Keep exceptions in regular text
      } else {
        return(paste0("<i>", word, "</i>"))  # Italicize other words
      }
    })
    # Combine parts back into a single string
    return(paste(formatted_parts, collapse = " "))
  })
} 

### okay this will work!!! make a stacked bar plot of each of these genome-wide, using these colors
# Retrotransposon colors (pastel reds and pinks)
retro_colors <- c(RLG="#F7B7B7", RLC="#F28E8E", RLX="#E67373")
# DNA transposon colors (cohesive blue gradient)
dna_colors <- c(DHH="#377EB8", DTA="#4B8EBA", DTC="#5DA9BD", DTH="#70C2BF", 
                DTM="#83D6C1", DTT="#96E2C3", DTS="#A8EEC5")
# Tandem repeat color (pastel lavender)
tandem_color <- c(TandemRepeat="#D4B4F4", TR='#D4B4F4')
# Combine all colors
te_colors <- c(retro_colors, dna_colors, tandem_color)


bardata=read.table('../transposable_elements/te_barplot_data.txt', header=T, sep='\t')
bardata$shortSpeciesLabel=shorttaxonnames[match(bardata$genome, names(shorttaxonnames))]
bardata$shortSpeciesLabel=factor(bardata$shortSpeciesLabel, levels=rev(shorttaxonnames))

# Apply the formatting to y-axis labels
formatted_labels <- format_labels(levels(bardata$shortSpeciesLabel), exceptions)
formatted_labels[4]='<i>&quot;A.&quot;</i> <i>burmanicus</i>'

bardata$supLegend <- factor(bardata$sup, levels = c(
  "DHH", "DTA", "DTC", "DTH", "DTM", "DTT", "RLC", "RLG", "RLX", "TandemRepeat"),
  labels = c("DHH", "DTA", "DTC", "DTH", "DTM", "DTT", "RLC", "RLG", "RLX", "TR")  # Change label to TR
)

tebar=ggplot(bardata[!bardata$genome%in%c('zmB735', 'tdacn2', 'tdacs2'),], aes(y=shortSpeciesLabel, x=bp, fill=supLegend)) + geom_col(position='fill') + scale_fill_manual(values=te_colors)+ylab('Genome')+xlab('Proportion of Repeat Base Pairs')+
  geom_vline(xintercept=c(0.25,0.5,0.75), color='snow', alpha=0.1)+labs(fill='TE superfamily')+
  theme(
    legend.position = "top",                # Place legend at bottom
    legend.title = element_text(hjust=0.5, size=9),            
    legend.text = element_text(size = 8),     # Adjust text size
    legend.key.height = unit(0.5, "cm"),       # Reduce height of keys
    legend.key.width = unit(0.5, "cm"),           # Adjust width of keys
    axis.text.y=element_markdown(size=9),
    legend.justification='center'
  ) +
  scale_y_discrete(labels=formatted_labels)+
  guides(fill = guide_legend(
    title.position = "top",                    # Moves legend title to top (optional)
    label.position = "bottom",                    # Places labels above boxes
    nrow = 1,                                  # Arrange in one row
    byrow = TRUE,                               # Fill row by row
    reverse=T
  ))+ylab('')


### distance to genes. 


combined_data_helixer=read.table('../transposable_elements/combined_tedist_helixer_genes.txt', header=T, sep='\t')
combined_data_helixer$ploidy=factor(combined_data_helixer$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))

### upstream
uph=ggplot(combined_data_helixer[combined_data_helixer$windowadj%in%1:200,], aes(group=genome_id, x=(windowadj-200)*100, y=mean, color=ploidy))+ geom_vline(xintercept=c(0:20*10*-100), color='whitesmoke', lty='dotted')+ geom_hline(yintercept=c(0.2,0.4,0.6,0.8), color='seashell2', lty='longdash') + geom_line() + scale_color_manual(values=ploidycolors)  + xlab('Base pairs away from TSS') + ylab('Mean TE base pairs\nin 100bp window') + theme(legend.position='NULL')
### downstream
downh=ggplot(combined_data_helixer[combined_data_helixer$windowadj%in%201:400,], aes(group=genome_id, x=(windowadj-200)*100, y=mean, color=ploidy))+ geom_vline(xintercept=c(0:20*10*100), color='whitesmoke', lty='dotted')+ geom_hline(yintercept=c(0.2,0.4,0.6,0.8), color='seashell2', lty='longdash') + geom_line() + scale_color_manual(values=ploidycolors)  + xlab('Base pairs away from TTS') + ylab('Mean TE base pairs\nin 100bp window') + theme(legend.position='NULL')

genedists=plot_grid(uph, downh, align='hv', axis='tb', ncol=2, labels=c('E', 'F'))



#### types of TEs no burst??

ages=read.table('../transposable_elements/ages_plot_data.txt', header=T, sep='\t')
ageplot=ggplot(ages[ages$genome!='zmB735',], aes(x=mya, y=(1-medianage)/6.5e-9/2/1e6, color=ploidy, size=copies)) + geom_point() + scale_color_manual(values=ploidycolors) + ylab('Median Insertion Time,\nLTR retrotransposon (Mya)')+xlab('Divergence Between\nParental Genomes (Mya)')+
                    theme(legend.position=c(0.55,0.8), legend.title = element_text(hjust=0.5, size=9), legend.text = element_text(size=8))+guides(color='none', fill='none')+ labs(color='Ploidy', size='Total LTR\nRetro. Copies')+
                    guides(size=guide_legend(ncol=1, byrow=T))


evenness=read.table('../transposable_elements/evenness10copies_plot_data.txt', header=T, sep='\t')
evennessplot=ggplot(evenness[evenness$genome!='zmB735',], aes(y=nfam, x=JN, color=ploidy)) + geom_point(aes(size=haploidRepeatSize/1e6)) + scale_color_manual(values=ploidycolors)+
                    xlab('Evenness (1 if all families equally sized)')+ylab('Number Families')+
                    theme(legend.position=c(0.1,0.8), legend.title = element_text(hjust=0.5, size=9), legend.text = element_text(size=8))+guides(color='none', fill='none')+ labs(color='Ploidy', size='Mb of Genomic\nRepeats')+
                    guides(size=guide_legend(ncol=1, byrow=T))



### 

toprowte=plot_grid(rpgg, teploidy, teploidyage, labels='AUTO', align='h', ncol=3)

#plot_grid(toprowte, genedists, ncol=1, rel_heights=c(0.7,1), align='hv')

rightbottom=plot_grid(evennessplot, ageplot, labels=c('G', 'H'), align='h', axis='tb',ncol=2)
rightte=plot_grid(genedists, rightbottom, align='v', axis='lr', ncol=1)

bottomrowte=plot_grid(tebar, rightte, labels=c('D', ''), ncol=2)#, align='h', axis='t')


pdf('../figures/figure5_tes-are-unpredictable.pdf', 14,10)

plot_grid(toprowte, bottomrowte, ncol=1, rel_heights=c(0.4,1), align='vh', axis='lb')

dev.off()


### species-specific figure 5!

#################################################
## generate subfigures highlighting each taxon

for(i in asize$V2){
pdf(paste0('../figures/fig5_by_sp/', i, '_fig5.pdf'), 14,10)

rpggsp=ggplot(asize, aes(x=haploidRepeatSize/(haploidAssemblySize-haploidNCount), y=medGGdist)) + 
  geom_point(aes( color=ploidy), size=3, color='gray') + 
  geom_point(data=asize[asize$V2==i,], aes( fill=ploidy), size=3, pch=21, color='black') + 
  scale_color_manual(values=ploidycolors) + scale_fill_manual(values=ploidycolors) + #geom_text(aes(label=V2, color=ploidy)) + 
  stat_smooth(data=asize[asize$ploidy!='Paleotetraploid',], method='lm', lty='longdash', color='gray', se=F)+ 
#  stat_smooth(data=asize, method='lm',color='gray', se=F) + 
  ylab('Median Syntenic\nGene-Gene Distance') + xlab('Repeat Proportion') + theme(legend.position='NULL')



teploidysp=ggplot(asize, aes(x=ploidy, y=repeatProp, group=ploidy, color=ploidy, fill=ploidy)) + 
#  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values=ploidycolors) + scale_fill_manual(values=ploidycolors) + ylim(0,1) + 
  ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=0.98) + 
  geom_hline(yintercept=c(median(asize$repeatProp[asize$ploidy=='Diploid'], na.rm=T)*c(1,2,2,3)), lty='dotted', color='darkgray')+
#  annotate("text", x = 4.35, y = median(asize$repeatProp[asize$ploidy=='Diploid']-0.07, na.rm=T), label = "Diploid\nMedian", vjust = -0.5) + 
#  annotate("text", x = 4.45, y = median(asize$repeatProp[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
#  annotate("text", x = 4.45, y = median(asize$repeatProp[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  geom_point(position = position_jitter(w=0.3, h=0,seed = 1), color='gray', fill='gray', size=3, pch=25)+ 
  geom_point(data=asize[asize$V2==i,], position = position_jitter(w=0.3, h=0,seed = 1), color='black', size=3, pch=25)+ 
  xlab('Ploidy') + ylab('Repeat Proportion') + theme(legend.position='NULL') +
  scale_x_discrete(labels = c("Diploid" = "Dip.", "Tetraploid" = "Tet.", "Hexaploid" = "Hex.", "Paleotetraploid" = "Paleo."))

teploidyagesp=ggplot(asize, aes(x=mya, y=haploidRepeatSize, color=ploidy)) + geom_vline(xintercept=c(2,4,6,8,10,12,14), color='gray90', lty='dotted', alpha=0.3) + geom_vline(xintercept=c(1,3,5,7,9,11,13), color='gray80', lty='dashed', alpha=0.3)+ 
  geom_point(size=3, color='gray') + 
  geom_point(data=asize[asize$V2==i,],pch=21,aes(fill=ploidy), color='black', size=3) + 
  scale_color_manual(values=ploidycolors)+scale_fill_manual(values=ploidycolors)+ xlab('Divergence between\nParental Subgenomes (Mya)') + 
#  stat_smooth(data=asizezt, method='lm', aes(group=NA), alpha=0.1, se=F, color='gray90')+
  stat_smooth(data=asize[asize$ploidy!='Paleotetraploid' & !asize$V2%in%lowQualAssemblies,], method='lm', lty='longdash', color='gray', se=F)+ 
  stat_smooth(data=asizezt, method='lm',color='gray', se=F) + 
  ylab('Repeat Base Pairs')  + theme(legend.position='NULL')

tebarsp=ggplot(bardata[!bardata$genome%in%c('zmB735', 'tdacn2', 'tdacs2'),], aes(y=shortSpeciesLabel, x=bp, fill=supLegend)) + 
geom_col(position='fill', alpha=0.2) + 
geom_col(data=bardata[bardata$genome==i,], position='fill') + 
scale_fill_manual(values=te_colors)+ylab('Genome')+xlab('Proportion of Repeat Base Pairs')+
  geom_vline(xintercept=c(0.25,0.5,0.75), color='snow', alpha=0.1)+labs(fill='TE superfamily')+
  theme(
    legend.position = "top",                # Place legend at bottom
    legend.title = element_text(hjust=0.5, size=9),            
    legend.text = element_text(size = 8),     # Adjust text size
    legend.key.height = unit(0.5, "cm"),       # Reduce height of keys
    legend.key.width = unit(0.5, "cm"),           # Adjust width of keys
    axis.text.y=element_markdown(size=9),
    legend.justification='center'
  ) +
  scale_y_discrete(labels=formatted_labels)+
  guides(fill = guide_legend(
    title.position = "top",                    # Moves legend title to top (optional)
    label.position = "bottom",                    # Places labels above boxes
    nrow = 1,                                  # Arrange in one row
    byrow = TRUE,                               # Fill row by row
    reverse=T
  ))+ylab('')

### upstream
uphsp=ggplot(combined_data_helixer[combined_data_helixer$windowadj%in%1:200,], aes(group=genome_id, x=(windowadj-200)*100, y=mean, color=ploidy))+ geom_vline(xintercept=c(0:20*10*-100), color='whitesmoke', lty='dotted')+ geom_hline(yintercept=c(0.2,0.4,0.6,0.8), color='seashell2', lty='longdash') + 
		geom_line(color='gray') + 
		geom_line(lwd=2, data=combined_data_helixer[combined_data_helixer$windowadj%in%1:200 & combined_data_helixer$genome_id==i,])+
		scale_color_manual(values=ploidycolors)  + xlab('Base pairs away from TSS') + ylab('Mean TE base pairs\nin 100bp window') + theme(legend.position='NULL')
### downstream
downhsp=ggplot(combined_data_helixer[combined_data_helixer$windowadj%in%201:400,], aes(group=genome_id, x=(windowadj-200)*100, y=mean, color=ploidy))+ geom_vline(xintercept=c(0:20*10*100), color='whitesmoke', lty='dotted')+ geom_hline(yintercept=c(0.2,0.4,0.6,0.8), color='seashell2', lty='longdash') + 
		geom_line(color='gray') + 
		geom_line(lwd=2, data=combined_data_helixer[combined_data_helixer$windowadj%in%201:400 & combined_data_helixer$genome_id==i,])+
	    scale_color_manual(values=ploidycolors)  + xlab('Base pairs away from TTS') + ylab('Mean TE base pairs\nin 100bp window') + theme(legend.position='NULL')

genedistssp=plot_grid(uphsp, downhsp, align='hv', axis='tb', ncol=2, labels=c('E', 'F'))


ageplotsp=ggplot(ages[ages$genome!='zmB735',], aes(x=mya, y=(1-medianage)/6.5e-9/2/1e6, color=ploidy, size=copies)) + 
                    geom_point(color='gray') + 
                    geom_point(data=ages[ages$genome==i,], aes(fill=ploidy), pch=21, color='black') + 
                    scale_color_manual(values=ploidycolors) + scale_fill_manual(values=ploidycolors)+ ylab('Median Insertion Time,\nLTR retrotransposon (Mya)')+xlab('Divergence Between\nParental Genomes (Mya)')+
                    theme(legend.position=c(0.55,0.8), legend.title = element_text(hjust=0.5, size=9), legend.text = element_text(size=8))+guides(color='none', fill='none')+ labs(color='Ploidy', size='Total LTR\nRetro. Copies')+
                    guides(size=guide_legend(ncol=1, byrow=T))


evennessplotsp=ggplot(evenness[evenness$genome!='zmB735',], aes(y=nfam, x=JN, color=ploidy)) + 
                    geom_point(aes(size=haploidRepeatSize/1e6), color='gray') + 
                    geom_point(data=evenness[evenness$genome==i,], aes(size=haploidRepeatSize/1e6, fill=ploidy), color='black', pch=21) + 
                    scale_color_manual(values=ploidycolors)+ scale_fill_manual(values=ploidycolors)+
                    xlab('Evenness (1 if all families equally sized)')+ylab('Number Families')+
                    theme(legend.position=c(0.1,0.8), legend.title = element_text(hjust=0.5, size=9), legend.text = element_text(size=8))+guides(color='none', fill='none')+ labs(color='Ploidy', size='Mb of Genomic\nRepeats')+
                    guides(size=guide_legend(ncol=1, byrow=T))



### 

toprowtesp=plot_grid(rpggsp, teploidysp, teploidyagesp, labels='AUTO', align='h', ncol=3)

#plot_grid(toprowte, genedists, ncol=1, rel_heights=c(0.7,1), align='hv')

rightbottomsp=plot_grid(evennessplotsp, ageplotsp, labels=c('G', 'H'), align='h', axis='tb',ncol=2)
righttesp=plot_grid(genedistssp, rightbottomsp, align='v', axis='lr', ncol=1)

bottomrowtesp=plot_grid(tebarsp, righttesp, labels=c('D', ''), ncol=2)#, align='h', axis='t')


figure5sp=plot_grid(toprowtesp, bottomrowtesp, ncol=1, rel_heights=c(0.4,1), align='vh', axis='lb')


print(figure5sp)
dev.off()

}





