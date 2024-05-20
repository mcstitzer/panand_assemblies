

asize=fread('../general_summaries/panand_assembly_sizes.txt', header=T, quote="", fill=T)
asize$ploidy=factor(asize$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))

## count rearrangments
diploids=c('cserru', 'irugos', 'sbicol', 'ppanic', 'ttrian', 'crefra', 'avirgi', 'smicro', 'rtuber')
tetraploids=c('snutan', 'hconto', 'ccitra', 'achine', 'sscopa', 'etrips',  'vcuspi')
hexaploids=c('udigit', 'agerar', 'hcompr', 'blagur')

countRearrangements <- function(filepath, color_palette=muted_colors, minBlock=20) {
  # Load data
  data <- read.table(filepath, header = TRUE)
  data <- data[data$gene != 'interanchor', ]
  
  # Reduce to blocks and calculate stats
  data <- data %>%
    group_by(blockIndex) %>%
    mutate(blockLength = n()) %>%
    group_by(queryChr) %>%
    mutate(freqStrand = names(which.max(table(strand))),
           maxChr = max(queryStart),
           freqRef = names(which.max(table(refChr))))
  
  # Filter data based on block length
  data <- data[data$blockLength > minBlock, ]
  data$refChr <- factor(data$refChr, levels = c(paste0('Chr0', 1:9), 'Chr10'))
  
  temp=data %>% filter(blockLength>minBlock) %>% group_by(queryChr)%>% summarize(n=length(unique(refChr)))
  sum(temp$n>1)
}

for(i in diploids){
  print(countRearrangements(Sys.glob(paste0('../syntenic_anchors/anchors/', i, '-Pv-*'))))
}

for(i in tetraploids){
  print(countRearrangements(Sys.glob(paste0('../syntenic_anchors/anchors//', i, '-Pv-*'))))
}
for(i in hexaploids){
  print(countRearrangements(Sys.glob(paste0('../syntenic_anchors/anchors/', i, '-Pv-*')), minBlock=50))
}                              
for(i in c('zmB735')){
  print(countRearrangements(Sys.glob(paste0('../syntenic_anchors/anchors/', i, '-Pv-*')), minBlock=50))
}  

lowQualAssemblies=c('telega', 'atenui', 'rrottb')

asize$transloc=sapply(asize$V2, function(x) countRearrangements(Sys.glob(paste0('../syntenic_anchors/anchors/', x, '-Pv-*')), minBlock=30))
asize$transloc=ifelse(asize$haploid, asize$transloc, asize$transloc/2) ## if haploid, this is true number, if allelic, don't count each allele
asize %>% group_by(ploidy) %>% summarize(mean(transloc, na.rm=T))


ggplot(asize[!asize$V2%in%lowQualAssemblies], aes(x=ploidy, y=transloc, group=ploidy, color=ploidy)) + geom_boxplot(outlier.shape=NA) +geom_point(position = position_jitter(seed = 1),  size=2, aes(shape=haploid))+ scale_color_manual(values=ploidycolors) + theme(legend.position='none')+ ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=10)

asize$scaledTransloc=ifelse(asize$ploidy%in%c('Paleotetraploid', 'Tetraploid'), asize$transloc/2, ifelse(asize$ploidy=='Hexaploid', asize$transloc/3, asize$transloc))
#### haploid or not scaled by triangle/circle
ggplot(asize[!asize$V2%in%lowQualAssemblies,], aes(x=ploidy, y=scaledTransloc, group=ploidy, color=ploidy)) + geom_boxplot(outlier.shape=NA) +geom_point(position = position_jitter(seed = 1),  size=2, aes(shape=haploid))+ scale_color_manual(values=ploidycolors) + theme(legend.position='none')+ ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=6) + ylab('Translocations per haploid equivalent') + xlab('Ploidy')


## final set
ggplot(asize[!asize$V2%in%lowQualAssemblies,], aes(x=ploidy, y=scaledTransloc, group=ploidy, color=ploidy)) + geom_boxplot(outlier.shape=NA) +geom_point(position = position_jitter(seed = 1),  size=3)+ scale_color_manual(values=ploidycolors) + theme(legend.position='none')+ ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=6) + ylab('Translocations per haploid equivalent') + xlab('Ploidy')
## with bad assemblies
ggplot(asize, aes(x=ploidy, y=scaledTransloc, group=ploidy, color=ploidy)) + geom_boxplot(outlier.shape=NA) +geom_point(position = position_jitter(seed = 1),  size=3)+ scale_color_manual(values=ploidycolors) + theme(legend.position='none')+ ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=6) + ylab('Translocations per haploid equivalent') + xlab('Ploidy')


ggplot(asize[!asize$V2%in%lowQualAssemblies,], aes(x=ploidy, y=scaledTransloc, group=ploidy, color=ploidy)) + geom_boxplot(outlier.shape=NA) +geom_point(position = position_jitter(seed = 1),  size=2, aes(shape=haploid))+ geom_text(aes(label=V2), position = position_jitter(seed = 1)) + scale_color_manual(values=ploidycolors) + theme(legend.position='none')+ ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=6) + ylab('Translocations per haploid equivalent') + xlab('Ploidy')
ggplot(asize, aes(x=ploidy, y=scaledTransloc, group=ploidy, color=ploidy)) + geom_boxplot(outlier.shape=NA) +geom_point(position = position_jitter(seed = 1),  size=2, aes(shape=haploid))+ geom_text(aes(label=V2), position = position_jitter(seed = 1)) + scale_color_manual(values=ploidycolors) + theme(legend.position='none')+ ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=6) + ylab('Translocations per haploid equivalent') + xlab('Ploidy')


## no clear relatioonship with parental age?? fair if it's all cytotype age...
ggplot(asize[!asize$V2%in%lowQualAssemblies,], aes(x=mya, y=scaledTransloc, color=ploidy)) + geom_point() + scale_color_manual(values=ploidycolors)

## translocations are related to chromosome count, reassuringly!!!
chrCount=read.table('../general_summaries/panand_chr_counts.txt', header=T)
asize$chrCount=chrCount$haploidchr[match(asize$V2, chrCount$sixlettercode)]
ggplot(asize[!asize$V2%in%lowQualAssemblies,], aes(x=chrCount, y=scaledTransloc, color=ploidy)) + geom_jitter() + scale_color_manual(values=ploidycolors) + xlab('Haploid Chromosome Number') + ylab('Translocations per haploid equivalent')
ggplot(asize, aes(x=chrCount, y=scaledTransloc, color=ploidy)) + geom_jitter() + scale_color_manual(values=ploidycolors)+ xlab('Haploid Chromosome Number') + ylab('Translocations per haploid equivalent')
ggplot(asize, aes(x=chrCount, y=scaledTransloc, color=ploidy)) + geom_jitter() + geom_text(aes(label=V2)) + scale_color_manual(values=ploidycolors)+ xlab('Haploid Chromosome Number') + ylab('Translocations per haploid equivalent')


ggplot(asize, aes(x=chrCount, y=transloc, color=ploidy)) + geom_jitter() + scale_color_manual(values=ploidycolors)+ xlab('Haploid Chromosome Number') + ylab('Translocations')


### are rearrangements happening at the same position relative to paspalum?
countRearrangementsPasp <- function(filepath, query='queryName', color_palette=muted_colors, minBlock=20) {
  # Load data
  data <- read.table(filepath, header = TRUE)
  data <- data[data$gene != 'interanchor', ]
  
  # Reduce to blocks and calculate stats
  data <- data %>%
    group_by(blockIndex) %>%
    mutate(blockLength = n(), blockEnd=max(referenceEnd)) %>%
    group_by(queryChr) %>%
    mutate(freqStrand = names(which.max(table(strand))),
           maxChr = max(queryStart),
           freqRef = names(which.max(table(refChr))))
  
  # Filter data based on block length
  data <- data[data$blockLength > minBlock, ]
  data$refChr <- factor(data$refChr, levels = c(paste0('Chr0', 1:9), 'Chr10'))
  
  temp=data %>% filter(blockLength>minBlock) %>% group_by(refChr, queryChr, blockIndex, blockLength, blockEndWindow=round(blockEnd,-6)) %>% dplyr::summarize(ick=n())
  temp$query=query
   return(temp[,c('refChr', 'blockLength', 'blockEndWindow', 'query')])
}

for(i in c('zmB735')){
  print(countRearrangementsPasp(Sys.glob(paste0('../syntenic_anchors/anchors/', i, '-Pv-*')), minBlock=50, query=i))
}  

## i think this is interesting but not helpful, since these aren't just the regions that 
paspBreaks=do.call(rbind, lapply(asize$V2, function(x) countRearrangementsPasp(Sys.glob(paste0('../syntenic_anchors/anchors/', x, '-Pv-*')), minBlock=30, query=x)))
paspBreaks$ploidy=asize$ploidy[match(paspBreaks$query, asize$V2)]
ggplot(paspBreaks, aes(x=blockEndWindow, group=ploidy, fill=ploidy, color=ploidy)) + geom_histogram(binwidth=1e6) + facet_wrap(~refChr, nrow=1) + scale_color_manual(values=ploidycolors) + scale_fill_manual(values=ploidycolors)

