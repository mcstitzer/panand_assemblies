library(data.table)
library(ggrepel)
library(dplyr)
library(rtracklayer)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

all=read.table('../panand_sp_ploidy.txt', header=F)

asize=fread('../general_summaries/panand_assembly_sizes.txt', header=T, quote="", fill=T)
asize$ploidy=factor(asize$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))

## count rearrangments
diploids=c('cserru', 'irugos', 'sbicol', 'ppanic', 'ttrian', 'crefra', 'avirgi', 'smicro', 'rtuber')
tetraploids=c('snutan', 'hconto', 'ccitra', 'achine', 'sscopa', 'etrips',  'vcuspi')
hexaploids=c('udigit', 'agerar', 'hcompr', 'blagur')
paleotetraploids=c("tdacn1", "tdacs1", "zdgigi", "zdmomo", "zluxur", "zmhuet", "zTIL18", "zTIL25", "zTIL01", "zTIL11", "znicar", "zmB735")

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
  # Step 4: Identify fusion points
  fusion_data <- data %>%
    group_by(queryChr) %>%
    arrange(queryChr, queryStart) %>%
    mutate(
      ref_transition = refChr != lag(refChr, default = refChr[1]), # Detect transitions in refChr
      fusion_point = ifelse(ref_transition, queryStart, NA)
    ) %>%
    filter(ref_transition) %>%
    select(queryChr, fusion_point, refChr) %>%
    tidyr::drop_na(fusion_point)
  
  # Step 5: Summarize counts of unique reference chromosomes
  rearrangement_counts <- fusion_data %>%
    group_by(queryChr) %>%
    summarize(
      n_fusions = n(), # Count the number of fusion points
      unique_refChrs = length(unique(refChr)), # Count unique reference chromosomes
      .groups = "drop"
    )
  
  # Combine results
  result <- list(
    fusion_data = fusion_data,
    rearrangement_summary = rearrangement_counts
  )
  
  return(sum(rearrangement_counts$n_fusions))
  
#   temp=data %>% filter(blockLength>minBlock) %>% group_by(queryChr)%>% summarize(n=length(unique(refChr)))
#   sum(temp$n>1)
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

lowQualAssemblies=c('telega', 'atenui', 'rrottb', 'ccitra')

asize$transloc=sapply(asize$V2, function(x) countRearrangements(Sys.glob(paste0('../syntenic_anchors/anchors/', x, '-Pv-*')), minBlock=30))
asize$transloc=ifelse(asize$haploid, asize$transloc, asize$transloc/2) ## if haploid, this is true number, if allelic, don't count each allele
asize %>% group_by(ploidy) %>% summarize(translocMean=mean(transloc, na.rm=T), scaledTranslocMean=mean(scaledTransloc, na.rm=T))
kruskal.test(transloc~ploidy, data=asize) ## pvalue=0.0001115 --> sig diff between group
pairwise.wilcox.test(asize$transloc, asize$ploidy,
                     p.adjust.method = "BH")
kruskal.test(transloc/chrCount~ploidy, data=asize) ## pvalue=0.0001932 --> sig diff between group
pairwise.wilcox.test(asize$transloc/asize$chrCount, asize$ploidy,
                     p.adjust.method = "BH")


ggplot(asize[!asize$V2%in%lowQualAssemblies], aes(x=ploidy, y=transloc, group=ploidy, color=ploidy)) + geom_boxplot(outlier.shape=NA) +geom_point(position = position_jitter(seed = 1),  size=2, aes(shape=haploid))+ scale_color_manual(values=ploidycolors) + theme(legend.position='none')+ ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=10)

asize$scaledTransloc=ifelse(asize$ploidy%in%c('Paleotetraploid', 'Tetraploid'), asize$transloc/2, ifelse(asize$ploidy=='Hexaploid', asize$transloc/3, asize$transloc))
#### haploid or not scaled by triangle/circle
ggplot(asize[!asize$V2%in%lowQualAssemblies,], aes(x=ploidy, y=scaledTransloc, group=ploidy, color=ploidy)) + geom_boxplot(outlier.shape=NA) +geom_point(position = position_jitter(seed = 1),  size=2, aes(shape=haploid))+ scale_color_manual(values=ploidycolors) + theme(legend.position='none')+ ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=6) + ylab('Translocations per haploid equivalent') + xlab('Ploidy')


## final set
ggplot(asize[!asize$V2%in%lowQualAssemblies,], aes(x=ploidy, y=scaledTransloc, group=ploidy, color=ploidy)) + geom_boxplot(outlier.shape=NA) +geom_point(position = position_jitter(seed = 1),  size=3)+ scale_color_manual(values=ploidycolors) + theme(legend.position='none')+ ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=6) + ylab('Rearrangements per haploid equivalent') + xlab('Ploidy')
## with bad assemblies
ggplot(asize, aes(x=ploidy, y=scaledTransloc, group=ploidy, color=ploidy)) + geom_boxplot(outlier.shape=NA) +geom_point(position = position_jitter(seed = 1),  size=3)+ scale_color_manual(values=ploidycolors) + theme(legend.position='none')+ ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=6) + ylab('Rearrangments per haploid equivalent') + xlab('Ploidy')

## and unscaled
## final set
ggplot(asize[!asize$V2%in%lowQualAssemblies,], aes(x=ploidy, y=transloc, group=ploidy, color=ploidy)) + geom_boxplot(outlier.shape=NA) +geom_point(position = position_jitter(seed = 1),  size=3)+ scale_color_manual(values=ploidycolors) + theme(legend.position='none')+ ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=9.5) + ylab('Rearrangements') + xlab('Ploidy')
## with bad assemblies
ggplot(asize, aes(x=ploidy, y=transloc, group=ploidy, color=ploidy)) + geom_boxplot(outlier.shape=NA) +geom_point(position = position_jitter(seed = 1),  size=3)+ scale_color_manual(values=ploidycolors) + theme(legend.position='none')+ ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=9.5) + ylab('Rearrangements') + xlab('Ploidy')


## numbers for paper
asize %>% filter(!V2 %in% lowQualAssemblies) %>% group_by(ploidy) %>% summarize(meanrearr=mean(transloc, na.rm=T), meanscaled=mean(scaledTransloc, na.rm=T), meanchrscaled=mean(transloc/chrCount, na.rm=T))
## sloppy way to get tripsacum and zea
asize %>% filter(!V2 %in% lowQualAssemblies) %>% group_by(V2%in%c('tdacs1', 'tdacn1'), substr(V2,1,1)=='z') %>% summarize(meanrearr=mean(transloc, na.rm=T), meanscaled=mean(scaledTransloc, na.rm=T), meanchrscaled=mean(transloc/chrCount, na.rm=T))


ggplot(asize[!asize$V2%in%lowQualAssemblies,], aes(x=ploidy, y=scaledTransloc, group=ploidy, color=ploidy)) + geom_boxplot(outlier.shape=NA) +geom_point(position = position_jitter(seed = 1),  size=2, aes(shape=haploid))+ geom_text(aes(label=V2), position = position_jitter(seed = 1)) + scale_color_manual(values=ploidycolors) + theme(legend.position='none')+ ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=6) + ylab('Translocations per haploid equivalent') + xlab('Ploidy')
ggplot(asize, aes(x=ploidy, y=scaledTransloc, group=ploidy, color=ploidy)) + geom_boxplot(outlier.shape=NA) +geom_point(position = position_jitter(seed = 1),  size=2, aes(shape=haploid))+ geom_text(aes(label=V2), position = position_jitter(seed = 1)) + scale_color_manual(values=ploidycolors) + theme(legend.position='none')+ ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=6) + ylab('Translocations per haploid equivalent') + xlab('Ploidy')


## no clear relatioonship with parental age?? fair if it's all cytotype age...
asize$mya=het$mya[match(asize$V2, het$genome)]
ggplot(asize[!asize$V2%in%lowQualAssemblies,], aes(x=mya, y=scaledTransloc, color=ploidy)) + geom_point() + scale_color_manual(values=ploidycolors)

## translocations are related to chromosome count, reassuringly!!!
chrCount=read.table('../general_summaries/panand_chr_counts.txt', header=T)
asize$chrCount=chrCount$haploidchr[match(asize$V2, chrCount$sixlettercode)]
ggplot(asize[!asize$V2%in%lowQualAssemblies,], aes(x=chrCount, y=scaledTransloc, color=ploidy)) + geom_jitter() + scale_color_manual(values=ploidycolors) + xlab('Haploid Chromosome Number') + ylab('Translocations per haploid equivalent')
ggplot(asize, aes(x=chrCount, y=scaledTransloc, color=ploidy)) + geom_jitter() + scale_color_manual(values=ploidycolors)+ xlab('Haploid Chromosome Number') + ylab('Translocations per haploid equivalent')
ggplot(asize, aes(x=chrCount, y=scaledTransloc, color=ploidy)) + geom_jitter() + geom_text(aes(label=V2)) + scale_color_manual(values=ploidycolors)+ xlab('Haploid Chromosome Number') + ylab('Translocations per haploid equivalent')


ggplot(asize, aes(x=chrCount, y=scaledTransloc, color=ploidy)) + geom_vline(xintercept=c(9,11:19,21:29), color='gray90', lty='dotted', alpha=0.3) + geom_vline(xintercept=c(10,20,30), color='gray80', lty='dashed', alpha=0.3)+ geom_jitter(size=4, width = 0.1, alpha=1) + scale_color_manual(values=ploidycolors)+ xlab('Haploid Chromosome Number') + ylab('Rearrangements per Diploid Equivalent') 

summary(lm(data=asize, scaledTransloc~chrCount))

## get medians for zea and tripsacum genera separately
asizezt=asize
asizezt$V2zt=ifelse(asize$V2 %in% c('tdacn1', 'tdacs1'), 'trips', asize$V2)
asizezt$V2zt=ifelse(asize$V2 %in% c("zdgigi", "zdmomo", "zluxur", "zmhuet", "zTIL18", "zTIL25", "zTIL01", "zTIL11", "znicar", "zmB735"), 'zea', asize$V2)
asizezt=asizezt%>% group_by(V2zt, ploidy) %>% summarize(chrCount=median(chrCount), scaledTransloc=median(scaledTransloc), mya=median(mya))


## with regression line
fig_chrrearr=ggplot(asize, aes(x=chrCount, y=scaledTransloc, color=ploidy)) + geom_vline(xintercept=c(9,11:19,21:29), color='gray90', lty='dotted', alpha=0.3) + geom_vline(xintercept=c(10,20,30), color='gray80', lty='dashed', alpha=0.3)+ 
  geom_point(size=4, position=position_jitter(seed=9,width = 0.1), alpha=0.8)+ #, pch=1, stroke=3) + 
  scale_color_manual(values=ploidycolors)+ 
  xlab('Haploid Chromosome Number') + ylab('Rearrangements per\nDiploid Equivalent') + 
  stat_smooth(data=asizezt, method='lm', aes(group=NA), se=F, color='gray80')
fig_chrrearr 

cor.test(asizezt$scaledTransloc, asizezt$chrCount)

### compare rearrangemnets to age
fig_agerearr=ggplot(asize, aes(x=mya, y=scaledTransloc, color=ploidy)) + geom_vline(xintercept=c(2,4,6,8,10,12,14), color='gray90', lty='dotted', alpha=0.3) + geom_vline(xintercept=c(1,3,5,7,9,11,13), color='gray80', lty='dashed', alpha=0.3)+ 
  geom_point(size=4, position=position_jitter(seed=9,width = 0.1), alpha=0.8)+#, pch=1, stroke=3) + 
  scale_color_manual(values=ploidycolors)+ xlab('Divergence between\nParental Subgenomes (Mya)') + 
#  stat_smooth(data=asizezt, method='lm', aes(group=NA), alpha=0.1, se=F, color='gray90')+
  ylab('Rearrangements per\nDiploid Equivalent') 
fig_agerearr

cor.test(asizezt$scaledTransloc, asizezt$mya)


##### laxy and putting this in the other document, count_anchors_and_plot.R


## add text label for toby
ggplot(asize, aes(x=chrCount, y=scaledTransloc, color=ploidy)) + geom_vline(xintercept=c(9,11:19,21:29), color='gray90', lty='dotted', alpha=0.3) + geom_vline(xintercept=c(10,20,30), color='gray80', lty='dashed', alpha=0.3)+ geom_point(size=4, position=position_jitter(seed=9,width = 0.1), alpha=0.8, pch=1, stroke=3) + scale_color_manual(values=ploidycolors)+ xlab('Haploid Chromosome Number') + ylab('Rearrangements per Diploid Equivalent') + stat_smooth(method='lm', aes(group=NA), se=F, color='gray80') + geom_text_repel(aes(label=V2))


## not scaled to diploid
ggplot(asize, aes(x=chrCount, y=transloc, color=ploidy)) + geom_vline(xintercept=c(9,11:19,21:29), color='gray90', lty='dotted', alpha=0.3) + geom_vline(xintercept=c(10,20,30), color='gray80', lty='dashed', alpha=0.3)+ geom_jitter(size=4, width = 0.1, alpha=1) + scale_color_manual(values=ploidycolors)+ xlab('Haploid Chromosome Number') + ylab('Rearrangements per Diploid Equivalent') 


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






## toby's suggestion - count proportion of assembly covered by syntenic blocks
syntenicCoverage <- function(filepath, color_palette=muted_colors, minBlock=20) {
  # Load data
  data <- read.table(filepath, header = TRUE)
  data <- data[data$gene != 'interanchor', ]
  
  # Reduce to blocks and calculate stats
  data <- data %>%
    group_by(blockIndex) %>%
    dplyr::summarize(blockLength =max(queryEnd)-min(queryStart), nAnchors=n()) 
  
  # Filter data based on block length
  data <- data[data$nAnchors > minBlock, ]
#  data$refChr <- factor(data$refChr, levels = c(paste0('Chr0', 1:9), 'Chr10'))
  
  return(sum(data$blockLength))
  
}

refCoverage <- function(filepath, color_palette=muted_colors, minBlock=20) {
  # Load data
  data <- read.table(filepath, header = TRUE)
  data <- data[data$gene != 'interanchor', ]
  
  # Reduce to blocks and calculate stats
  data <- data %>%
    group_by(blockIndex, refChr) %>%
    dplyr::summarize(referenceEnd=max(referenceEnd), referenceStart=min(referenceStart), nAnchors=n()) 
  
  # Filter data based on block length
  data <- data[data$nAnchors > minBlock, ]
  #  data$refChr <- factor(data$refChr, levels = c(paste0('Chr0', 1:9), 'Chr10'))
#  data=data[data$refChr%in%c(paste0('Chr0', 1:9), 'Chr10'),]
  gr=reduce(GRanges(seqnames=data$refChr, IRanges(start=data$referenceStart, end=data$referenceEnd)))
  
  return(sum(width(gr)))
  
}


refCoverageGenes <- function(filepath, color_palette=muted_colors, minBlock=20) {
  # Load data
  data <- read.table(filepath, header = TRUE)
  data <- data[data$gene != 'interanchor', ]
  
  # Reduce to blocks and calculate stats
  data <- data %>%
    group_by(gene, refChr) %>%
    dplyr::summarize(geneCount=n()) 
  
  # Filter data based on block length
#  data <- data[data$nAnchors > minBlock, ]
  #  data$refChr <- factor(data$refChr, levels = c(paste0('Chr0', 1:9), 'Chr10'))
  #  data=data[data$refChr%in%c(paste0('Chr0', 1:9), 'Chr10'),]

  return(nrow(data))
  
}



asize$bpCovered30=sapply(asize$V2, function(x) syntenicCoverage(Sys.glob(paste0('../syntenic_anchors/anchors/', x, '-Pv-*')), minBlock=30))
asize$bpCovered=sapply(asize$V2, function(x) syntenicCoverage(Sys.glob(paste0('../syntenic_anchors/anchors/', x, '-Pv-*')), minBlock=2))

asize$refCovered=sapply(asize$V2, function(x) refCoverage(Sys.glob(paste0('../syntenic_anchors/anchors/', x, '-Pv-*')), minBlock=2))

asize$refGenes=sapply(asize$V2, function(x) refCoverageGenes(Sys.glob(paste0('../syntenic_anchors/anchors/', x, '-Pv-*'))))


paspsize=651047655
summary(asize$refCovered/paspsize)

pa=read.table('../syntenic_anchors/anchors/pvagin-Pv-2', header=T)
pa=pa[pa$gene!='interanchor',]
paspanchors=nrow(pa)

summary(asize$refGenes/paspanchors)


for(i in diploids){
  print(syntenicCoverage(Sys.glob(paste0('../syntenic_anchors/anchors/', i, '-Pv-*')), minBlock = 3))
}



#tes=import.gff3(paste0('trash/', all$V1[all$V2==genotype], '_EDTAandTandemRepeat.gff3'))

#### also now count knobs
countKnobsCent <- function(filepath, tefilepath, genome='', distance=10000,color_palette=muted_colors, minBlock=20) {
  # Load data
  anchors <- read.table(filepath, header = TRUE)
  anchors <- anchors[anchors$gene != 'interanchor', ]
  
  # Reduce to blocks and calculate stats
  anchors <- anchors %>%
    group_by(blockIndex) %>%
    mutate(blockLength = n()) %>%
    group_by(queryChr) %>%
    mutate(freqStrand = names(which.max(table(strand))),
           maxChr = max(queryStart),
           freqRef = names(which.max(table(refChr))))
  
  # Filter data based on block length
  anchors <- anchors[anchors$blockLength > minBlock, ]
 # Identify fusions and calculate query positions
  fusion_data <- anchors %>%
    filter(blockLength > minBlock) %>%
    group_by(queryChr) %>%
    arrange(queryChr, queryStart) %>% # Ensure sorted order by queryStart
    mutate(
      ref_transition = refChr != lag(refChr, default = refChr[1]), # Detect transitions in refChr
      fusion_point = ifelse(ref_transition, queryStart, NA) # Record query position of the transition
    ) %>%
    filter(ref_transition) %>% # Retain only rows with transitions
    select(queryChr, fusion_point, refChr)%>% tidyr::drop_na(fusion_point)
#     summarize(
#       queryChr = dplyr::first(queryChr),
#       fusion_positions = paste(na.omit(fusion_point), collapse = ", "),
#       involved_refChrs = paste(unique(refChr), collapse = ", ")
#     ) %>% tidyr::unnest(fusion_positions)
#   
  # Print fusion information
  print(fusion_data)

 tes=import.gff(tefilepath)
 knob180=tes[grepl('179bp_repeat|180bp_repeat', tes$Target),]
 if(length(knob180)>0){knob180$repeattype='knob180'}
 knobtr1=tes[grepl('357bp_repeat|358bp_repeat', tes$Target),]
 if(length(knobtr1)>0){knobtr1$repeattype='knobtr1'}
 cent=tes[grepl('154bp_repeat|155bp_repeat|156bp_repeat', tes$Target),]
 if(length(cent)>0){cent$repeattype='cent'}
 tands=c(knob180,knobtr1,cent)
 print(table(tands$repeattype))
 ## only get from chromsomes, as there are those giant contigs of knob that we don't have in the assembly chr to see if they're syntenic...
 print(sum(width(reduce(knob180[seqnames(knob180)%in%paste0('chr',1:10),], min.gapwidth=10))))
# print(sum(width(knob180)))

if(genome=='zmB735'){
print(sum(width(tes[grepl('knob180', tes$Classification),])))
tands=tes[tes$type%in%c('centromeric_repeat', 'knob'),]
tands$repeattype=tands$type
print(head(seqlevels(tands)))
print(class(seqlevels(tands)))
tands_df <- as.data.frame(tands)
# Step 2: Remove the "B73_" prefix from seqnames
tands_df$seqnames <- gsub("^B73_", "", tands_df$seqnames)
# Step 3: Convert the data frame back to GRanges
tands <- makeGRangesFromDataFrame(tands_df, keep.extra.columns = TRUE)


}

breakpoints_gr=GRanges(seqnames=fusion_data$queryChr, ranges=IRanges(start=fusion_data$fusion_point, end=fusion_data$fusion_point))
if(length(tands)>0){
overlaps=findOverlaps(breakpoints_gr, tands, maxgap=distance, ignore.strand=T)

  # Summarize results
  overlap_results <- as.data.frame(overlaps) %>%
    mutate(
      queryChr = as.character(seqnames(breakpoints_gr))[queryHits],
      refChr = fusion_data$refChr[queryHits],
      fusion_point = start(breakpoints_gr)[queryHits],
      repeat_info = mcols(tands)$repeattype[subjectHits]
    ) %>%
    group_by(queryChr, fusion_point, repeat_info, refChr) %>%
    summarize(count=n())
  
  print(overlap_results)
     # Add counts to fusion_data
    fusion_data <- fusion_data %>%
      left_join(overlap_results %>%
                  group_by(queryChr, fusion_point) %>%
                  summarize(total_repeats = sum(count), .groups = "drop"),
                by = c("queryChr", "fusion_point"))
    
    # Replace NA counts with 0
    fusion_data$total_repeats[is.na(fusion_data$total_repeats)] <- 0
  } else {
    fusion_data$total_repeats <- 0
  }
if(genome!=''){fusion_data$genome=genome}
return(fusion_data)



}

## /Users/mcs368/Downloads/repeatmask_tandems/Td-FL_9056069_6-REFERENCE-PanAnd-2.0a_EDTAandTandemRepeat.gff3

for(i in paleotetraploids){
print(i)
countKnobsCent(paste0('../syntenic_anchors/anchors/',i, '-Pv-', all$V3[all$V2==i]*2), paste0('~/Downloads/repeatmask_tandems/',all$V1[all$V2==i], '_EDTAandTandemRepeat.gff3'), genome=i, distance=1e6, minBlock=30)
}

fused_list=lapply(paleotetraploids, function(i){
countKnobsCent(paste0('../syntenic_anchors/anchors/',i, '-Pv-', all$V3[all$V2==i]*2), paste0('~/Downloads/repeatmask_tandems/',all$V1[all$V2==i], '_EDTAandTandemRepeat.gff3'), genome=i, distance=1e6, minBlock=30)
})

fused=do.call(rbind, fused_list)


## get unique copies
fused[fused$total_repeats>0,] %>% group_by(queryChr, refChr) %>% summarize(mean=mean(total_repeats))