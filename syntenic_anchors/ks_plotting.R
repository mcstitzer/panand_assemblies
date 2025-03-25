library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(stringr)




synt=read.table('~/Downloads/sharedSyntenicAnchors.txt', header=F)


b=read.table('../syntenic_anchors/gene_anchorTable_AnchorWave_Paspalum.txt', header=T)
b=b[,-which(colnames(b)%in%c('svirid',"pprate", "osativ", "eophiu", "bdista", "agerjg", 'tdactm', 'tzopol'))]
b=b[,-which(colnames(b)%in%c('tdacs2', 'tdacn2'))]



table(rowSums(b[,-1]>0))
bb=b[rowSums(b[,-1]>0)>=32,]

sum(synt$V1==str_split_fixed(bb$gene, '\\.',2)[,1])


## now get positions in the genome for thsese
all=read.table('../panand_sp_ploidy.txt')
all$V3[all$V2=='rtuber']=2
all$V3[all$V2=='telega']=2

anchors=lapply(all$V2, function(x) {
  a=read.table(paste0('anchors/',x, '-Pv-', 2*all$V3[all$V2==x]), header=T)
  a$genome=x
  return(a)
})

ab=Reduce(function(...) merge(..., all=T), anchors)

b=ab[ab$gene!='interanchor',] ## keep only genes
###
### okay so b is now positions in each queyr genome
bs=b[b$gene%in%bb$gene,]


asize=fread('../general_summaries/panand_assembly_sizes.txt', header=T, quote="", fill=T)
asize$syntAnchors=colSums(bb[,-1])[match(asize$V2, names(colSums(bb[,-1])))]
asize$syntAnchorsCount=colSums(bb[,-1]>0)[match(asize$V2, names(colSums(bb[,-1]>0)))]
asize$syntAnchors=ifelse(asize$haploid, asize$syntAnchors, asize$syntAnchors/2)
asize$ploidy=factor(asize$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))
asize$diploidEquivalent=ifelse(asize$ploidy=='Diploid', asize$syntAnchors, ifelse(asize$ploidy%in%c('Tetraploid', 'Paleotetraploid'), asize$syntAnchors/2, ifelse(asize$ploidy=='Hexaploid', asize$syntAnchors/3, NA)))


#### whoops now making separate files for self-to-self of each!!!


ks=fread('ks_self/blagur_matches.txt')
ks$gene=str_split_fixed(ks$V3, '_',3)[,2]
kss=ks[ks$gene%in%synt$V1,]
kss$queryChr=str_split_fixed(kss$V3, '_', 3)[,3]
kss$queryChr[substr(kss$queryChr,1,1)%in%c('a', 's')]=str_extract(kss$V3, "(scaf_[0-9]+|alt-scaf_[0-9]+)")[substr(kss$queryChr,1,1)%in%c('a', 's')]
kss$queryChr[substr(kss$queryChr,1,3)%in%c('ctg')]=str_extract(kss$V3, "(ctg_[0-9]+)")[substr(kss$queryChr,1,3)%in%c('ctg')]
kss$queryStart <- as.integer(str_match(kss$V3, "_(\\d+)-(\\d+)$")[, 2])
kss$queryEnd <- as.integer(str_match(kss$V3, "_(\\d+)-(\\d+)$")[, 3])
#outFiltered$dnds=ks$V19[match(paste(outFiltered$genome, outFiltered$first_gene), paste(ks$genome, ks$gene))]
kss$refGenome=substr(kss$V4,1,6)
kss$refChr=str_split_fixed(kss$V4, '_', 3)[,3]
kss$refChr[substr(kss$refChr,1,1)%in%c('a', 's')]=str_extract(kss$V4, "(scaf_[0-9]+|alt-scaf_[0-9]+)")[substr(kss$queryChr,1,1)%in%c('a', 's')]
kss$refChr[substr(kss$refChr,1,3)%in%c('ctg')]=str_extract(kss$V4, "(ctg_[0-9]+)")[substr(kss$queryChr,1,3)%in%c('ctg')]
kss$refStart <- as.integer(str_match(kss$V4, "_(\\d+)-(\\d+)$")[, 2])
kss$refEnd <- as.integer(str_match(kss$V4, "_(\\d+)-(\\d+)$")[, 3])

kss$paspChr=substr(kss$gene,6,7)


#### now do for each 
merged_results <- list()

# Loop through each file, process it, and append it to the merged_results list
for (file in all$V2) {
  # Read in the ks file
  ks <- fread(paste0('../syntenic_anchors/ks_self/', file, '_matches.txt'))
  
  # Extract gene information
  ks$gene <- str_split_fixed(ks$V3, '_', 3)[, 2]
  
  # Filter based on 'synt' gene list
  kss <- ks[ks$gene %in% synt$V1, ]
  
  # Process queryChr
  kss$queryChr <- str_split_fixed(kss$V3, '_', 4)[, 3]
  kss$queryChr[substr(kss$queryChr, 1, 1) %in% c('a', 's')] <- str_extract(kss$V3, "(scaf_[0-9]+|alt-scaf_[0-9]+)")[substr(kss$queryChr, 1, 1) %in% c('a', 's')]
  kss$queryChr[substr(kss$queryChr, 1, 3) %in% c('ctg')] <- str_extract(kss$V3, "(ctg_[0-9]+)")[substr(kss$queryChr, 1, 3) %in% c('ctg')]
  
  # Extract queryStart and queryEnd
  kss$queryStart <- as.integer(str_match(kss$V3, "_(\\d+)-(\\d+)$")[, 2])
  kss$queryEnd <- as.integer(str_match(kss$V3, "_(\\d+)-(\\d+)$")[, 3])
  
  # Process reference information
  kss$refGenome <- substr(kss$V4, 1, 6)
  kss$refChr <- str_split_fixed(kss$V4, '_', 4)[, 3]
  kss$refChr[substr(kss$refChr, 1, 1) %in% c('a', 's')] <- str_extract(kss$V4, "(scaf_[0-9]+|alt-scaf_[0-9]+)")[substr(kss$refChr, 1, 1) %in% c('a', 's')]
  kss$refChr[substr(kss$refChr, 1, 3) %in% c('ctg')] <- str_extract(kss$V4, "(ctg_[0-9]+)")[substr(kss$refChr, 1, 3) %in% c('ctg')]
  
  # Extract refStart and refEnd
  kss$refStart <- as.integer(str_match(kss$V4, "_(\\d+)-(\\d+)$")[, 2])
  kss$refEnd <- as.integer(str_match(kss$V4, "_(\\d+)-(\\d+)$")[, 3])
  
  # Add paspChr column
  kss$paspChr <- substr(kss$gene, 6, 7)
  
  # Add a column for the genome from the file name (assuming 'blagur' is part of the file name)
  kss$genome <- file
  
  # Append the result to the merged_results list
  merged_results[[file]] <- kss
}

# Combine all processed data frames into one
final_merged_data <- rbindlist(merged_results)
final_merged_data$ploidy=asize$ploidy[match(final_merged_data$genome, asize$V2)]
final_merged_data$speciesLabel=asize$speciesLabel[match(final_merged_data$genome, asize$V2)]

final_merged_data %>% group_by(queryChr, refChr, paspChr, genome) %>%summarize(ks=median(V17), n=n()) %>% filter(n>10) %>% ggplot() + geom_histogram(aes(x=ks), binwidth=0.001) + facet_wrap(~paspChr, ncol=1) 


pdf('~/Downloads/ks_self.pdf',10,20)
final_merged_data %>% group_by(queryChr, refChr, paspChr, genome,ploidy, speciesLabel) %>%summarize(ks=median(V17), n=length(unique(gene))) %>% filter(n>10) %>% ggplot()+ geom_vline(xintercept=seq(from=0,to=0.14,by=0.01), color='gray', lty='dashed', alpha=0.2) + geom_histogram(aes(x=ks, color=ploidy, fill=ploidy), binwidth=0.0001) + facet_wrap(~speciesLabel, ncol=1, scales='free_y') + scale_color_manual(values=ploidycolors)+ scale_fill_manual(values=ploidycolors) 
final_merged_data %>% group_by(queryChr, refChr, paspChr, genome,ploidy, speciesLabel) %>%summarize(ks=median(V17), n=length(unique(gene))) %>% filter(n>10) %>% ggplot() + geom_vline(xintercept=seq(from=0,to=0.14,by=0.01), color='gray', lty='dashed', alpha=0.2)+ geom_histogram(aes(x=ks, color=ploidy, fill=ploidy), binwidth=0.001) + facet_wrap(~speciesLabel, ncol=1, scales='free_y') + scale_color_manual(values=ploidycolors)+ scale_fill_manual(values=ploidycolors)

final_merged_data %>% group_by(queryChr, refChr, paspChr, genome,ploidy, speciesLabel) %>%summarize(ks=median(V17), n=length(unique(gene))) %>% filter(n>10) %>% ggplot()+ geom_vline(xintercept=seq(from=0,to=0.14,by=0.01), color='gray', lty='dashed', alpha=0.2) + geom_histogram(aes(x=ks, color=ploidy, fill=ploidy, weight=n), binwidth=0.0001) + facet_wrap(~speciesLabel, ncol=1, scales='free_y') + scale_color_manual(values=ploidycolors)+ scale_fill_manual(values=ploidycolors)
final_merged_data %>% group_by(queryChr, refChr, paspChr, genome,ploidy, speciesLabel) %>%summarize(ks=median(V17), n=length(unique(gene))) %>% filter(n>10) %>% ggplot()+ geom_vline(xintercept=seq(from=0,to=0.14,by=0.01), color='gray', lty='dashed', alpha=0.2) + geom_histogram(aes(x=ks, color=ploidy, fill=ploidy, weight=n), binwidth=0.001) + facet_wrap(~speciesLabel, ncol=1, scales='free_y') + scale_color_manual(values=ploidycolors)+ scale_fill_manual(values=ploidycolors)


final_merged_data %>% group_by(queryChr, refChr, paspChr, genome,ploidy, speciesLabel) %>%summarize(ks=median(V17), n=length(unique(gene))) %>% filter(n>100) %>% ggplot()+ geom_vline(xintercept=seq(from=0,to=0.14,by=0.01), color='gray', lty='dashed', alpha=0.2) + geom_histogram(aes(x=ks, color=ploidy, fill=ploidy), binwidth=0.0001) + facet_wrap(~speciesLabel, ncol=1, scales='free_y') + scale_color_manual(values=ploidycolors)+ scale_fill_manual(values=ploidycolors) 
final_merged_data %>% group_by(queryChr, refChr, paspChr, genome,ploidy, speciesLabel) %>%summarize(ks=median(V17), n=length(unique(gene))) %>% filter(n>100) %>% ggplot() + geom_vline(xintercept=seq(from=0,to=0.14,by=0.01), color='gray', lty='dashed', alpha=0.2)+ geom_histogram(aes(x=ks, color=ploidy, fill=ploidy), binwidth=0.001) + facet_wrap(~speciesLabel, ncol=1, scales='free_y') + scale_color_manual(values=ploidycolors)+ scale_fill_manual(values=ploidycolors)

final_merged_data %>% group_by(queryChr, refChr, paspChr, genome,ploidy, speciesLabel) %>%summarize(ks=median(V17), n=length(unique(gene))) %>% filter(n>100) %>% ggplot()+ geom_vline(xintercept=seq(from=0,to=0.14,by=0.01), color='gray', lty='dashed', alpha=0.2) + geom_histogram(aes(x=ks, color=ploidy, fill=ploidy, weight=n), binwidth=0.0001) + facet_wrap(~speciesLabel, ncol=1, scales='free_y') + scale_color_manual(values=ploidycolors)+ scale_fill_manual(values=ploidycolors)
final_merged_data %>% group_by(queryChr, refChr, paspChr, genome,ploidy, speciesLabel) %>%summarize(ks=median(V17), n=length(unique(gene))) %>% filter(n>100) %>% ggplot()+ geom_vline(xintercept=seq(from=0,to=0.14,by=0.01), color='gray', lty='dashed', alpha=0.2) + geom_histogram(aes(x=ks, color=ploidy, fill=ploidy, weight=n), binwidth=0.001) + facet_wrap(~speciesLabel, ncol=1, scales='free_y') + scale_color_manual(values=ploidycolors)+ scale_fill_manual(values=ploidycolors)

dev.off()


pdf('~/Downloads/ks_self.mean.pdf',10,20)
final_merged_data %>% group_by(queryChr, refChr, paspChr, genome,ploidy, speciesLabel) %>%summarize(ks=mean(V17), n=length(unique(gene))) %>% filter(n>10) %>% ggplot()+ geom_vline(xintercept=seq(from=0,to=0.14,by=0.01), color='gray', lty='dashed', alpha=0.2) + geom_histogram(aes(x=ks, color=ploidy, fill=ploidy), binwidth=0.0001) + facet_wrap(~speciesLabel, ncol=1, scales='free_y') + scale_color_manual(values=ploidycolors)+ scale_fill_manual(values=ploidycolors) 
final_merged_data %>% group_by(queryChr, refChr, paspChr, genome,ploidy, speciesLabel) %>%summarize(ks=mean(V17), n=length(unique(gene))) %>% filter(n>10) %>% ggplot() + geom_vline(xintercept=seq(from=0,to=0.14,by=0.01), color='gray', lty='dashed', alpha=0.2)+ geom_histogram(aes(x=ks, color=ploidy, fill=ploidy), binwidth=0.001) + facet_wrap(~speciesLabel, ncol=1, scales='free_y') + scale_color_manual(values=ploidycolors)+ scale_fill_manual(values=ploidycolors)

final_merged_data %>% group_by(queryChr, refChr, paspChr, genome,ploidy, speciesLabel) %>%summarize(ks=mean(V17), n=length(unique(gene))) %>% filter(n>10) %>% ggplot()+ geom_vline(xintercept=seq(from=0,to=0.14,by=0.01), color='gray', lty='dashed', alpha=0.2) + geom_histogram(aes(x=ks, color=ploidy, fill=ploidy, weight=n), binwidth=0.0001) + facet_wrap(~speciesLabel, ncol=1, scales='free_y') + scale_color_manual(values=ploidycolors)+ scale_fill_manual(values=ploidycolors)
final_merged_data %>% group_by(queryChr, refChr, paspChr, genome,ploidy, speciesLabel) %>%summarize(ks=mean(V17), n=length(unique(gene))) %>% filter(n>10) %>% ggplot()+ geom_vline(xintercept=seq(from=0,to=0.14,by=0.01), color='gray', lty='dashed', alpha=0.2) + geom_histogram(aes(x=ks, color=ploidy, fill=ploidy, weight=n), binwidth=0.001) + facet_wrap(~speciesLabel, ncol=1, scales='free_y') + scale_color_manual(values=ploidycolors)+ scale_fill_manual(values=ploidycolors)

dev.off()

mmm=final_merged_data %>% group_by(queryChr, refChr, paspChr, genome,ploidy, speciesLabel) %>%summarize(ks=median(V17), n=n()) %>% filter(n>10)


## use this to look at ks and dnds for each!
ks=fread('../general_summaries/ks_to_look_for_mixtures.txt', header=T, quote='')
#ks=data.frame(ks)[ks$ks>0.001,]
#mks=ks[ks$ploidy %in%c('Tetraploid', 'Paleotetraploid', 'Hexaploid'),] %>% group_by(genome, ploidy, species, haploid) %>% dplyr::summarize(median=median(ks, na.rm=T), nonallelic=median(ks[ks>0.005], na.rm=T))
#mks$medianNonAllelicCorr=ifelse(mks$nonallelic-mks$median>0.01, mks$nonallelic, mks$median)



mu=6.5e-9
mks$mya=mks$medianNonAllelicCorr/2/mu/1e6
ks$mya=ks$ks/2/mu/1e6
ggplot(data.frame(ks)[ks$ploidy %in% c('Tetraploid', 'Paleotetraploid', 'Hexaploid') & !ks$genome%in%c('zdmomo', 'zdgigi', 'zdnica'),], aes(x=mya, color=ploidy, group=genome)) + geom_density() + scale_color_manual(values=ploidycolors) + facet_wrap(~ploidy, ncol=1, scales='free_y') + geom_vline(data=mks[!mks$genome %in% c('zdmomo', 'zdgigi', 'zdnica'),], aes(xintercept=mya))


asize$medianKs=mks$medianNonAllelicCorr[match(asize$V2, mks$genome)]
asize$mya=mks$mya[match(asize$V2, mks$genome)]
asize$mya[asize$ploidy=='Diploid']=0
asize$medianKs[asize$ploidy=='Diploid']=0



tppp=reshape2::melt(bb, id.vars=c('genome', 'Count'))
tppp$genome=factor(tppp$genome, levels=names(taxonnames))
tppp$ploidy=gs$ploidy[match(tppp$genome, gs$V2)]
tppp$phyleticcol=paste0(tppp$ploidy, tppp$variable)
## make singletons "monophyletic"
tppp$phyleticcol[tppp$Count==1]=paste0(tppp$ploidy[tppp$Count==1], 'Monophyletic')
tppp$species=taxonnames[match(tppp$genome, names(taxonnames))]
tppp$species=factor(tppp$species, levels=taxonnames)
tppp=tppp[!is.na(tppp$genome),]

## copy bar plots, no label
cpb=ggplot(tppp[tppp$doubledCount%in% 1:6 & !is.na(tppp$species),], aes(x=doubledCount, y=value, group=ploidy, fill=ploidy)) + geom_hline(yintercept=0, color='snow2', linetype='dotted') + geom_vline(xintercept=c(1,3,5), color='snow2', linetype='dotted') + geom_vline(xintercept=c(2,4,6), color='snow3', linetype='dotted') + geom_histogram(stat='identity', position='stack') + 
  facet_wrap(~speciesLabel, ncol=1, strip.position='left', drop=F) + scale_fill_manual(values=ploidycolors) +   theme( strip.background = element_blank(), strip.text.y.left = element_blank(), panel.spacing = unit(3, "pt"), axis.text.y=element_blank(), axis.text.x=element_text(size=9))+ theme(legend.position = "none", plot.margin = margin(l = -15)) + ylab('')+ xlab('Syntenic Gene\nCopy Number')+ scale_y_continuous(n.breaks = 3)
## auto/allo pies, no label
aap=tppp %>% group_by(speciesLabel, variable) %>% summarize(n=n(), count=sum(value)) %>% filter(variable!='NotApplicable') %>% mutate(pct = count/sum(count)*100, width=sum(count))%>%
  ggplot(aes(x=width/2, y=pct, fill=variable, width=width)) + geom_bar(stat='identity', position='fill') + coord_polar(theta='y') + facet_wrap(~speciesLabel, ncol=1, strip.position='left')+   theme( strip.background = element_blank(), strip.text.y.left = element_blank(), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + theme(axis.ticks=element_blank(), axis.text=element_blank(), panel.grid=element_blank(), panel.border=element_blank()) + scale_fill_manual(values=c('#5F4B8BFF', '#E69A8DFF')) + ylab('Proportion\nDuplicates\nMonophyletic') + xlab('')
## ks distributions, no label
#ksp=ggplot(ks[!is.na(ks$speciesLabel) ,], aes(x=V17,  color=ploidy, fill=ploidy))  + geom_vline(xintercept=c(0.05,0.15,0.25), color='snow2', linetype='dotted') + geom_vline(xintercept=c(0.1,0.2,0.3), color='snow3', linetype='dotted')+ geom_density() + scale_color_manual(values=ploidycolors) + scale_fill_manual(values=ploidycolors) +  facet_wrap(~speciesLabel, ncol=1, strip.position='left', scales='free_y', drop=F) + theme( strip.background = element_blank(), strip.text.y.left = element_blank(), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + ylab('')+ xlab('Ks Between\nDuplicates') + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + geom_vline(aes(xintercept=median), color='snow2', alpha=0.6)
## wait add mohamed's idea of histograms!
ksp=ggplot(ks[!is.na(ks$speciesLabel) & ks$V17>=0.001 & !is.na(ks$V17),], aes(x=V17,  color=ploidy, fill=ploidy))  + geom_vline(xintercept=c(0.05,0.15,0.25), color='snow2', linetype='dotted') + geom_vline(xintercept=c(0.1,0.2,0.3), color='snow3', linetype='dotted')+ geom_histogram(aes(y = after_stat(density)), binwidth=0.005, alpha=0.7, color=NA) + geom_density(alpha=1, fill=NA) + scale_color_manual(values=ploidycolors) + scale_fill_manual(values=ploidycolors) +  facet_wrap(~speciesLabel, ncol=1, strip.position='left', scales='free_y', drop=F) + theme( strip.background = element_blank(), strip.text.y.left = element_blank(), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + ylab('')+ xlab('Ks Between\nDuplicates') + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + geom_vline(aes(xintercept=median), color='dimgray', alpha=0.5) # + geom_vline(aes(xintercept=mtm1), color='dimgray', alpha=0.5)

