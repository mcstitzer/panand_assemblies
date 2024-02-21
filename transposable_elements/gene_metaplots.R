library(ggtree) ## /programs/R-4.2.1-r9/bin/R
library(ggplot2) 
library(cowplot)
theme_set(theme_cowplot())
library(rtracklayer)
library(reshape2)
library(RColorBrewer)
library(tidypaleo) ## facet species names in italics!!!!!
library(tidyverse)
library(viridisLite)

all=read.table('../panand_sp_ploidy.txt')
all=all[!all$V2 %in% c('pprate', 'tdactm', 'tzopol', 'osativ', 'bdista', 'agerjg', 'svirid', 'eophiu', 'tdacs2', 'tdacn2'),]

genomecountlist=vector(mode = "list", length = length(all$V2))
names(genomecountlist)=all$V2

repeatlengthslist=vector(mode = "list", length = length(all$V2))
names(repeatlengthslist)=all$V2


gfftypestokeep=c("Gypsy_LTR_retrotransposon", 
"LTR_retrotransposon", "Copia_LTR_retrotransposon", "CACTA_TIR_transposon", 
"helitron", "hAT_TIR_transposon", "PIF_Harbinger_TIR_transposon", 
"Tc1_Mariner_TIR_transposon", "Mutator_TIR_transposon", 'tandem_repeat') 

for( genotype in all$V2){
##import gff3
a=import.gff3(paste0('../trash/repeatmask_tandems/', all$V1[all$V2==genotype], '_EDTATandemRepeat.gff3'))
a=a[a$type %in% gfftypestokeep,]
a$genome=genotype
# ## get each fam
# temp=data.frame(a) #%>% group_by(fam=gsub('_LTR', '', gsub('_INT', '', Name)), Classification) %>% dplyr::filter(!is.na(as.numeric(Identity))) 
# if('Note' %in% colnames(temp)){temp$Note=''}
# if(genotype=='zmB735'){temp$Parent=''
#                       temp$ltr_identity=NA
#                       temp=temp[,colnames(temp) %in% colnames(genomecountlist[[1]])]}
#   if(genotype=='znicar'){ ## these guys edta is weird - chr not named with chr!!
#                       temp$seqnames[temp$seqnames %in% 1:10]=paste0('chr', temp$seqnames[temp$seqnames %in% 1:10])}
# temp$genome=genotype

genomecountlist[[genotype]]=a
}

# genomecount=do.call(c, genomecountlist)
# genomecount$Identity=as.numeric(genomecount$Identity)

## add a superfamily based on matchign up to the classification field - note most of these NAs are relics from the B73 annotation being included!!!!!
# genomecount$sup=c(NA, 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', 'DHH', NA,NA,NA,NA,NA,NA,'RLC', 'RLG', 'RLG', 'RLX', 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', NA,NA,NA)[match(genomecount$Classification, c("Cent/CentC", "DNA/DTA", "DNA/DTC", "DNA/DTH", "DNA/DTM", "DNA/DTT", 
# "DNA/Helitron", "knob/knob180", "knob/TR-1", "LINE/L1", "LINE/RTE", 
# "LINE/unknown", "Low_complexity", "LTR/Copia", "LTR/CRM", "LTR/Gypsy", 
# "LTR/unknown", "MITE/DTA", "MITE/DTC", "MITE/DTH", "MITE/DTM", 
# "MITE/DTT", "rDNA/spacer", "Simple_repeat", "subtelomere/4-12-1"))]

## tandem repeat maybe want to do by length class???
#genomecount$sup[genomecount$source=='TRASH']=genomecount$Name[genomecount$source=='TRASH']
# genomecount$sup[genomecount$source=='TRASH']='TandemRepeat'



syn=fread('~/transfer/anchorPositions_AnchorWave_Paspalum.txt.gz')
syn=syn[!syn$genome %in% c('pprate', 'svirid', 'tdacn2', 'tdacs2', 'tdactm', 'tzopol', 'bdista', 'agerjg', 'eophiu', 'osativ'),]

syngr=GRanges(seqnames=syn$queryChr, IRanges(start=syn$queryStart, end=syn$queryEnd), strand=syn$strand)
mcols(syngr)$genome=syn$genome
mcols(syngr)$quickgene=syn$quickgene

syngr3=syngr
start(syngr3)=end(syngr)
## get upstream
flankspace=10000
geneflanks=promoters(syngr, upstream=flankspace, downstream=1000)
geneflanks3=promoters(syngr3, upstream=1000, downstream=flankspace)


#geneflanksranges=slide_ranges(geneflanks, width=100, step=10)
## easier to not overlap tiles...
geneflanksranges=tile_ranges(geneflanks, width=100)
#geneflanksranges=slide_ranges(geneflanks, width=100, step=10)
geneflanksranges3=tile_ranges(geneflanks3, width=100)

## NOW NEED TO PUT IN 3' END - could concatenate the granges together? need to set up naming of windows carefully then
## or keep as separate side-by-side plots?
geneflanksranges3$genome=geneflanks3$genome[geneflanksranges3$partition]
geneflanksranges3$ogstrand=strand(geneflanks3)[geneflanksranges3$partition]

geneflanksranges$genome=geneflanks$genome[geneflanksranges$partition]
geneflanksranges$ogstrand=strand(geneflanks)[geneflanksranges$partition]
# geneflanksranges$window=rep(1:191,length(geneflanks))
# metaplot=data.frame(window=1:191)
nwindows=names(table(table(geneflanksranges$partition)))
geneflanksranges$window=rep(1:nwindows,length(geneflanks))
geneflanksranges3$window=rep((as.numeric(nwindows)+1):(as.numeric(nwindows)*2),length(geneflanks))

metaplot=data.frame(window=1:(as.numeric(nwindows)*2))
metaplot75=data.frame(window=1:(as.numeric(nwindows)*2))
metaplot25=data.frame(window=1:(as.numeric(nwindows)*2))
metaplotmean=data.frame(window=1:(as.numeric(nwindows)*2))

### combine them here!!!
geneflanksranges=c(geneflanksranges, geneflanksranges3)

pdf('~/transfer/try_te_metaplot.pdf',12,8)

for(genome in all$V2){
  tesg=genomecountlist[[genome]]
 tes=unstrand(tesg) ### why is this so stupid
 gf=geneflanksranges[geneflanksranges$genome==genome,]
  if(genome=='znicar'){seqnames(tes)[seqnames(tes) %in% 1:10]=paste0('chr', seqnames(tes)[seqnames(tes) %in% 1:10])}
  if(genome=='sbicol'){seqlevels(gf)[seqlevels(gf) %in% 1:9]=paste0('Chr0', seqlevels(gf)[seqlevels(gf) %in% 1:9])
                      seqlevels(gf)[seqlevels(gf)=='10']='Chr10'}

 tewindow=join_overlap_intersect(unstrand(gf), tes)      ## cut tes at boundaries of ranges

posplot=data.frame(tewindow[tewindow$ogstrand=='+',])[,c('window', 'width', 'partition')]
  posplot=posplot %>% complete(partition, window, fill=list(width=0))

  ## now ask which is the closest/furthest TE for each subgenome (gene copy)
  posplot$quickgene=geneflanks$quickgene[posplot$partition]
  posplothigh=posplot %>% group_by(quickgene,window) %>% summarize(partition[which.max(width)])
print(  ggplot(posplot, aes(x=window, y=width, group=window)) + geom_boxplot(outlier.shape=NA) + ggtitle(genome) )
# print( ggplot(posplot, aes(x=window, y=width, group=partition)) + geom_line(alpha=0.01) + ggtitle(genome) )

  metaplot[,genome]=(posplot %>% group_by(window) %>% summarize(medianTE=median(width)))$medianTE
  metaplot75[,genome]=(posplot %>% group_by(window) %>% summarize(percTE=quantile(width, 0.75)))$percTE
  metaplot25[,genome]=(posplot %>% group_by(window) %>% summarize(percTE=quantile(width, 0.25)))$percTE
  metaplotmean[,genome]=(posplot %>% group_by(window) %>% summarize(meanTE=median(width)))$meanTE

  
}

metaplotmelt=melt(metaplot, id.vars='window')
metaplotmeanmelt=melt(metaplotmean, id.vars='window')
metaplot75melt=melt(metaplot75, id.vars='window')
metaplot25melt=melt(metaplot25, id.vars='window')

gs=read.table('~/transfer/panand_assembly_sizes.txt', header=T, sep='\t')
## to color by te prooportion
gs$teprop=gs$haploidRepeatSize/(gs$haploidAssemblySize-gs$haploidNCount)

metaplotmelt$ploidy=gs$ploidy[match(metaplotmelt$variable, gs$V2)]
metaplotmelt$ploidy=factor(metaplotmelt$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))
metaplot25melt$ploidy=gs$ploidy[match(metaplot25melt$variable, gs$V2)]
metaplot25melt$ploidy=factor(metaplot25melt$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))
metaplot75melt$ploidy=gs$ploidy[match(metaplot75melt$variable, gs$V2)]
metaplot75melt$ploidy=factor(metaplot75melt$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))
metaplotmeanmelt$ploidy=gs$ploidy[match(metaplotmeanmelt$variable, gs$V2)]
metaplotmeanmelt$ploidy=factor(metaplotmeanmelt$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))

ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')


ggplot(metaplotmelt, aes(group=variable, x=window, y=value, color=ploidy))+ geom_vline(xintercept=(1:(((flankspace*2)+2000)/1000)*10), color='whitesmoke') + geom_line() + scale_color_manual(values=ploidycolors)  + xlab('Base pairs away from TSS/TTS') + ylab('Median TE base pairs in 100bp window') + geom_vline(xintercept=c(flankspace/100, (flankspace/100)+10, (flankspace/100)+20), lty=c('dashed', 'solid', 'dashed')) + scale_x_continuous(breaks=c(1, flankspace/100, (flankspace/100)+10, (flankspace/100)+20, ((flankspace*2)+2000)/100), labels=c(paste0('-', flankspace), 'TranslationSS', 'X', 'TranslationTermS', paste0('+', flankspace))) 
ggplot(metaplot75melt, aes(group=variable, x=window, y=value, color=ploidy))+ geom_vline(xintercept=(1:(((flankspace*2)+2000)/1000)*10), color='whitesmoke') + geom_line() + scale_color_manual(values=ploidycolors)  + xlab('Base pairs away from TSS/TTS') + ylab('75th percentile TE base pairs in 100bp window') + geom_vline(xintercept=c(flankspace/100, (flankspace/100)+10, (flankspace/100)+20), lty=c('dashed', 'solid', 'dashed'))+ scale_x_continuous(breaks=c(1, flankspace/100, (flankspace/100)+10, (flankspace/100)+20, ((flankspace*2)+2000)/100), labels=c(paste0('-', flankspace), 'TranslationSS', 'X', 'TranslationTermS', paste0('+', flankspace)))
ggplot(metaplot25melt, aes(group=variable, x=window, y=value, color=ploidy))+ geom_vline(xintercept=(1:(((flankspace*2)+2000)/1000)*10), color='whitesmoke') + geom_line() + scale_color_manual(values=ploidycolors)  + xlab('Base pairs away from TSS/TTS') + ylab('25th percentile TE base pairs in 100bp window') + geom_vline(xintercept=c(flankspace/100, (flankspace/100)+10, (flankspace/100)+20), lty=c('dashed', 'solid', 'dashed'))+ scale_x_continuous(breaks=c(1, flankspace/100, (flankspace/100)+10, (flankspace/100)+20, ((flankspace*2)+2000)/100), labels=c(paste0('-', flankspace), 'TranslationSS', 'X', 'TranslationTermS', paste0('+', flankspace)))
ggplot(metaplotmeanmelt, aes(group=variable, x=window, y=value, color=ploidy)) + geom_vline(xintercept=(1:(((flankspace*2)+2000)/1000)*10), color='whitesmoke')+ geom_line() + scale_color_manual(values=ploidycolors)  + xlab('Base pairs away from TSS/TTS') + ylab('Mean TE base pairs in 100bp window') + geom_vline(xintercept=c(flankspace/100, (flankspace/100)+10, (flankspace/100)+20), lty=c('dashed', 'solid', 'dashed'))+ scale_x_continuous(breaks=c(1, flankspace/100, (flankspace/100)+10, (flankspace/100)+20, ((flankspace*2)+2000)/100), labels=c(paste0('-', flankspace), 'TranslationSS', 'X', 'TranslationTermS', paste0('+', flankspace)))

metaplotmelt$teprop=gs$teprop[match(metaplotmelt$variable, gs$V2)]
ggplot(metaplotmelt, aes(group=variable, x=window, y=value, color=teprop)) + geom_vline(xintercept=(1:(((flankspace*2)+2000)/1000)*10), color='whitesmoke') + geom_line() + scale_color_viridis_c(option='inferno')  + xlab('window index relative to TranslationSS (dashed)') + ylab('Median TEs in 100bp window') + geom_vline(xintercept=c(flankspace/100, (flankspace/100)+10, (flankspace/100)+20), lty=c('dashed', 'solid', 'dashed'))+ scale_x_continuous(breaks=c(1, flankspace/100, (flankspace/100)+10, (flankspace/100)+20, ((flankspace*2)+2000)/100), labels=c(paste0('-', flankspace), 'TranslationSS', 'X', 'TranslationTermS', paste0('+', flankspace)))

## also want to do this with ks - time since polyploidy!!@! 
## younger polyploids might not have had TE invasions yet
## if i get subgenomes, same question!!!!!

ks=fread('~/transfer/ks_to_look_for_mixtures.txt', sep='\t', quote='')
ks=ks %>% group_by(genome) %>% filter(ks>0.0001) %>% summarize(ks=median(ks))
metaplotmelt$ks=ks$ks[match(metaplotmelt$variable, ks$genome)]
ggplot(metaplotmelt, aes(group=variable, x=window, y=value, color=ks))+ geom_vline(xintercept=(1:(((flankspace*2)+2000)/1000)*10), color='whitesmoke') + geom_line() + scale_color_viridis_c(option='viridis')  + xlab('Base pairs away from TSS/TTS') + ylab('Median TE base pairs in 100bp window') + geom_vline(xintercept=c(flankspace/100, (flankspace/100)+10, (flankspace/100)+20), lty=c('dashed', 'solid', 'dashed'))+ scale_x_continuous(breaks=c(1, flankspace/100, (flankspace/100)+10, (flankspace/100)+20, ((flankspace*2)+2000)/100), labels=c(paste0('-', flankspace), 'TranslationSS', 'X', 'TranslationTermS', paste0('+', flankspace)))
ggplot(metaplotmelt, aes(group=variable, x=window, y=value, color=ks)) + geom_vline(xintercept=(1:(((flankspace*2)+2000)/1000)*10), color='whitesmoke')+ geom_line() + scale_color_viridis_c(limits = c(0, 0.06), oob = scales::squish, option='viridis')  + xlab('window index relative to TranslationSS (dashed)') + ylab('Median TEs in 100bp window') + geom_vline(xintercept=c(flankspace/100, (flankspace/100)+10, (flankspace/100)+20), lty=c('dashed', 'solid', 'dashed'))+ scale_x_continuous(breaks=c(1, flankspace/100, (flankspace/100)+10, (flankspace/100)+20, ((flankspace*2)+2000)/100), labels=c(paste0('-', flankspace), 'TranslationSS', 'X', 'TranslationTermS', paste0('+', flankspace)))

metaplotmelt$genomesize=(gs$haploidAssemblySize-gs$haploidNCount)[match(metaplotmelt$variable, gs$V2)]/1e6

ggplot(metaplotmelt, aes(group=variable, x=window, y=value, color=genomesize)) + geom_vline(xintercept=(1:(((flankspace*2)+2000)/1000)*10), color='whitesmoke')+ geom_line() + scale_color_viridis_c( option='magma')  + xlab('Base pairs away from TSS/TTS') + ylab('Median TE base pairs in 100bp window') + geom_vline(xintercept=c(flankspace/100, (flankspace/100)+10, (flankspace/100)+20), lty=c('dashed', 'solid', 'dashed'))+ scale_x_continuous(breaks=c(1, flankspace/100, (flankspace/100)+10, (flankspace/100)+20, ((flankspace*2)+2000)/100), labels=c(paste0('-', flankspace), 'TranslationSS', 'X', 'TranslationTermS', paste0('+', flankspace)))

dev.off()
