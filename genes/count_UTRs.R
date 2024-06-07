library(rtracklayer)
library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())



ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')


setwd('/Users/mcs368/Documents/GitHub/panand_assemblies/genes')
asize=fread('../general_summaries/panand_assembly_sizes.txt', header=T, quote="", fill=T)
asize$ploidy=factor(asize$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))

## ploidy groups
diploids=c('cserru', 'irugos', 'sbicol', 'ppanic', 'ttrian', 'crefra', 'avirgi', 'smicro', 'rtuber')
tetraploids=c('snutan', 'hconto', 'ccitra', 'achine', 'sscopa', 'etrips',  'vcuspi')
hexaploids=c('udigit', 'agerar', 'hcompr', 'blagur')

## function to generate per-geneID (ignoring transcript) count of what has a UTR
countUTRs <- function(filepath, id='') {
  # Load data
  data <- import.gff3(filepath)
  utrcount=data.frame(data) %>% group_by(gene=substr(ID,1,15)) %>% summarize(fiveutr=sum(type=='five_prime_UTR'), threeutr=sum(type=='three_prime_UTR'))
  
  utrcount$genome=id
  return(utrcount)
}

## file from google sheets that has gff3 name
genepath=read.table('panand_gene_annotations.txt', header=F, sep='\t')

asize=asize[asize$V2!='zluxur',]
asize$propFiveUTR=sapply(asize$V2, function(x) {
                          print(x)
                          blah=countUTRs(paste0('~/Downloads/Annotations_v4.1/', genepath$V4[genepath$V2==x]))
                         return(sum(blah$fiveutr==0)/nrow(blah))})

asize$propThreeUTR=sapply(asize$V2, function(x) {
  print(x)
  blah=countUTRs(paste0('~/Downloads/Annotations_v4.1/', genepath$V4[genepath$V2==x]))
  sum(blah$threeutr==0)/nrow(blah)})


asize$propAnyUTR=sapply(asize$V2, function(x) {
  print(x)
  blah=countUTRs(paste0('~/Downloads/Annotations_v4.1/', genepath$V4[genepath$V2==x]))
  sum(blah$fiveutr==0 | blah$threeutr==0)/nrow(blah)})


asize$nGenes=sapply(asize$V2, function(x) {
  print(x)
  blah=countUTRs(paste0('~/Downloads/Annotations_v4.1/', genepath$V4[genepath$V2==x]))
  nrow(blah)})

## remove maize and sorghum
asize=asize[!asize$V2 %in% c('zmB735', 'sbicol'),]

ggplot(asize, aes(x=propFiveUTR, y=propThreeUTR, color=ploidy, shape=haploid))+ geom_point() + scale_color_manual(values=ploidycolors) + geom_text(aes(label=V2))

ggplot(asize, aes(x=propAnyUTR, y=nGenes, color=ploidy, shape=haploid))+ geom_point() + scale_color_manual(values=ploidycolors)
ggplot(asize, aes(x=propAnyUTR, y=nGenes, color=ploidy, shape=haploid))+ geom_point() + scale_color_manual(values=ploidycolors) + geom_text(aes(label=V2))
## scale by haploid
ggplot(asize, aes(x=propAnyUTR, y=ifelse(haploid, nGenes, nGenes/2), color=ploidy, shape=haploid))+ geom_point(size=3) + scale_color_manual(values=ploidycolors) + ylab('haploid equivalent gene count') + xlab('proportion genes with any UTR on the gene')
ggplot(asize, aes(x=propAnyUTR, y=ifelse(haploid, nGenes, nGenes/2), color=ploidy, shape=haploid))+ geom_point(size=3) + scale_color_manual(values=ploidycolors) + geom_text(aes(label=V2)) + ylab('haploid equivalent gene count')+ xlab('proportion genes with any UTR on the gene')

ggplot(asize[asize$ploidy!='Paleotetraploid',], aes(x=propAnyUTR, y=ifelse(haploid, nGenes, nGenes/2), color=ploidy, shape=haploid))+ geom_point(size=3) + scale_color_manual(values=ploidycolors) + ylab('haploid equivalent gene count') + xlab('proportion genes with any UTR on the gene')
ggplot(asize[asize$ploidy!='Paleotetraploid',], aes(x=propAnyUTR, y=ifelse(haploid, nGenes, nGenes/2), color=ploidy, shape=haploid))+ geom_point(size=3) + scale_color_manual(values=ploidycolors) + geom_text(aes(label=V2)) + ylab('haploid equivalent gene count')+ xlab('proportion genes with any UTR on the gene')

## diploid equivalents

mediandip=median(ifelse(asize$haploid, asize$nGenes, asize$nGenes/2)[asize$ploidy=='Diploid'])
ggplot(asize[asize$ploidy!='Paleotetraploid',], aes(x=propAnyUTR, y=ifelse(haploid, nGenes, nGenes/2), color=ploidy, shape=haploid))+ geom_point(size=2) + scale_color_manual(values=ploidycolors) + ylab('haploid equivalent gene count') + xlab('proportion genes with any UTR on the gene') + geom_hline(yintercept=mediandip*1:3, color='gray')
ggplot(asize[asize$ploidy!='Paleotetraploid',], aes(x=propAnyUTR, y=ifelse(haploid, nGenes, nGenes/2), color=ploidy, shape=haploid))+ geom_point(size=2) + scale_color_manual(values=ploidycolors) + geom_text(aes(label=V2)) + ylab('haploid equivalent gene count')+ xlab('proportion genes with any UTR on the gene')+ geom_hline(yintercept=mediandip*1:3, color='gray')



## sloppy and counting here because they're not in the supplemental table :(
#asize$rnaseqlibs=c(4,7,5,8,9,5,7,6,7,7,6,7,7,5,8,7,6,7,9,8,8,8,7,7,7,5,5,7,2,2,2,6,0,1,NA)
rnaseqlibs=c(achine=4, agerar=7, atenui=5, avirgi=8,blagur=9,
             cserru=5, ccitra=6 , crefra=7, etrips=6, hcompr=7, hconto=6,
             irugos=7, ppanic=7, rrottb=5, rtuber=8, smicro=7, sscopa=6,
             snutan=7, telega=9, ttrian=8, tdacn1=9, tdacs1=8, udigit=7, vcuspi=5,
             zdgigi=10, zdmomo=7, zmhuet=2, zTIL18=2, zTIL25=1, zTIL01=6, zTIL11=0, 
             znicar=1)
asize$rnaseqlibs=rnaseqlibs[match(asize$V2, names(rnaseqlibs))]




## scale by haploid
ggplot(asize, aes(x=propAnyUTR, y=ifelse(haploid, nGenes, nGenes/2), color=ploidy, shape=haploid, size=rnaseqlibs))+ geom_point() + scale_color_manual(values=ploidycolors) + ylab('haploid equivalent gene count') + xlab('proportion genes with any UTR on the gene')
ggplot(asize, aes(x=propAnyUTR, y=ifelse(haploid, nGenes, nGenes/2), color=ploidy, shape=haploid, size=rnaseqlibs))+ geom_point() + scale_color_manual(values=ploidycolors) + geom_text(aes(label=V2)) + ylab('haploid equivalent gene count')+ xlab('proportion genes with any UTR on the gene')

ggplot(asize[asize$ploidy!='Paleotetraploid',], aes(x=propAnyUTR, y=ifelse(haploid, nGenes, nGenes/2), color=ploidy, shape=haploid, size=rnaseqlibs))+ geom_point() + scale_color_manual(values=ploidycolors) + ylab('haploid equivalent gene count') + xlab('proportion genes with any UTR on the gene')
ggplot(asize[asize$ploidy!='Paleotetraploid',], aes(x=propAnyUTR, y=ifelse(haploid, nGenes, nGenes/2), color=ploidy, shape=haploid, size=rnaseqlibs))+ geom_point() + scale_color_manual(values=ploidycolors) + geom_text(aes(label=V2)) + ylab('haploid equivalent gene count')+ xlab('proportion genes with any UTR on the gene')


summary(lm(asize$propAnyUTR~asize$rnaseqlibs+asize$haploid+asize$nGenes+asize$ploidy))
summary(lm(asize$propAnyUTR~asize$rnaseqlibs+asize$haploid*asize$nGenes))



## what if it's homeologous divergence 

## now add ks to these values
ks=fread('../general_summaries/ks_to_look_for_mixtures.txt', header=T, quote='')

ks=data.frame(ks)[ks$ks>0.001,]

mks=ks[ks$ploidy %in%c('Tetraploid', 'Paleotetraploid', 'Hexaploid'),] %>% group_by(genome, ploidy, species, haploid) %>% dplyr::summarize(median=median(ks, na.rm=T), nonallelic=median(ks[ks>0.005], na.rm=T))

mks$medianNonAllelicCorr=ifelse(mks$nonallelic-mks$median>0.01, mks$nonallelic, mks$median)

mu=6.5e-9
mks$mya=mks$medianNonAllelicCorr/2/mu/1e6
ks$mya=ks$ks/2/mu/1e6


asize$medianKs=mks$medianNonAllelicCorr[match(asize$V2, mks$genome)]
asize$mya=mks$mya[match(asize$V2, mks$genome)]
asize$mya[asize$ploidy=='Diploid']=0
asize$medianKs[asize$ploidy=='Diploid']=0



## scale by haploid, add age
ggplot(asize, aes(x=propAnyUTR, y=mya, color=ploidy, shape=haploid, size=rnaseqlibs))+ geom_point() + scale_color_manual(values=ploidycolors) + ylab('mya') + xlab('proportion genes with any UTR on the gene')
ggplot(asize, aes(x=propAnyUTR, y=mya, color=ploidy, shape=haploid, size=rnaseqlibs))+ geom_point() + scale_color_manual(values=ploidycolors) + geom_text(aes(label=V2)) + ylab('mya')+ xlab('proportion genes with any UTR on the gene')

ggplot(asize[asize$ploidy!='Paleotetraploid',], aes(x=propAnyUTR, y=ifelse(haploid, nGenes, nGenes/2), color=ploidy, shape=haploid, size=rnaseqlibs))+ geom_point() + scale_color_manual(values=ploidycolors) + ylab('haploid equivalent gene count') + xlab('proportion genes with any UTR on the gene')
ggplot(asize[asize$ploidy!='Paleotetraploid',], aes(x=propAnyUTR, y=ifelse(haploid, nGenes, nGenes/2), color=ploidy, shape=haploid, size=rnaseqlibs))+ geom_point() + scale_color_manual(values=ploidycolors) + geom_text(aes(label=V2)) + ylab('haploid equivalent gene count')+ xlab('proportion genes with any UTR on the gene')


summary(lm(asize$propAnyUTR~asize$rnaseqlibs+asize$haploid+asize$nGenes+asize$ploidy + asize$mya))


summary(lm(asize$propAnyUTR~asize$rnaseqlibs*asize$mya+asize$haploid+asize$nGenes+asize$ploidy))

summary(lm(asize$propAnyUTR~asize$rnaseqlibs+asize$haploid*asize$nGenes))




### zong yan generated 
utrlen=read.csv('~/Downloads/Annotations_v4.1_UTR_medians.csv', header=T)
utrlen$V2=genepath$V2[match(utrlen$Filename, genepath$V4)]
asize$fiveUTRMedianLen=utrlen$Five_Prime_UTR_Median[match(asize$V2, utrlen$V2)]
asize$threeUTRMedianLen=utrlen$Three_Prime_UTR_Median[match(asize$V2, utrlen$V2)]

ggplot(asize, aes(x=fiveUTRMedianLen, y=threeUTRMedianLen, color=ploidy, shape=haploid, size=rnaseqlibs))+ geom_point() + scale_color_manual(values=ploidycolors) + ylab('median threeUTRLen') + xlab('median fiveUTRLen')
ggplot(asize, aes(x=fiveUTRMedianLen, y=threeUTRMedianLen, color=ploidy, shape=haploid, size=rnaseqlibs))+ geom_point() + scale_color_manual(values=ploidycolors) + geom_text(aes(label=V2)) + ylab('median threeUTRLen')+ xlab('median fiveUTRLen')


