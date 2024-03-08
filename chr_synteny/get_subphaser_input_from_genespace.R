## or at least try!!
library(GENESPACE)
library(ggplot2)
library(cowplot)

## arun previously ran genespace!
# gpar <- init_genespace(
#   wd = "/work/mash-covid/genespace/GeneSpace_run6",
#   path2mcscanx = "/opt/MCScanX/")
# gpar <- run_genespace(gsParam = gpar)

gpar <- init_genespace(
  wd = "/local/workdir/mcs368/panand_htt/genespace/GeneSpace_run6", 
  path2mcscanx = "/programs/MCScanX/")

## okay it's smart enough to find what's there??
#out <- run_genespace(gpar, overwrite = F)

## this will load in gsParam object
load('results/gsParams.rda')

## omfg i have to change all these paths because they're fucking hard coded. oh and these are all considered haploid wtf i hate this
gsParam$paths=lapply(gsParam$paths, function(x) gsub('/work/mash-covid/genespace/', '/local/workdir/mcs368/panand_htt/genespace/', x))
gsParam$synteny$combBed=gsub('/work/mash-covid/genespace/', '/local/workdir/mcs368/panand_htt/genespace/', gsParam$synteny$combBed)
gsParam$synteny$SpeciesIDs=gsub('/work/mash-covid/genespace/', '/local/workdir/mcs368/panand_htt/genespace/', gsParam$synteny$SpeciesIDs)
gsParam$synteny$SequenceIDs=gsub('/work/mash-covid/genespace/', '/local/workdir/mcs368/panand_htt/genespace/', gsParam$synteny$SequenceIDs)
gsParam$synteny$ogs=gsub('/work/mash-covid/genespace/', '/local/workdir/mcs368/panand_htt/genespace/', gsParam$synteny$ogs)
gsParam$synteny$blast$queryBlast=sapply(gsParam$synteny$blast$queryBlast, function(x) gsub('/work/mash-covid/genespace/', '/local/workdir/mcs368/panand_htt/genespace/', x))
gsParam$synteny$blast$targetBlast=sapply(gsParam$synteny$blast$targetBlast, function(x) gsub('/work/mash-covid/genespace/', '/local/workdir/mcs368/panand_htt/genespace/', x))
gsParam$synteny$blast$allBlast=sapply(gsParam$synteny$blast$allBlast, function(x) gsub('/work/mash-covid/genespace/', '/local/workdir/mcs368/panand_htt/genespace/', x))
gsParam$synteny$blast$synHits=sapply(gsParam$synteny$blast$synHits, function(x) gsub('/work/mash-covid/genespace/', '/local/workdir/mcs368/panand_htt/genespace/', x))

setwd('../genespace/')

## this will load in gsParam object
load('results/gsParams.rda')

setwd('GeneSpace_run6')

## this will load in gsParam object
load('results/gsParams.rda')


                                     ## omfg i have to change all these paths because they're fucking hard coded. oh and these are all considered haploid wtf i hate this
gsParam$paths=lapply(gsParam$paths, function(x) gsub('/work/mash-covid/genespace/', '/local/workdir/mcs368/panand_htt/genespace/', x))
gsParam$synteny$combBed=gsub('/work/mash-covid/genespace/', '/local/workdir/mcs368/panand_htt/genespace/', gsParam$synteny$combBed)
gsParam$synteny$SpeciesIDs=gsub('/work/mash-covid/genespace/', '/local/workdir/mcs368/panand_htt/genespace/', gsParam$synteny$SpeciesIDs)
gsParam$synteny$SequenceIDs=gsub('/work/mash-covid/genespace/', '/local/workdir/mcs368/panand_htt/genespace/', gsParam$synteny$SequenceIDs)
gsParam$synteny$ogs=gsub('/work/mash-covid/genespace/', '/local/workdir/mcs368/panand_htt/genespace/', gsParam$synteny$ogs)
gsParam$synteny$blast$queryBlast=sapply(gsParam$synteny$blast$queryBlast, function(x) gsub('/work/mash-covid/genespace/', '/local/workdir/mcs368/panand_htt/genespace/', x))
gsParam$synteny$blast$targetBlast=sapply(gsParam$synteny$blast$targetBlast, function(x) gsub('/work/mash-covid/genespace/', '/local/workdir/mcs368/panand_htt/genespace/', x))
gsParam$synteny$blast$allBlast=sapply(gsParam$synteny$blast$allBlast, function(x) gsub('/work/mash-covid/genespace/', '/local/workdir/mcs368/panand_htt/genespace/', x))
gsParam$synteny$blast$synHits=sapply(gsParam$synteny$blast$synHits, function(x) gsub('/work/mash-covid/genespace/', '/local/workdir/mcs368/panand_htt/genespace/', x))


## order bottom to top!!!
#rip=c('zB73v5', 'zluxur', 'tdacs1', 'tdacn1', 'udigit', 'vcuspi', 'hcompr', 'etrips', 'avirgi', 'achine', 'agerar', 'hconto', 'ppanic', 'sbicol', 'snutan', 'pvagin')
rip=c('pvagin', 'sbicol', 'ppanic', 'hconto', 'avirgi', 'etrips', 'udigit', 'tdacn1', 'tdacs1', 'znicar', 'zTIL25', 'zB73v5')
simpledips=c('pvagin', 'sbicol', 'avirgi', 'tdacn1', 'tdacs1', 'znicar', 'zTIL25', 'zB73v5')
otherpoly=c('pvagin', 'sbicol',  'hconto', 'hcompr','udigit')

invchr <- data.frame(
  genome = c("avirgi", "avirgi", "avirgi", "avirgi", rep(c('tdacs1', 'tdacn1'),each=6), 'znicar', 'zTIL25', 'zB73v5', 'znicar', 'zTIL25', 'zB73v5'), 
  chr = c('chr5','chr1','chr10','chr7', rep(paste0('chr', c(10,3,2,6,9,14)),2), 'chr1', 'chr1', 'chr1', 'chr6', 'chr6', 'chr6'))

invchr <- data.frame(
  genome = c("avirgi", "avirgi", "avirgi", "avirgi", rep(c('tdacs1', 'tdacn1'),each=6), 'zTIL11', 'zTIL01', 'zTIL18', 'zdmomo', 'zdgigi','znicar', 'zmhuet', 'zTIL25', 'zB73v5', 'znicar', 'zmhuet', 'zTIL25', 'zB73v5','zTIL11', 'zTIL01', 'zTIL18', 'zdmomo', 'zdgigi'), 
  chr = c('chr5','chr1','chr10','chr7', rep(paste0('chr', c(10,3,2,6,9,14)),2), 'chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr6', 'chr6', 'chr6', 'chr6', 'chr6', 'chr6', 'chr6', 'chr6', 'chr6'))



## subgenomes
sg=read.table('~/transfer/panand_all_subgenomes_10052023.txt', header=F)
sgb73=sg[grepl('zB73v5', sg$V1),]
sgb73$queryChr=gsub('zB73v5\\.','', sgb73$V1)
sgb73=sgb73[!is.na(sgb73$V2) & !is.na(sgb73$V3) & !is.na(sgb73$queryChr),]
## this was a bad idea: sgb73$referenceChr=gsub('zB73v5\\.','', sgb73$V1) ## jsut do both so i can use it for both :)
library(rtracklayer)
sgr1=reduce(GRanges(seqnames=sgb73$queryChr[sgb73$V13=='maize1' & !is.na(sgb73$V13)], IRanges(start=sgb73$V2[sgb73$V13=='maize1' & !is.na(sgb73$V13)], end=sgb73$V3[sgb73$V13=='maize1' & !is.na(sgb73$V13)])), min.gapwidth=3e6)
sgr2=reduce(GRanges(seqnames=sgb73$queryChr[sgb73$V13=='maize2' & !is.na(sgb73$V13)], IRanges(start=sgb73$V2[sgb73$V13=='maize2' & !is.na(sgb73$V13)], end=sgb73$V3[sgb73$V13=='maize2' & !is.na(sgb73$V13)])), min.gapwidth=3e6)
sgr1$subgenome='maize1' 
sgr2$subgenome='maize2'
sgrs=c(sgr1,sgr2)            ## hahah it's fine they don't have overlapping chromosomes :)       
sgrs=sgrs[width(sgrs)>5e5,]
subgenome=data.frame(genome='zB73v5', chr=seqnames(sgrs), start=start(sgrs), end=end(sgrs), color=ifelse(sgrs$subgenome=='maize1','blue','red')) ## how embarassing i forgot james's colors

ls()
head(subgenome)
tail(subgenome)
tetraploid=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=c('pvagin', 'snutan', 'hconto', 'ccitra', 'achine', 'sscopa', 'etrips',  'vcuspi'), #all$V2[all$ploidy=='Tetraploid']), #'atenui', 'rrottb',
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=5e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr, xlabel='',
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar','zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zdmomo', 'zdgigi', 'tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes
                    )

head(tetraploid$plotData$sourceData$chromosomes
)
ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))

tetraploid=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=c('pvagin', 'snutan', 'hconto', 'ccitra', 'achine', 'sscopa', 'etrips',  'vcuspi'), #all$V2[all$ploidy=='Tetraploid']), #'atenui', 'rrottb',
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=5e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr, xlabel='',
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar','zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zdmomo', 'zdgigi', 'tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                   braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes
                    )

head(tetraploid$plotData$sourceData$chromosomes)
head(tetraploid$plotData$sourceData)
head(tetraploid$plotData)
head(tetraploid$plotData$sourceData)
head(tetraploid$blks)
head(tetraploid$blks[tetraploid$blks$genome1=='pvagin',])
blk
blk=tetraploid$blks
head(blk
)
table(blk$genome1)
table(blk$genome2)
table(blk$genome2[blk$genome1=='pvagin'])
table(blk$genome1[blk$genome2=='pvagin'])
head(blk)
blkp=blk[blk$genome1=='pvagin',]
head(blkp)
sum(blkp$genome2=='pvagin')
blkp %>% group_by(genome2, chr2, refChr) %>% summarize(n=n())
library(dplyr)
blkp %>% group_by(genome2, chr2, refChr) %>% summarize(n=n())
blkp %>% group_by(genome2, chr2, refChr) %>% summarize(n=n()) %>% pivot_wider(names_from=refChr, values_from=chr2)
library(tidyr)
blkp %>% group_by(genome2, chr2, refChr) %>% summarize(n=n()) %>% pivot_wider(names_from=refChr, values_from=chr2)
blkp %>% group_by(genome2, chr2, refChr) %>% summarize(n=n()) %>% pivot_wider(names_from=n, values_from=chr2)
blkp %>% group_by(genome2, chr2, refChr) %>% summarize(n=n()) %>% pivot_wider(names_from=chr2, values_from=n)
blkp %>% group_by(genome2, chr2, refChr) %>% summarize(n=n()) %>% pivot_wider(names_from=chr2, values_from=n, values_fill=0)
blkpp=blkp %>% group_by(genome2, chr2, refChr) %>% summarize(n=n()) %>% pivot_wider(names_from=chr2, values_from=n, values_fill=0)
rowSums(blkpp[,-c(1,2)])
head(blkpp)
blkp %>% group_by(genome2, chr2, refChr) %>% summarize(n=n()) %>% pivot_wider(names_from=chr2, values_from=n, values_fill=0)
head(blkp)
blkp %>% filter_by(nHits2>100) %>% group_by(genome2, chr2, refChr) %>% summarize(n=n()) %>% pivot_wider(names_from=chr2, values_from=n, values_fill=0)
blkp %>% filter(nHits2>100) %>% group_by(genome2, chr2, refChr) %>% summarize(n=n()) %>% pivot_wider(names_from=chr2, values_from=n, values_fill=0)
blkpp=blkp %>% filter(nHits2>100) %>% group_by(genome2, chr2, refChr) %>% summarize(n=n()) %>% pivot_wider(names_from=chr2, values_from=n, values_fill=0)
rowSums(blkpp[,-c(1,2)])
data.frame(blkpp$genome2, blkpp$refChr, rowSums(blkpp[,-c(1,2)]))
blkp %>% filter(nHits2>100) %>% group_by(genome2, chr2, refChr) %>% summarize(n=n())
blkpp=blkp %>% filter(nHits2>100) %>% group_by(genome2, refChr) %>% summarize(val=paste0(chr2, collapse=',')) %>% pivot_wider(names_from=refChr, values_from=val)
head(blkpp)
blkpp
tail(blkpp)
tail(blkpp,20)
blkpp=blkp %>% filter(nHits2>100) %>% group_by(genome2, refChr) %>% summarize(val=paste0(unique(chr2), collapse=',')) %>% pivot_wider(names_from=refChr, values_from=val)
tail(blkpp,20)
write.table(blkpp, '~/transfer/contig_try_separate_subgenomes.txt', quote=F, sep='\t', row.names=F, col.names=T)
blkpp[blkpp$genome2=='tdacs1',]
dput(blkpp[blkpp$genome2=='tdacs1',])
dput(blkpp[blkpp$genome2=='blagur',])
dput(blkpp[blkpp$genome2=='vcuspi',])
dput(blkpp[blkpp$genome2=='sscopa',])
blkpp[blkpp$genome2=='sscopa',]
dput(blkpp[blkpp$genome2=='sscopa',])
dput(blkpp[blkpp$genome2=='agerar',])

                                     
