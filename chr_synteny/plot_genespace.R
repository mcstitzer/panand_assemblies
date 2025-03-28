## need /programs/R-4.3.2/bin/R
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

## for local!
gpar <- init_genespace(
  wd = "/home/mcs368/GeneSpace_run6")
## this will load in gsParam object
load('/home/mcs368/GeneSpace_run6/results/gsParams.rda')

## omfg i have to change all these paths because they're fucking hard coded. oh and these are all considered haploid wtf i hate this
gsParam$paths=lapply(gsParam$paths, function(x) gsub('/work/mash-covid/genespace/', '/home/mcs368/', x))
gsParam$synteny$combBed=gsub('/work/mash-covid/genespace/', '/home/mcs368/', gsParam$synteny$combBed)
gsParam$synteny$SpeciesIDs=gsub('/work/mash-covid/genespace/', '/home/mcs368/', gsParam$synteny$SpeciesIDs)
gsParam$synteny$SequenceIDs=gsub('/work/mash-covid/genespace/', '/home/mcs368/', gsParam$synteny$SequenceIDs)
gsParam$synteny$ogs=gsub('/work/mash-covid/genespace/', '/home/mcs368/', gsParam$synteny$ogs)
gsParam$synteny$blast$queryBlast=sapply(gsParam$synteny$blast$queryBlast, function(x) gsub('/work/mash-covid/genespace/', '/home/mcs368/', x))
gsParam$synteny$blast$targetBlast=sapply(gsParam$synteny$blast$targetBlast, function(x) gsub('/work/mash-covid/genespace/', '/home/mcs368/', x))
gsParam$synteny$blast$allBlast=sapply(gsParam$synteny$blast$allBlast, function(x) gsub('/work/mash-covid/genespace/', '/home/mcs368/', x))
gsParam$synteny$blast$synHits=sapply(gsParam$synteny$blast$synHits, function(x) gsub('/work/mash-covid/genespace/', '/home/mcs368/', x))
                       

## select a few for riparian

## order bottom to top!!!
#rip=c('zB73v5', 'zluxur', 'tdacs1', 'tdacn1', 'udigit', 'vcuspi', 'hcompr', 'etrips', 'avirgi', 'achine', 'agerar', 'hconto', 'ppanic', 'sbicol', 'snutan', 'pvagin')
rip=c('pvagin', 'sbicol', 'ppanic', 'hconto', 'avirgi', 'etrips', 'udigit', 'tdacn1', 'tdacs1', 'znicar', 'zTIL25', 'zB73v5')
simpledips=c('pvagin', 'sbicol', 'avirgi', 'tdacn1', 'tdacs1', 'znicar', 'zTIL25', 'zB73v5')
otherpoly=c('pvagin', 'sbicol',  'hconto', 'hcompr','udigit')
                                      
ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))


muted_colors <- c("#b34064", "#459abf", "#68b488", "#b3ac40", "#8d4cba", 
                  "#bf9140", "#ae459a", "#99aabf", "#409f90", "#405973")
cust_colors <- function(n = 10){
  cols <- muted_colors
  pal <- colorRampPalette(cols)
  return(pal(n))
}


                                     
pdf('~/transfer/rip_plot.pdf',10,6)
um=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=rip,
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=10e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes, palette=cust_colors
                    )

##if need to invert specific chr
invchr <- data.frame(
  genome = c("avirgi", "avirgi", "avirgi", "avirgi", rep(c('tdacs1', 'tdacn1'),each=6), 'znicar', 'zTIL25', 'zB73v5', 'znicar', 'zTIL25', 'zB73v5'), 
  chr = c('chr5','chr1','chr10','chr7', rep(paste0('chr', c(10,3,2,6,9,14)),2), 'chr1', 'chr1', 'chr1', 'chr6', 'chr6', 'chr6'))
excuse=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=simpledips,
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=10e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr,
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes, palette=cust_colors
                    )
me=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=rip,
                     useOrder=T, ## keep chr position info there!!!
                     minChrLen2plot=100, ## since we're using chr size, we're only doing 10 Mb scafs
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes, palette=cust_colors
                    )
## wtf does this store these these stupid way jsut let me plot

invchr <- data.frame(
  genome = c("avirgi", "avirgi", "avirgi", "avirgi", rep(c('tdacs1', 'tdacn1'),each=6), 'zTIL11', 'zTIL01', 'zTIL18', 'zdmomo', 'zdgigi','znicar', 'zmhuet', 'zTIL25', 'zB73v5', 'znicar', 'zmhuet', 'zTIL25', 'zB73v5','zTIL11', 'zTIL01', 'zTIL18', 'zdmomo', 'zdgigi'), 
  chr = c('chr5','chr1','chr10','chr7', rep(paste0('chr', c(10,3,2,6,9,14)),2), 'chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr6', 'chr6', 'chr6', 'chr6', 'chr6', 'chr6', 'chr6', 'chr6', 'chr6'))
excuse=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=otherpoly,
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=10e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr,
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zdmomo', 'zdgigi', 'tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes, palette=cust_colors
                    )


diploid=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=c('pvagin', 'cserru', 'irugos', 'sbicol', 'ppanic', 'ttrian', 'crefra', 'avirgi', 'smicro', 'rtuber'),#all$V2[all$ploidy=='Diploid']),'telega',
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=5e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr, xlabel='',
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zdmomo', 'zdgigi','tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes, palette=cust_colors
                    )
tetraploid=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=c('pvagin', 'snutan', 'hconto', 'ccitra', 'achine', 'sscopa', 'etrips',  'vcuspi'), #all$V2[all$ploidy=='Tetraploid']), #'atenui', 'rrottb',
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=5e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr, xlabel='',
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar','zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zdmomo', 'zdgigi', 'tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes, palette=cust_colors
                    )
hexaploid=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=c('pvagin', 'blagur', 'agerar', 'hcompr', 'udigit'), #all$V2[all$ploidy=='Hexaploid']),
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=5e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr, xlabel='',
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zdmomo', 'zdgigi','tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes, palette=cust_colors
                    )
paleotetraploid=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=c('pvagin', 'tdacn1', 'tdacs1', 'zdgigi', 'zdmomo', 'znicar', 'zmhuet', 'zTIL25', 'zTIL18', 'zTIL01', 'zTIL11',
 'zB73v5'), #all$V2[all$ploidy=='Paleotetraploid' & all$V2!='zluxur']),
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=15e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr, xlabel='',
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zmhuet', 'zdmomo', 'zdgigi','tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes, palette=cust_colors
                    )
gerardi=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=c('pvagin', 'avirgi', 'achine', 'agerar', 'sscopa', 'smicro'),
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=5e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr, xlabel='',
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zdmomo', 'zdgigi','tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes, palette=cust_colors
                    )     
scoparium=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=c('pvagin', 'avirgi', 'sscopa', 'smicro'),
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=5e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr, xlabel='',
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zdmomo', 'zdgigi','tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes, palette=cust_colors
                    )   
##for real, add scalePlotWidth, to make so diploids, tetraploids, paleotetraploids all same x scale even when chr are bigger!



dev.off()

### ugh this didn't work - maybe i have to remove maize??
roi_sg1=data.frame(genome=c(rep('tdacn1',10), rep('tdacs1',10),  rep('pvagin',10)),
                   chr=c(paste0('chr',c(1,3,2,7,14,8,11,18,12,4)),paste0('chr',c(1,3,2,7,14,8,11,18,12,4)),
                         paste0('Chr0',1:9),'Chr10'))

sg1=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=c('pvagin', 'tdacn1', 'tdacs1', 
 'zB73v5'), #all$V2[all$ploidy=='Paleotetraploid' & all$V2!='zluxur']),
                     useOrder=F, ## keep chr position info there!!!
                     highlightBed = roi_sg1, ## only plot chr in roi!! 
                    backgroundColor = NULL, 
                     minChrLen2plot=15e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr, xlabel='',
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zmhuet', 'zdmomo', 'zdgigi','tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes, palette=cust_colors
                    )
roi_sg2=data.frame(genome=c(rep('tdacn1',10), rep('tdacs1',10), rep('zB73v5', 10), rep('pvagin',10)),
                   chr=c(paste0('chr',c(10,5,6,9,5,15,17,2,16,13)),paste0('chr',c(10,5,6,9,5,15,17,2,16,13)),
                         paste0('chr',1:10), paste0('Chr0',1:9),'Chr10'))

sg2=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=c('pvagin', 'tdacn1', 'tdacs1', 
 'zB73v5'), #all$V2[all$ploidy=='Paleotetraploid' & all$V2!='zluxur']),
                     useOrder=F, ## keep chr position info there!!!
                     highlightBed = roi_sg2, 
                    backgroundColor = NULL, 
                     minChrLen2plot=15e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr, xlabel='',
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zmhuet', 'zdmomo', 'zdgigi','tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes, palette=cust_colors
                    )


## above is not working :( i think b/c maize has homology to the absent trip crh
                                     ## plot separately, backgroundColor=NULL might give me all maize chr, highlihgting only each subgenome link.

roi_sg1simple=data.frame(genome=c(rep('tdacs1',10)),
                   chr=c(paste0('chr',c(1,3,2,7,14,8,11,18,12,4))), color='red')

sg1simplea=plot_riparian(gsParam=gsParam, refGenome='tdacs1', forceRecalcBlocks=F, genomeIDs=c('tdacs1', 
 'zB73v5'), #all$V2[all$ploidy=='Paleotetraploid' & all$V2!='zluxur']),
                     useOrder=F, ## keep chr position info there!!!
                     highlightBed = roi_sg1simple, ## only plot chr in roi!! 
                    backgroundColor = NULL, 
                        customRefChrOrder=c(paste0('chr',c(1,3,2,7,14,8,11,18,12,4)), paste0('chr',c(10,5,6,9,15,17,16,13))),
                     minChrLen2plot=15e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr, xlabel='',
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zmhuet', 'zdmomo', 'zdgigi','tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes, palette=cust_colors
                    )
roi_sg2simple=data.frame(genome=c(rep('tdacs1',8)),
                   chr=c(paste0('chr',c(10,5,6,9,15,17,16,13))), color='blue')

sg2simplea=plot_riparian(gsParam=gsParam, refGenome='tdacs1', forceRecalcBlocks=F, genomeIDs=c('tdacs1', 
 'zB73v5'), #all$V2[all$ploidy=='Paleotetraploid' & all$V2!='zluxur']),
                     useOrder=F, ## keep chr position info there!!!
                     highlightBed = roi_sg2simple, ## only plot chr in roi!! 
                    backgroundColor = NULL, 
                        customRefChrOrder=c(paste0('chr',c(1,3,2,7,14,8,11,18,12,4)), paste0('chr',c(10,5,6,9,15,17,16,13))),
                     minChrLen2plot=15e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr, xlabel='',
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zmhuet', 'zdmomo', 'zdgigi','tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes, palette=cust_colors
                    )

                                     
## way to add labels?? this might get me species names
p2 <- p1 + 
  scale_y_discrete(limits = c("chicky", "people", "mice", "duckBilled", "dangerNoodle"))+ ## bottom to top!
  theme(axis.text.y = element_text(size = 16))
                                     

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

maizesubgenomes=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=c('pvagin', 'tdacs1', 'zB73v5'), #all$V2[all$ploidy=='Paleotetraploid' & all$V2!='zluxur']),
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=15e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr, xlabel='',
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zmhuet', 'zdmomo', 'zdgigi','tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .5, chrFill = "lightgrey", addThemes = ggthemes,
                     highlightBed=subgenome, backgroundColor='snow2'#,
#                     customRefChrOrder=paste0('chr', c(1,17,8,5,2,14,9,10,7,12,13,3,6,16,4,18,11,15)) ## can set up to make maize go  1:10 for all our friends...
#                     customRefChrOrder=paste0('chr', c(1,10,4,3,5,2,6,7,9,14,8,15,11,17,18,12,16,13)) ## can set up to make maize go  1:10 for all our friends...
                    )
### it does not like generating this plot :D lots of invalid pointers - could get somebody who understands r to help me?
                                     
## to switch out my own taxon names, use this!! in glab
#scale_y_continuous(breaks = (glab$y1+glab$y2)/2, labels = glab$genome, expand = c(0.01, 0.01), name = NULL)
# maybe this wil update?!?!?!
hexaploid$plotData$ggplotObj + scale_y_continuous(labels = c('P. vaginatum', 'A. gerardi', 'B. laguroides', 'H. compressa', 'U. digitatum'))

pdf(paste0('~/transfer/rip_plot.', Sys.Date(), '.pdf'), 10,8)
## to switch out my own taxon names, use this!! in glab
#scale_y_continuous(breaks = (glab$y1+glab$y2)/2, labels = glab$genome, expand = c(0.01, 0.01), name = NULL)
# maybe this wil update?!?!?! - have to separately load ggplot2 lib
hex=hexaploid$plotData$ggplotObj + scale_y_continuous(labels = c('P. vaginatum', 'A. gerardi', 'B. laguroides', 'H. compressa', 'U. digitatum')) + ylab('')
tet=tetraploid$plotData$ggplotObj + scale_y_continuous(labels = c('P. vaginatum', 'S. nutans', 'H. contortus', 'C. citratus', 'A. chinensis', 'S. scoparium', 'E. tripsacoides', 'V. cuspidata')) + ylab('')
dip=diploid$plotData$ggplotObj + scale_y_continuous(labels = c('P. vaginatum', 'C. serrulatus', 'I. rugosum', 'S. bicolor', 'P. paniceum', 'T. triandra', 'C. refracta', 'A. virginicus', 'S. microstachyum', 'R. tuberosum')) + ylab('')
ptet=paleotetraploid$plotData$ggplotObj + scale_y_continuous(labels = c('P. vaginatum', 'T. dactyloides KS', 'T. dactyloides FL', 'Z. diploperennis Gigi', 'Z. diploperennis Momo', 'Z. nicaraguensis', 'Z. mays ssp. huehuetengensis', 'Z. mays ssp. mexicana TIL25', 'Z. mays ssp. mexicana TIL18', 'Z. mays ssp. parviglumis TIL01', 'Z. mays ssp. parviglumis TIL11', 'Z. mays ssp. mays B73v5')) + ylab('')
msubg=maizesubgenomes$plotData$ggplotObj + scale_y_continuous(labels = c('T. dactyloides FL', 'Z. mays ssp. mays B73v5')) + ylab('')
                                     
#plot_grid(dip, tet, ptet, hex, align='hv', labels=c('a Diploid', 'b Tetraploid', 'c Paleotetraploid', 'd Hexaploid'), ncol=2)
plot_grid(plot_grid(diploid$plotData$ggplotObj, tetraploid$plotData$ggplotObj,  hex, align='hv', labels=c('a Diploid', 'b Tetraploid', 'c  Hexaploid'), 
ncol=3, rel_widths=c(0.5,1,1)),
paleotetraploid$plotData$ggplotObj, ncol=1, rel_heights=c(1,1.5), labels=c('', 'd Paleotetraploid'))

plot_grid(plot_grid(diploid$plotData$ggplotObj, tetraploid$plotData$ggplotObj,  hex, align='hv', labels=c('a Diploid', 'b Tetraploid', 'c  Hexaploid'), 
ncol=3, rel_widths=c(1,1,1)),
paleotetraploid$plotData$ggplotObj, ncol=1, rel_heights=c(1,1.5), labels=c('', 'd Paleotetraploid'))

plot_grid(plot_grid(diploid$plotData$ggplotObj, tetraploid$plotData$ggplotObj,  hex, align='hv', labels=c('a Diploid', 'b Tetraploid', 'c  Hexaploid'), 
ncol=3, rel_widths=c(1,1,1)),
paleotetraploid$plotData$ggplotObj, maizesubgenomes$plotData$ggplotObj, ncol=1, rel_heights=c(1,1.5,0.5), labels=c('', 'd Paleotetraploid', 'e Maize Subgenomes'))

                                     
dev.off()


## put knobs on chromsomes
ptk=paleotetraploid$plotData$sourceData$chromosomes

## get from te plotting
# tand=genomecount[genomecount$sup=='TandemRepeat' & genomecount$genome %in% c('pvagin', 'tdacn1', 'tdacs1', 'zdgigi', 'zdmomo', 'znicar', 'zmhuet', 'zTIL25', 'zTIL18', 'zTIL01', 'zTIL11',
 'zB73v5'),c(1,2,3,10,11,22,23)]
# write.table(tand, '~/transfer/tandem_repeats_panand_rm_from_genomecount.txt', row.names=F, col.names=T, sep='\t', quote=F)
tand=read.table('~/transfer/tandem_repeats_panand_rm_from_genomecount.txt', header=T)    
tand$rl=str_split_fixed(tand$Name,'_',4)[,3]
tand$class=NA
tand$class[tand$rl%in%c('180bp', '179bp')]='knob180'
tand$class[tand$rl%in%c('359.5bp', '359bp')]='knobTr1'
tand$class[tand$rl%in%c('156bp','155bp')]='centromere'       
tand=tand[!is.na(tand$class),]
gr180=unlist(reduce(split(GRanges(seqnames=tand$chr[tand$class=='knob180'], IRanges(start=tand$start[tand$class=='knob180'], end=tand$end[tand$class=='knob180'])), ~genome=tand$genome[tand$class=='knob180'])))                                    
grtr1=unlist(reduce(split(GRanges(seqnames=tand$chr[tand$class=='knobTr1'], IRanges(start=tand$start[tand$class=='knobTr1'], end=tand$end[tand$class=='knobTr1'])), ~tand$genome[tand$class=='knobTr1'])))                                    
grcent=unlist(reduce(split(GRanges(seqnames=tand$chr[tand$class=='centromere'], IRanges(start=tand$start[tand$class=='centromere'], end=tand$end[tand$class=='centromere'])), ~tand$genome[tand$class=='centromere'])))                                    
ptand=c(gr180[width(gr180)>1e4],grtr1[width(grtr1)>1e4],grcent[width(grcent)>1e4])
ptand$class=c(rep('knob180', sum(width(gr180)>1e4)), rep('knobTr1', sum(width(grtr1)>1e4)),rep('centromere', sum(width(grcent)>1e4)))
ptand$genome=names(ptand)
ptand=data.frame(ptand)
## put knobs in middle of y's for this genome
ptand$x1=ptk$x1[match(paste(ptand$seqnames, ptand$genome),paste(ptk$chr,ptk$genome))]
ptand$x2=ptk$x2[match(paste(ptand$seqnames, ptand$genome),paste(ptk$chr,ptk$genome))]
ptand$y1=ptk$y1[match(paste(ptand$seqnames, ptand$genome),paste(ptk$chr,ptk$genome))]
ptand$y2=ptk$y2[match(paste(ptand$seqnames, ptand$genome),paste(ptk$chr,ptk$genome))]
colnames(ptand)[1]='chr'                                      
ptand$xadjstart=ptand$x1+ptand$start
ptand$xadjend=ptand$x1+ptand$end
ptand$yadj=(ptand$y1+ptand$y2)/2
knobcolors=c('#86B049', '#9FE7F5', '#525B88')
names(knobcolors)=c('knob180', 'knobTr1', 'centromere')
                                      

pdf(paste0('~/transfer/riparian_fig2.', Sys.Date(), '.pdf'), 10,10)

plot_grid(plot_grid(diploid$plotData$ggplotObj+ xlab(''), tetraploid$plotData$ggplotObj+ xlab(''),  hex+ xlab(''), align='hv', labels=c('a Diploid', 'b Tetraploid', 'c  Hexaploid'), 
ncol=3, rel_widths=c(1,1,1)),
          paleotetraploid$plotData$ggplotObj + xlab('') + geom_rect(data=ptand[ptand$class!='knobTr1',], aes(xmin=xadjstart,xmax=xadjend,ymin=yadj-0.05,ymax=yadj+0.05, color=class), alpha=0.5) + scale_color_manual(values=knobcolors)+ theme(legend.position='none'),
          maizesubgenomes$plotData$ggplotObj + xlab(''), ncol=1, rel_heights=c(1,1.5,0.7), labels=c('', 'd Paleotetraploid', 'e Maize Subgenomes'))

#paleotetraploid$plotData$ggplotObj + xlab('') + geom_point(data=ptand, aes(x=xadj,y=yadj,color=class), pch='|') + scale_color_manual(values=knobcolors)                                  
#paleotetraploid$plotData$ggplotObj + xlab('') + geom_rect(data=ptand[ptand$class!='knobTr1',], aes(xmin=xadjstart,xmax=xadjend,ymin=yadj-0.05,ymax=yadj+0.05, color=class), alpha=0.7) + scale_color_manual(values=knobcolors)                           

dev.off()



## just avirgi to show paspalum identity
avirgi=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=c('pvagin',  'avirgi'),#all$V2[all$ploidy=='Diploid']),'telega',
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=5e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr, xlabel='',
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zdmomo', 'zdgigi','tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes, palette=cust_colors
                    )
avirgi
udigit=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=c('pvagin',  'udigit'),#all$V2[all$ploidy=='Diploid']),'telega',
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=5e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr, xlabel='',
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zdmomo', 'zdgigi','tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes, palette=cust_colors
                    )
udigit


blagur=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=c('pvagin',  'blagur'),#all$V2[all$ploidy=='Diploid']),'telega',
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=5e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr, xlabel='',
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zdmomo', 'zdgigi','tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes, palette=cust_colors
                    )
blagur
agerar=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=c('pvagin',  'agerar'),#all$V2[all$ploidy=='Diploid']),'telega',
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=5e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr, xlabel='',
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zdmomo', 'zdgigi','tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes, palette=cust_colors
                    )
agerar
                                      
## compare to dotplot
av=read.table('~/Downloads/avirgi-Pv-2', header=T)
av=av[av$gene!='interanchor',]
## reduce to 10 block lengths
av=av %>% group_by(blockIndex) %>% mutate(blockLength=n())    
av=av%>% group_by(queryChr)%>% mutate(freqStrand=names(which.max(table(strand))), maxChr=max(queryStart), freqRef=names(which.max(table(refChr))))


av=av[av$blockLength>10,]
av$refChr=factor(av$refChr, levels=c(paste0('Chr0', 1:9), 'Chr10'))
names(muted_colors)=c(paste0('Chr0', 1:9), 'Chr10')

## reverse strand if there is a 
                                      
library(dplyr)
av=av %>% arrange(freqRef, referenceStart, queryStart)
av$queryChr=factor(av$queryChr, levels=rev(av$queryChr[!duplicated(av$queryChr)]))
av$revQueryStart=av$queryStart
av$revQueryStart[av$freqStrand=='-']=abs(av$queryStart-av$maxChr)[av$freqStrand=='-']
theme_set(theme_cowplot())
ggplot(av[av$queryChr%in%paste0('chr', 1:10) & av$refChr%in%names(muted_colors),], aes(x=referenceStart/1e6, y=revQueryStart/1e6, color=refChr)) + geom_point() + facet_grid(queryChr~refChr, scales='free', space='free')  + scale_color_manual(values=muted_colors)+theme(legend.position='none')+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


## compare to dotplot
ud=read.table('~/Downloads/udigit-Pv-6', header=T)
ud=ud[ud$gene!='interanchor',]
## reduce to 10 block lengths
ud=ud %>% group_by(blockIndex) %>% mutate(blockLength=n())    
ud=ud%>% group_by(queryChr)%>% mutate(freqStrand=names(which.max(table(strand))), maxChr=max(queryStart), freqRef=names(which.max(table(refChr))))


ud=ud[ud$blockLength>10,]
ud$refChr=factor(ud$refChr, levels=c(paste0('Chr0', 1:9), 'Chr10'))
names(muted_colors)=c(paste0('Chr0', 1:9), 'Chr10')

## reverse strand if there is a 
                                      
library(dplyr)
ud=ud %>% arrange(freqRef, referenceStart, queryStart)
ud$queryChr=factor(ud$queryChr, levels=rev(ud$queryChr[!duplicated(ud$queryChr)]))
ud$revQueryStart=ud$queryStart
ud$revQueryStart[ud$freqStrand=='-']=abs(ud$queryStart-ud$maxChr)[ud$freqStrand=='-']
theme_set(theme_cowplot())
ggplot(ud[ud$refChr%in%names(muted_colors),], aes(x=referenceStart/1e6, y=revQueryStart/1e6, color=refChr)) + geom_point() + facet_grid(queryChr~refChr, scales='free', space='free')  + scale_color_manual(values=muted_colors)+theme(legend.position='none')+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


                                      
pdf('evoday_dotplots.pdf',14,12)
ggplot(av[av$queryChr%in%paste0('chr', 1:10) & av$refChr%in%names(muted_colors),], aes(x=referenceStart/1e6, y=revQueryStart/1e6, color=refChr)) + geom_point() + facet_grid(queryChr~refChr, scales='free', space='free')  + scale_color_manual(values=muted_colors)+theme(legend.position='none')+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(ud[ud$refChr%in%names(muted_colors),], aes(x=referenceStart/1e6, y=revQueryStart/1e6, color=refChr)) + geom_point() + facet_grid(queryChr~refChr, scales='free', space='free')  + scale_color_manual(values=muted_colors)+theme(legend.position='none')+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
                                  

## add lines
ggplot(av[av$queryChr%in%paste0('chr', 1:10) & av$refChr%in%names(muted_colors),], aes(x=referenceStart/1e6, y=revQueryStart/1e6, color=refChr)) + geom_point() + facet_grid(queryChr~refChr, scales='free', space='free')  + scale_color_manual(values=muted_colors)+theme(legend.position='none')+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + geom_hline(aes(yintercept=maxChr/1e6), lty='dashed', color='gray')
ggplot(ud[ud$refChr%in%names(muted_colors),], aes(x=referenceStart/1e6, y=revQueryStart/1e6, color=refChr)) + geom_point() + facet_grid(queryChr~refChr, scales='free', space='free')  + scale_color_manual(values=muted_colors)+theme(legend.position='none')+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + geom_hline(aes(yintercept=maxChr/1e6), lty='dashed', color='gray')

                                      
## tetra
process_anchors_to_dotplot('~/Downloads/vcuspi-Pv-4', minBlock=20, title='vcuspi')
process_anchors_to_dotplot('~/Downloads/etrips-Pv-4', minBlock=20, title='etrips')
process_anchors_to_dotplot('~/Downloads/hconto-Pv-4', minBlock=20, title='hconto')
process_anchors_to_dotplot('~/Downloads/sscopa-Pv-4', minBlock=20, title='sscopa')
process_anchors_to_dotplot('~/Downloads/achine-Pv-4', minBlock=20, title='achine')
process_anchors_to_dotplot('~/Downloads/ccitra-Pv-4', minBlock=20, title='ccitra')
process_anchors_to_dotplot('~/Downloads/snutan-Pv-4', minBlock=20, title='snutan')

##hex                                      
process_anchors_to_dotplot('~/Downloads/udigit-Pv-6', minBlock=50, title='udigit')
process_anchors_to_dotplot('~/Downloads/blagur-Pv-6', minBlock=50, title='blagur')
process_anchors_to_dotplot('~/Downloads/hcompr-Pv-6', minBlock=50, title='hcompr')
process_anchors_to_dotplot('~/Downloads/agerar-Pv-6', minBlock=50, title='agerar')                                      
dev.off()            

                                      
pdf('evoday_alluvial.pdf',14,4)
avirgi
udigit
blagur
agerar
dev.off()

pdf('evoday_alluvial_groups.pdf',8,8)
diploid
tetraploid
hexaploid
paleotetraploid
dev.off()
pdf('evoday_alluvial_groups.zt.pdf',16,8)
paleotetraploid
dev.off()
                                      
## generalize!
  # Set the theme and color palette
  theme_set(theme_cowplot())
 # Define color palette
  muted_colors <- c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf')
  names(muted_colors) <- c(paste0('Chr0', 1:9), 'Chr10')
  
process_anchors_to_dotplot <- function(filepath, color_palette=muted_colors, minBlock=10, title='') {
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

  # Reverse strand calculations
  data <- data %>%
    arrange(freqRef, referenceStart, queryStart)
  data$queryChr <- factor(data$queryChr, levels = rev(data$queryChr[!duplicated(data$queryChr)]))
  data$revQueryStart <- data$queryStart
  data$revQueryStart[data$freqStrand == '-'] <- abs(data$queryStart - data$maxChr)[data$freqStrand == '-']


  # Create the plot
  ggplot(data[ data$refChr %in% names(color_palette), ],
         aes(x = referenceStart / 1e6, y = revQueryStart / 1e6, color = refChr)) +
    geom_point() +
    facet_grid(queryChr ~ refChr, scales = 'free', space = 'free') +
    scale_color_manual(values = color_palette) +
    theme(legend.position = 'none') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    geom_hline(aes(yintercept=maxChr/1e6), lty='dashed', color='gray') +
    ggtitle(title)
}


# Example of usage
process_anchors_to_dotplot('~/Downloads/agerar-Pv-6', minBlock=100)

## count rearrangments
 diploids=c('cserru', 'irugos', 'sbicol', 'ppanic', 'ttrian', 'crefra', 'avirgi', 'smicro', 'rtuber')
 tetraploids=c('snutan', 'hconto', 'ccitra', 'achine', 'sscopa', 'etrips',  'vcuspi')
hexaploids=c('udigit', 'agerar', 'hcompr', 'blagur')

countRearrangements <- function(filepath, color_palette=muted_colors, minBlock=10) {
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
  print(countRearrangements(Sys.glob(paste0('~/Downloads/', i, '-Pv-*'))))
  }

 for(i in tetraploids){
  print(countRearrangements(Sys.glob(paste0('~/Downloads/', i, '-Pv-*'))))
  }
        for(i in hexaploids){
  print(countRearrangements(Sys.glob(paste0('~/Downloads/', i, '-Pv-*')), minBlock=50))
  }                              
        for(i in c('zmB735')){
  print(countRearrangements(Sys.glob(paste0('~/Downloads/', i, '-Pv-*')), minBlock=50))
  }  
