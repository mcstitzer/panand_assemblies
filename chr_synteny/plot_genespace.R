## need /programs/R-4.3.2/bin/R
library(GENESPACE)

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


## select a few for riparian

## order bottom to top!!!
#rip=c('zB73v5', 'zluxur', 'tdacs1', 'tdacn1', 'udigit', 'vcuspi', 'hcompr', 'etrips', 'avirgi', 'achine', 'agerar', 'hconto', 'ppanic', 'sbicol', 'snutan', 'pvagin')
rip=c('pvagin', 'sbicol', 'ppanic', 'hconto', 'avirgi', 'etrips', 'udigit', 'tdacn1', 'tdacs1', 'zdgigi', 'zTIL25', 'zB73v5')
simpledips=c('pvagin', 'sbicol', 'avirgi', 'tdacn1', 'tdacs1', 'zdgigi', 'zTIL25', 'zB73v5')
                                     
ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))
pdf('~/transfer/rip_plot.pdf',10,10)
um=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=rip,
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=10e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes
                    )

##if need to invert specific chr
invchr <- data.frame(
  genome = c("avirgi", "avirgi", "avirgi", "avirgi"), 
  chr = c('chr5','chr1','chr10','chr7'))
excuse=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=simpledips,
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=10e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr,
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes
                    )
me=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=rip,
                     useOrder=T, ## keep chr position info there!!!
                     minChrLen2plot=100, ## since we're using chr size, we're only doing 10 Mb scafs
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes
                    )
## wtf does this store these these stupid way jsut let me plot
dev.off()



