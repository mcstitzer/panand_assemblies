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


## select a few for riparian

## order bottom to top!!!
#rip=c('zB73v5', 'zluxur', 'tdacs1', 'tdacn1', 'udigit', 'vcuspi', 'hcompr', 'etrips', 'avirgi', 'achine', 'agerar', 'hconto', 'ppanic', 'sbicol', 'snutan', 'pvagin')
rip=c('pvagin', 'sbicol', 'ppanic', 'hconto', 'avirgi', 'etrips', 'udigit', 'tdacn1', 'tdacs1', 'znicar', 'zTIL25', 'zB73v5')
simpledips=c('pvagin', 'sbicol', 'avirgi', 'tdacn1', 'tdacs1', 'znicar', 'zTIL25', 'zB73v5')
otherpoly=c('pvagin', 'sbicol',  'hconto', 'hcompr','udigit')
                                      
ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))
pdf('~/transfer/rip_plot.pdf',10,6)
um=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=rip,
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=10e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes
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
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes
                    )
me=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=rip,
                     useOrder=T, ## keep chr position info there!!!
                     minChrLen2plot=100, ## since we're using chr size, we're only doing 10 Mb scafs
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes
                    )
## wtf does this store these these stupid way jsut let me plot

invchr <- data.frame(
  genome = c("avirgi", "avirgi", "avirgi", "avirgi", rep(c('tdacs1', 'tdacn1'),each=6), 'zTIL11', 'zTIL01', 'zTIL18', 'zdmomo', 'zdgigi','znicar', 'zTIL25', 'zB73v5', 'znicar', 'zTIL25', 'zB73v5','zTIL11', 'zTIL01', 'zTIL18', 'zdmomo', 'zdgigi'), 
  chr = c('chr5','chr1','chr10','chr7', rep(paste0('chr', c(10,3,2,6,9,14)),2), 'chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr6', 'chr6', 'chr6', 'chr6', 'chr6', 'chr6', 'chr6', 'chr6'))
excuse=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=otherpoly,
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=10e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr,
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zdmomo', 'zdgigi', 'tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes
                    )


diploid=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=c('pvagin', 'cserru', 'irugos', 'sbicol', 'ppanic', 'ttrian', 'crefra', 'avirgi', 'smicro', 'rtuber'),#all$V2
[all$ploidy=='Diploid']),'telega',
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=5e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr,
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zdmomo', 'zdgigi','tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes
                    )
tetraploid=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=c('pvagin', 'snutan', 'hconto', 'ccitra', 'achine', 'sscopa', 'etrips',  'vcuspi'), #all$V2[all$ploidy=='T
etraploid']), #'atenui', 'rrottb',
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=5e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr,
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar','zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zdmomo', 'zdgigi', 'tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes
                    )
hexaploid=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=c('pvagin', 'blagur', 'agerar', 'hcompr', 'udigit'), #all$V2[all$ploidy=='Hexaploid']),
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=5e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr,
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zdmomo', 'zdgigi','tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes
                    )
paleotetraploid=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=c('pvagin', 'tdacn1', 'tdacs1', 'zdgigi', 'zdmomo', 'znicar', 'zmhuet', 'zTIL25', 'zTIL18', 'zTIL01', 'zTIL11',
 'zB73v5'), #all$V2[all$ploidy=='Paleotetraploid' & all$V2!='zluxur']),
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=15e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr,
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zmhuet', 'zdmomo', 'zdgigi','tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes
                    )
gerardi=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=c('pvagin', 'avirgi', 'achine', 'agerar', 'sscopa', 'smicro'),
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=5e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr,
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zdmomo', 'zdgigi','tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes
                    )     
scoparium=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=c('pvagin', 'avirgi', 'sscopa', 'smicro'),
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=5e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr,
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zdmomo', 'zdgigi','tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes
                    )   
##for real, add scalePlotWidth, to make so diploids, tetraploids, paleotetraploids all same x scale even when chr are bigger!

dev.off()

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
#plot_grid(dip, tet, ptet, hex, align='hv', labels=c('a Diploid', 'b Tetraploid', 'c Paleotetraploid', 'd Hexaploid'), ncol=2)
plot_grid(plot_grid(diploid$plotData$ggplotObj, tetraploid$plotData$ggplotObj,  hex, align='hv', labels=c('a Diploid', 'b Tetraploid', 'c  Hexaploid'), 
ncol=3, rel_widths=c(0.5,1,1)),
paleotetraploid$plotData$ggplotObj, ncol=1, rel_heights=c(1,1.5), labels=c('', 'd Paleotetraploid'))

plot_grid(plot_grid(diploid$plotData$ggplotObj, tetraploid$plotData$ggplotObj,  hex, align='hv', labels=c('a Diploid', 'b Tetraploid', 'c  Hexaploid'), 
ncol=3, rel_widths=c(1,1,1)),
paleotetraploid$plotData$ggplotObj, ncol=1, rel_heights=c(1,1.5), labels=c('', 'd Paleotetraploid'))

dev.off()
