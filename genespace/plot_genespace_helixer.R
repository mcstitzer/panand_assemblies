library(GENESPACE)
library(cowplot)
theme_set(theme_cowplot())

setwd('genespace_noog')

genomeRepo <- "/workdir/mcs368/genespace/genomeRepo"
wd <- "/workdir/mcs368/genespace/"
path2mcscanx <- "/programs/MCScanX" ## directory that contains MCScanX_h


asize=read.table('../panand_sp_ploidy.txt')
#h=read.table('../general_summaries/panand_assembly_sizes.txt', header=T, sep='\t')
h=read.table('panand_assembly_sizes.txt', header=T, sep='\t')

#h$ploidy=asize$V3[match(h$V2, asize$V2)]*ifelse(h$haploid, 2,1) ## need to adjust for haploid assemblies!
asize$ploidy=asize$V3*ifelse(asize$V2%in%h$V2[h$haploid], 1,2)
asize$ploidy[asize$V2%in%c('tdacn1','tdacs1')]= 2

### even though step 8 is not finishing, enough is here to make riparian plots!!!!!!

  wd <- "/workdir/mcs368/genespace_noog" # no outgroup
  
  gpar <- init_genespace(
  wd = wd,
  path2mcscanx = path2mcscanx,
  genomeIDs=c('pvagin', asize$V2 ),
    ploidy= c(1, asize$ploidy),
#    outgroup='pvagin',
     dotplots = "never",
     orthofinderInBlk=F, ## since orthofinder gets the species tree wrong, don't make hogs!!
    nCores=64)
gsParam=gpar
                                    
ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))


muted_colors <- c("#b34064", "#459abf", "#68b488", "#b3ac40", "#8d4cba", 
                  "#bf9140", "#ae459a", "#99aabf", "#409f90", "#405973")
cust_colors <- function(n = 10){
  cols <- muted_colors
  pal <- colorRampPalette(cols)
  return(pal(n))
}

##if need to invert specific chr
invchr <- data.frame(
  genome = c("avirgi", "avirgi", "avirgi", "avirgi", rep(c('tdacs1', 'tdacn1'),each=6), 'znicar', 'zTIL25', 'zB73v5', 'znicar', 'zTIL25', 'zB73v5'), 
  chr = c('chr5','chr1','chr10','chr7', rep(paste0('chr', c(10,3,2,6,9,14)),2), 'chr1', 'chr1', 'chr1', 'chr6', 'chr6', 'chr6'))


### ready to plot!!!!




pdf(paste0('rip_plot.', Sys.Date(), '.pdf'), 10,8)

## pick four examples for each level (since i only have 4 hexaploids)
## is there a way i could trick the haploid assemblies to be allelic???? hmmmm

diploidsample=c('ppanic', 'crefra', 'avirgi', 'smicro')
tetrasample=c('snutan', 'hconto', 'achine', 'etrips')
paleosample=c('tdacs1', 'tdacn1', 'zluxur', 'zTIL01')

diploid=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=c('pvagin', diploidsample),#all$V2[all$ploidy=='Diploid']),'telega',
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=5e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr, xlabel='',
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zdmomo', 'zdgigi','tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes, palette=cust_colors
                    )
tetraploid=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=c('pvagin', tetrasample), #all$V2[all$ploidy=='Tetraploid']), #'atenui', 'rrottb',
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
paleotetraploid=plot_riparian(gsParam=gsParam, refGenome='pvagin', forceRecalcBlocks=F, genomeIDs=c('pvagin', paleosample),#'tdacn1', 'tdacs1', 'zdgigi', 'zdmomo', 'znicar', 'zluxur','zmhuet', 'zTIL25', 'zTIL18', 'zTIL01', 'zTIL11', 'zmB735'), #all$V2[all$ploidy=='Paleotetraploid' & all$V2!='zluxur']),
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=15e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr, xlabel='',
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zmhuet', 'zdmomo', 'zdgigi','tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes, palette=cust_colors
                    )

dip=diploid$plotData$ggplotObj + scale_y_continuous(labels = c('P. vaginatum', 'P. paniceum', 'C. refracta', 'A. virginicum', 'S. microstachyum')) + ylab('')

hex=hexaploid$plotData$ggplotObj + scale_y_continuous(labels = c('P. vaginatum', 'A. gerardi', 'B. laguroides', 'H. compressa', 'U. digitatum')) + ylab('')
tet=tetraploid$plotData$ggplotObj + scale_y_continuous(labels = c('P. vaginatum', 'S. nutans', 'H. contortus', 'A. chinensis','E. tripsacoides')) + ylab('')
ptet=paleotetraploid$plotData$ggplotObj + scale_y_continuous(labels = c('P. vaginatum', 'T. dactyloides KS', 'T. dactyloides FL', 'Z. luxurians', 'Z. m. ssp. parviglumis TIL01')) + ylab('')# 'Z. diploperennis Gigi', 'Z. diploperennis Momo', 'Z. nicaraguensis', 'Z. mays ssp. huehuetengensis', 'Z. mays ssp. mexicana TIL25', 'Z. mays ssp. mexicana TIL18', 'Z. mays ssp. parviglumis TIL01', 'Z. mays ssp. parviglumis TIL11', 'Z. mays ssp. mays B73v5')) + ylab('')
dev.off()

plot_data <- ggplot_build(diploid$plotData$ggplotObj)$data[[1]]
# Find the maximum chromosome length
max_chr_length <- max(plot_data$x, na.rm = TRUE)


dipwidth=max(ggplot_build(diploid$plotData$ggplotObj)$data[[1]]$x, na.rm=T)-min(ggplot_build(diploid$plotData$ggplotObj)$data[[1]]$x, na.rm=T)
tetwidth=max(ggplot_build(tetraploid$plotData$ggplotObj)$data[[1]]$x, na.rm=T)-min(ggplot_build(tetraploid$plotData$ggplotObj)$data[[1]]$x, na.rm=T)
hexwidth=max(ggplot_build(hexaploid$plotData$ggplotObj)$data[[1]]$x, na.rm=T)-min(ggplot_build(hexaploid$plotData$ggplotObj)$data[[1]]$x, na.rm=T)
paleowidth=max(ggplot_build(paleotetraploid$plotData$ggplotObj)$data[[1]]$x, na.rm=T)-min(ggplot_build(paleotetraploid$plotData$ggplotObj)$data[[1]]$x, na.rm=T)

pdf(paste0('rip_plot.scaled.', Sys.Date(), '.pdf'), 30,4)

plot_grid(diploid$plotData$ggplotObj, tetraploid$plotData$ggplotObj, hexaploid$plotData$ggplotObj, paleotetraploid$plotData$ggplotObj, ncol=4, rel_widths=c(dipwidth, tetwidth, hexwidth, paleowidth))
plot_grid(diploid$plotData$ggplotObj, tetraploid$plotData$ggplotObj, hexaploid$plotData$ggplotObj, paleotetraploid$plotData$ggplotObj, ncol=4, rel_widths=c(dipwidth*2, tetwidth, hexwidth, paleowidth))
margin=-80
riparians=plot_grid(dip+theme(axis.text = element_text(color = "gray30"), axis.line = element_blank(), axis.text.y=element_text(hjust=0, vjust=-1, margin = margin(r = margin-20))), tet+theme(axis.text = element_text(color = "gray30"), axis.line = element_blank(), axis.text.y=element_text(hjust=0,vjust=-1, margin = margin(r =margin))), hex+theme(axis.text = element_text(color = "gray30"), axis.line = element_blank(), axis.text.y=element_text(hjust=0,vjust=-1, margin = margin(r = margin))), ptet+theme(axis.text = element_text(color = "gray30"), axis.line = element_blank(), axis.text.y=element_text(hjust=0,vjust=-1, margin = margin(r = margin-80))), ncol=4, rel_widths=c(dipwidth, tetwidth, hexwidth, paleowidth),
  labels='AUTO')
riparians

dev.off()



### try again without paspalum????
####### for the most part this didn't work!!!


pdf(paste0('rip_plot.nopasp.', Sys.Date(), '.pdf'), 10,8)

## pick four examples for each level (since i only have 4 hexaploids)
## is there a way i could trick the haploid assemblies to be allelic???? hmmmm

diploidsample=c('ppanic', 'crefra', 'avirgi', 'smicro')
tetrasample=c('snutan', 'hconto', 'achine', 'etrips')
paleosample=c('tdacs1', 'tdacn1', 'zluxur', 'zTIL01')

onecol=plot_riparian(gsParam=gsParam, refGenome=diploidsample[1], forceRecalcBlocks=F, genomeIDs=c(diploidsample, tetrasample, c('blagur', 'agerar', 'hcompr', 'udigit'), paleosample),#all$V2[all$ploidy=='Diploid']),'telega',
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=5e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr, xlabel='',
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zdmomo', 'zdgigi','tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes, palette=cust_colors
                    )
tetraploid=plot_riparian(gsParam=gsParam, refGenome=tetrasample[1], forceRecalcBlocks=F, genomeIDs=c(tetrasample), #all$V2[all$ploidy=='Tetraploid']), #'atenui', 'rrottb',
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=5e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr, xlabel='',
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar','zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zdmomo', 'zdgigi', 'tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes, palette=cust_colors
                    )
hexaploid=plot_riparian(gsParam=gsParam, refGenome='blagur', forceRecalcBlocks=F, genomeIDs=c( 'blagur', 'agerar', 'hcompr', 'udigit'), #all$V2[all$ploidy=='Hexaploid']),
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=5e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr, xlabel='',
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zdmomo', 'zdgigi','tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes, palette=cust_colors
                    )
paleotetraploid=plot_riparian(gsParam=gsParam, refGenome=paleosample[1], forceRecalcBlocks=F, genomeIDs=c( paleosample),#'tdacn1', 'tdacs1', 'zdgigi', 'zdmomo', 'znicar', 'zluxur','zmhuet', 'zTIL25', 'zTIL18', 'zTIL01', 'zTIL11', 'zmB735'), #all$V2[all$ploidy=='Paleotetraploid' & all$V2!='zluxur']),
                     useOrder=F, ## keep chr position info there!!!
                     minChrLen2plot=15e6, ## since we're using chr size, we're only doing 10 Mb scafs
                     invertTheseChrs = invchr, xlabel='',
                     chrLabFontSize = 7, labelTheseGenomes = c('zB73v5', 'zTIL25', 'znicar', 'zTIL11', 'zTIL01', 'zTIL18', 'znicar', 'zmhuet', 'zdmomo', 'zdgigi','tdacs1', 'tdacn1','avirgi','sbicol','pvagin'),
                     braidAlpha = .75, chrFill = "lightgrey", addThemes = ggthemes, palette=cust_colors
                    )

dip=diploid$plotData$ggplotObj + scale_y_continuous(labels = c('P. paniceum', 'C. refracta', 'A. virginicum', 'S. microstachyum')) + ylab('')

hex=hexaploid$plotData$ggplotObj + scale_y_continuous(labels = c('A. gerardi', 'B. laguroides', 'H. compressa', 'U. digitatum')) + ylab('')
tet=tetraploid$plotData$ggplotObj + scale_y_continuous(labels = c('S. nutans', 'H. contortus', 'A. chinensis','E. tripsacoides')) + ylab('')
ptet=paleotetraploid$plotData$ggplotObj + scale_y_continuous(labels = c( 'T. dactyloides KS', 'T. dactyloides FL', 'Z. luxurians', 'Z. m. ssp. parviglumis TIL01')) + ylab('')# 'Z. diploperennis Gigi', 'Z. diploperennis Momo', 'Z. nicaraguensis', 'Z. mays ssp. huehuetengensis', 'Z. mays ssp. mexicana TIL25', 'Z. mays ssp. mexicana TIL18', 'Z. mays ssp. parviglumis TIL01', 'Z. mays ssp. parviglumis TIL11', 'Z. mays ssp. mays B73v5')) + ylab('')
dev.off()



pdf(paste0('rip_plot.onecol.nopasp.', Sys.Date(), '.pdf'), 6,10)
onecol
onecol$plotData$ggplotObj + scale_y_continuous(labels=c('P. paniceum', 'C. refracta', 'A. virginicum', 'S. microstachyum','S. nutans', 'H. contortus', 'A. chinensis','E. tripsacoides', 'A. gerardi', 'B. laguroides', 'H. compressa', 'U. digitatum', 'T. dactyloides KS', 'T. dactyloides FL', 'Z. luxurians', 'Z. m. ssp. parviglumis TIL01'))

onecol$plotData$ggplotObj + scale_y_continuous(labels=c('P. paniceum', 'C. refracta', 'A. virginicum', 'S. microstachyum','S. nutans', 'H. contortus', 'A. chinensis','E. tripsacoides', 'A. gerardi', 'B. laguroides', 'H. compressa', 'U. digitatum', 'T. dactyloides KS', 'T. dactyloides FL', 'Z. luxurians', 'Z. m. ssp. parviglumis TIL01')) +theme(axis.text = element_text(color = "gray30"), axis.line = element_blank(), axis.text.y=element_text(hjust=0, vjust=-1, margin = margin(r = margin-20)))
dev.off()

plot_data <- ggplot_build(diploid$plotData$ggplotObj)$data[[1]]
# Find the maximum chromosome length
max_chr_length <- max(plot_data$x, na.rm = TRUE)


dipwidth=max(ggplot_build(diploid$plotData$ggplotObj)$data[[1]]$x, na.rm=T)-min(ggplot_build(diploid$plotData$ggplotObj)$data[[1]]$x, na.rm=T)
tetwidth=max(ggplot_build(tetraploid$plotData$ggplotObj)$data[[1]]$x, na.rm=T)-min(ggplot_build(tetraploid$plotData$ggplotObj)$data[[1]]$x, na.rm=T)
hexwidth=max(ggplot_build(hexaploid$plotData$ggplotObj)$data[[1]]$x, na.rm=T)-min(ggplot_build(hexaploid$plotData$ggplotObj)$data[[1]]$x, na.rm=T)
paleowidth=max(ggplot_build(paleotetraploid$plotData$ggplotObj)$data[[1]]$x, na.rm=T)-min(ggplot_build(paleotetraploid$plotData$ggplotObj)$data[[1]]$x, na.rm=T)

pdf(paste0('rip_plot.scaled.nopasp.', Sys.Date(), '.pdf'), 30,4)

plot_grid(diploid$plotData$ggplotObj, tetraploid$plotData$ggplotObj, hexaploid$plotData$ggplotObj, paleotetraploid$plotData$ggplotObj, ncol=4, rel_widths=c(dipwidth, tetwidth, hexwidth, paleowidth))
plot_grid(diploid$plotData$ggplotObj, tetraploid$plotData$ggplotObj, hexaploid$plotData$ggplotObj, paleotetraploid$plotData$ggplotObj, ncol=4, rel_widths=c(dipwidth*2, tetwidth, hexwidth, paleowidth))
margin=-80
plot_grid(dip+theme(axis.text = element_text(color = "gray30"), axis.line = element_blank(), axis.text.y=element_text(hjust=0, vjust=-1, margin = margin(r = margin-20))), tet+theme(axis.text = element_text(color = "gray30"), axis.line = element_blank(), axis.text.y=element_text(hjust=0,vjust=-1, margin = margin(r =margin))), hex+theme(axis.text = element_text(color = "gray30"), axis.line = element_blank(), axis.text.y=element_text(hjust=0,vjust=-1, margin = margin(r = margin))), ptet+theme(axis.text = element_text(color = "gray30"), axis.line = element_blank(), axis.text.y=element_text(hjust=0,vjust=-1, margin = margin(r = margin-80))), ncol=4, rel_widths=c(dipwidth, tetwidth, hexwidth, paleowidth))


dev.off()













