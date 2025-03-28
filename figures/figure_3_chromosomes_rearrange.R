library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(rtracklayer)



### chromosome number and age rearrangments
### ../syntenic_anchors/count_rearrangements.R
fig_chrrearr
fig_agerearr


chrplots=plot_grid(fig_chrrearr + theme(legend.position='NULL'), 
                   fig_agerearr+theme(legend.position='NULL'), 
                   ncol=2, align='h', axis='tb', labels=c('G', 'H'))



## riparians from ../genespace/plot_genespace_helixer.R
## have to make them on the server, so open them
#readRDS('~/Downloads/riparians.RDS')
riparians ## these have ABCD already


## dotplots
## from ../syntenic_anchors/count_anchors_and_plot.R

## avirgi
av
## blagur
bl



leftside=plot_grid(av, chrplots, ncol=1, rel_heights=c(1,0.7), labels=c('E',''), align='v', axis='l')

bottompart=plot_grid(leftside, bl, ncol=2, labels=c('', 'F'), align='h', axis='tb')

## on mine
saveRDS(bottompart, 'fig3bottom.RDS')

### on cbsu where riparians that are 1gb are loaded
bottompart=readRDS('fig3bottom.RDS')

pdf('figure3_chromosomes-rearrange.pdf', 12,10)
plot_grid(riparians, bottompart, ncol=1, rel_heights=c(0.3,1))
dev.off()



## testing sizing
pdf('~/Downloads/test_chr.pdf',12,10)
plot_grid(NULL, bottompart, ncol=1, rel_heights=c(0.3,1))
dev.off()


