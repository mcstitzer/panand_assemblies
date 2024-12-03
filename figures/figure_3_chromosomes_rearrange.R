library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(rtracklayer)



### chromosome number and age rearrangments
### ../syntenic_anchors/count_rearrangements.R
fig_chrrearr
fig_agerearr


chrplots=plot_grid(fig_chrrearr + theme(legend.position='NULL'), fig_agerearr+theme(legend.position='NULL'), ncol=2, align='v', axis='tb', labels=c('G', 'H'))



## riparians from ../genespace/plot_genespace_helixer.R
riparians ## these have ABCD already


## dotplots
## from ../syntenic_anchors/count_anchors_and_plot.R

## avirgi
av
## blagur
bl



leftside=plot_grid(av, chrplots, ncol=1, rel_heights=c(1,0.7), labels=c('E',''), align='v', axis='l')

bottompart=plot_grid(leftside, bl, ncol=2, labels=c('', 'F'), align='h', axis='tb')

plot_grid(riparians, bottompart, ncol=1, rel_heights=c(0.3,1))