library(ggtree)  ## /programs/R-4.2.1-r9/bin/R
library(ggplot2) ## either scinet or cbsu
library(cowplot)
theme_set(theme_cowplot())
library(ape)
library(treeio)
library(phytools)
library(reshape2)
library(dplyr)
library(ggridges)



filenames=list.files('../scinet_trees/', pattern='RAxML_bipartitionsBranchLabels.')

all=read.table('../panand_sp_ploidy.txt')                                


b=sapply(1:length(filenames), function(i){ ## ep2 is not there

awt=read.raxml(paste0('../scinet_trees/',filenames[i]))
awt=as.phylo(awt)
awt$tip.label=gsub('_R_', '', awt$tip.label)
  ## add to make sure outgroup is there, and is monophyletic - otherwise skip this tree
  if(any(substr(as.phylo(awt)$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate'))){
    if(is.monophyletic(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate')])){
  awt=root(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate')])

## now keep only six digit code
awt$tip.label=substr(awt$tip.label,1,6)
awt$tip.label[grepl('Pavag', awt$tip.label)]='paspal'
awt=drop.tip(awt, awt$tip.label[grepl('pprate', awt$tip.label)])
awt=drop.tip(awt, awt$tip.label[grepl('pvagin', awt$tip.label)])
## drop non=panand sp
awt=drop.tip(awt, awt$tip.label[grepl('eophiu', awt$tip.label)])
awt=drop.tip(awt, awt$tip.label[grepl('agerjg', awt$tip.label)])
return(awt)
}
}
                                              })
                                              
                                              
d=do.call("c",b)

write.tree(d, paste0('paspalum_anchors_aster.', Sys.Date(), '.tre'))

