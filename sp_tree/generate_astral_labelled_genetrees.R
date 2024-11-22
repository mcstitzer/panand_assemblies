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



filenames=list.files('../gene_trees/trees/', pattern='RAxML_bipartitionsBranchLabels.*0')

all=read.table('../panand_sp_ploidy.txt')                                


b=lapply(c(1:length(filenames))[-40], function(i){ ## ep2 is not there

print(i)
awto=read.raxml(paste0('../gene_trees/trees/',filenames[i]))
awt=as.phylo(awto)
awt$node.label <- awto@data$bootstrap
awt$tip.label=gsub('_R_', '', awt$tip.label)
  ## add to make sure outgroup is there, and is monophyletic - otherwise skip this tree
  if(any(substr(as.phylo(awt)$tip.label,1,6) %in% c('pvagin'))){
#    if(is.monophyletic(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('pvagin')])){
  awt=root(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('pvagin')], resolve.root=T)

## now keep only six digit code
awt$tip.label=substr(awt$tip.label,1,6)

return(awt)
}
#}
                                              })
                                              
                                              
d=do.call("c",b)

write.tree(d, paste0('paspalum_anchors_aster.', Sys.Date(), '.tre'))

