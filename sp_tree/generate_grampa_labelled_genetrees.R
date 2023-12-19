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
library(stringr)


filenames=list.files('../scinet_trees/', pattern='RAxML_bipartitionsBranchLabels.')

all=read.table('../panand_sp_ploidy.txt')                                


b=lapply(1:length(filenames), function(i){ ## ep2 is not there

awto=read.raxml(paste0('../scinet_trees/',filenames[i]))
awt=as.phylo(awto)
awt$node.label <- awto@data$bootstrap
awt$tip.label=gsub('_R_', '', awt$tip.label)
  ## add to make sure outgroup is there, and is monophyletic - otherwise skip this tree
  if(any(substr(as.phylo(awt)$tip.label,1,5) %in% c('Pavag'))){
    if(is.monophyletic(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,5) %in% c('Pavag')])){
  awt=root(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,5) %in% c('Pavag')])

## now keep only six digit code
awt$tip.label=substr(awt$tip.label,1,6)
awt$tip.label[grepl('Pavag', awt$tip.label)]='paspal'
awt=drop.tip(awt, awt$tip.label[grepl('pprate', awt$tip.label)])
awt=drop.tip(awt, awt$tip.label[grepl('pvagin', awt$tip.label)])
## drop non=panand sp
awt=drop.tip(awt, awt$tip.label[grepl('eophiu', awt$tip.label)])
awt=drop.tip(awt, awt$tip.label[grepl('agerjg', awt$tip.label)])
awt=drop.tip(awt, awt$tip.label[grepl('tdacs2', awt$tip.label)])
awt=drop.tip(awt, awt$tip.label[grepl('tdacn2', awt$tip.label)])
awt=drop.tip(awt, awt$tip.label[grepl('tdactm', awt$tip.label)])
awt=drop.tip(awt, awt$tip.label[grepl('tzopol', awt$tip.label)])
awt=drop.tip(awt, awt$tip.label[grepl('osativ', awt$tip.label)])
awt=drop.tip(awt, awt$tip.label[grepl('bdista', awt$tip.label)])
awt=drop.tip(awt, awt$tip.label[grepl('svirid', awt$tip.label)])

## grampa needs blank_species, so use gene name, _, species
awt$tip.label=paste0(str_split_fixed(filenames[i], '\\.',2)[,2], 1:length(awt$tip.label), '_', awt$tip.label) ## not sure, addign in 1:length(awt$tip.label), to see if unique first part is necesary??
awt$edge.length=NULL
awt$node.label=NULL      
return(awt)
}
}
                                              })
                                              
                                              
d=do.call("c",b)




write.tree(d, paste0('paspalum_anchors_forgrampa.', Sys.Date(), '.tre'))
write.tree(d[1:100], paste0('paspalum_anchors_forgrampa.test.', Sys.Date(), '.tre'))
