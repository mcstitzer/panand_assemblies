library(ggtree)  ## /programs/R-4.2.1-r9/bin/R
library(ggplot2) ## on cbsu
library(cowplot)
theme_set(theme_cowplot())
library(ape)
library(treeio)
library(phytools)
library(reshape2)
library(dplyr)
library(ggridges)



filenames=list.files('.', pattern='RAxML_bipartitionsBranchLabels.')
                                       
outdf=data.frame(filenames=filenames)
all=read.table('../panand_sp_ploidy.txt')                                
#for(sp in c('vcuspi', 'rtuber', 'blagur', 'hcompr', 'udigit', 'telega')){
for(sp in all$V2){

  outdf[,paste0(sp, 'Count')]=NA
  outdf[,paste0(sp, 'Monophyletic')]=NA
#brlenlist=list()

for(i in (1:length(filenames))){ ## ep2 is not there

awt=read.raxml(paste0('',filenames[i]))
#tryCatch(("reroot" (
#awt=root(awt, awt$tip.label[substr(awt$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate')]))), error = function(err) print("Outgroups already at the root"))
awt=as.phylo(awt)
awt$tip.label=gsub('_R_', '', awt$tip.label)
  ## add to make sure outgroup is there, and is monophyletic - otherwise skip this tree
  if(any(substr(as.phylo(awt)$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate'))){
    if(is.monophyletic(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate')])){
awt=root(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate')])

vc=sum(substr(as.phylo(awt)$tip.label,1,6) %in% sp)
outdf[i, paste0(sp, 'Count')]=vc
if(vc>1){
mrcanode=getMRCA(as.phylo(awt), as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c(sp)])
getDescendants(as.phylo(awt), mrcanode)
tips=as.phylo(awt)$tip.label[getDescendants(as.phylo(awt), mrcanode)]
tips=tips[!is.na(tips)]
outdf[i, paste0(sp, 'Monophyletic')]=all(substr(tips,1,6)==sp)
  ## get terminal branch length for each copy
  ## first get the node numbers of the tips
nodes<-sapply(as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c(sp)],function(x,y) which(y==x),y=as.phylo(awt)$tip.label)
## then get the edge lengths for those nodes
edge.lengths<-setNames(awt$edge.length[sapply(nodes,function(x,y) which(y==x),y=awt$edge[,2])],names(nodes))
#brlenlist[i]=edge.lengths
}}
}
}
                                              }
                                              
                                              


## plot out each species, with barplot of Count,  and barplot of monophyletic count (or split bars of monophyly on first barplot!)
melt(outdf[,-1])                                              
## 80 columns, take and combine
outdflist=lapply(1:40*2, function(x){ 
  temp=outdf[,c(x, x+1)]
  sp=substr(colnames(temp),1,6)[1]
  colnames(temp)=gsub(sp, '', colnames(temp))
  temp$genome=sp
  return(temp)
  })
tp=do.call(rbind, outdflist)
##dplyr is being stupid, or i am, so just make columsn
tp$Polyphyletic=!tp$Monophyletic
tp$NotApplicable=is.na(tp$Monophyletic)
tpp=tp[!is.na(tp$Count),] %>% group_by(Count, genome) %>% dplyr::summarize(Monophyletic=sum(Monophyletic, na.rm=T), Polyphyletic=sum(Polyphyletic, na.rm=T), NotApplicable=sum(NotApplicable))

pdf(paste0('~/transfer/genetree_counting.', Sys.Date(), '.pdf'), 35,5)
ggplot(tp, aes(x=genome,  group=Count, fill=Monophyletic)) + geom_histogram(stat='count', position='dodge') + scale_fill_manual(values=c('slateblue', 'gray'))
ggplot(tp, aes(x=factor(Count), y=genome,  group=Count, fill=Monophyletic)) + stat_binline() + scale_fill_manual(values=c('slateblue', 'gray'))
ggplot(tp[!is.na(tp$Count),], aes(x=factor(Count),  group=Count, fill=Monophyletic)) + geom_histogram(stat='count', position='dodge') + facet_wrap(~genome, nrow=1) + scale_fill_manual(values=c('slateblue', 'forestgreen'))
ggplot(melt(tpp[!is.na(tpp$Count),], id.vars=c('genome', 'Count')), aes(x=factor(Count), y=value, group=value, fill=variable)) + geom_histogram(stat='identity', position='stack') + facet_wrap(~genome, nrow=1) + scale_fill_manual(values=c('slateblue', 'forestgreen', 'gray'))

                                              
dev.off()   
 
