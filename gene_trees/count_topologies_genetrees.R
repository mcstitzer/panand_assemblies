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



ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')
ploidycolorsmonophyletic=sapply(c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf'), function(x) colorspace::lighten(x, amount=0.7))
names(ploidycolorsmonophyletic)=paste0(c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid'), 'Monophyletic')
ploidycolorspolyphyletic=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolorspolyphyletic)=paste0(c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid'), 'Polyphyletic')

ploidyphyletic=c(ploidycolorsmonophyletic, ploidycolorspolyphyletic)
                                
taxonnames=c("Zea mays subsp. parviglumis TIL11", "Zea mays subsp. mays B73v5", "Zea mays subsp. parviglumis TIL01", "Zea mays subsp. mexicana TIL25", "Zea mays subsp. mexicana TIL18", "Zea mays subsp. huehuetengensis", 
"Zea luxurians", "Zea nicaraguensis", "Zea diploperennis Momo", "Zea diploperennis Gigi", "Tripsacum zoloptense", "Tripsacum dactyloides Southern Hap1", "Tripsacum dactyloides Southern Hap2", 
"Tripsacum dactyloides Northern Hap2", "Tripsacum dactyloides Northern Hap1", "Tripsacum dactyloides tetraploid", "Urelytrum digitatum", "Vossia cuspidata", "Rhytachne rottboellioides", "Rottboellia tuberculosa", 
"Hemarthria compressa", "Elionurus tripsacoides", "Schizachyrium scoparium", "Schizachyrium microstachyum", "Andropogon virginius", "Andropogon chinensis", "Andropogon gerardi", 
"Cymbopogon refractus", "Cymbopogon citratus", "Heteropogon contortus", "Themeda triandra", "Bothriochloa laguroides", "Pogonatherum paniceum", "Sorghum bicolor", 
"Ischaemum rugosum", "Sorghastrum nutans", "Andropogon tenuifolius", "Thelopogon elegans", "Chrysopogon serrulatus", "Paspalum vaginatum")
names(taxonnames)=c("zTIL11", "zmB735", "zTIL01", "zTIL25", "zTIL18", "zmhuet", 
"zluxur", "znicar", "zdmomo", "zdgigi", "tzopol", "tdacs1", "tdacs2", 
"tdacn2", "tdacn1", "tdactm", "udigit", "vcuspi", "rrottb", "rtuber", 
"hcompr", "etrips", "sscopa", "smicro", "avirgi", "achine", "agerar", 
"crefra", "ccitra", "hconto", "ttrian", "blagur", "ppanic", "sbicol", 
"irugos", "snutan", "atenui", "telega", "cserru", "pvagin")

                                
all=read.table('../panand_sp_ploidy.txt', header=F)
all=all[!all$V2 %in% c('tdactm', 'agerjg', 'bdista', 'eophiu', 'osativ', 'svirid', 'tdactn', 'tdacts'),]
all$boxplotx=all$V3*2
all$polyploid=all$V3>1
all$trip=all$V2 %in% c('tdacn1', 'tdacn2', 'tdacs1', 'tdacs2', 'tdactm', 'zdgigi', 'zdmomo', 'zluxur', 'zmB735', 'zmhuet', 'znicar', 'zTIL01', 'zTIL11', 'zTIL18', 'zTIL25')         


all$boxplotx[all$boxplotx==2]='Diploid'
all$boxplotx[all$boxplotx==4 & !all$trip]='Tetraploid'
#all$boxplotx[all$boxplotx==4 & all$trip]='Paleotetraploid'
all$boxplotx[all$trip]='Paleotetraploid'
all$boxplotx[all$boxplotx==6]='Hexaploid'
all$boxplotx[all$boxplotx==8]='Octaploid'
all$boxplotx=factor(all$boxplotx, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid', 'Octaploid'))                                                   

tppp=melt(tpp[!is.na(tpp$Count) & !tpp$genome %in% c('tdactm', 'tzopol', 'osativ', 'pprate'),], id.vars=c('genome', 'Count'))
tppp$genome=as.factor(tppp$genome, levels=names(taxonnames))
tppp$ploidy=all$boxplotx[match(tppp$genome, all$V2)]
tppp$phyleticcol=paste0(tppp$ploidy, tppp$variable)
## make singletons "monophyletic"
tppp$phyleticcol[tppp$Count==1]=paste0(tppp$ploidy[tppp$Count==1], 'Monophyletic')
tppp$species=taxonnames[match(tppp$genome, names(taxonnames))]
tppp$species=factor(tppp$species, levels=taxonnames)

                               
#pdf(paste0('~/transfer/genetree_synteny.', Sys.Date(), '.pdf'), 5,10)
pdf(paste0('genetree_synteny.', Sys.Date(), '.pdf'), 5,12)

ggplot(tppp, aes(x=factor(Count), y=value, group=phyleticcol, fill=phyleticcol)) + geom_histogram(stat='identity', position='stack') + facet_wrap(~species, ncol=1, strip.position='left') + scale_fill_manual(values=ploidyphyletic) +   theme(strip.placement = "outside",  strip.background = element_blank(), strip.text.y.left = element_text(angle=0), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + ylab('') + xlab('Syntenic Gene\nCopy Number')+ scale_y_continuous(n.breaks = 2)
ggplot(tppp[tppp$Count!=0,], aes(x=factor(Count), y=value, group=phyleticcol, fill=phyleticcol)) + geom_histogram(stat='identity', position='stack') + facet_wrap(~species, ncol=1, strip.position='left') + scale_fill_manual(values=ploidyphyletic) +   theme(strip.placement = "outside",  strip.background = element_blank(), strip.text.y.left = element_text(angle=0), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + ylab('')+ xlab('Syntenic Gene\nCopy Number')+ scale_y_continuous(n.breaks = 2)

dev.off()
                                              
