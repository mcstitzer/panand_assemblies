library(phangorn) ### in here! /workdir/mcs368/panand_htt/scinet_trees
library(ape)
library(ggtree)  ## /programs/R-4.2.1-r9/bin/R
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(ape)
library(treeio)
library(phytools)
library(reshape2)
library(dplyr)
library(ggridges)
library(stringr)
##source('../sp_tree/ggdensitree.R')


all=read.table('../panand_sp_ploidy.txt')                 


filenames=list.files('.', pattern='RAxML_bipartitionsBranchLabels.')
outdf=data.frame(filenames=filenames)
outtrees=list()
#ts=lapply(filenames, function(x){as.phylo(read.raxml(x))})
#tsc=do.call("c",ts)
## pick trees based on bootstrap support across nodes

for(x in 1:length(filenames)){
ts1=read.raxml(filenames[x])
outdf$minBS[x]=min(ts1@data$bootstrap, na.rm=T)
outdf$meanBS[x]=mean(ts1@data$bootstrap, na.rm=T)
outdf$ntipsorig[x]=length(ts1@phylo$tip.label)

## drop tips at random so one per species
ts2=drop.tip(as.phylo(ts1), ts1@phylo$tip.label[duplicated(str_split_fixed(ts1@phylo$tip.label, '_', 2)[,1])])
#awt=root(ts2, as.phylo(ts2)$tip.label[substr(as.phylo(ts2)$tip.label,1,6) %in% c('osativ', 'bdista')])
awt=ts2
awt$tip.label=substr(awt$tip.label,1,6)
awt=drop.tip(awt, c('Pavag0', 'Pavag1', 'agerjg', 'eophiu', 'svirid', 'bdista', 'osativ', 'tdacn2', 'tdacs2', 'tzopol', 'tdactm'))
outdf$ntipssingle[x]=length(awt$tip.label)
outdf$ntaxa[x]=length(unique(awt$tip.label))
outtrees[[x]]=awt
           }

## use backbone of astral tree


## plot densitree

rsp=read.tree('~/transfer/paspalum_anchors_ASTRALout.2023-08-23.tre')
rsp$tip.label[rsp$tip.label=='paspal']='pvagin'

#a=do.call('c', outtrees[which(outdf$minBS>40 & outdf$meanBS>70 & outdf$ntaxa==40)]) ## removed outgroups, so not 43, it's 40!!!!
## run twice - tighter with 33 trees, looser with 958 trees
a=do.call('c', outtrees[which(outdf$minBS>30 & outdf$meanBS>70 & outdf$ntaxa==36)])
a=do.call('c', outtrees[which(outdf$meanBS>50  & outdf$ntaxa==36)])

a=root(a, 'pvagin')
t2=drop.tip(rsp, c('svirid', 'bdista', 'osativ', 'tdacs2', 'tdacn2', 'tdactm', 'tzopol'))
## rotate nodes to put zea-(dactyloides/zapolotense) closer
#t2=rotateNodes(t2, MRCA(t2, c('tdacs1', 'tzopol')))
t2=rotateConstr(t2, rev(t2$tip.label)) ## plots bottom to top
b=do.call('c', c(t2, a)) ## to put the tip order i want up first!!

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

bn=do.call('c', sapply(b, function(x){ x$tip.label=taxonnames[x$tip.label]
          return(x)}))
an=do.call('c', lapply(a, function(x){ x$tip.label=taxonnames[x$tip.label]
          return(x)}))


pdf('~/transfer/densitree_sub.pdf',20,20)
densiTree(a)
densiTree(a, scaleX=T)
densiTree(a, consensus=rsp$tip.label)
densiTree(a, consensus=rsp$tip.label, scaleX=T)
densiTree(a, consensus=rsp$tip.label, type="phylogram")
densiTree(a, consensus=rsp$tip.label, scaleX=T, type="phylogram")
#densityTree(b,type="cladogram",nodes="intermediate") ## failing on something now that i don't understand 
#densityTree(b,type="cladogram",nodes="intermediate", fix.depth=T)
#densityTree(b,use.edge.length=FALSE,type="cladogram",nodes="centered", alpha=0.03)
#densityTree(b,use.edge.length=FALSE,type="cladogram",nodes="centered", alpha=0.03, use.gradient=T)
densityTree(b,use.edge.length=FALSE,type="cladogram",nodes="centered", alpha=0.03, use.gradient=T, compute.consensus=F)
densityTree(b,color='springgreen4',use.edge.length=FALSE,type="cladogram",nodes="centered", alpha=0.03, compute.consensus=F)
densityTree(b,color='springgreen4',use.edge.length=FALSE,type="cladogram",nodes="centered", alpha=0.1, use.gradient=T, compute.consensus=F)
densityTree(b,color='springgreen4',use.edge.length=FALSE,type="cladogram",nodes="centered", alpha=0.1, compute.consensus=F)
densityTree(b,color='springgreen4',use.edge.length=FALSE,type="cladogram",nodes="centered", alpha=0.1, use.gradient=T, compute.consensus=F, fix.depth=T)
densityTree(b,color='springgreen4',use.edge.length=FALSE,type="cladogram",nodes="centered", alpha=0.1, compute.consensus=F, fix.depth=T)

densityTree(bn,use.edge.length=FALSE,type="cladogram",nodes="centered", alpha=0.03, use.gradient=T, compute.consensus=F, fix.depth=T)

# for(i in a){
#            print(ggtree(i)  + geom_tiplab())
#            }

dev.off()

anc=do.call('c', lapply(an, function(x){ chronos(x, lambda=0)
          }))
write.tree(anc, 'temp.tre')
anc=read.tree('temp.tre')

ancl=lapply(a, function(x) force.ultrametric(x, method='extend'))
anc=do.call('c', ancl)

## remove really long trees???
longtrees=which(sapply(anc, function(x) max(nodeHeights(x)[,2]))<0.5)
ancs=anc[longtrees]            
## or scale branches? 
    tree$edge.length <- tree$edge.length / max(nodeHeights(tree)[,2]) * new.tree.length
ans=do.call('c', lapply(anc, function(x){ x$edge.length=x$edge.length/max(nodeHeights(x)[,2])
         return(force.ultrametric(x, method='extend')) }))
##ant <- force.ultrametric(a, method = "extend")

#write.tree(d, paste0('paspalum_anchors_aster.', Sys.Date(), '.tre'))
pdf('~/transfer/densitree_new.pdf',20,20)
#ggtree(anc) + geom_tiplab()
#ggdensitree(an[1:100], layout="rectangular", lwd=1,tip.order=rev(taxonnames), align.tips=T, color="lightblue", alpha=.3,) + geom_tiplab(cex=1)
#ggdensitree(anc[1:200], layout="rectangular",tip.order=rev(taxonnames), align.tips=T, color="lightblue", alpha=.3,) + geom_tiplab(cex=1)
#ggdensitree(ans[1:100], layout="rectangular", lwd=1,tip.order=rev(taxonnames), align.tips=T, color="lightblue", alpha=.3,) + geom_tiplab(cex=1)
#ggdensitree(ans[1:100], layout='fan', lwd=1,tip.order=rev(taxonnames), align.tips=T, color="lightblue", alpha=.3,) + geom_tiplab(cex=1)

options(ignore.negative.edge=TRUE)
ggdensitree(anc[1:200], layout="rectangular",tip.order=names(rev(taxonnames)), align.tips=T, color="lightblue", alpha=.3,) + geom_tiplab(cex=1)
ggdensitree(anc[1:200], layout="rectangular",tip.order=names(rev(taxonnames)), align.tips=F, color="lightblue", alpha=.3,) + geom_tiplab(cex=1)
ggdensitree(ans[1:200], layout="rectangular",tip.order=names(rev(taxonnames)), align.tips=T, color="lightblue", alpha=.3,) + geom_tiplab(cex=1)
ggdensitree(ancs[1:200], layout="rectangular",tip.order=names(rev(taxonnames)), align.tips=T, color="lightblue", alpha=.3,) + geom_tiplab(cex=1)
ggdensitree(ancs[1:200], layout="slanted",lwd=1,tip.order=names(rev(taxonnames)), align.tips=T, color="lightblue", alpha=.3,) + geom_tiplab(cex=1)

ggtree(anc[[1]])
ggtree(ans[[1]])
dev.off()


#densityTree(b,use.edge.length=FALSE,type="cladogram",nodes="centered", alpha=0.03, use.gradient=T, compute.consensus=F)


## k need to get base r densitytree from phytools into ggplot format - i'm giving up on this ggdensitree stuff

# define a function that emits the desired plot
p1 <- function() {
  par(
    mar = c(3, 3, 1, 1),
    mgp = c(2, 1, 0)
  )
  densityTree(b,use.edge.length=FALSE,type="cladogram",nodes="centered", alpha=0.03, use.gradient=T, compute.consensus=F)

}

densit=ggdraw(p1)
