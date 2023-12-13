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
library(tidypaleo) ## facet species names in italics!!!!!


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
"Hemarthria compressa", "Elionurus tripsacoides", "Schizachyrium scoparium", "Schizachyrium microstachyum", "Andropogon virginicus", "Andropogon chinensis", "Andropogon gerardi", 
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
all$V3[all$V2=='ccitra']=2
all$V3[all$V2 %in%c('telega', 'rtuber')]=1 ## flow shows this shouldn't be doubled - it's a diploid!!!
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

tppp=melt(tpp[!is.na(tpp$Count) & !tpp$genome %in% c('tdactm', 'tzopol', 'osativ', 'pprate', 'tdacs2', 'tdacn2'),], id.vars=c('genome', 'Count'))
tppp$genome=factor(tppp$genome, levels=names(taxonnames))
tppp$ploidy=all$boxplotx[match(tppp$genome, all$V2)]
tppp$phyleticcol=paste0(tppp$ploidy, tppp$variable)
## make singletons "monophyletic"
tppp$phyleticcol[tppp$Count==1]=paste0(tppp$ploidy[tppp$Count==1], 'Monophyletic')
tppp$species=taxonnames[match(tppp$genome, names(taxonnames))]
tppp$species=factor(tppp$species, levels=taxonnames)

                                
## double up those haploids
gc=read.table('../panand_sp_genecopyguess.txt', header=T, sep='\t')
tg=merge(all, gc, by.x='V2', by.y='X6.letter.code', all=T)
tg$homolog.state[substr(tg$V2,1,4)=='tdac']='haploid' ## make sure you've removed tdactm!!!!
tg$homolog.state[tg$V2 %in% c('sbicol', 'zdmomo', 'zdgigi', 'znicar', 'telega', 'rrottb')]='haploid' ## gigi and momo are clearly not haploid, but do it anyways!!
tg$homolog.state[tg$V2%in%c('telega')]='allelic'
tg$homolog.state[tg$V2%in%c('rtuber')]='haploid'
# tg=tg[tg$genome!='pprate',]
# tg$V3[tg$genome=='ccitra']=2
tppp$haploid=(tg$homolog.state=='haploid')[match(tppp$genome, tg$V2)]
tppp$doubledCount=ifelse(tppp$haploid, tppp$Count*2, tppp$Count)

tppp$speciesLabel=ifelse(tppp$haploid, paste0(tppp$species, '*'), as.character(tppp$species))
tppp$speciesLabel=tppp$species
levels(tppp$speciesLabel)[levels(tppp$speciesLabel) %in% tppp$species[tppp$haploid]]=paste0(levels(tppp$speciesLabel)[levels(tppp$speciesLabel) %in% tppp$species[tppp$haploid]], '*')

tppp$linetype=NA
tppp$linetype[tppp$doubledCount%in%1:6]=rep(c('dotted', 'dashed'),3)[tppp$doubledCount[tppp$doubledCount%in%1:6]]

## switch it back, not thinking this through!
tppp$variable[tppp$Count==1]='NotApplicable'                                

## assembly size, since we don't have flow for everybody
asize=read.table('../panand_assembly_sizes.txt', header=F)
asize=asize[asize$V2 %in% tppp$genome,]
asize$haploid=(tg$homolog.state=='haploid')[match(asize$V2, tg$V2)]

asize$doubledAssembly=ifelse(asize$haploid, asize$V3*2, asize$V3)
asize$species=taxonnames[match(asize$V2, names(taxonnames))]
asize$species=factor(asize$species, levels=taxonnames)

asize$speciesLabel=ifelse(asize$haploid, paste0(asize$species, '*'), as.character(asize$species))
asize$speciesLabel=asize$species
levels(asize$speciesLabel)[levels(asize$speciesLabel) %in% asize$species[asize$haploid]]=paste0(levels(asize$speciesLabel)[levels(asize$speciesLabel) %in% asize$species[asize$haploid]], '*')
asize$ploidy=all$boxplotx[match(asize$V2, all$V2)]

## add flow, for supp
flow=read.table('../panand_flow_cyt.txt', header=T, sep='\t') 
asize$flow=flow[,2][match(asize$V2, flow[,1])]
cor.test(asize$doubledAssembly/2, asize$flow, use='complete.obs')

pdf(paste0('~/transfer/supp_flow_assembly.', Sys.Date(), '.pdf'), 4,4)
## "haploid" assembly size
ggplot(asize, aes(x=doubledAssembly/1e9/2, y=flow/1000, color=ploidy)) +  scale_color_manual(values=ploidycolors)  + geom_point() + ylab('Genome Size, Flow Cytometry (Gb)')+ xlab('Haploid Assembly\nSize (Gb)') 
                                
dev.off()

## also write means in text
asize %>% group_by(ploidy) %>% summarize(flow=mean(flow, na.rm=T), gs=mean(doubledAssembly/2, na.rm=T))
## MAKE SURE I HAVE THIS
asize$haploidAssemblySize=asize$doubledAssembly/2
write.table(asize[,c('V2', 'haploid', 'species', 'speciesLabel', 'ploidy', 'flow', 'haploidAssemblySize')], '~/transfer/panand_assembly_sizes.txt', sep='\t', quote=F, row.names=F, col.names=T)


## okay now try ks!
#ks=read.table('~/transfer/paspalum_omega.txt', header=F)
#ks=ks[substr(ks$V3,1,3)=='_Pa',]
# paspdups=unique(ks$V3[substr(ks$V3,1,3)=='_Pa' & substr(ks$V4,1,3)=='_Pa'])
# ks=ks[!ks$V3 %in% paspdups,]
ks=fread('../orthoFinderOG_omega_skh/*.omega')
ks=ks[substr(ks$V3,1,8)==substr(ks$V4,1,8),] ## witin species
## this is temp for now!!
taxon=c("_Ab00001", "_Ac00001", "_Ag00001", "_Av00001", "_Bl00001", 
"_Cc00001", "_Cr00001", "_Cs00001", "_Et00001", "_Hc00001", "_Hp00001", 
"_Ir00001", "_Pi00001", "_Rr00001", "_Rt00001", 
 "_Sm00001", "_Sn00001", 
"_Sobic.0", "_Sobic.K", "_Ss00002", "_Td00001", "_Td00002", 
"_Te00001", "_Tt00001", "_Ud00001", "_Vc00001", "_Zd00001", 
"_Zd00003", "_Zh00001", "_Zm00001", "_Zn00001", "_Zv00001", "_Zv00002", 
"_Zx00002", "_Zx00003")
names(taxon)=c("atenui", "achine", "agerar",  "avirgi", "blagur",  
"ccitra", "crefra", "cserru","etrips", "hconto","hcompr",  "irugos",
"ppanic", "rrottb", "rtuber", "smicro",  "snutan",'sbicol', 'sbicol', "sscopa", 
 "tdacs1", "tdacn1", "telega","ttrian", "udigit", 
"vcuspi", "zdgigi", "zdmomo","zmhuet","zmB735","znicar", "zTIL01", "zTIL11", "zTIL18", "zTIL25")
ks$genome=names(taxon)[match(substr(ks$V4,1,8), taxon)]   
ks$genome[ks$genome=='zmB735']='zluxur' ##for now, to make the plot look okay
                                
ks$haploid=(tg$homolog.state=='haploid')[match(ks$genome, tg$V2)]

ks$species=taxonnames[match(ks$genome, names(taxonnames))]
ks$species=factor(ks$species, levels=taxonnames)

ks$speciesLabel=ifelse(ks$haploid, paste0(ks$species, '*'), as.character(ks$species))
ks$speciesLabel=ks$species
levels(ks$speciesLabel)[levels(ks$speciesLabel) %in% ks$species[ks$haploid]]=paste0(levels(ks$speciesLabel)[levels(ks$speciesLabel) %in% ks$species[ks$haploid]], '*')
ks$ploidy=all$boxplotx[match(ks$genome, all$V2)]
ks=ks[ks$V17<0.3,]
## don't want to plot those with 0's! 
dupnums=ks %>% group_by(genome) %>% summarize(n=n())
ks$V17[ks$genome %in% dupnums$genome[dupnums$n<5000]]=NA      


                                
pdf(paste0('~/transfer/genetree_synteny.', Sys.Date(), '.pdf'), 5,12)
#pdf(paste0('genetree_synteny.', Sys.Date(), '.pdf'), 5,12)

ggplot(tppp, aes(x=factor(Count), y=value, group=phyleticcol, fill=phyleticcol)) + geom_histogram(stat='identity', position='stack') + facet_wrap(~species, ncol=1, strip.position='left') + scale_fill_manual(values=ploidyphyletic) +   theme(strip.placement = "outside",  strip.background = element_blank(), strip.text.y.left = element_text(angle=0), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + ylab('') + xlab('Syntenic Gene\nCopy Number')+ scale_y_continuous(n.breaks = 2)
ggplot(tppp[tppp$Count!=0,], aes(x=factor(Count), y=value, group=phyleticcol, fill=phyleticcol)) + geom_histogram(stat='identity', position='stack') + facet_wrap(~species, ncol=1, strip.position='left') + scale_fill_manual(values=ploidyphyletic) +   theme(strip.placement = "outside",  strip.background = element_blank(), strip.text.y.left = element_text(angle=0), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + ylab('')+ xlab('Syntenic Gene\nCopy Number')+ scale_y_continuous(n.breaks = 2)

ggplot(tppp, aes(x=factor(doubledCount), y=value, group=phyleticcol, fill=phyleticcol)) + geom_histogram(stat='identity', position='stack') + facet_wrap(~species, ncol=1, strip.position='left') + scale_fill_manual(values=ploidyphyletic) +   theme(strip.placement = "outside",  strip.background = element_blank(), strip.text.y.left = element_text(angle=0), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + ylab('') + xlab('Syntenic Gene\nCopy Number')+ scale_y_continuous(n.breaks = 2)
ggplot(tppp[tppp$doubledCount!=0,], aes(x=factor(doubledCount), y=value, group=phyleticcol, fill=phyleticcol)) + geom_histogram(stat='identity', position='stack') + facet_wrap(~species, ncol=1, strip.position='left') + scale_fill_manual(values=ploidyphyletic) +   theme(strip.placement = "outside",  strip.background = element_blank(), strip.text.y.left = element_text(angle=0), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + ylab('')+ xlab('Syntenic Gene\nCopy Number')+ scale_y_continuous(n.breaks = 2)


## remove auto/allo distinction
ggplot(tppp, aes(x=factor(doubledCount, levels=c(0:6)), y=value, group=ploidy, fill=ploidy)) + geom_histogram(stat='identity', position='stack') + facet_wrap(~speciesLabel, ncol=1, strip.position='left') + scale_fill_manual(values=ploidycolors) +   theme(strip.placement = "outside",  strip.background = element_blank(), strip.text.y.left = element_text(angle=0), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + ylab('') + xlab('Syntenic Gene\nCopy Number')+ scale_y_continuous(n.breaks = 2)
ggplot(tppp[tppp$doubledCount%in% 1:6 & !is.na(tppp$species),], aes(x=doubledCount, y=value, group=ploidy, fill=ploidy)) + geom_vline(xintercept=c(1,3,5), color='snow2', linetype='dotted') + geom_vline(xintercept=c(2,4,6), color='snow3', linetype='dotted') + geom_histogram(stat='identity', position='stack') + 
        facet_wrap(~speciesLabel, ncol=1, strip.position='left', labeller=purrr::partial(label_species, dont_italicize=c('subsp.', 'TIL11', 'TIL01', 'TIL25', 'TIL18', 'Momo', 'Gigi', 'Southern Hap1', 'Northern Hap1', '\\*'))) + scale_fill_manual(values=ploidycolors) +   theme(strip.placement = "outside",  strip.background = element_blank(), strip.text.y.left = element_text(angle=0), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + ylab('')+ xlab('Syntenic Gene\nCopy Number')+ scale_y_continuous(n.breaks = 2)


## pies for auto/allo distinction
tppp %>% group_by(speciesLabel, variable) %>% summarize(n=n(), count=sum(value)) %>% filter(variable!='NotApplicable') %>% mutate(pct = count/sum(count)*100, width=sum(count))%>%
                                ggplot(aes(x=width/2, y=pct, fill=variable, width=width)) + geom_bar(stat='identity', position='fill') + coord_polar(theta='y') + facet_wrap(~speciesLabel, ncol=1, strip.position='left')+   theme(strip.placement = "outside",  strip.background = element_blank(), strip.text.y.left = element_text(angle=0), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + theme(axis.ticks=element_blank(), axis.text=element_blank(), panel.grid=element_blank(), panel.border=element_blank()) + scale_fill_manual(values=c('#5F4B8BFF', '#E69A8DFF'))

## assembly size
ggplot(asize, aes(x=doubledAssembly/1e6, y=1, color=ploidy)) + geom_segment(aes(y=1,yend=1, x=0, xend=doubledAssembly/1e6)) + geom_vline(xintercept=c(1,3,5,7,9)*1000, color='snow2', linetype='dotted') + geom_vline(xintercept=c(2,4,6,8,10)*1000, color='snow3', linetype='dotted') + scale_color_manual(values=ploidycolors)  + geom_point(size=4)+ facet_wrap(~speciesLabel, ncol=1, strip.position='left') + theme(strip.placement = "outside",  strip.background = element_blank(), strip.text.y.left = element_text(angle=0), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + ylab('')+ xlab('Assembly Size (Mb)') + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())

## "haploid" assembly size
ggplot(asize, aes(x=doubledAssembly/1e9/2, y=1, color=ploidy)) + geom_segment(aes(y=1,yend=1, x=0, xend=doubledAssembly/1e9/2)) + geom_vline(xintercept=c(2,4), color='snow2', linetype='dotted') + geom_vline(xintercept=c(1,3,5), color='snow3', linetype='dotted') + scale_color_manual(values=ploidycolors)  + geom_point(size=4)+ facet_wrap(~speciesLabel, ncol=1, strip.position='left', labeller=purrr::partial(label_species, dont_italicize=c('subsp.', 'TIL11', 'TIL01', 'TIL25', 'TIL18', 'Momo', 'Gigi', 'Southern Hap1', 'Northern Hap1', '\\*'))) + theme(strip.placement = "outside",  strip.background = element_blank(), strip.text.y.left = element_text(angle=0), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + ylab('')+ xlab('Haploid Assembly\nSize (Gb)') + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())


## ks
ggplot(ks, aes(x=V17,  color=ploidy, fill=ploidy)) + geom_density() + scale_color_manual(values=ploidycolors) + scale_fill_manual(values=ploidycolors) +  facet_wrap(~speciesLabel, ncol=1, strip.position='left', labeller=purrr::partial(label_species, dont_italicize=c('subsp.', 'TIL11', 'TIL01', 'TIL25', 'TIL18', 'Momo', 'Gigi', 'Southern Hap1', 'Northern Hap1', '\\*'))) + theme(strip.placement = "outside",  strip.background = element_blank(), strip.text.y.left = element_text(angle=0), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + ylab('')+ xlab('Ks between duplicates') + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
 
                                
dev.off()



         library(gridGraphics)                                     
pdf(paste0('~/transfer/genetree_synteny.fig1combo.', Sys.Date(), '.pdf'), 9,10)
## "haploid" assembly size
hgs=ggplot(asize, aes(x=doubledAssembly/1e9/2, y=1, color=ploidy)) + geom_segment(aes(y=1,yend=1, x=0, xend=doubledAssembly/1e9/2)) + geom_vline(xintercept=c(2,4), color='snow2', linetype='dotted') + geom_vline(xintercept=c(1,3,5), color='snow3', linetype='dotted') + scale_color_manual(values=ploidycolors)  + geom_point(size=4)+ facet_wrap(~speciesLabel, ncol=1, strip.position='left', labeller=purrr::partial(label_species, dont_italicize=c('subsp.', 'TIL11', 'TIL01', 'TIL25', 'TIL18', 'Momo', 'Gigi', 'Southern Hap1', 'Northern Hap1', '\\*'))) + theme(strip.placement = "outside",  strip.background = element_blank(), strip.text.y.left = element_text(angle=0), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + ylab('')+ xlab('Haploid Assembly\nSize (Gb)') + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
## copy bar plots, no label
cpb=ggplot(tppp[tppp$doubledCount%in% 1:6 & !is.na(tppp$species),], aes(x=doubledCount, y=value, group=ploidy, fill=ploidy)) + geom_hline(yintercept=1, color='snow2', linetype='dotted') + geom_vline(xintercept=c(1,3,5), color='snow2', linetype='dotted') + geom_vline(xintercept=c(2,4,6), color='snow3', linetype='dotted') + geom_histogram(stat='identity', position='stack') + 
        facet_wrap(~speciesLabel, ncol=1, strip.position='left') + scale_fill_manual(values=ploidycolors) +   theme( strip.background = element_blank(), strip.text.y.left = element_blank(), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + ylab('')+ xlab('Syntenic Gene\nCopy Number')+ scale_y_continuous(n.breaks = 2)
## auto/allo pies, no label
aap=tppp %>% group_by(speciesLabel, variable) %>% summarize(n=n(), count=sum(value)) %>% filter(variable!='NotApplicable') %>% mutate(pct = count/sum(count)*100, width=sum(count))%>%
                                ggplot(aes(x=width/2, y=pct, fill=variable, width=width)) + geom_bar(stat='identity', position='fill') + coord_polar(theta='y') + facet_wrap(~speciesLabel, ncol=1, strip.position='left')+   theme( strip.background = element_blank(), strip.text.y.left = element_blank(), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + theme(axis.ticks=element_blank(), axis.text=element_blank(), panel.grid=element_blank(), panel.border=element_blank()) + scale_fill_manual(values=c('#5F4B8BFF', '#E69A8DFF')) + ylab('Proportion\nDuplicates\nMonophyletic') + xlab('')
## ks distributions, no label
ksp=ggplot(ks[!is.na(ks$speciesLabel),], aes(x=V17,  color=ploidy, fill=ploidy)) + geom_density() + scale_color_manual(values=ploidycolors) + scale_fill_manual(values=ploidycolors) +  facet_wrap(~speciesLabel, ncol=1, strip.position='left', scales='free_y') + theme( strip.background = element_blank(), strip.text.y.left = element_blank(), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + ylab('')+ xlab('Ks Between\nDuplicates') + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
                                
## densit isn't in this file - done in sp_tree densitree :( will combine in a true figure code when I'm happy!!
plot_grid( hgs, cpb,ksp, aap,  align='hv',axis='tb', ncol=4, rel_widths=c(0.7,0.3,0.18,0.2), labels=c('b', 'c', 'd', 'e'))
#plot_grid(densit, hgs, cpb, aap, align='hv',axis='tb', ncol=4, rel_widths=c(0.5,0.7,0.3,0.18), labels=c('a','b', 'c', 'd'))

dev.off()



                                
