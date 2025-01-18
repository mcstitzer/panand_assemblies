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
library(ggh4x) ## facet strips spanning groups (subtribe)


## need to have run in same environment "../sp_tree/plot_sp_tree_pies.R" to have chronopienotip object here

## also to load all the proper mapping packages, "../general_summaries/plot_map_panand.R"
### to get range map working correctly

## also for heterozygosity... "../general_summaries/plot_heterozygosity.R"



ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')


## assembly size, since we don't have flow for everybody
asize=read.table('../general_summaries/panand_assembly_sizes.txt', header=T, sep='\t')
asize$species[3]='\"Andropogon\" burmanicus'
asize$haploid[asize$V2=='zmB735']=T
                                
                 
### reordering here!
taxonnames=c("Zea mays ssp. parviglumis TIL11", "Zea mays ssp. parviglumis TIL01", "Zea mays ssp. mays B73v5", "Zea mays ssp. mexicana TIL25", "Zea mays ssp. mexicana TIL18", "Zea mays ssp. huehuetenangensis", 
"Zea luxurians", "Zea nicaraguensis", "Zea diploperennis Momo", "Zea diploperennis Gigi", "Tripsacum zoloptense", "Tripsacum dactyloides FL", "Tripsacum dactyloides Southern Hap2", 
"Tripsacum dactyloides Northern Hap2", "Tripsacum dactyloides KS", "Tripsacum dactyloides tetraploid", "Urelytrum digitatum", "Vossia cuspidata", "Rhytachne rottboellioides", "Rottboellia tuberculosa", 
"Hemarthria compressa", "Elionurus tripsacoides", "Schizachyrium scoparium", "Schizachyrium microstachyum", "Anatherum virginicum", "Andropogon chinensis", "Andropogon gerardi", 
"Cymbopogon refractus", "Cymbopogon citratus", "Heteropogon contortus", "Themeda triandra", "Bothriochloa laguroides", "Pogonatherum paniceum", "Sorghum bicolor", 
"Ischaemum rugosum", "Sorghastrum nutans", '"Andropogon" burmanicus', "Thelepogon elegans", "Chrysopogon serrulatus", "Paspalum vaginatum")
names(taxonnames)=c("zTIL11", "zTIL01", "zmB735", "zTIL25", "zTIL18", "zmhuet", 
"zluxur", "znicar", "zdmomo", "zdgigi", "tzopol", "tdacs1", "tdacs2", 
"tdacn2", "tdacn1", "tdactm", "udigit", "vcuspi", "rrottb", "rtuber", 
"hcompr", "etrips", "sscopa", "smicro", "avirgi", "achine", "agerar", 
"crefra", "ccitra", "hconto", "ttrian", "blagur", "ppanic", "sbicol", 
"irugos", "snutan", "atenui", "telega", "cserru", "pvagin")



                                                                                 
## subtribes
subtribenames=c("Tripsacinae", "Tripsacinae", "Tripsacinae", "Tripsacinae", "Tripsacinae", "Tripsacinae", 
"Tripsacinae", "Tripsacinae", "Tripsacinae", "Tripsacinae", "Tripsacinae", "Tripsacinae", "Tripsacinae", 
"Tripsacinae", "Tripsacinae", "Tripsacinae", "Rhytachninae", "Rhytachninae", "Rhytachninae", "Ratzeburgiinae", 
"Ratzeburgiinae", "Ratzeburgiinae", "Andropogoninae", "Andropogoninae", "Andropogoninae", "Andropogoninae", "Andropogoninae", 
"Anthistiriinae", "Anthistiriinae", "Anthistiriinae", "Anthistiriinae", "Anthistiriinae", "Germainiinae", "Sorghinae", 
"Ischaeminae", "Apludinae", 'incertae sedis', "incertae sedis", "Chrysopogoninae", "Paspalum vaginatum")
names(subtribenames)=c("zTIL11", "zmB735", "zTIL01", "zTIL25", "zTIL18", "zmhuet", 
"zluxur", "znicar", "zdmomo", "zdgigi", "tzopol", "tdacs1", "tdacs2", 
"tdacn2", "tdacn1", "tdactm", "udigit", "vcuspi", "rrottb", "rtuber", 
"hcompr", "etrips", "sscopa", "smicro", "avirgi", "achine", "agerar", 
"crefra", "ccitra", "hconto", "ttrian", "blagur", "ppanic", "sbicol", 
"irugos", "snutan", "atenui", "telega", "cserru", "pvagin")

subtribenamesSimple=c("Tripsacinae", "Tripsacinae", "Tripsacinae", "Tripsacinae", "Tripsacinae", "Tripsacinae", 
"Tripsacinae", "Tripsacinae", "Tripsacinae", "Tripsacinae", "Tripsacinae", "Tripsacinae", "Tripsacinae", 
"Tripsacinae", "Tripsacinae", "Tripsacinae", "Rhytachninae", "Rhytachninae", "Rhytachninae", "Ratzeburgiinae", 
"Ratzeburgiinae", "Ratzeburgiinae", "Andropogoninae", "Andropogoninae", "Andropogoninae", "Andropogoninae", "Andropogoninae", 
"Anthistiriinae", "Anthistiriinae", "Anthistiriinae", "Anthistiriinae", "Anthistiriinae", "1", "2", 
"3", "4", '5', "6", "7")#, "Paspalum vaginatum")
## 1 Germainiinae
## 2 Sorghinae
## 3 Ischaeminae
## 4 Apludinae
## 5 incertae sedis
## 6 incertae sedis
## 7 Chrysopogoninae
                                
names(subtribenamesSimple)=c("zTIL11", "zmB735", "zTIL01", "zTIL25", "zTIL18", "zmhuet", 
"zluxur", "znicar", "zdmomo", "zdgigi", "tzopol", "tdacs1", "tdacs2", 
"tdacn2", "tdacn1", "tdactm", "udigit", "vcuspi", "rrottb", "rtuber", 
"hcompr", "etrips", "sscopa", "smicro", "avirgi", "achine", "agerar", 
"crefra", "ccitra", "hconto", "ttrian", "blagur", "ppanic", "sbicol", 
"irugos", "snutan", "atenui", "telega", "cserru")#, "pvagin")

                                
subtribedf=data.frame(subtribe=subtribenamesSimple, genome=names(subtribenamesSimple))
subtribedf$haploid=asize$haploid[match(subtribedf$genome, asize$V2)]
subtribedf$species=taxonnames[match(subtribedf$genome, names(taxonnames))]
subtribedf$species=factor(subtribedf$species, levels=taxonnames)
subtribedf$speciesLabel=ifelse(subtribedf$haploid, paste0(subtribedf$species, '*'), as.character(subtribedf$species))
subtribedf$speciesLabel=subtribedf$species
levels(subtribedf$speciesLabel)[levels(subtribedf$speciesLabel) %in% subtribedf$species[subtribedf$haploid]]=paste0(levels(subtribedf$speciesLabel)[levels(subtribedf$speciesLabel) %in% subtribedf$species[subtribedf$haploid]], '*')





## the tree
rsp=read.tree('../sp_tree/paspalum_anchors_aster.2024-11-22.ASTRALPRO3OUT.support.tre')
rsp=drop.tip(rsp, c('pvagin', 'tdacn2', 'tdacs2'))
#rsp=rotateConstr(rsp, rev(names(taxonnames) ))

#taxonnames=taxonnames[rsp$tip.label]



#densit=ggtree(force.ultrametric(rsp, method='extend'))
p2 <- tibble(ymin = c(1,which(!duplicated(subtribedf$subtribe[subtribedf$genome %in% asize$V2]))[-1])-1, ymax = c(which(!duplicated(subtribedf$subtribe[subtribedf$genome %in% asize$V2]))[-1],sum(subtribedf$genome %in% asize$V2)+1)-1, fill = unique(subtribedf$subtribe)) %>%  ggplot() +
  geom_rect(aes(xmin = 0.1, xmax = 0.9, ymin = ymin+0.1, ymax = ymax), color='black', fill='snow2') +
  geom_text(aes(x = .5, y = (ymin  + ymax) / 2, label = fill), angle = 90, size=2.4, fontface = "bold") +scale_y_reverse(breaks = seq(1, 10), expand = expansion(mult = c(0, 0))) +scale_x_continuous(breaks = c(0), expand = expansion(mult = c(0, 0))) +guides(fill = FALSE) +theme_void()


tip_order <- chronopienotip$data %>%
  filter(isTip) %>%          # Keep only tip labels
  arrange(y) %>%             # Sort by y-coordinate (top to bottom)
  pull(label)                # Extract the label column
tip_order=tip_order[!tip_order%in%c('pvagin', 'tdacn2', 'tdacs2')]

asize$species=factor(asize$species, levels=rev(taxonnames[tip_order]))
asize$speciesLabel=ifelse(asize$haploid, paste0(asize$species, '*'), as.character(asize$species))
asize$speciesLabel=asize$species
levels(asize$speciesLabel)[levels(asize$speciesLabel) %in% asize$species[asize$haploid]]=paste0(levels(asize$speciesLabel)[levels(asize$speciesLabel) %in% asize$species[asize$haploid]], '*')
## quotes are hard to parse!!


### plot the haploid genome size, but nee
## and fake a legend
legend_data <- data.frame(
  shape = c("Genome", "TE + TR"),
  fill='black',
  color='black',
  x = c(1,1),
  y = c(100, 100) ## otuside of ylim??
)

hgs=ggplot(asize, aes(x=(haploidAssemblySize-haploidNCount)/1e9, y=1, color=ploidy)) + geom_vline(xintercept=c(2,4), color='snow2', linetype='dotted') + geom_vline(xintercept=c(1,3,5), color='snow3', linetype='dotted') + 
geom_segment(aes(y=1,yend=1, x=0, xend=haploidAssemblySize/1e9)) + scale_color_manual(values=ploidycolors) + scale_fill_manual(values=ploidycolors) + geom_point(aes(x=haploidAssemblySize/1e9), color='snow3', size=2)+ geom_point(size=4)+ geom_point(aes(x=haploidRepeatSize/1e9, bg=ploidy, y=1.6), shape=25, size=3)+ facet_wrap(~speciesLabel, ncol=1, strip.position='left', labeller=purrr::partial(label_species, dont_italicize=c('subsp.', 'ssp.', 'TIL11', 'TIL01', 'TIL25', 'TIL18', 'Momo', 'Gigi', 'Southern Hap1', 'Northern Hap1', 'FL', 'KS',  '\\*', '\\"', 'B73v5'))) + theme(strip.placement = "outside",  strip.background = element_blank(), strip.text.y.left = element_text(angle=0, hjust=1), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ #theme(legend.position = "none") + 
ylab('')+ xlab('Haploid Size (Gb)') + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + ylim(0,2) +
geom_point(data=legend_data, aes(x=x, y=y, shape=shape), size=3, fill='black', color='black')+
scale_shape_manual(values=c('Genome'=21, 'TE + TR'=25))+
theme(legend.position=c(0.71,0.3), legend.text = element_text(size=10)) +
  labs(color='Ploidy', pch='') +guides(color='none', fill='none')



## the tree
## from ../sp_tree/plot_sp_tree_pies.R
# chronopie <- ggtree(chronogram) +
# geom_tiplab()+
# #  geom_text2(aes(label = node), subset = .(isTip == FALSE)) +
#   geom_inset(pies, x = "node", height=0.1, width=0.1)+
#   theme_tree2()+ scale_x_continuous(expand = expansion(mult = 1.5), labels = abs)
chronopienotip


## the map
## from plot_map_panand.R
rangemap

## the heterozygosity
## from ../general_summaries/plot_heterozygosity.R
hethomeo

## the ploidy/flow
ploidyflow=ggplot(asize, aes(x=c(haploidAssemblySize-haploidNCount)/1e9, y=flow/1000, color=ploidy)) +  scale_color_manual(values=ploidycolors)  + geom_point(size=3) + ylab('Genome Size,\nFlow Cytometry (Gb)')+ xlab('Haploid Assembly\nSize (Gb)') 




##### put it all together

maxage=max(chronogram$edge.length[chronogram$edge.length<25]) ## not the paspalum branch
tree_plot <- chronopienotip+  scale_x_continuous(breaks=c( maxage-20+4,maxage-10+4, maxage+4), labels=c('20','10','0')) + xlab('Divergence (Mya)')

treesubtribegs=plot_grid(tree_plot, p2, NULL, hgs, align='h', axis='tb', ncol=4, rel_widths=c(0.1,0.02,-0.01,0.3), labels=c('A', '', 'B', '', ''))

hetflow=plot_grid(hethomeo, ploidyflow + theme(legend.position='NULL'), ncol=2, align='hv', axis='tb', labels=c('D', 'E'))
maphetflow=plot_grid(rangemap, hetflow, ncol=1, align='v', labels=c('C',''), rel_heights=c(1,0.7))

figure1=plot_grid(treesubtribegs, maphetflow, align='hv', axis='tb', ncol=2, rel_widths=c(1,0.8), labels='')


### rethinking because these figures are small!!

## zuerst, try to put points of polyploidy timing
maxage=max(chronogram$edge.length[chronogram$edge.length<25]) ## not the paspalum branch
tree_data=fortify(chronogram)
point_data=merge(het, tree_data, by.x='genome', by.y='label')
perennialmedian=median(point_data$mya[point_data$genome%in%c('tdacn1', 'tdacs1', 'znicar', 'zluxur', 'zdgigi', 'zdmomo')])
perennialy=26.4 ## right above nicaraguensis
point_data=point_data[,c('genome', 'mya', 'y', 'ploidy')]
point_data$mya[point_data$ploidy%in%c('Diploid', 'Paleotetraploid')]=NA
point_data=rbind(point_data, data.frame(genome='tripsacinae', mya=perennialmedian, y=perennialy, ploidy='Paleotetraploid'))


offsettip=3.4 ## offset for the plot margin
tree_plot <- chronopienotip+  
#scale_x_continuous(breaks=c( maxage-20+4,maxage-10+4, maxage+4), labels=c('20','10','0')) + 
scale_x_continuous(breaks=c(maxage-20+offsettip,maxage-15+offsettip,maxage-10+offsettip,maxage-5+offsettip,maxage+offsettip), labels=c('20','15', '10', '5', '0'))+ ## extra 3.5 because of plot margin
xlab('Divergence (Mya)') + 
geom_point(data=point_data, aes(x=maxage-mya+offsettip, y=y), shape=18, size=4.4, color='black')+ ## don't mess up fill legend
geom_point(data=point_data, aes(x=maxage-mya+offsettip, y=y, color=factor(ploidy,levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid') )),shape=18, size=4)+ scale_color_manual(values=ploidycolors) +theme(legend.box.spacing = unit(-1, "cm"), legend.margin = margin(-0.5, 0, 0, 0, unit='cm'), axis.title.x=element_text(size=11), legend.text = element_text(size=10))+ labs(color='', fill='')#+ theme(legend.box = "horizontal")



treesubtribegs=plot_grid(tree_plot, p2, NULL, hgs, align='h', axis='tb', ncol=4, rel_widths=c(0.1,0.01,-0.008,0.2), labels=c('A', '', 'B', '', ''))

maphet=plot_grid(NULL, ggdraw(rangemap), hethomeo, ncol=3, align='v', labels=c('','C','D'), rel_widths=c(-0.15,1,0.5), label_x=c(NA, 0.17, 0.01), label_y=c(NA, 1,    1))

figure1=plot_grid(treesubtribegs, maphet, align='hv', axis='tb', ncol=1, rel_heights=c(1,0.4), labels='')


#plot_grid(densit, p2, NULL, hgs, cpb,NULL,aap, ksp, align='hv',axis='tb', ncol=8, rel_widths=c(0.2,0.05,-0.05,0.7,0.2,-0.03,0.2,0.3), labels=c('a', '','b','', 'c', 'd','', 'e')),
#  bottomlayer,newbottom, ncol=1, align='hv', rel_heights=c(0.7,0.2,0.2))

#pdf('../figures/figure1_intro-to-the-andropogoneae.pdf',18,8)

pdf('../figures/figure1_intro-to-the-andropogoneae.pdf',10,12)
figure1
dev.off()


### subset figure 1 by species!


for(i in asize$V2){
#ksplot
ksdat=final_merged_data %>% filter(ploidy!='Diploid', !is.na(shortspeciesLabel)) %>% group_by(queryChr, refChr, paspChr, genome,ploidy, speciesLabel, shortspeciesLabel) %>%summarize(ks=mean(V17), n=length(unique(gene))) %>% filter(n>10) 
pdf(paste0('../figures/fig1_by_sp/', i, '_fig1.pdf'), 18,8)

## gosh this dang facet order is tricky!! rerun above

hgssp=ggplot(asize, aes(x=(haploidAssemblySize-haploidNCount)/1e9, y=1)) + geom_vline(xintercept=c(2,4), color='snow2', linetype='dotted') + geom_vline(xintercept=c(1,3,5), color='snow3', linetype='dotted') + 
geom_segment(aes(y=1,yend=1, x=0, xend=haploidAssemblySize/1e9), color='gray') + 
scale_color_manual(values=ploidycolors) + scale_fill_manual(values=ploidycolors) + 
geom_point(aes(x=haploidAssemblySize/1e9), color='snow3', size=2)+ 
geom_point(size=4, color='gray')+ 
geom_point(aes(x=haploidRepeatSize/1e9, y=1.6), shape=25, size=3, color='gray', bg='gray')+ 
## these are the colored ones!
geom_segment(data=asize[asize$V2==i,], aes(y=1,yend=1, x=0, xend=haploidAssemblySize/1e9, color=ploidy)) + 
geom_point(data=asize[asize$V2==i,], aes(x=haploidAssemblySize/1e9), color='snow3', size=2)+ geom_point(data=asize[asize$V2==i,], size=4, aes(color=ploidy))+ 
geom_point(data=asize[asize$V2==i,], aes(x=haploidRepeatSize/1e9, bg=ploidy, y=1.6), shape=25, size=3)+ 
facet_wrap(~speciesLabel, ncol=1, strip.position='left', labeller=purrr::partial(label_species, dont_italicize=c('subsp.', 'ssp.', 'TIL11', 'TIL01', 'TIL25', 'TIL18', 'Momo', 'Gigi', 'Southern Hap1', 'Northern Hap1', 'FL', 'KS',  '\\*', '\\"', 'B73v5'))) + 
theme(strip.placement = "outside",  strip.background = element_blank(), strip.text.y.left = element_text(angle=0, hjust=1), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + ylab('')+ xlab('Haploid Size (Gb)') + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + ylim(0,2)



rangemapsp=ggplot() +
  geom_sf(data = ne_countries(type = "countries",returnclass = "sf"), fill = "gray90", color = "gray80", lwd = 0.5) +
  geom_spatraster(data = range_andro_do100, aes(fill = layer), alpha = 0.9) +
  scale_fill_gradient(low = "white", high = "#238b45", na.value = NA) +
  new_scale_fill() +
  #  geom_point(data=asize, aes(x=longloose, y=latloose), color='black', size=2.2) +
  geom_point(data=asize, aes(x=longloose, y=latloose), pch=21, fill='white',color='black',size=2) + 
    geom_point(data=asize[asize$V2==i,], aes(x=longloose, y=latloose, fill=ploidy), pch=21, color='black',size=2) + 
  scale_fill_manual(values=ploidycolors) + theme(legend.position='NULL')+
  xlab('')+ylab('')

hethomeosp=ggplot(het[!is.na(het$type),], aes(x=het, y=ks, pch=ifelse(type=='autopolyploid', 'autopolyploid', ifelse(ploidy=='Diploid', 'diploid', 'allopolyploid'))))+geom_abline(slope=1, intercept=0, lwd=5, color='gray90') +
  geom_point(size=3, alpha=0.8, color='gray')+
  geom_point(data=het[!is.na(het$type) & het$genome==i,], size=3, alpha=0.8, pch=21, aes(fill=ploidy), color='black')+
scale_color_manual(values=ploidycolors)+scale_fill_manual(values=ploidycolors)+theme(legend.position=c(0.61,0.9), legend.text = element_text(size=8)) +
  ylab('Ks between\nduplicates') + xlab('Heterozygosity') + labs(color='Ploidy', pch='') +guides(color='none', fill='none')+
  theme() # Remove y-axis and other clutter

ploidyflowsp=ggplot(asize, aes(x=c(haploidAssemblySize-haploidNCount)/1e9, y=flow/1000)) +  scale_color_manual(values=ploidycolors)  + scale_fill_manual(values=ploidycolors)  +
geom_point(size=3, color='gray')+
geom_point(data=asize[asize$V2==i,], aes(fill=ploidy),pch=21, color='black', size=3) + 
ylab('Genome Size,\nFlow Cytometry (Gb)')+ xlab('Haploid Assembly\nSize (Gb)') 

## don't change the tree...
tree_plot <- chronopienotip+  scale_x_continuous(breaks=c( maxage-20+4,maxage-10+4, maxage+4), labels=c('20','10','0')) + xlab('Divergence (Mya)')

treesubtribegssp=plot_grid(tree_plot, p2, NULL, hgssp, align='h', axis='tb', ncol=4, rel_widths=c(0.1,0.02,-0.01,0.3), labels=c('A', '', 'B', '', ''))

hetflowsp=plot_grid(hethomeosp, ploidyflowsp + theme(legend.position='NULL'), ncol=2, align='hv', axis='tb', labels=c('D', 'E'))
maphetflowsp=plot_grid(rangemapsp, hetflowsp, ncol=1, align='v', labels=c('C',''), rel_heights=c(1,0.7))

figure1sp=plot_grid(treesubtribegssp, maphetflowsp, align='hv', axis='tb', ncol=2, rel_widths=c(1,0.8), labels='')

print(figure1sp)

dev.off()

}











