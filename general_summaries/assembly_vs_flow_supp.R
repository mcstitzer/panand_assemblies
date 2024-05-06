library(ggplot2) ## on macbook
library(cowplot)
theme_set(theme_cowplot())
library(data.table)


setwd('~/Documents/GitHub/panand_assemblies/general_summaries/')

all=read.table('../panand_sp_ploidy.txt', header=F)
all=all[!all$V2 %in% c('tdactm', 'agerjg', 'bdista', 'eophiu', 'osativ', 'svirid', 'tdactn', 'tdacts'),]
## changed these in file as of Dec 22!
all$V3[all$V2=='ccitra']=2
all$V3[all$V2 %in%c('telega', 'rtuber')]=1 ## flow shows this shouldn't be doubled - it's a diploid!!!
all$boxplotx=all$V3*2
all$polyploid=all$V3>1
all$trip=all$V2 %in% c('tdacn1', 'tdacn2', 'tdacs1', 'tdacs2', 'tdactm', 'zdgigi', 'zdmomo', 'zluxur', 'zmB735', 'zmhuet', 'znicar', 'zTIL01', 'zTIL11', 'zTIL18', 'zTIL25')         


ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')
ploidycolorsmonophyletic=sapply(c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf'), function(x) colorspace::lighten(x, amount=0.7))
names(ploidycolorsmonophyletic)=paste0(c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid'), 'Monophyletic')
ploidycolorspolyphyletic=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolorspolyphyletic)=paste0(c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid'), 'Polyphyletic')

ploidyphyletic=c(ploidycolorsmonophyletic, ploidycolorspolyphyletic)

taxonnames=c("Zea mays ssp. parviglumis TIL11", "Zea mays ssp. mays B73v5", "Zea mays ssp. parviglumis TIL01", "Zea mays ssp. mexicana TIL25", "Zea mays ssp. mexicana TIL18", "Zea mays ssp. huehuetengensis", 
             "Zea luxurians", "Zea nicaraguensis", "Zea diploperennis Momo", "Zea diploperennis Gigi", "Tripsacum zoloptense", "Tripsacum dactyloides FL", "Tripsacum dactyloides Southern Hap2", 
             "Tripsacum dactyloides Northern Hap2", "Tripsacum dactyloides KS", "Tripsacum dactyloides tetraploid", "Urelytrum digitatum", "Vossia cuspidata", "Rhytachne rottboellioides", "Rottboellia tuberculosa", 
             "Hemarthria compressa", "Elionurus tripsacoides", "Schizachyrium scoparium", "Schizachyrium microstachyum", "Anatherum virginicum", "Andropogon chinensis", "Andropogon gerardi", 
             "Cymbopogon refractus", "Cymbopogon citratus", "Heteropogon contortus", "Themeda triandra", "Bothriochloa laguroides", "Pogonatherum paniceum", "Sorghum bicolor", 
             "Ischaemum rugosum", "Sorghastrum nutans", '"Andropogon" burmanicus', "Thelepogon elegans", "Chrysopogon serrulatus", "Paspalum vaginatum")
names(taxonnames)=c("zTIL11", "zmB735", "zTIL01", "zTIL25", "zTIL18", "zmhuet", 
                    "zluxur", "znicar", "zdmomo", "zdgigi", "tzopol", "tdacs1", "tdacs2", 
                    "tdacn2", "tdacn1", "tdactm", "udigit", "vcuspi", "rrottb", "rtuber", 
                    "hcompr", "etrips", "sscopa", "smicro", "avirgi", "achine", "agerar", 
                    "crefra", "ccitra", "hconto", "ttrian", "blagur", "ppanic", "sbicol", 
                    "irugos", "snutan", "atenui", "telega", "cserru", "pvagin")

## assembly size, since we don't have flow for everybody
asize=fread('../general_summaries/panand_assembly_sizes.txt', header=T, quote="", fill=T)
#asize=asize[asize$V2 %in% tppp$genome,]
#asize$haploid=gs$haploid[match(asize$V2, gs$V2)]
#asize$haploid[asize$V2=='zmB735']=T

#asize$doubledAssembly=ifelse(asize$haploid, (asize$haploidAssemblySize-asize$haploidNCount)*2, asize$haploidAssemblySize-asize$haploidNCount)
## becaue flow cyt measueres in diploid
asize$doubledAssembly=asize$haploidAssemblySize*2
asize$doubledAssemblyNoN=(asize$haploidAssemblySize-asize$haploidNCount) *2

## there are 1.5 Gb of Ns in one of these assemblies screamemoji
#asize$nCountDoubled=ifelse(asize$haploid, asize$haploidNCount*2, asize$haploidNCount)
#asize$nCount=asize$rawNCount

## also check GC content because it's cool - these genomes are so big i get integer overflows adding?????
#asize$gc=(asize$V5/1e6+asize$V6/1e6)/(asize$V4/1e6+asize$V5/1e6+asize$V6/1e6+asize$V7/1e6)
## nm it's kinda boring 43-47% 

asize$species=taxonnames[match(asize$V2, names(taxonnames))]
asize$species=factor(asize$species, levels=taxonnames)

asize$speciesLabel=ifelse(asize$haploid, paste0(asize$species, '*'), as.character(asize$species))
asize$speciesLabel=asize$species
levels(asize$speciesLabel)[levels(asize$speciesLabel) %in% asize$species[asize$haploid]]=paste0(levels(asize$speciesLabel)[levels(asize$speciesLabel) %in% asize$species[asize$haploid]], '*')
asize$ploidy=gs$ploidy[match(asize$V2, gs$V2)]

## add flow, for supp
flow=read.table('panand_flow_cyt.txt', header=T, sep='\t') 
asize$flow=flow[,2][match(asize$V2, flow[,1])]
cor.test((asize$doubledAssembly)/2, asize$flow, use='complete.obs')
#cor.test((asize$doubledAssembly-asize$nCountDoubled)/2, asize$flow, use='complete.obs')


pdf(paste0('~/Downloads/supp_flow_assembly.', Sys.Date(), '.pdf'), 4,4)
## "haploid" assembly size
ggplot(asize, aes(x=doubledAssembly/1e9/2, y=flow/1000, color=ploidy)) +  scale_color_manual(values=ploidycolors)  + geom_point(size=3) + ylab('Genome Size, Flow Cytometry (Gb)')+ xlab('Haploid Assembly\nSize (Gb)') 
#ggplot(asize, aes(x=(doubledAssembly-nCountDoubled)/1e9/2, y=flow/1000, color=ploidy)) +  scale_color_manual(values=ploidycolors)  + geom_point() + ylab('Genome Size, Flow Cytometry (Gb)')+ xlab('Haploid Assembly\n Not N Size (Gb)') 
ggplot(asize, aes(x=doubledAssembly/1e9/2, y=flow/1000, color=ploidy)) +  scale_color_manual(values=ploidycolors)  + geom_text(aes(label=speciesLabel)) + ylab('Genome Size, Flow Cytometry (Gb)')+ xlab('Haploid Assembly\nSize (Gb)') 

## for talks, plot in gray first
ggplot(asize, aes(x=doubledAssembly/1e9/2, y=flow/1000)) +  scale_color_manual(values=ploidycolors)  + geom_point(size=3) + ylab('Genome Size, Flow Cytometry (Gb)')+ xlab('Haploid Assembly\nSize (Gb)') 
ggplot(asize, aes(x=doubledAssembly/1e9/2, y=flow/1000, color=ploidy)) +  scale_color_manual(values=ploidycolors)  + geom_point(size=3) + ylab('Genome Size, Flow Cytometry (Gb)')+ xlab('Haploid Assembly\nSize (Gb)')  + theme(legend.position='none')


dev.off()
