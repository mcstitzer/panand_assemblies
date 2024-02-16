library(rtracklayer)
library(stringr)
library(plyranges)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggpubr)
theme_set(theme_cowplot())

all=read.table('../panand_sp_ploidy.txt')
all=all[!all$V2 %in% c('pprate', 'tdactm', 'tzopol', 'osativ', 'bdista', 'agerjg', 'svirid', 'eophiu'),]

yw=lapply(all$V2, function(x) {
 a=import.gff(paste0(x, '.ywminiprot.gff'))
 a=a[a$type=='mRNA',]
  ## eventually deal with intersections (rank by score)
# a=
 a$genome=x
 return(a)
 })


at=do.call(c, yw)
at$name=str_split_fixed(at$Target, ' ', 3)[,1]


tegroup=read.table('yuan_wessler_2011.supnames.txt', header=T)
at$sup=tegroup$superfamily[match(at$name, tegroup$fasta_entry)]

te=data.frame(at) %>% dplyr::group_by(sup, genome) %>% dplyr::summarize(n=n()) %>% pivot_wider(names_from=sup, values_from=n, values_fill=0)


gs=read.table('panand_assembly_sizes.txt', header=T, sep='\t')

teg=merge(gs, te, by.x='V2', by.y='genome')

ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')

syntcount=c(achine = 30037, agerar = 62404, avirgi = 14549, blagur = 79949, 
ccitra = 52037, crefra = 14712, cserru = 23638, etrips = 26922, 
hcompr = 65110, hconto = 54530, irugos = 16094, ppanic = 14777, 
rrottb = 27766, sbicol = 15053, smicro = 14200, sscopa = 42060, 
tdacn1 = 18980, tdacn2 = 18784, tdacs1 = 19151, tdacs2 = 19130, 
telega = 28664, ttrian = 23530, udigit = 40017, vcuspi = 41068, 
zTIL01 = 17394, zTIL11 = 17057, zTIL18 = 17311, zTIL25 = 17338, 
zdgigi = 22451, zdmomo = 18246, zluxur = 17052, zmB735 = 17032, 
zmhuet = 17500, znicar = 18980, atenui = 49008, rtuber = 14986, 
snutan = 35298) ## don't worry about b73 name, since we don't want to compare to the b73 gene annotation anyways
teg$synteniccount=syntcount[match(teg$V2, names(syntcount))]
teg$doubledsyntenic=teg$synteniccount
teg$doubledsyntenic[teg$haploid]=teg$doubledsyntenic[teg$haploid]*2

teg$doubledCacta=teg$CMC
teg$doubledCacta[teg$haploid]=teg$doubledCacta[teg$haploid]*2
teg$doubledMutator=teg$MULE
teg$doubledMutator[teg$haploid]=teg$doubledMutator[teg$haploid]*2
teg$doubledPifHarbinger=teg$`PIF/Harbinger`
teg$doubledPifHarbinger[teg$haploid]=teg$doubledPifHarbinger[teg$haploid]*2
teg$doubledTc1Mariner=teg$`Tc1/mariner`
teg$doubledTc1Mariner[teg$haploid]=teg$doubledTc1Mariner[teg$haploid]*2
teg$doubledHat=teg$hAT
teg$doubledHat[teg$haploid]=teg$doubledHat[teg$haploid]*2


teg %>% group_by(ploidy) %>% summarize_at(vars(doubledsyntenic:doubledHat), median) %>% data.frame
teg$ploidy=factor(teg$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))

## compare to earl grey
fil=Sys.glob('../earlgrey/earlGreyOutputs/*/*/*.high*.txt')
teg$eg=NA
teg$egtot=NA
teg$ltr=NA
teg$dna=NA
teg$line=NA
teg$satellite=NA
teg$rc=NA
teg$unclass=NA
teg$sine=NA
for(i in 1:length(fil)){
 b=read.table(fil[i], sep='\t', header=T)
teg$eg[teg$V2==substr(fil[i],29,34)]=sum(b$proportion)
teg$egtot[teg$V2==substr(fil[i],29,34)]=sum(b$cov)
teg$ltr[teg$V2==substr(fil[i],29,34)]=b$cov[b$tclassif=='LTR']
teg$dna[teg$V2==substr(fil[i],29,34)]=b$cov[b$tclassif=='DNA']
teg$line[teg$V2==substr(fil[i],29,34)]=b$cov[b$tclassif=='LINE']
teg$satellite[teg$V2==substr(fil[i],29,34)]=b$cov[b$tclassif=='Other (Simple Repeat, Microsatellite, RNA)']
teg$rc[teg$V2==substr(fil[i],29,34)]=b$cov[b$tclassif=='Rolling Circle']
teg$unclass[teg$V2==substr(fil[i],29,34)]=b$cov[b$tclassif=='Unclassified']
teg$sine[teg$V2==substr(fil[i],29,34)]=ifelse('SINE' %in% b$tclassif, b$cov[b$tclassif=='SINE'],0)
}

teg$doubledegtot=teg$egtot
teg$doubledegtot[teg$haploid]=teg$doubledegtot[teg$haploid]*2

teg$doubledltr=teg$ltr
teg$doubledltr[teg$haploid]=teg$doubledltr[teg$haploid]*2
teg$doubleddna=teg$dna
teg$doubleddna[teg$haploid]=teg$doubleddna[teg$haploid]*2
teg$doubledline=teg$line
teg$doubledline[teg$haploid]=teg$doubledline[teg$haploid]*2
teg$doubledsatellite=teg$satellite
teg$doubledsatellite[teg$haploid]=teg$doubledsatellite[teg$haploid]*2
teg$doubledrc=teg$rc
teg$doubledrc[teg$haploid]=teg$doubledrc[teg$haploid]*2
teg$doubledunclass=teg$unclass
teg$doubledunclass[teg$haploid]=teg$doubledunclass[teg$haploid]*2
teg$doubledsine=teg$sine
teg$doubledsine[teg$haploid]=teg$doubledsine[teg$haploid]*2




pdf('te_autonomous_counts.pdf',8,8)
#for(i in c('doubledsyntenic', 'doubledCacta', 'doubledMutator', 'doubledPifHarbinger', 'doubledTc1Mariner', 'doubledHat')){


cacta= ggplot(teg, aes(x=ploidy, y=doubledCacta/2, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=110) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('CACTA Protein Copies') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
mutator= ggplot(teg, aes(x=ploidy, y=doubledMutator/2, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=1500) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('Mutator Protein Copies') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
pif= ggplot(teg, aes(x=ploidy, y=doubledPifHarbinger/2, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=550) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('Pif/Harbinger Protein Copies') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
tc1= ggplot(teg, aes(x=ploidy, y=doubledTc1Mariner/2, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=250) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('Tc1/Mariner Protein Copies') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
hat= ggplot(teg, aes(x=ploidy, y=doubledHat/2, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=260) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('hAT Protein Copies') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
cacta
mutator
pif
tc1
hat

cacta + geom_text(aes(label=V2))
mutator+ geom_text(aes(label=V2))
pif+ geom_text(aes(label=V2))
tc1+ geom_text(aes(label=V2))
hat+ geom_text(aes(label=V2))

ggplot(teg, aes(x=ploidy, y=(haploidAssemblySize-haploidNCount)/1e6, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=6000) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('Haploid Assembly Size (Gbp)') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(teg, aes(x=ploidy, y=doubledsyntenic/2, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=41000) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('Haploid Syntenic Copy Number') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(teg, aes(x=ploidy, y=(haploidAssemblySize-haploidNCount)/1e6, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=6000) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('Haploid Assembly Size (Gbp)') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ geom_text(aes(label=V2))
ggplot(teg, aes(x=ploidy, y=doubledsyntenic/2, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=41000) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('Haploid Syntenic Copy Number') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ geom_text(aes(label=V2))


ggplot(teg, aes(x=(haploidAssemblySize-haploidNCount)/1e6, y=doubledCacta/2, color=ploidy)) + geom_point()  + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Haploid Assembly Size (Gbp)') + ylab('CACTA Protein Copies') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(teg, aes(x=(haploidAssemblySize-haploidNCount)/1e6, y=doubledMutator/2, color=ploidy)) + geom_point()  + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Haploid Assembly Size (Gbp)') + ylab('Mutator Protein Copies') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(teg, aes(x=(haploidAssemblySize-haploidNCount)/1e6, y=doubledPifHarbinger/2, color=ploidy)) + geom_point()  + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Haploid Assembly Size (Gbp)') + ylab('Pif/Harbinger Protein Copies') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(teg, aes(x=(haploidAssemblySize-haploidNCount)/1e6, y=doubledTc1Mariner/2, color=ploidy)) + geom_point()  + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Haploid Assembly Size (Gbp)') + ylab('Tc1/Mariner Protein Copies') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(teg, aes(x=(haploidAssemblySize-haploidNCount)/1e6, y=doubledHat/2, color=ploidy)) + geom_point() + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Haploid Assembly Size (Gbp)') + ylab('hAT Protein Copies') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
## earlgrey vs edta
## this isn't right, but to be fair put N's into edta denominator
# ggplot(teg, aes(x=haploidRepeatSize/(haploidAssemblySize-haploidNCount), y=eg, color=ploidy)) + geom_point() + 
#                                       scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('EDTA Repeat Proportion') + ylab('EarlGrey/RepeatModeler2 Repeat Proportion') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
## this is right, recalculated it
ggplot(teg, aes(x=haploidRepeatSize/(haploidAssemblySize-haploidNCount), y=(doubledegtot/2)/(haploidAssemblySize-haploidNCount), color=ploidy)) + geom_point() + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('EDTA Repeat Proportion') + ylab('EarlGrey/RepeatModeler2 Repeat Proportion') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggplot(teg, aes(x=ploidy, y=(doubledegtot/2)/(haploidAssemblySize-haploidNCount), color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=1) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('EarlGrey/RepeatModeler2 Repeat Proportion') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(teg, aes(x=ploidy, y=(doubledegtot/2)/(haploidAssemblySize-haploidNCount), color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=1) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('EarlGrey/RepeatModeler2 Repeat Proportion') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + geom_text(aes(label=V2))

 dev.off()
