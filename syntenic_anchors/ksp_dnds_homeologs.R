library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())


ksp=fread('cat anchorAln_omega_v2/Pavag*.fasta')
head(ksp)
ksp$genome=substr(ksp$V3,1,6)                       
ksp$gene=str_split_fixed(ksp$V4, '_', 3)[,2]

kspall=ksp
ksp=ksp[substr(ksp$V3,1,6)=='pvagin' | substr(ksp$V4,1,6)=='pvagin',] ## to paspalum

taxonnames=c("Z. mays ssp. parviglumis TIL11", "Z. mays ssp. mays B73v5", "Z. mays ssp. parviglumis TIL01", "Z. mays ssp. mexicana TIL25", "Z. mays ssp. mexicana TIL18", "Z. mays ssp. huehuetengensis", 
"Z. luxurians", "Z. nicaraguensis", "Z. diploperennis Momo", "Z. diploperennis Gigi", "T. zoloptense", "T. dactyloides FL", "T. dactyloides Southern Hap2", 
"T. dactyloides Northern Hap2", "T. dactyloides KS", "T. dactyloides tetraploid", "U. digitatum", "V. cuspidata", "R. rottboellioides", "R. tuberculosa", 
"H. compressa", "E. tripsacoides", "S. scoparium", "S. microstachyum", "A. virginicum", "A. chinensis", "A. gerardi", 
"C. refractus", "C. citratus", "H. contortus", "T. triandra", "B. laguroides", "P. paniceum", "S. bicolor", 
"I. rugosum", "S. nutans", '"A." burmanicus', "T. elegans", "C. serrulatus", "P. vaginatum")
names(taxonnames)=c("zTIL11", "zmB735", "zTIL01", "zTIL25", "zTIL18", "zmhuet", 
"zluxur", "znicar", "zdmomo", "zdgigi", "tzopol", "tdacs1", "tdacs2", 
"tdacn2", "tdacn1", "tdactm", "udigit", "vcuspi", "rrottb", "rtuber", 
"hcompr", "etrips", "sscopa", "smicro", "avirgi", "achine", "agerar", 
"crefra", "ccitra", "hconto", "ttrian", "blagur", "ppanic", "sbicol", 
"irugos", "snutan", "atenui", "telega", "cserru", "pvagin")

ksp=ksp[ksp$genome %in% names(taxonnames),]

ksp=ksp%>% group_by(genome, gene) %>% mutate(is_min=V19==min(V19), is_max=V19==max(V19))


ksp %>% group_by(genome, gene) %>% summarize(maxmin=max(V19)-min(V19), V19=median(V19, na.rm=T)) %>% ungroup() %>% group_by(genome) %>% summarize(mean=mean(maxmin, na.rm=T), median=median(maxmin, na.rm=T), mediandnds=median(V19,na.rm=T)) %>% data.frame %>% dput()


## get spatial position along chromosomes
a$string=paste0(a$genome, '_', a$quickgene, '_', a$queryChr, '_', a$queryStart, '-', a$queryEnd)
ak=merge(a, ksp, by.x=c('genome', 'quickgene', 'string'), by.y=c('genome', 'gene', 'V3'), all.y=T)

## tyry for a more interesting genome!!! but there are diffs here for this one!!!
ak %>% group_by( refChr, queryChr) %>% filter(genome=='udigit') %>% summarize(median=median(V19, na.rm=T), mean=mean(V19, na.rm=T), n=n())  %>% filter(n>100) %>% print(n=39)


aks=ak %>% group_by( refChr, queryChr, genome) %>% summarize(median=median(V19, na.rm=T), mean=mean(V19, na.rm=T), n=n())

ksf=merge(sgfrac, aks[aks$refChr=='Chr01',], by=c('queryChr', 'genome'))


ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')

pdf('dnds_fractionation.pdf', 8,8)
ggplot(ksf[ksf$n>100,], aes(x=mya, y=mean, color=ploidy, label=genome, group=genome)) + geom_point(size=3) + geom_line() + geom_text()+ scale_color_manual(values=ploidycolors)
ggplot(ksf[ksf$n>100,], aes(x=mya, y=median, color=ploidy, label=genome, group=genome)) + geom_point(size=3) + geom_line() + geom_text()+ scale_color_manual(values=ploidycolors)
ggplot(ksf[ksf$n>100,], aes(y=medFrac, x=mean, color=ploidy, label=genome, group=genome)) + geom_point(size=3) + geom_line() + geom_text()+ scale_color_manual(values=ploidycolors)
ggplot(ksf[ksf$n>100,], aes(y=medFrac, x=median, color=ploidy, label=genome, group=genome)) + geom_point(size=3) + geom_line() + geom_text()+ scale_color_manual(values=ploidycolors)

ggplot(ksf[ksf$n>100,], aes(x=mya, y=mean, color=ploidy, label=genome, group=genome)) + geom_point(size=3) + geom_line() +  scale_color_manual(values=ploidycolors)
ggplot(ksf[ksf$n>100,], aes(x=mya, y=median, color=ploidy, label=genome, group=genome)) + geom_point(size=3) + geom_line() +  scale_color_manual(values=ploidycolors)
ggplot(ksf[ksf$n>100,], aes(y=medFrac, x=mean, color=ploidy, label=genome, group=genome)) + geom_point(size=3) + geom_line() + scale_color_manual(values=ploidycolors)
ggplot(ksf[ksf$n>100,], aes(y=medFrac, x=median, color=ploidy, label=genome, group=genome)) + geom_point(size=3) + geom_line() + scale_color_manual(values=ploidycolors)



dev.off()



## plotting on my computer from dput data frame 

pdf('dnds_along_pavag01.pdf',14,28)
ggplot(ksp[grepl('Pavag01',ksp$gene),], aes(x=gene, y=V19, color=ifelse(is_min, 'blue', ifelse(is_max, 'red', 'gray')))) + geom_point() + facet_wrap(~genome, ncol=1)+ylim(0,1) + scale_color_manual(values=c('blue', 'gray', 'red'))
ggplot(kspmm[grepl('Pavag01', kspmm$gene),], aes(x=gene, y=maxmin)) + geom_point() + facet_wrap(~genome, ncol=1) + ylim(0,1)
ggplot(kspmm[grepl('Pavag01', kspmm$gene),], aes(x=gene, y=maxmin)) + stat_smooth(se=F) + facet_wrap(~genome, ncol=1) + ylim(0,1)
ggplot(kspmm[grepl('Pavag01', kspmm$gene),], aes(x=gene, y=maxmin)) + stat_smooth(se=F) + facet_wrap(~genome, ncol=1) + ylim(0,0.1)
ggplot(kspmm[grepl('Pavag01', kspmm$gene),], aes(x=gene, y=maxmin)) + geom_point() + facet_wrap(~genome, ncol=1) + ylim(0,0.1)
dev.off()



