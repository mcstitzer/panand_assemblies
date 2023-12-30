library(ggtree) ## /programs/R-4.2.1-r9/bin/R
library(ggplot2) 
library(cowplot)
theme_set(theme_cowplot())
library(rtracklayer)
library(reshape2)
library(RColorBrewer)
library(tidypaleo) ## facet species names in italics!!!!!
library(data.table)

gs=read.table('~/transfer/panand_assembly_sizes.txt', header=T, sep='\t')
## this drops sorghum, since it's not our annotation
g=read.csv('../supplement_table_annotations.csv')
g$six=c('atenui', 'achine', 'agerar', 'avirgi', 'blagur', 'ccitra', 'crefra', 'cserru', 'etrips', 'hcompr', 'hconto', 'irugos', 'ppanic', 'rrottb', 'rtuber', 'smicro', 'snutan', 'sscopa', 'tdacs1','tdacs2', 'tdacn1', 'tdacn2', 'telega', 'ttrian', 'udigit', 'vcuspi', 'zdgigi', 'zdmomo', 'zmhuet', 'znicar', 'zTIL01', 'zTIL11', 'zTIL18', 'zTIL25', 'zluxur')
gg=merge(gs, g, by.x='V2', by.y='six')
gg$genecount=as.numeric(gsub(',','',gg$Genes))
gg$meangenelength=as.numeric(gsub(',','',gg$Avg.Gene.length..bp.))
gg$totalcds=as.numeric(gsub(',','',gg$Total.CDS.region..bp.))
gg$meancdslength=as.numeric(gsub(',','',gg$Avg.CDS.length..bp.))
gg$meanintronlength=gg$meangenelength-gg$meancdslength


ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')


gg$ploidy=factor(gg$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))

## come on you have to adjust for haploid assemblies!@!
gg$doubledgenecount=gg$genecount
gg$doubledgenecount[gg$haploid]=gg$doubledgenecount[gg$haploid]*2
gg$doubledtotalcds=gg$totalcds
gg$doubledtotalcds[gg$haploid]=gg$doubledtotalcds[gg$haploid]*2
#### then reduce them back down in plot, because people are used to seeing haploid values

## genome size not correlated to cds length!!
cor.test(gg$haploidAssemblySize-gg$haploidNCount, gg$meancdslength)
## but yes to intron size!!
cor.test(gg$haploidAssemblySize-gg$haploidNCount, gg$meanintronlength)

## gettin real sloppy here - this is from aum, which is from counting genes from gene table for fractionated/resistant genes
#> dput(colSums(a[,-1]))
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
gg$synteniccount=syntcount[match(gg$V2, names(syntcount))]
gg$doubledsyntenic=gg$synteniccount
gg$doubledsyntenic[gg$haploid]=gg$doubledsyntenic[gg$haploid]*2



pdf(paste0('~/transfer/gene_ploidy.', Sys.Date(), '.pdf'), 6,6)
ggplot(gg, aes(x=ploidy, y=doubledgenecount/2, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=120000) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('Annotated Genes') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(gg, aes(x=ploidy, y=meangenelength, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=5000) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('Average Gene Length (bp)') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(gg, aes(x=ploidy, y=doubledtotalcds/1e6/2, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=210000000/1e6) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('Total CDS Length (Mbp)') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(gg, aes(x=ploidy, y=meancdslength, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=1300) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('Average CDS Length (bp)') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(gg, aes(x=ploidy, y=meangenelength-meancdslength, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=3700) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('Average Intron Length (bp)') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(gg, aes(x=meangenelength, y=meangenelength-meancdslength, color=ploidy))  + geom_point()     +      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Average Gene length (bp)') + ylab('Average Intron Length (bp)') 

dev.off()

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

ksp=fread('cat ../anchorAln_omega_v2/Pavag*.fasta')
ksp=ksp[substr(ksp$V3,1,6)=='pvagin' | substr(ksp$V4,1,6)=='pvagin',] ## to paspalum
ksp$genome=substr(ksp$V3,1,6)                       

ksp$haploid=gs$haploid[match(ksp$genome, gs$V2)]
ksp$species=taxonnames[match(ksp$genome, names(taxonnames))]
ksp$species=factor(ksp$species, levels=taxonnames)
ksp$speciesLabel=ifelse(ksp$haploid, paste0(ksp$species, '*'), as.character(ksp$species))
ksp$speciesLabel=ksp$species
levels(ksp$speciesLabel)[levels(ksp$speciesLabel) %in% ksp$species[ksp$haploid]]=paste0(levels(ksp$speciesLabel)[levels(ksp$speciesLabel) %in% ksp$species[ksp$haploid]], '*')
ksp$ploidy=gs$ploidy[match(ksp$genome, gs$V2)]
ksp=ksp[!ksp$genome %in% c('bdista', 'eophiu', 'osativ', 'svirid', 'tdacn2', 'tdacs2', 'tdactm', 'tzopol', 'pvagin'),]
ksp=ksp[ksp$V17<0.4,] ## scary is this a good decision i think it's okay, this is 30mya


## keep only dups within species
ksd=fread('cat ../anchorAln_omega_v2/Pavag*.fasta')
ksd=ksd[substr(ksd$V3,1,6)==substr(ksd$V4,1,6),] ## witin species
ksd$genome=substr(ksd$V3,1,6)                       

ksp$intraspecificdup=ksp$V3 %in% c(ksd$V3, ksd$V4)

## for gene text of paper
## numbers in text about proportion lost come from this - take mean tetraploid / (mean diploid * 2 ) or subtract
gg %>% group_by(ploidy) %>% summarize(mean(doubledsyntenic/2, na.rm=T), mean(doubledgenecount/2, na.rm=T))

numbergenes=ggplot(gg, aes(x=ploidy, y=doubledgenecount/2, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=120000) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('Annotated Genes') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

numbersyngenes=ggplot(gg, aes(x=ploidy, y=doubledsyntenic/2, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=50000) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('Syntenic Gene Regions') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
meanintronlength=ggplot(gg, aes(x=ploidy, y=meangenelength-meancdslength, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=3700) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('Average Intron Length (bp)') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

dnds=ggplot(ksp[!is.na(ksp$speciesLabel),], aes(x=speciesLabel, y=V19, color=ploidy, group=speciesLabel)) + geom_hline(yintercept=c(0,0.2,0.4,0.6,0.8,1), linetype='dotted', color='snow3', alpha=0.5)+ geom_hline(yintercept=c(0.1,0.3,0.5,0.7,0.9), linetype='dotted', color='snow2', alpha=0.5) + geom_violin(trim=F, adjust=0.4) + stat_summary(fun.data="mean_sdl", geom="pointrange") + scale_color_manual(values=ploidycolors)+   theme( axis.text.x=element_text(size=9)) + theme(legend.position='none',axis.text.x = element_text(angle = 30, hjust = 1, vjust = 0.5)) + ylab('dN/dS to Paspalum') + xlab('Species') + ylim(0,1) #+ geom_sina(alpha=0.1)
dndsIntra=ggplot(ksp[!is.na(ksp$speciesLabel) & ksp$intraspecificdup,], aes(x=speciesLabel, y=V19, color=ploidy, group=speciesLabel)) + geom_hline(yintercept=c(0,0.2,0.4,0.6,0.8,1), linetype='dotted', color='snow3', alpha=0.5)+ geom_hline(yintercept=c(0.1,0.3,0.5,0.7,0.9), linetype='dotted', color='snow2', alpha=0.5) + geom_violin(trim=F, adjust=0.4) + stat_summary(fun.data="mean_sdl", geom="pointrange") + scale_color_manual(values=ploidycolors)+   theme( axis.text.x=element_text(size=9)) + theme(legend.position='none',axis.text.x = element_text(angle = 30, hjust = 1, vjust = 0.5)) + ylab('dN/dS to Paspalum') + xlab('Species') + ylim(0,1) + scale_x_discrete(guide = guide_axis(angle = 30)) #+ geom_sina(alpha=0.1)

legend <- get_legend(numbergenes) 

pdf(paste0('~/transfer/gene_panand_fig3.', Sys.Date(), '.pdf'), 15,15)

                                        
## these objects are all from plot_te_summaries_along_chr.R AAAHHHHHHH
plot_grid(plot_grid(numbergenes + theme(legend.position='NULL'), numbersyngenes+ theme(legend.position='NULL'), meanintronlength+ theme(legend.position='NULL'),
          legend, align='hv', nrow=1,rel_widths=c(1,1,1,0.4),labels=c('a', 'b', 'c', '')),
          dndsIntra,nrow=2, align='hv', rel_heights=c(0.5,0.5), labels=c('', 'd'))
                                                                                                                                                                                                                                                                                                                              
dev.off()

## get dnds to paspalum for each gene duplicate
ksd$dndsp1=ksp$V19[match(ksd$V3, ksp$V3)]
ksd$dndsp2=ksp$V19[match(ksd$V4, ksp$V3)]
ksd$deltadnds=abs(ksd$dndsp1-ksd$dndsp2)

cor.test(ksd$deltadnds, ksd$V17, use='complete')

ksd$genome=substr(ksd$V3,1,6)                       

ksd$haploid=gs$haploid[match(ksd$genome, gs$V2)]
ksd$species=taxonnames[match(ksd$genome, names(taxonnames))]
ksd$species=factor(ksd$species, levels=taxonnames)
ksd$speciesLabel=ifelse(ksd$haploid, paste0(ksd$species, '*'), as.character(ksd$species))
ksd$speciesLabel=ksd$species
levels(ksd$speciesLabel)[levels(ksd$speciesLabel) %in% ksd$species[ksd$haploid]]=paste0(levels(ksd$speciesLabel)[levels(ksd$speciesLabel) %in% ksd$species[ksd$haploid]], '*')
ksd$ploidy=gs$ploidy[match(ksd$genome, gs$V2)]
ksd=ksd[!ksd$genome %in% c('bdista', 'eophiu', 'osativ', 'svirid', 'tdacn2', 'tdacs2', 'tdactm', 'tzopol', 'pvagin'),]
ksd=ksd[ksd$V17<0.3,] ## scary is this a good decision i think it's okay, it gets rid of grass wgd



png(paste0('~/transfer/dnds_vs_ks.', Sys.Date(), '.png'), 15,15, unit='in', res=300)

                                        
ggplot(ksd[ksd$V17>0.001 & (ksd$dndsp1<1|ksd$dndsp2<1),], aes(x=V17, y=deltadnds, color=ploidy)) +geom_point(alpha=0.2)+ scale_color_manual(values=ploidycolors)+   theme( axis.text.x=element_text(size=9)) + theme(legend.position='none',axis.text.x = element_text(angle = 30, hjust = 1, vjust = 0.5)) + ylab('Delta dN/dS to Paspalum') + xlab('Ks to Intraspecific Duplicate') + ylim(0,1)  #+ geom_sina(alpha=0.1)
                                                                                                                                                                                                                                                                                                                        
dev.off()
png(paste0('~/transfer/dnds_vs_ks.facet.', Sys.Date(), '.png'), 15,15, unit='in', res=300)

ggplot(ksd[ksd$V17>0.001 & (ksd$dndsp1<1|ksd$dndsp2<1),], aes(x=V17, y=deltadnds, color=ploidy)) +geom_point(alpha=0.2)+ scale_color_manual(values=ploidycolors)+  facet_wrap(~genome) + theme( axis.text.x=element_text(size=9)) + theme(legend.position='none',axis.text.x = element_text(angle = 30, hjust = 1, vjust = 0.5)) + ylab('Delta dN/dS to Paspalum') + xlab('Ks to Intraspecific Duplicate') + ylim(0,1) + stat_smooth(method='lm', se=F) #+ geom_sina(alpha=0.1)
dev.off()


pdf(paste0('~/transfer/dnds_vs_ks.', Sys.Date(), '.pdf'), 15,15)
ggplot(ksd[ksd$V17>0.001 & (ksd$dndsp1<1|ksd$dndsp2<1) & !is.na(ksd$V17) & !is.na(ksd$deltadnds) & !is.na(ksd$ploidy),] %>% group_by(genome, ploidy) %>% summarize(V17=median(V17), upperks=quantile(V17,0.75), lowerks=quantile(V17,0.25), deltadnds=median(deltadnds), upperdnds=quantile(deltadnds,0.75), lowerdnds=quantile(deltadnds,0.25)), aes(x=V17, y=deltadnds, color=ploidy, ymin=deltadnds-lowerdnds, ymax=deltadnds+upperdnds, xmin=V17-lowerks, xmax=V17+upperks )) +geom_point(alpha=0.2)+ geom_pointrange() + geom_errorbarh(height=0) + scale_color_manual(values=ploidycolors) + geom_errorbarh()
dev.off()
