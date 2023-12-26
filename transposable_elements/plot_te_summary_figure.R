library(ggtree) ## /programs/R-4.2.1-r9/bin/R
library(ggplot2) 
library(cowplot)
theme_set(theme_cowplot())
library(rtracklayer)
library(reshape2)
library(RColorBrewer)
library(tidypaleo) ## facet species names in italics!!!!!


all=read.table('../panand_sp_ploidy.txt')
all=all[!all$V2 %in% c('pprate', 'tdactm', 'tzopol', 'osativ', 'bdista', 'agerjg', 'svirid', 'eophiu', 'tdacs2', 'tdacn2'),]

genomecountlist=vector(mode = "list", length = length(all$V2))
names(genomecountlist)=all$V2

repeatlengthslist=vector(mode = "list", length = length(all$V2))
names(repeatlengthslist)=all$V2


gfftypestokeep=c("Gypsy_LTR_retrotransposon", 
"LTR_retrotransposon", "Copia_LTR_retrotransposon", "CACTA_TIR_transposon", 
"helitron", "hAT_TIR_transposon", "PIF_Harbinger_TIR_transposon", 
"Tc1_Mariner_TIR_transposon", "Mutator_TIR_transposon", 'tandem_repeat') 

for( genotype in all$V2){
##import gff3
a=import.gff3(paste0('../trash/repeatmask_tandems/', all$V1[all$V2==genotype], '_EDTATandemRepeat.gff3'))
all[all$V2==genotype,gfftypestokeep]=sapply(gfftypestokeep, function (x) sum(width(reduce(a[a$type==x,]))))
all[all$V2==genotype,'tebp']=sum(width(reduce(a[a$type %in% gfftypestokeep & a$type!='tandem_repeat',])))

## get each fam
temp=data.frame(a) #%>% group_by(fam=gsub('_LTR', '', gsub('_INT', '', Name)), Classification) %>% dplyr::filter(!is.na(as.numeric(Identity))) 
if('Note' %in% colnames(temp)){temp$Note=''}
if(genotype=='zmB735'){temp$Parent=''
                      temp$ltr_identity=NA
                      temp=temp[,colnames(temp) %in% colnames(genomecountlist[[1]])]}
  if(genotype=='znicar'){ ## these guys edta is weird - chr not named with chr!!
                      temp$seqnames[temp$seqnames %in% 1:10]=paste0('chr', temp$seqnames[temp$seqnames %in% 1:10])}
temp$genome=genotype

genomecountlist[[genotype]]=temp
repeatlengthslist[[genotype]]=sum(width(reduce(a, ignore.strand=T)))
}

genomecount=do.call(rbind, genomecountlist)
genomecount$Identity=as.numeric(genomecount$Identity)

## add a superfamily based on matchign up to the classification field - note most of these NAs are relics from the B73 annotation being included!!!!!
genomecount$sup=c(NA, 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', 'DHH', NA,NA,NA,NA,NA,NA,'RLC', 'RLG', 'RLG', 'RLX', 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', NA,NA,NA)[match(genomecount$Classification, c("Cent/CentC", "DNA/DTA", "DNA/DTC", "DNA/DTH", "DNA/DTM", "DNA/DTT", 
"DNA/Helitron", "knob/knob180", "knob/TR-1", "LINE/L1", "LINE/RTE", 
"LINE/unknown", "Low_complexity", "LTR/Copia", "LTR/CRM", "LTR/Gypsy", 
"LTR/unknown", "MITE/DTA", "MITE/DTC", "MITE/DTH", "MITE/DTM", 
"MITE/DTT", "rDNA/spacer", "Simple_repeat", "subtelomere/4-12-1"))]

## tandem repeat maybe want to do by length class???
#genomecount$sup[genomecount$source=='TRASH']=genomecount$Name[genomecount$source=='TRASH']
genomecount$sup[genomecount$source=='TRASH']='TandemRepeat'



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


## get statistics for the paper!
gs=read.table('~/transfer/panand_assembly_sizes.txt', header=T, sep='\t')
cor.test(gs$haploidAssemblySize-gs$haploidNCount, gs$haploidRepeatSize)
gs$polyploid=gs$ploidy!='Diploid'

gs$species=taxonnames[match(gs$V2, names(taxonnames))]
gs$species=factor(gs$species, levels=taxonnames)
gs$speciesLabel=ifelse(gs$haploid, paste0(gs$species, '*'), as.character(gs$species))
gs$speciesLabel=gs$species
levels(gs$speciesLabel)[levels(gs$speciesLabel) %in% gs$species[gs$haploid]]=paste0(levels(gs$speciesLabel)[levels(gs$speciesLabel) %in% gs$species[gs$haploid]], '*')


am=melt(all[,-which(colnames(all)%in%c('V1', 'V3'))], id.vars=c('V2'))
am=am[!am$variable %in% c('target_site_duplication', 'repeat_region', 'long_terminal_repeat'),]
am$variable=as.character(am$variable)

amm=merge(am,gs, by='V2')
amm$doubledValue=amm$value
amm$doubledValue[amm$haploid]=amm$doubledValue[amm$haploid]*2
## okay get rid of all those stupid Ns in the denominator
amm$prop=amm$doubledValue/((amm$haploidAssemblySize-amm$haploidNCount)*2) ## proportion of genome for each, can change below in plotting code
                                            
oldnames=c("CACTA_TIR_transposon", "Copia_LTR_retrotransposon", "Gypsy_LTR_retrotransposon", 
"hAT_TIR_transposon", "helitron", "LTR_retrotransposon", "Mutator_TIR_transposon", 
"PIF_Harbinger_TIR_transposon", "tandem_repeat", "Tc1_Mariner_TIR_transposon")
names(oldnames)=c('DTC', 'RLC', 'RLG', 'DTA', 'DHH', 'RLX', 'DTM', 'DTH', 'TR', 'DTT')
amm$sup=factor(names(oldnames)[match(amm$variable, oldnames)], levels=c('DHH', 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', 'RLC', 'RLG', 'RLX', 'TR'))

ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')

                                            
pdf('~/transfer/te_sup_panand_bubble.pdf',12,8)
ggplot(amm[amm$variable!='tebp',], aes(x=prop, y=1, size=prop, color=ploidy)) +geom_hline(yintercept=1, linetype='dotted', alpha=0.5) + geom_vline(xintercept=0, linetype='dotted', color='snow3', alpha=0.5) + geom_point() + facet_grid(speciesLabel~sup, scales='free_x', space='free_x', switch='y', labeller=purrr::partial(label_species, dont_italicize=c('subsp.', 'ssp.', 'TIL11', 'TIL01', 'TIL25', 'TIL18', 'Momo', 'Gigi', 'Southern Hap1', 'Northern Hap1', 'FL', 'KS',  '\\*', '\\"')))+ scale_color_manual(values=ploidycolors)+   theme( strip.background = element_blank(),  panel.spacing = unit(3, "pt"), axis.text.y=element_blank(), axis.text.x=element_text(size=9)) + theme(legend.position='none', strip.text.y.left = element_text(angle=0), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.text.x = element_text(angle = 30)) + ylab('') + xlab('Genome Proportion') + scale_x_continuous(n.breaks = 3) 

ggplot(amm[amm$variable!='tebp',], aes(x=prop, y=1, size=prop, color=ploidy)) +geom_hline(yintercept=1, linetype='dotted', alpha=0.5) + geom_vline(xintercept=0, linetype='dotted', color='snow3', alpha=0.5) + geom_point() + facet_grid(speciesLabel~factor(sup, levels=c('RLG', 'RLC', 'RLX', 'DHH', 'DTM', 'DTC', 'DTT','DTA', 'DTH', 'TR')), scales='free_x', space='free_x', switch='y', labeller=purrr::partial(label_species, dont_italicize=c('subsp.', 'ssp.', 'TIL11', 'TIL01', 'TIL25', 'TIL18', 'Momo', 'Gigi', 'Southern Hap1', 'Northern Hap1', 'FL', 'KS',  '\\*', '\\"')))+ scale_color_manual(values=ploidycolors)+   theme( strip.background = element_blank(),  panel.spacing = unit(3, "pt"), axis.text.y=element_blank(), axis.text.x=element_text(size=9)) + theme(legend.position='none', strip.text.y.left = element_text(angle=0), strip.text.x = element_text(angle=90), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.text.x = element_text(angle = 30)) + ylab('') + xlab('Genome Proportion') + scale_x_continuous(n.breaks = 3) 
ggplot(amm[amm$variable!='tebp',], aes(x=prop, y=1, size=prop)) +geom_hline(yintercept=1, linetype='dotted',  aes(color=ploidy)) + geom_vline(xintercept=0, linetype='dotted', color='snow3', alpha=0.5) + geom_point() + facet_grid(speciesLabel~factor(sup, levels=c('RLG', 'RLC', 'RLX', 'DHH', 'DTM', 'DTC', 'DTT','DTA', 'DTH', 'TR')), scales='free_x', space='free_x', switch='y', labeller=purrr::partial(label_species, dont_italicize=c('subsp.', 'ssp.', 'TIL11', 'TIL01', 'TIL25', 'TIL18', 'Momo', 'Gigi', 'Southern Hap1', 'Northern Hap1', 'FL', 'KS',  '\\*', '\\"')))+ scale_color_manual(values=ploidycolors)+   theme( strip.background = element_blank(),  panel.spacing = unit(3, "pt"), axis.text.y=element_blank(), axis.text.x=element_text(size=9)) + theme(legend.position='none', strip.text.y.left = element_text(angle=0), strip.text.x = element_text(angle=90), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.text.x = element_text(angle = 30)) + ylab('') + xlab('Genome Proportion') + scale_x_continuous(n.breaks = 3) 

dev.off()



## ignore for now, or grab the plot from over there
# ## subtribes?? - highilight that the bottom bumpiness is subtribe?

# subtribenamesSimple=c("Tripsacinae", "Tripsacinae", "Tripsacinae", "Tripsacinae", "Tripsacinae", "Tripsacinae", 
# "Tripsacinae", "Tripsacinae", "Tripsacinae", "Tripsacinae", "Tripsacinae", "Tripsacinae", "Tripsacinae", 
# "Tripsacinae", "Tripsacinae", "Tripsacinae", "Rhytachninae", "Rhytachninae", "Rhytachninae", "Ratzeburgiinae", 
# "Ratzeburgiinae", "Ratzeburgiinae", "Andropogoninae", "Andropogoninae", "Andropogoninae", "Andropogoninae", "Andropogoninae", 
# "Anthistiriinae", "Anthistiriinae", "Anthistiriinae", "Anthistiriinae", "Anthistiriinae", "1", "2", 
# "3", "4", '5', "6", "7")#, "Paspalum vaginatum")
# ## 1 Germainiinae
# ## 2 Sorghinae
# ## 3 Ischaeminae
# ## 4 Apludinae
# ## 5 incertae sedis
# ## 6 incertae sedis
# ## 7 Chrysopogoninae
                                
# names(subtribenamesSimple)=c("zTIL11", "zmB735", "zTIL01", "zTIL25", "zTIL18", "zmhuet", 
# "zluxur", "znicar", "zdmomo", "zdgigi", "tzopol", "tdacs1", "tdacs2", 
# "tdacn2", "tdacn1", "tdactm", "udigit", "vcuspi", "rrottb", "rtuber", 
# "hcompr", "etrips", "sscopa", "smicro", "avirgi", "achine", "agerar", 
# "crefra", "ccitra", "hconto", "ttrian", "blagur", "ppanic", "sbicol", 
# "irugos", "snutan", "atenui", "telega", "cserru")#, "pvagin")

                                
# subtribedf=data.frame(subtribe=subtribenamesSimple, genome=names(subtribenamesSimple))
# subtribedf$haploid=gs$haploid[match(subtribedf$genome, gs$V2)]
# subtribedf$species=taxonnames[match(subtribedf$genome, names(taxonnames))]
# subtribedf$species=factor(subtribedf$species, levels=taxonnames)
# subtribedf$speciesLabel=ifelse(subtribedf$haploid, paste0(subtribedf$species, '*'), as.character(subtribedf$species))
# subtribedf$speciesLabel=subtribedf$species
# levels(subtribedf$speciesLabel)[levels(subtribedf$speciesLabel) %in% subtribedf$species[subtribedf$haploid]]=paste0(levels(subtribedf$speciesLabel)[levels(subtribedf$speciesLabel) %in% subtribedf$species[subtribedf$haploid]], '*')
# p2 <- tibble(ymin = c(1,which(!duplicated(subtribedf$subtribe[subtribedf$genome %in% asize$V2]))[-1])-1, ymax = c(which(!duplicated(subtribedf$subtribe[subtribedf$genome %in% asize$V2]))[-1],sum(subtribedf$genome %in% asize$V2)+1)-1, fill = unique(subtribedf$subtribe)) %>%  ggplot() +
#   geom_rect(aes(xmin = 0.1, xmax = 0.9, ymin = ymin+0.1, ymax = ymax), color='black', fill='snow2') +
#   geom_text(aes(x = .5, y = (ymin  + ymax) / 2, label = fill), angle = 90, size=2, fontface = "bold") +scale_y_reverse(breaks = seq(1, 10), expand = expansion(mult = c(0, 0))) +scale_x_continuous(breaks = c(0), expand = expansion(mult = c(0, 0))) +guides(fill = FALSE) +theme_void()

amm$ploidy=factor(amm$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))


## get gene info to compare to helitrons
g=read.csv('../supplement_table_annotations.csv')
g$six=c('atenui', 'achine', 'agerar', 'avirgi', 'blagur', 'ccitra', 'crefra', 'cserru', 'etrips', 'hcompr', 'hconto', 'irugos', 'ppanic', 'rrottb', 'rtuber', 'smicro', 'snutan', 'sscopa', 'tdacs1','tdacs2', 'tdacn1', 'tdacn2', 'telega', 'ttrian', 'udigit', 'vcuspi', 'zdgigi', 'zdmomo', 'zmhuet', 'znicar', 'zTIL01', 'zTIL11', 'zTIL18', 'zTIL25', 'zluxur')
gg=merge(amm, g, by.x='V2', by.y='six')
gg$genecount=as.numeric(gsub(',','',gg$Genes))
gg$meangenelength=as.numeric(gsub(',','',gg$Avg.Gene.length..bp.))
gg$totalcds=as.numeric(gsub(',','',gg$Total.CDS.region..bp.))
gg$meancdslength=as.numeric(gsub(',','',gg$Avg.CDS.length..bp.))
gg$meanintronlength=gg$meangenelength-gg$meancdslength
gg$ploidy=factor(gg$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))
## come on you have to adjust for haploid assemblies!@!
gg$doubledgenecount=gg$genecount
gg$doubledgenecount[gg$haploid]=gg$doubledgenecount[gg$haploid]*2
gg$doubledtotalcds=gg$totalcds
gg$doubledtotalcds[gg$haploid]=gg$doubledtotalcds[gg$haploid]*2
#### then reduce them back down in plot, because people are used to seeing haploid values

### gettin real sloppy here - this is from aum, which is from counting genes from gene table for fractionated/resistant genes
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
                                            
dtt=ggplot(amm[amm$sup=='DTT' & !is.na(amm$sup),], aes(x=ploidy, y=doubledValue/2, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=3e8) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('DTT base pairs') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dttfam=ggplot(nfam, aes(x=ploidy, y=DTT, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=24000) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('DTT Family Count') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
helitron=ggplot(amm[amm$sup=='DHH' & !is.na(amm$sup),], aes(x=ploidy, y=doubledValue/2, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=6e8) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('DHH base pairs') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
helgene=ggplot(gg[gg$sup=='DHH',], aes(x=doubledValue/2, y=doubledgenecount/2, color=ploidy)) +  geom_point() + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('DHH base pairs') + ylab('Total Haploid Gene Count') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
helsynt=ggplot(gg[gg$sup=='DHH',], aes(x=doubledValue/2, y=doubledsyntenic/2, color=ploidy)) +  geom_point() + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('DHH base pairs') + ylab('Total Haploid\nSyntenic Gene Count') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
helnonsynt=ggplot(gg[gg$sup=='DHH',], aes(x=doubledValue/2, y=(doubledgenecount-doubledsyntenic)/2, color=ploidy)) +  geom_point() + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('DHH base pairs') + ylab('Total Haploid\nNonsyntenic Gene Count') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

helgeneprop=ggplot(gg[gg$sup=='DHH',], aes(x=prop, y=doubledgenecount/2, color=ploidy)) +  geom_point() + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('DHH base pairs') + ylab('Total Haploid Gene Count') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
helsyntprop=ggplot(gg[gg$sup=='DHH',], aes(x=prop, y=doubledsyntenic/2, color=ploidy)) +  geom_point() + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('DHH base pairs') + ylab('Total Haploid\nSyntenic Gene Count') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

summary(lm(gg$doubledgenecount[gg$sup=='DHH']~gg$doubledValue[gg$sup=='DHH']))
summary(lm(gg$doubledgenecount[gg$sup=='DHH' & gg$ploidy!='Paleotetraploid']~gg$doubledValue[gg$sup=='DHH' & gg$ploidy!='Paleotetraploid']))
summary(lm(gg$doubledsyntenic[gg$sup=='DHH']~gg$doubledValue[gg$sup=='DHH']))
summary(lm(gg$doubledsyntenic[gg$sup=='DHH' & gg$ploidy!='Paleotetraploid']~gg$doubledValue[gg$sup=='DHH' & gg$ploidy!='Paleotetraploid']))
summary(lm(gg$doubledgenecount[gg$sup=='DHH']~gg$doubledValue[gg$sup=='DHH'] + gg$doubledsyntenic[gg$sup=='DHH']))


syntannot=ggplot(gg[gg$sup=='DHH',], aes(x=doubledgenecount/2, y=doubledsyntenic/2, color=ploidy)) +  geom_point() + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Total Haploid Gene Count') + ylab('Total Haploid\nSyntenic Gene Count') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

pdf(paste0('~/transfer/te_panand_fig5.', Sys.Date(), '.pdf'), 15,15)


superfams=ggplot(amm[amm$variable!='tebp',], aes(x=prop, y=1, size=prop, color=ploidy)) +geom_hline(yintercept=1, linetype='dotted', alpha=0.5) + geom_vline(xintercept=0, linetype='dotted', color='snow3', alpha=0.5) + geom_point() + facet_grid(speciesLabel~factor(sup, levels=c('RLG', 'RLC', 'RLX', 'DHH', 'DTM', 'DTC', 'DTT','DTA', 'DTH', 'TR')), scales='free_x', space='free_x', switch='y', labeller=purrr::partial(label_species, dont_italicize=c('subsp.', 'ssp.', 'TIL11', 'TIL01', 'TIL25', 'TIL18', 'Momo', 'Gigi', 'Southern Hap1', 'Northern Hap1', 'FL', 'KS',  '\\*', '\\"')))+ scale_color_manual(values=ploidycolors)+   theme( strip.background = element_blank(),  panel.spacing = unit(3, "pt"), axis.text.y=element_blank(), axis.text.x=element_text(size=9)) + theme(legend.position='none', strip.text.y.left = element_text(angle=0), strip.text.x = element_text(angle=90), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.text.x = element_text(angle = 30)) + ylab('') + xlab('Genome Proportion') + scale_x_continuous(n.breaks = 3) 

bottomlayer=plot_grid(dtt + theme(legend.position='NULL'), dttfam + theme(legend.position='NULL'), helitron + theme(legend.position='NULL'), helnonsynt + theme(legend.position='NULL') + geom_smooth(method = "lm", se = FALSE, alpha=0.3),
                      align='hv', nrow=1,rel_widths=c(1,1,1,1),labels=c('f', 'g', 'h', 'i'))
                                            
## these objects are all from plot_te_summaries_along_chr.R AAAHHHHHHH
plot_grid(plot_grid(bp1 + theme(legend.position='NULL'), pr2+ theme(legend.position='NULL'), fa3+ theme(legend.position='NULL'), fa4+ theme(legend.position='NULL'),
          legend, align='hv', nrow=1,rel_widths=c(1,1,1,1,0.4),labels=c('a', 'b', 'c', 'd', '')),
          superfams, bottomlayer, nrow=3, align='hv', rel_heights=c(0.5,1,0.5), labels=c('', 'e', ''))
                                                                                                                                                                                                                                                                                                                              
dev.off()
                                            

amm$ltrretro=amm$variable %in% c('Copia_LTR_retrotransposon', 'Gypsy_LTR_retrotransposon', 'LTR_retrotransposon')
amm %>% group_by(V2, ltrretro)%>% mutate(prop=doubledValue/(haploidAssemblySize*2)) %>% summarize(a=sum(prop)) %>% filter(ltrretro==T) %>% arrange(-a)
amm %>% group_by(V2, ltrretro)%>% mutate(prop=doubledValue/(haploidAssemblySize*2)) %>% summarize(a=sum(prop)) %>% filter(ltrretro==T) %>% arrange(a)

                                            
## do summaries like fig for each TE superfamily!!!
pdf(paste0('~/transfer/te_panand_superfambox.', Sys.Date(), '.pdf'), 15,12)
for(sup in unique(amm$variable)){

bp1=ggplot(amm[amm$variable==sup,], aes(x=ploidy, y=value, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=2e9) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('TE base pairs') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
pr2=ggplot(amm[amm$variable==sup,], aes(x=ploidy, y=prop*100, color=ploidy)) + geom_boxplot(outlier.shape=NA) +geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=50) + 
                                      scale_color_manual(values=ploidycolors) + xlab('Ploidy') + ylab('TE proportion')+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# fa3=ggplot(nfam, aes(x=ploidy, y=totFam, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=210000) + 
#                                       scale_color_manual(values=ploidycolors) + xlab('Ploidy') + ylab('TE family count')+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# fa4=ggplot(nfam100, aes(x=ploidy, y=totFam, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=23000) + 
#                                       scale_color_manual(values=ploidycolors) + xlab('Ploidy') + ylab('TE family count (>100 copies)')+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#ag5=ggplot(gs, aes(x=ploidy, y=meanage, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid",label.y=0.89) + 
#                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('LTR sequence identity') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

print(plot_grid(bp1 + theme(legend.position='NULL')  + ggtitle(sup), pr2+ theme(legend.position='NULL') + ggtitle(sup),
          legend, align='hv', nrow=1))
  }

## dtt beauty full glory
ggplot(amm[amm$variable=='Tc1_Mariner_TIR_transposon',], aes(x=ploidy, y=value, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=4e8) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('TE base pairs') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

dev.off()



pdf(paste0('~/transfer/te_panand_tc1mariner.', Sys.Date(), '.pdf'), 6,6)
ggplot(amm[amm$variable=='Tc1_Mariner_TIR_transposon',], aes(x=ploidy, y=value, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=4e8) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('TE base pairs') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle('DTT')
ggplot(amm[amm$variable=='PIF_Harbinger_TIR_transposon',], aes(x=ploidy, y=value, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=2e8) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('TE base pairs') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle('DTH')
ggplot(amm[amm$variable=='hAT_TIR_transposon',], aes(x=ploidy, y=value, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=1e8) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('TE base pairs') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle('DTA')
ggplot(amm[amm$variable=='Mutator_TIR_transposon',], aes(x=ploidy, y=value, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=1e8) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('TE base pairs') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle('DTM')
ggplot(amm[amm$variable=='helitron',], aes(x=ploidy, y=value, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=6e8) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('TE base pairs') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle('DHH')
ggplot(amm[amm$variable=='CACTA_TIR_transposon',], aes(x=ploidy, y=value, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=4e8) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('TE base pairs') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle('DTC')

dev.off()
    
