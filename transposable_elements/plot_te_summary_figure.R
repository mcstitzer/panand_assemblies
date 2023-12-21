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
cor.test(gs$haploidAssemblySize, gs$haploidRepeatSize)
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

oldnames=c("CACTA_TIR_transposon", "Copia_LTR_retrotransposon", "Gypsy_LTR_retrotransposon", 
"hAT_TIR_transposon", "helitron", "LTR_retrotransposon", "Mutator_TIR_transposon", 
"PIF_Harbinger_TIR_transposon", "tandem_repeat", "Tc1_Mariner_TIR_transposon")
names(oldnames)=c('DTC', 'RLC', 'RLG', 'DTA', 'DHH', 'RLX', 'DTM', 'DTH', 'TR', 'DTT')
amm$sup=factor(names(oldnames)[match(amm$variable, oldnames)], levels=c('DHH', 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', 'RLC', 'RLG', 'RLX', 'TR'))

ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')

                                            
pdf('~/transfer/te_sup_panand_bubble.pdf',12,8)
ggplot(amm[amm$variable!='tebp',], aes(x=doubledValue/(haploidAssemblySize*2), y=1, size=doubledValue/(haploidAssemblySize*2), color=ploidy)) +geom_hline(yintercept=1, linetype='dotted', alpha=0.5) + geom_vline(xintercept=0, linetype='dotted', color='snow3', alpha=0.5) + geom_point() + facet_grid(speciesLabel~sup, scales='free_x', space='free_x', switch='y', labeller=purrr::partial(label_species, dont_italicize=c('subsp.', 'ssp.', 'TIL11', 'TIL01', 'TIL25', 'TIL18', 'Momo', 'Gigi', 'Southern Hap1', 'Northern Hap1', 'FL', 'KS',  '\\*', '\\"')))+ scale_color_manual(values=ploidycolors)+   theme( strip.background = element_blank(),  panel.spacing = unit(3, "pt"), axis.text.y=element_blank(), axis.text.x=element_text(size=9)) + theme(legend.position='none', strip.text.y.left = element_text(angle=0), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.text.x = element_text(angle = 30)) + ylab('') + xlab('Genome Proportion') + scale_x_continuous(n.breaks = 3) 

ggplot(amm[amm$variable!='tebp',], aes(x=doubledValue/(haploidAssemblySize*2), y=1, size=doubledValue/(haploidAssemblySize*2), color=ploidy)) +geom_hline(yintercept=1, linetype='dotted', alpha=0.5) + geom_vline(xintercept=0, linetype='dotted', color='snow3', alpha=0.5) + geom_point() + facet_grid(speciesLabel~factor(sup, levels=c('RLG', 'RLC', 'RLX', 'DHH', 'DTM', 'DTC', 'DTT','DTA', 'DTH', 'TR')), scales='free_x', space='free_x', switch='y', labeller=purrr::partial(label_species, dont_italicize=c('subsp.', 'ssp.', 'TIL11', 'TIL01', 'TIL25', 'TIL18', 'Momo', 'Gigi', 'Southern Hap1', 'Northern Hap1', 'FL', 'KS',  '\\*', '\\"')))+ scale_color_manual(values=ploidycolors)+   theme( strip.background = element_blank(),  panel.spacing = unit(3, "pt"), axis.text.y=element_blank(), axis.text.x=element_text(size=9)) + theme(legend.position='none', strip.text.y.left = element_text(angle=0), strip.text.x = element_text(angle=90), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.text.x = element_text(angle = 30)) + ylab('') + xlab('Genome Proportion') + scale_x_continuous(n.breaks = 3) 
ggplot(amm[amm$variable!='tebp',], aes(x=doubledValue/(haploidAssemblySize*2), y=1, size=doubledValue/(haploidAssemblySize*2))) +geom_hline(yintercept=1, linetype='dotted',  aes(color=ploidy)) + geom_vline(xintercept=0, linetype='dotted', color='snow3', alpha=0.5) + geom_point() + facet_grid(speciesLabel~factor(sup, levels=c('RLG', 'RLC', 'RLX', 'DHH', 'DTM', 'DTC', 'DTT','DTA', 'DTH', 'TR')), scales='free_x', space='free_x', switch='y', labeller=purrr::partial(label_species, dont_italicize=c('subsp.', 'ssp.', 'TIL11', 'TIL01', 'TIL25', 'TIL18', 'Momo', 'Gigi', 'Southern Hap1', 'Northern Hap1', 'FL', 'KS',  '\\*', '\\"')))+ scale_color_manual(values=ploidycolors)+   theme( strip.background = element_blank(),  panel.spacing = unit(3, "pt"), axis.text.y=element_blank(), axis.text.x=element_text(size=9)) + theme(legend.position='none', strip.text.y.left = element_text(angle=0), strip.text.x = element_text(angle=90), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.text.x = element_text(angle = 30)) + ylab('') + xlab('Genome Proportion') + scale_x_continuous(n.breaks = 3) 

dev.off()


