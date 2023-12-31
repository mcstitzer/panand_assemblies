library(rtracklayer)
library(dplyr)
library(data.table)
library(stringr)
library(plyranges)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

all=read.table('../panand_sp_ploidy.txt')
all=all[!all$V2 %in% c('pprate', 'tdactm', 'tzopol', 'osativ', 'bdista', 'agerjg', 'svirid', 'eophiu'),]

genomecountlist=vector(mode = "list", length = length(all$V2))
names(genomecountlist)=all$V2

repeatlengthslist=vector(mode = "list", length = length(all$V2))
names(repeatlengthslist)=all$V2


for( genotype in all$V2){
##import gff3
a=import.gff3(paste0('../trash/repeatmask_tandems/', all$V1[all$V2==genotype], '_EDTATandemRepeat.gff3'))
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


## get the largest LTR retro families by genome - structural copies

genomecount$collapsedfam=gsub('_INT','',gsub('_LTR', '', genomecount$Name))
bigfam=genomecount %>% filter(sup%in%c('RLC', 'RLG', 'RLX') & Method=='structural' & type=='LTR_retrotransposon') %>%group_by(genome, collapsedfam, sup) %>% summarize(ncopies=n(), meanID=mean(as.numeric(ltr_identity), na.rm=T)) # %>% slice_max(ncopies, n=10)
bigfam=bigfam%>% group_by(genome) %>%  mutate(rank=order(order(ncopies, decreasing=T)))
### wtaf there are no structural LTR retros annotated to superfamily in these fucking genomes
## what am i supposed to do with this
#### how the fuck is this doing anything?

genomecount$haploid=gs$haploid[match(genomecount$genome, gs$V2)]
genomecount$species=taxonnames[match(genomecount$genome, names(taxonnames))]
genomecount$species=factor(genomecount$species, levels=taxonnames)
genomecount$speciesLabel=ifelse(genomecount$haploid, paste0(genomecount$species, '*'), as.character(genomecount$species))
genomecount$speciesLabel=genomecount$species
levels(genomecount$speciesLabel)[levels(genomecount$speciesLabel) %in% genomecount$species[genomecount$haploid]]=paste0(levels(genomecount$speciesLabel)[levels(genomecount$speciesLabel) %in% genomecount$species[genomecount$haploid]], '*')

genomecount$ploidy=gs$ploidy[match(genomecount$genome, gs$V2)]
ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')
genomecount$ploidy=factor(genomecount$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))

genomecount$rank=bigfam$rank[match(paste(genomecount$genome, genomecount$collapsedfam), paste(bigfam$genome, bigfam$collapsedfam))]



## ugh get rid of these stupid repeated entries for ltr retros - should be doing this for plotting etc
genomecount=genomecount[genomecount$type %in% c("tandem_repeat", "Gypsy_LTR_retrotransposon", "LTR_retrotransposon", "Copia_LTR_retrotransposon", "CACTA_TIR_transposon", "helitron", 
"hAT_TIR_transposon", "PIF_Harbinger_TIR_transposon", "Tc1_Mariner_TIR_transposon", 
"Mutator_TIR_transposon", "L1_LINE_retrotransposon",
"RTE_LINE_retrotransposon", "LINE_element"),]

## sup assignemnt is there in homologgy copies!
sups=genomecount %>% filter(Method=='homology') %>% group_by(genome, collapsedfam, sup) %>% summarize(n=n())
genomecount$realsup=genomecount$sup
genomecount$realsup[genomecount$Method=='structural' & genomecount$type=='LTR_retrotransposon']=sups$sup[match(paste(genomecount$genome, genomecount$collapsedfam)[genomecount$Method=='structural' & genomecount$type=='LTR_retrotransposon'], paste(sups$genome, sups$collapsedfam))]

## want to put on reverse scale, match to ks distribution
## ks in black, each te sup in color
#### also want plot of relative size of each of top 10 families
bigfams=bigfam %>% filter(rank<11) %>% group_by(genome) %>% mutate(prop=ncopies/sum(ncopies))
bigfams$haploid=gs$haploid[match(bigfams$genome, gs$V2)]
bigfams$species=taxonnames[match(bigfams$genome, names(taxonnames))]
bigfams$species=factor(bigfams$species, levels=taxonnames)
bigfams$speciesLabel=ifelse(bigfams$haploid, paste0(bigfams$species, '*'), as.character(bigfams$species))
bigfams$speciesLabel=bigfams$species
levels(bigfams$speciesLabel)[levels(bigfams$speciesLabel) %in% bigfams$species[bigfams$haploid]]=paste0(levels(bigfams$speciesLabel)[levels(bigfams$speciesLabel) %in% bigfams$species[bigfams$haploid]], '*')
bigfams$ploidy=gs$ploidy[match(bigfams$genome, gs$V2)]
bigfams$ploidy=factor(bigfams$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))
bigfams$rank=factor(bigfams$rank,levels=rev(1:10))

mks=read.table('~/transfer/median_ks.txt', header=T)
genomecount$medianks=mks$median[match(genomecount$genome, mks$genome)]


## okay, let's get bp of the top 100? 20? families in each genome, nomatter what they are!!
bigfamall=genomecount %>% filter(sup!='TandemRepeat') %>% group_by(genome, collapsedfam, sup) %>% summarize(ncopies=n(), bp=sum(width))
bigfamall=bigfamall[!(bigfamall$genome=='znicar' & bigfamall$collapsedfam=='TE_00015576'),] # this is still knob contamination like osed and sela in MTEC :( the assembly is bad enough it didn't pass tandem repeat finding wth trash... # %>% slice_max(ncopies, n=10)
bigfamall=bigfamall%>% group_by(genome) %>%  mutate(cnrank=order(order(ncopies, decreasing=T)), bprank=order(order(bp, decreasing=T))) 
bigfamall$haploid=gs$haploid[match(bigfamall$genome, gs$V2)]
bigfamall$species=taxonnames[match(bigfamall$genome, names(taxonnames))]
bigfamall$species=factor(bigfamall$species, levels=taxonnames)
bigfamall$speciesLabel=ifelse(bigfamall$haploid, paste0(bigfamall$species, '*'), as.character(bigfamall$species))
bigfamall$speciesLabel=bigfamall$species
levels(bigfamall$speciesLabel)[levels(bigfamall$speciesLabel) %in% bigfamall$species[bigfamall$haploid]]=paste0(levels(bigfamall$speciesLabel)[levels(bigfamall$speciesLabel) %in% bigfamall$species[bigfamall$haploid]], '*')
bigfamall$ploidy=gs$ploidy[match(bigfamall$genome, gs$V2)]
bigfamall$ploidy=factor(bigfamall$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))
totrepeat=read.table('../transposable_elements/total_repeat_bp.txt', header=T)
bigfamall$totrepeat=totrepeat$repeatbp[match(bigfamall$genome, totrepeat$genome)]


bigfamsallcn=bigfamall %>% filter(cnrank<31) %>% group_by(genome) %>% mutate(prop=ncopies/sum(ncopies), n=n())
bigfamsallbp=bigfamall %>% filter(bprank<31) %>% group_by(genome) %>% arrange(-bp) %>% mutate(prop=bp/sum(bp), n=n())
bigfamsallbpfilt=bigfamall %>% filter(bprank<501) %>% group_by(genome) %>% arrange(-bp) %>% mutate(cumbp=cumsum(as.numeric(bp)), cumprop=cumbp/totrepeat) %>% filter(cumprop<0.15) %>% mutate(prop=bp/sum(bp), n=n())
## so keep fams that get you to 10%? 30% of the genomic repeats by bp
bigfamsallcn$cnrank=factor(bigfamsallcn$cnrank,levels=rev(1:31))
bigfamsallbp$bprank=factor(bigfamsallbp$bprank,levels=rev(1:31))
bigfamsallbpfilt$bprank=factor(bigfamsallbpfilt$bprank,levels=rev(1:max(bigfamsallbpfilt$bprank)))
#levels(bigfamsallbpfilt$speciesLabel)=rev(levels(bigfamsallbpfilt$speciesLabel))

keep=c("zTIL11", "zTIL01", "zTIL25", "zTIL18", "zmhuet", 
"zluxur", "znicar", "zdmomo", "zdgigi", "tdacs1",  
"tdacn1",  "udigit", "vcuspi", "rrottb", "rtuber", 
"hcompr", "etrips", "sscopa", "smicro", "avirgi", "achine", "agerar", 
"crefra", "ccitra", "hconto", "ttrian", "blagur", "ppanic", "sbicol", 
"irugos", "snutan", "atenui", "telega", "cserru")

tea=ggplot(genomecount[genomecount$genome%in% keep &genomecount$Method=='structural' & genomecount$type=='LTR_retrotransposon' & paste(genomecount$genome, genomecount$collapsedfam)%in% paste(bigfam$genome, bigfam$collapsedfam)[bigfam$rank %in% 1:10],], 
           aes(x=1-as.numeric(ltr_identity), group=rank, alpha=rank,color=ploidy))  + geom_vline(xintercept=1-c(0.85,0.95), color='snow2', linetype='dotted') + geom_vline(xintercept=1-c(0.8,0.9), color='snow3', linetype='dotted') + geom_density() + scale_color_manual(values=ploidycolors) +  facet_wrap(~speciesLabel, ncol=1, strip.position='left', scales='free_y') + theme( strip.background = element_blank(), strip.text.y.left = element_text(angle=0), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + ylab('')+ xlab('LTR Divergence (Gene Duplicate\nDivergence in Gray Bar)') + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +xlim(0,0.2) + geom_vline(aes(xintercept=medianks), color='dimgray', alpha=0.5)


pdf('~/transfer/te_age_try2.pdf', 10,14)
tea

ggplot(bigfams, aes(x=speciesLabel, y=prop)) + geom_bar(width=0.7,stat='summary', fun='sum',aes(color=ploidy), fill=NA)+ geom_col(width=0.7, aes(fill=rank))  + coord_flip() + scale_color_manual(values=ploidycolors)  + theme( strip.background = element_blank(), strip.text.y.left = element_text(angle=0), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + ylab('')+ xlab('LTR Divergence') + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + scale_fill_manual(values= hcl.colors(n=10, palette='Mint')) # + geom_vline(aes(xintercept=mtm1), color='dimgray', alpha=0.5)

ggplot(bigfamsallcn, aes(x=speciesLabel, y=prop))+ geom_bar(width=0.7,stat='summary', fun='sum',aes(color=ploidy), fill=NA) + geom_col(width=0.7, aes(fill=cnrank)) + coord_flip() + scale_color_manual(values=ploidycolors)  + theme( strip.background = element_blank(), strip.text.y.left = element_text(angle=0), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + ylab('')+ xlab('LTR Divergence') + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + scale_fill_manual(values= hcl.colors(n=31, palette='Mako'))+ geom_text(aes(label=n),  colour = "black", y=-0.03) # + geom_vline(aes(xintercept=mtm1), color='dimgray', alpha=0.5)
ggplot(bigfamsallbp, aes(x=speciesLabel, y=prop))+ geom_bar(width=0.7,stat='summary', fun='sum',aes(color=ploidy), fill=NA) + geom_col(width=0.7, aes(fill=bprank)) + coord_flip() + scale_color_manual(values=ploidycolors)  + theme( strip.background = element_blank(), strip.text.y.left = element_text(angle=0), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + ylab('')+ xlab('LTR Divergence') + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + scale_fill_manual(values= hcl.colors(n=31, palette='Mako')) + geom_text(aes(label=n), colour = "black", y=-0.03)# + geom_vline(aes(xintercept=mtm1), color='dimgray', alpha=0.5)
ggplot(bigfamsallbpfilt, aes(x=speciesLabel, y=prop))+ geom_bar(width=0.7,stat='summary', fun='sum',aes(color=ploidy), fill=NA) + geom_col(width=0.7, aes(fill=bprank)) + coord_flip() + scale_color_manual(values=ploidycolors)  + theme( strip.background = element_blank(), strip.text.y.left = element_text(angle=0), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + xlab('')+ ylab('Proportion abundance of TE families\nin top 15% of repeat space') + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + scale_fill_manual(values= hcl.colors(n=max(as.numeric(bigfamsallbpfilt$n)), palette='Mako')) + geom_text(aes(label=n), colour = "black", y=-0.03)  # + geom_vline(aes(xintercept=mtm1), color='dimgray', alpha=0.5)


plot_grid(
## superfamilies are  mess, ignore theme!
# tea=ggplot(genomecount[genomecount$Method=='structural' & genomecount$type=='LTR_retrotransposon' & paste(genomecount$genome, genomecount$collapsedfam)%in% paste(bigfam$genome, bigfam$collapsedfam)[bigfam$rank %in% 1:10],], 
#            aes(x=as.numeric(ltr_identity), group=rank, fill=realsup, color=ploidy))  + geom_vline(xintercept=c(0.75,0.85,0.95), color='snow2', linetype='dotted') + geom_vline(xintercept=c(0.7,0.8,0.9), color='snow3', linetype='dotted') + geom_density(alpha=0.1) + scale_color_manual(values=ploidycolors) + scale_fill_manual(values=ploidycolors)+  facet_wrap(~speciesLabel, ncol=1, strip.position='left', scales='free_y', drop=F) + theme( strip.background = element_blank(), strip.text.y.left = element_text(angle=0), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + ylab('')+ xlab('LTR Divergence') + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +xlim(0.9,1)  # + geom_vline(aes(xintercept=mtm1), color='dimgray', alpha=0.5)
# tea

dev.off()

pdf('~/transfer/te_age_try3.pdf', 10,14)

famabund=ggplot(bigfamsallbpfilt[bigfamsallbpfilt$genome %in% keep,], aes(x=fct_rev(speciesLabel), y=prop))+ geom_bar(width=0.7,stat='summary', fun='sum',aes(color=ploidy), fill=NA) + geom_col(width=0.7, aes(fill=bprank)) + coord_flip() + scale_color_manual(values=ploidycolors)  + theme( strip.background = element_blank(), strip.text.y.left = element_text(angle=0), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + xlab('')+ ylab('Proportion abundance of TE families\nin top 15% of repeat space') + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + scale_fill_manual(values= hcl.colors(n=max(as.numeric(bigfamsallbpfilt$n)), palette='Mako')) + geom_text(aes(label=n), colour = "black", y=-0.03)  # + geom_vline(aes(xintercept=mtm1), color='dimgray', alpha=0.5)
plot_grid(tea, famabund, rel_widths=c(1,0.7), nrow=1, align='hv', axis='tb')

  dev.off()






  

