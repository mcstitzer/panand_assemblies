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

genomecount$rank=bigfam$rank[match(genomecount$genome, bigfam$genome)]


pdf('~/transfer/te_age_try2.pdf', 10,14)
tea=ggplot(genomecount[genomecount$Method=='structural' & genomecount$type=='LTR_retrotransposon' & paste(genomecount$genome, genomecount$collapsedfam)%in% paste(bigfam$genome, bigfam$collapsedfam)[bigfam$rank %in% 1:10],], 
           aes(x=as.numeric(ltr_identity), alpha=rank, group=rank, color=ploidy))  + geom_vline(xintercept=c(0.75,0.85,0.95), color='snow2', linetype='dotted') + geom_vline(xintercept=c(0.7,0.8,0.9), color='snow3', linetype='dotted') + geom_density() + scale_color_manual(values=ploidycolors) +  facet_wrap(~speciesLabel, ncol=1, strip.position='left', scales='free_y', drop=F) + theme( strip.background = element_blank(), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + ylab('')+ xlab('LTR Divergence') + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +xlim(0.85,1) # + geom_vline(aes(xintercept=mtm1), color='dimgray', alpha=0.5)
tea







