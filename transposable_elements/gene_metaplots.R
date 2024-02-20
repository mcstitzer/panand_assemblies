library(ggtree) ## /programs/R-4.2.1-r9/bin/R
library(ggplot2) 
library(cowplot)
theme_set(theme_cowplot())
library(rtracklayer)
library(reshape2)
library(RColorBrewer)
library(tidypaleo) ## facet species names in italics!!!!!
library(tidyverse)

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
a=a[a$type %in% gfftypestokeep,]
a$genome=genotype
# ## get each fam
# temp=data.frame(a) #%>% group_by(fam=gsub('_LTR', '', gsub('_INT', '', Name)), Classification) %>% dplyr::filter(!is.na(as.numeric(Identity))) 
# if('Note' %in% colnames(temp)){temp$Note=''}
# if(genotype=='zmB735'){temp$Parent=''
#                       temp$ltr_identity=NA
#                       temp=temp[,colnames(temp) %in% colnames(genomecountlist[[1]])]}
#   if(genotype=='znicar'){ ## these guys edta is weird - chr not named with chr!!
#                       temp$seqnames[temp$seqnames %in% 1:10]=paste0('chr', temp$seqnames[temp$seqnames %in% 1:10])}
# temp$genome=genotype

genomecountlist[[genotype]]=a
}

# genomecount=do.call(c, genomecountlist)
# genomecount$Identity=as.numeric(genomecount$Identity)

## add a superfamily based on matchign up to the classification field - note most of these NAs are relics from the B73 annotation being included!!!!!
# genomecount$sup=c(NA, 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', 'DHH', NA,NA,NA,NA,NA,NA,'RLC', 'RLG', 'RLG', 'RLX', 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', NA,NA,NA)[match(genomecount$Classification, c("Cent/CentC", "DNA/DTA", "DNA/DTC", "DNA/DTH", "DNA/DTM", "DNA/DTT", 
# "DNA/Helitron", "knob/knob180", "knob/TR-1", "LINE/L1", "LINE/RTE", 
# "LINE/unknown", "Low_complexity", "LTR/Copia", "LTR/CRM", "LTR/Gypsy", 
# "LTR/unknown", "MITE/DTA", "MITE/DTC", "MITE/DTH", "MITE/DTM", 
# "MITE/DTT", "rDNA/spacer", "Simple_repeat", "subtelomere/4-12-1"))]

## tandem repeat maybe want to do by length class???
#genomecount$sup[genomecount$source=='TRASH']=genomecount$Name[genomecount$source=='TRASH']
# genomecount$sup[genomecount$source=='TRASH']='TandemRepeat'



syn=fread('~/transfer/anchorPositions_AnchorWave_Paspalum.txt.gz')
syn=syn[!syn$genome %in% c('pprate', 'svirid', 'tdacn2', 'tdacs2', 'tdactm', 'tzopol', 'bdista', 'agerjg', 'eophiu', 'osativ'),]

syngr=GRanges(seqnames=syn$queryChr, IRanges(start=syn$queryStart, end=syn$queryEnd), strand=syn$strand)
mcols(syngr)$genome=syn$genome
mcols(syngr)$quickgene=syn$quickgene

## get upstream
geneflanks=promoters(syngr, upstream=2000, downstream=0)


geneflanksranges=slide_ranges(geneflanks, width=100, step=10)
geneflanksranges$genome=geneflanks$genome[geneflanksranges$partition]
geneflanksranges$ogstrand=strand(geneflanks)[geneflanksranges$partition]
geneflanksranges$window=rep(1:191,length(geneflanks))


pdf('~/transfer/try_te_metaplot.pdf',12,8)

for(genome in all$V2){
  tesg=genomecountlist[[genome]]
 tes=unstrand(tesg) ### why is this so stupid
 gf=geneflanksranges[geneflanksranges$genome==genome,]
 tewindow=join_overlap_intersect(unstrand(gf), tes)      ## cut tes at boundaries of ranges

posplot=data.frame(tewindow[tewindow$ogstrand=='+',])[,c('window', 'width', 'partition')]
  posplot=posplot %>% complete(partition, window, fill=list(width=0))

  ggplot(posplot, aes(x=window, y=width, group=window)) + geom_boxplot(outlier.shape=NA) + ggtitle(genome)
ggplot(posplot, aes(x=window, y=width, group=partition)) + geom_line(alpha=0.01) + ggtitle(genome)

}


dev.off()
