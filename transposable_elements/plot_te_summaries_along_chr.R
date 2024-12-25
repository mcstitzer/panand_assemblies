library(rtracklayer)
library(dplyr)
library(data.table)
library(stringr)
library(plyranges)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(ggtext)

## xm01 /workdir/mcs368/panand_assemblies/repeats/
all=read.table('../../panand_sp_ploidy.txt')
all=all[!all$V2 %in% c('pprate', 'tdactm', 'tzopol', 'osativ', 'bdista', 'agerjg', 'svirid', 'eophiu'),]

all=all[all$V2!='sbicol',] ## dummy me lost it

genomecountlist=vector(mode = "list", length = length(all$V2))
names(genomecountlist)=all$V2

repeatlengthslist=vector(mode = "list", length = length(all$V2))
names(repeatlengthslist)=all$V2


for( genotype in all$V2){
##import gff3
if(genotype=='snutan'){
a=import.gff3(paste0('trash/', all$V1[all$V2==genotype], '_EDTAandTandemRepeat.negRemoved.gff3'))
}else{
a=import.gff3(paste0('trash/', all$V1[all$V2==genotype], '_EDTAandTandemRepeat.gff3'))
}
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

#write.table(data.frame(genome=names(repeatlengthslist), repeatbp=unlist(repeatlengthslist)),'total_repeat_bp.txt', row.names=F, col.names=T, sep='\t')

## add a superfamily based on matchign up to the classification field - note most of these NAs are relics from the B73 annotation being included!!!!!
genomecount$sup=c(NA, 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', 'DHH', NA,NA,NA,NA,NA,NA,'RLC', 'RLG', 'RLG', 'RLX', 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', NA,NA,NA)[match(genomecount$Classification, c("Cent/CentC", "DNA/DTA", "DNA/DTC", "DNA/DTH", "DNA/DTM", "DNA/DTT", 
"DNA/Helitron", "knob/knob180", "knob/TR-1", "LINE/L1", "LINE/RTE", 
"LINE/unknown", "Low_complexity", "LTR/Copia", "LTR/CRM", "LTR/Gypsy", 
"LTR/unknown", "MITE/DTA", "MITE/DTC", "MITE/DTH", "MITE/DTM", 
"MITE/DTT", "rDNA/spacer", "Simple_repeat", "subtelomere/4-12-1"))]

## tandem repeat maybe want to do by length class???
#genomecount$sup[genomecount$source=='TRASH']=genomecount$Name[genomecount$source=='TRASH']
genomecount$sup[genomecount$source=='TRASH']='TandemRepeat'

genomecount$sup[genomecount$source=='RepeatMasker']='TandemRepeat'
genomecount$Classification[genomecount$source=='RepeatMasker']='TandemRepeat'

## there's still knob in this one znicar family ...
genomecount=genomecount %>% group_by(genome, fam=gsub('_LTR', '', gsub('_INT', '', Name))) %>% filter(!(fam == "TE_00015576" & genome == "znicar")) %>% mutate(fam=fam)

################## skip to the end for now
#  agc=unlist(reduce(split(as_granges(genomecount, keep_mcols=T), ~c(Name,genome)))) ## merge bookended or overlapping TRs if they're the same repeat consensus


## simple count per genome
gcn=genomecount %>% dplyr::group_by(genome, Classification) %>% dplyr::summarize(nfam=length(unique(Name)),nfamsize=median(table(Name)),maxfamsize=max(table(Name)), nfam10=sum(table(Name)>10),nfam100=sum(table(Name)>100),nfam100size=median(table(Name)[table(Name)>100]), nfam1000=sum(table(Name)>1000), ncopy=dplyr::n(), nbp=sum(width))
gcng=dcast(gcn, genome~Classification, value.var='nfam')
## remove na that are nam specific cols
gcng=gcng[,!is.na(gcng[1,])]
nfam=data.frame(genome=gcng$genome, DTA=gcng$`DNA/DTA`+gcng$`MITE/DTA`, DTC=gcng$`DNA/DTC`+gcng$`MITE/DTC`, DTH=gcng$`DNA/DTH`+gcng$`MITE/DTH`,
                  DTM=gcng$`DNA/DTM`+gcng$`MITE/DTM`, DTT=gcng$`DNA/DTT`+gcng$`MITE/DTT`, DHH=gcng$`DNA/Helitron`, 
                  RLC=gcng$`LTR/Copia`, RLG=gcng$`LTR/Gypsy`, RLX=gcng$`LTR/unknown`)
## double haploids
nfam[nfam$genome %in% gs$V2[gs$haploid],-1]=nfam[nfam$genome %in% gs$V2[gs$haploid],-1]*2

gcng=dcast(gcn, genome~Classification, value.var='nfam100')
## remove na that are nam specific cols
gcng=gcng[,!is.na(gcng[1,])]
nfam100=data.frame(genome=gcng$genome, DTA=gcng$`DNA/DTA`+gcng$`MITE/DTA`, DTC=gcng$`DNA/DTC`+gcng$`MITE/DTC`, DTH=gcng$`DNA/DTH`+gcng$`MITE/DTH`,
                  DTM=gcng$`DNA/DTM`+gcng$`MITE/DTM`, DTT=gcng$`DNA/DTT`+gcng$`MITE/DTT`, DHH=gcng$`DNA/Helitron`, 
                  RLC=gcng$`LTR/Copia`, RLG=gcng$`LTR/Gypsy`, RLX=gcng$`LTR/unknown`)
## double haploids
nfam100[nfam100$genome %in% gs$V2[gs$haploid],-1]=nfam100[nfam100$genome %in% gs$V2[gs$haploid],-1]*2


## ncopy
gcng=dcast(gcn, genome~Classification, value.var='ncopy')
## remove na that are nam specific cols
gcng=gcng[,!is.na(gcng[1,])]
ncopy=data.frame(genome=gcng$genome, DTA=gcng$`DNA/DTA`+gcng$`MITE/DTA`, DTC=gcng$`DNA/DTC`+gcng$`MITE/DTC`, DTH=gcng$`DNA/DTH`+gcng$`MITE/DTH`,
                  DTM=gcng$`DNA/DTM`+gcng$`MITE/DTM`, DTT=gcng$`DNA/DTT`+gcng$`MITE/DTT`, DHH=gcng$`DNA/Helitron`, 
                  RLC=gcng$`LTR/Copia`, RLG=gcng$`LTR/Gypsy`, RLX=gcng$`LTR/unknown`)
## double haploids
ncopy[ncopy$genome %in% gs$V2[gs$haploid],-1]=ncopy[ncopy$genome %in% gs$V2[gs$haploid],-1]*2


##nbp
gcng=dcast(gcn, genome~Classification, value.var='nbp')
## remove na that are nam specific cols
gcng=gcng[,!is.na(gcng[1,])]
nbp=data.frame(genome=gcng$genome, DTA=gcng$`DNA/DTA`+gcng$`MITE/DTA`, DTC=gcng$`DNA/DTC`+gcng$`MITE/DTC`, DTH=gcng$`DNA/DTH`+gcng$`MITE/DTH`,
                  DTM=gcng$`DNA/DTM`+gcng$`MITE/DTM`, DTT=gcng$`DNA/DTT`+gcng$`MITE/DTT`, DHH=gcng$`DNA/Helitron`, 
                  RLC=gcng$`LTR/Copia`, RLG=gcng$`LTR/Gypsy`, RLX=gcng$`LTR/unknown`)
## double haploids
nbp[nbp$genome %in% gs$V2[gs$haploid],-1]=nbp[nbp$genome %in% gs$V2[gs$haploid],-1]*2

## get statistics for the paper!
gs=read.table('~/transfer/panand_assembly_sizes.txt', header=T, sep='\t')
cor.test(gs$haploidAssemblySize-gs$haploidNCount, gs$haploidRepeatSize)
cor.test(amm$haploidAssemblySize[amm$variable=='tebp']-amm$haploidNCount[amm$variable=='tebp'], amm$doubledValue[amm$variable=='tebp']/2)
gs$polyploid=gs$ploidy!='Diploid'
## this is TE only
t.test(amm$doubledValue[amm$variable=='tebp' & amm$polyploid]/2, amm$doubledValue[amm$variable=='tebp' & !amm$polyploid]/2)
## this contains TR and TE
t.test((gs$haploidAssemblySize-gs$haploidNCount)[gs$polyploid], (gs$haploidAssemblySize-gs$haploidNCount)[!gs$polyploid])

t.test(gs$haploidRepeatSize[gs$polyploid]/(gs$haploidAssemblySize-gs$haploidNCount)[gs$polyploid], gs$haploidRepeatSize[!gs$polyploid]/(gs$haploidAssemblySize--gs$haploidNCount)[!gs$polyploid])
t.test(gs$haploidRepeatSize[gs$polyploid & gs$ploidy!='Paleotetraploid']/(gs$haploidAssemblySize-gs$haploidNCount)[gs$polyploid & gs$ploidy!='Paleotetraploid'], gs$haploidRepeatSize[!gs$polyploid & gs$ploidy!='Paleotetraploid']/(gs$haploidAssemblySize-gs$haploidNCount)[!gs$polyploid & gs$ploidy!='Paleotetraploid'])

t.test(rowSums(nfam[,2:10])[nfam$genome%in%gs$V2[gs$polyploid]], rowSums(nfam[,2:10])[nfam$genome%in%gs$V2[!gs$polyploid]])
t.test(rowSums(ncopy[,2:10])[nfam$genome%in%gs$V2[gs$polyploid]], rowSums(ncopy[,2:10])[ncopy$genome%in%gs$V2[!gs$polyploid]])
t.test(rowSums(nbp[,2:10])[nfam$genome%in%gs$V2[gs$polyploid]], rowSums(nbp[,2:10])[nbp$genome%in%gs$V2[!gs$polyploid]])
t.test(rowSums(nfam[,2:10])[nfam$genome%in%gs$V2[gs$polyploid& gs$ploidy!='Paleotetraploid']], rowSums(nfam[,2:10])[nfam$genome%in%gs$V2[!gs$polyploid& gs$ploidy!='Paleotetraploid']])
t.test(rowSums(ncopy[,2:10])[nfam$genome%in%gs$V2[gs$polyploid& gs$ploidy!='Paleotetraploid']], rowSums(ncopy[,2:10])[ncopy$genome%in%gs$V2[!gs$polyploid& gs$ploidy!='Paleotetraploid']])
t.test(rowSums(nbp[,2:10])[nfam$genome%in%gs$V2[gs$polyploid& gs$ploidy!='Paleotetraploid']], rowSums(nbp[,2:10])[nbp$genome%in%gs$V2[!gs$polyploid& gs$ploidy!='Paleotetraploid']])

t.test(rowSums(nfam100[,2:10])[nfam100$genome%in%gs$V2[gs$polyploid]], rowSums(nfam100[,2:10])[nfam100$genome%in%gs$V2[!gs$polyploid]])
t.test(rowSums(nfam100[,2:10])[nfam100$genome%in%gs$V2[gs$polyploid& gs$ploidy!='Paleotetraploid']], rowSums(nfam100[,2:10])[nfam100$genome%in%gs$V2[!gs$polyploid& gs$ploidy!='Paleotetraploid']])

t.test(rowSums(nfam100[,2:10])[nfam100$genome%in%gs$V2[ gs$ploidy=='Tetraploid']], rowSums(nfam100[,2:10])[nfam100$genome%in%gs$V2[gs$ploidy=='Paleotetraploid']])

## make summaries of age
ages=genomecount %>% filter(Method=='structural') %>% dplyr::group_by(genome, Method) %>% dplyr::filter(grepl('LTR', type)) %>% dplyr::summarize(copies=dplyr::n(), meanage=mean(as.numeric(ltr_identity), na.rm=T), meanage90=mean(as.numeric(ltr_identity)[as.numeric(ltr_identity)>0.9], na.rm=T), nIdentical=sum(as.numeric(ltr_identity)==1, na.rm=T), nYoung=sum(as.numeric(ltr_identity)>0.99, na.rm=T)) %>% data.frame()                                  
gs$meanage=ages$meanage[match(gs$V2, ages$genome)]
gs$meanage90=ages$meanage90[match(gs$V2, ages$genome)]
gs$nIdentical=ages$nIdentical[match(gs$V2, ages$genome)]
gs$nYoung=ages$nYoung[match(gs$V2, ages$genome)]
gs$propIdentical=gs$nIdentical/ages$copies[match(gs$V2, ages$genome)]
gs$propYoung=gs$nYoung/ages$copies[match(gs$V2, ages$genome)]

solos=genomecount %>% dplyr::group_by(genome, Method) %>% dplyr::filter(grepl('LTR', type)) %>% dplyr::summarize(copies=dplyr::n(), nInt=sum(grepl('_INT', Name) & width>1000 & width<50000, na.rm=T), nLTR=sum(grepl('_LTR', Name) & width>100 & width<3000, na.rm=T), nStructural=sum(!is.na(ltr_identity), na.rm=T)) %>% data.frame()                                  
solos$nStructural[solos$Method=='homology']=solos$nStructural[solos$Method=='structural']
solos$nIntact=solos$nStructural+solos$nInt
solos$ratioLTRINT=(solos$nLTR-2*solos$nInt)/solos$nIntact
solos$SIratio=(solos$nLTR-(solos$nInt*2))/solos$nInt
solos$ltrStructRatio=solos$nLTR/solos$nStructural
gs$soloProp=solos$ratioLTRINT[match(gs$V2, solos$genome)]
gs$SIratio=solos$SIratio[match(gs$V2, solos$genome)]
gs$ltrStructRatio=solos$ltrStructRatio[match(gs$V2, solos$genome)]

nfam$ploidy=gs$ploidy[match(nfam$genome, gs$V2)]
nfam100$ploidy=gs$ploidy[match(nfam100$genome, gs$V2)]
nfam$ploidy[nfam$genome=='zmB735']='Paleotetraploid'
nfam100$ploidy[nfam100$genome=='zmB735']='Paleotetraploid'
nfam=nfam[!nfam$genome%in%c('tdacn2', 'tdacs2', 'zmB735'),] ## got to leave out b73 because it's annotated sooooo differently
nfam100=nfam100[!nfam100$genome%in%c('tdacn2', 'tdacs2', 'zmB735'),]
nfam$totFam=rowSums(nfam[,2:10])
nfam100$totFam=rowSums(nfam100[,2:10])

ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')

gs$ploidy=factor(gs$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))
nfam$ploidy=factor(nfam$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))
nfam100$ploidy=factor(nfam100$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))

amm$ploidy=factor(amm$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))


pdf(paste0('~/transfer/te_panand_fig.', Sys.Date(), '.pdf'), 15,4)
                                      
bp1=ggplot(amm[amm$variable=='tebp',], aes(x=ploidy, y=doubledValue/2, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=4e9) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('TE base pairs') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
pr2=ggplot(amm[amm$variable=='tebp',], aes(x=ploidy, y=prop*100, color=ploidy)) + geom_boxplot(outlier.shape=NA) +geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=100) + 
                                      scale_color_manual(values=ploidycolors) + xlab('Ploidy') + ylab('TE proportion')+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
fa3=ggplot(nfam, aes(x=ploidy, y=totFam, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=210000) + 
                                      scale_color_manual(values=ploidycolors) + xlab('Ploidy') + ylab('TE family count')+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
fa4=ggplot(nfam100, aes(x=ploidy, y=totFam, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=23000) + 
                                      scale_color_manual(values=ploidycolors) + xlab('Ploidy') + ylab('TE family count (>100 copies)')+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ag5=ggplot(gs, aes(x=ploidy, y=meanage, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid",label.y=1) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('LTR sequence identity') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ag590=ggplot(gs, aes(x=ploidy, y=meanage90, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid",label.y=0.8) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('LTR sequence identity') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
id6=ggplot(gs, aes(x=ploidy, y=propIdentical, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) +ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid",label.y=0.2) +  
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('Proportion LTRs Identical') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
tebpident=ggplot(gs, aes(x=propIdentical,y=haploidRepeatSize,  color=ploidy)) + geom_point()  + #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + #ggpubr::stat_compare_means(comparisons=lapply(1:3, function(x) c(combn(unique(gr$group)[1:3], 2)[1,x], combn(unique(gr$group)[1:3], 2)[2,x])),label = 'p.signif', show.legend = F) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + ylab('TE base pairs') + xlab('Proportion LTRs Identical') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
tepropident=ggplot(gs, aes( x=propIdentical,y=haploidRepeatSize/haploidAssemblySize*100, color=ploidy)) + geom_point()  + #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + #ggpubr::stat_compare_means(comparisons=lapply(1:3, function(x) c(combn(unique(gr$group)[1:3], 2)[1,x], combn(unique(gr$group)[1:3], 2)[2,x])),label = 'p.signif', show.legend = F) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + ylab('TE proportion') + xlab('Proportion LTRs Identical') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
tesolo=ggplot(gs, aes( x=ploidy,y=soloProp, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge())  + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid",label.y=9) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + ylab('LTR:Intact Ratio') + xlab('Ploidy') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

legend <- get_legend(bp1)                                      
plot_grid(bp1 + theme(legend.position='NULL'), pr2+ theme(legend.position='NULL'), fa3+ theme(legend.position='NULL'), fa4+ theme(legend.position='NULL'),
          legend, nrow=1,labels=c('a', 'b', 'c', 'd', ''))
plot_grid(tebpident + theme(legend.position='NULL'), tepropident+ theme(legend.position='NULL'), ag5+ theme(legend.position='NULL'), id6+ theme(legend.position='NULL'),
          legend, nrow=1,labels=c('a', 'b', 'c', 'd', ''))
plot_grid(bp1 + theme(legend.position='NULL'), pr2+ theme(legend.position='NULL'), fa3+ theme(legend.position='NULL'), ag5+ theme(legend.position='NULL'), tepropident+ theme(legend.position='NULL'),
          legend, align='hv', nrow=1,labels=c('a', 'b', 'c', 'd', 'e', ''))
                                                                                                                                                                                                                                                                                                                              
dev.off()

pdf(paste0('~/transfer/te_panand_fig5.', Sys.Date(), '.pdf'), 15,4)
plot_grid(bp1 + theme(legend.position='NULL'), pr2+ theme(legend.position='NULL'), fa3+ theme(legend.position='NULL'), fa4+ theme(legend.position='NULL'), ag5+ theme(legend.position='NULL'),
          legend, align='hv', nrow=1,labels=c('a', 'b', 'c', 'd', 'e', ''))
                                                                                                                                                                                                                                                                                                                              
dev.off()


#### read in genes!
genecountlist=vector(mode = "list", length = length(all$V2))
names(genecountlist)=all$V2

for( i in all$V2){
##import gff3
  if(i!='zluxur'){
    if(!i %in% c('sbicol', 'zmB735')){
a=import.gff3(Sys.glob(paste0('../genes/', all$V1[all$V2==i], '*.2.gff3')))
      }
    if(i=='sbicol'){
      a=import.gff3('../genes/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.55.gff3')
      mcols(a)=mcols(a)[,colnames(mcols(a)) %in% colnames(genecountlist[[1]])]
      mcols(a)$canonical_transcript=NA
    }
    if(i=='zmB735'){
      a=import.gff3('../genes/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3')
      mcols(a)=mcols(a)[,colnames(mcols(a)) %in% colnames(genecountlist[[1]])]
      }
## get each fam
a$genome=i
a$Name=NULL ## remove Name from liftoff genes
genecountlist[[i]]=data.frame(a)
}
}

genes=do.call(rbind, genecountlist)
genes=genes[genes$type=='gene',]
    
## simple counts of repeats
pdf(paste0('~/transfer/chromosomes_panand.', Sys.Date(), '.pdf'), 20, 10)
for(genome in unique(genomecount$genome)){
for(i in unique(genomecount$seqnames[genomecount$end>50e6 & genomecount$genome==genome])){
print(ggplot(genomecount[genomecount$seqnames==i & genomecount$genome==genome,], aes(x=start, fill=factor(sup))) + 
geom_histogram(binwidth=1e6, position='stack') + ggtitle(paste0(genome, i)))}
}
dev.off()

### these edta annotations look really bad at centromeres - they're either calling dtm or rlc in these regions
### fixed that with trash!!!!
### instead, weight by amount of sequence to see if it's actually that bad??

pdf(paste0('~/transfer/chromosomes_panand.mbscaled.genes.', Sys.Date(), '.pdf'), 20, 10)
for(genome in all$V2){
for(i in unique(genomecount$seqnames[genomecount$end>50e6 & genomecount$genome==genome])){
print(plot_grid(ggplot(genomecount[genomecount$seqnames==i & genomecount$genome==genome,], aes(x=start, fill=factor(sup), weight=width)) + 
geom_histogram(binwidth=1e6, position='stack') + ggtitle(paste0(genome, i))+ theme(legend.position = "none"), 
               ggplot(genes[genes$seqnames==i & genes$genome==genome,], aes(x=start)) + geom_histogram(binwidth=1e6),
               align='hv', ncol=1, rel_heights=c(1,0.2)))}
}
dev.off()



pdf(paste0('~/transfer/chromosomes_panand.TRmbscaled.genes.', Sys.Date(), '.pdf'), 20, 10)
for(genome in all$V2){
for(i in unique(genomecount$seqnames[genomecount$end>50e6 & genomecount$genome==genome])){
print(plot_grid(ggplot(genomecount[genomecount$seqnames==i & genomecount$genome==genome & genomecount$sup=='TandemRepeat',], aes(x=start, fill=factor(Name), weight=width)) + 
geom_histogram(binwidth=1e6, position='stack') + ggtitle(paste0(genome, i)), 
               ggplot(genes[genes$seqnames==i & genes$genome==genome,], aes(x=start)) + geom_histogram(binwidth=1e6),
               align='hv', ncol=1, rel_heights=c(1,0.2)))}
}
dev.off()

pdf(paste0('~/transfer/chromosomes_panand.DTTmbscaled.genes.', Sys.Date(), '.pdf'), 20, 10)
for(genome in all$V2){
for(i in unique(genomecount$seqnames[genomecount$end>50e6 & genomecount$genome==genome])){
print(plot_grid(ggplot(genomecount[genomecount$seqnames==i & genomecount$genome==genome & genomecount$sup=='DTT',], aes(x=start, weight=width)) + 
geom_histogram(binwidth=1e6, position='stack') + ggtitle(paste0(genome, i)), 
               ggplot(genes[genes$seqnames==i & genes$genome==genome,], aes(x=start)) + geom_histogram(binwidth=1e6),
               align='hv', ncol=1, rel_heights=c(1,0.2)))}
}
dev.off()

ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')
asize=fread('../panand_assembly_sizes.txt', header=T, quote="", fill=T)
## factor so correct order
asize$ploidy=factor(asize$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))

ages=genomecount %>% filter(Method=='structural') %>% dplyr::group_by(genome, Method) %>% dplyr::filter(grepl('LTR', type)) %>% dplyr::summarize(copies=dplyr::n(), medianage=median(as.numeric(ltr_identity), na.rm=T), meanage=mean(as.numeric(ltr_identity), na.rm=T), meanage90=mean(as.numeric(ltr_identity)[as.numeric(ltr_identity)>0.9], na.rm=T), nIdentical=sum(as.numeric(ltr_identity)==1, na.rm=T), nYoung=sum(as.numeric(ltr_identity)>0.99, na.rm=T)) %>% data.frame()  
ages$propIdentical=ages$nIdentical/ages$copies
ages$propYoung=ages$nYoung/ages$copies ## young here is 1% divergence :)
ages$mya=asize$mya[match(ages$genome, asize$V2)]
ages$ploidy=asize$ploidy[match(ages$genome, asize$V2)]

evenness=genomecount %>% group_by(genome, sup, fam=gsub('_LTR', '', gsub('_INT', '', Name))) %>% filter(!(fam == "TE_00015576" & genome == "znicar"), !is.na(fam)) %>% summarize(nCopies=n(), bp=sum(width), avgwidth=mean(width)) 
evennesstoplot=evenness %>% group_by(genome) %>% mutate(total_bp=sum(bp), p_i=bp/total_bp, totalcopies=sum(nCopies), p_iN=nCopies/totalcopies)%>% summarize(nfam=dplyr::n(), H=-sum(p_i*log(p_i)), S=dplyr::n(), HN=-sum(p_iN*log(p_iN))) %>% mutate(J=H/log(S), JN=HN/log(S)) %>% ungroup()
evennesstoplot4copies=evenness %>% filter(nCopies>3)%>% group_by(genome) %>% mutate(total_bp=sum(bp), p_i=bp/total_bp, totalcopies=sum(nCopies), p_iN=nCopies/totalcopies)%>% summarize(nfam=dplyr::n(), H=-sum(p_i*log(p_i)), S=dplyr::n(), HN=-sum(p_iN*log(p_iN))) %>% mutate(J=H/log(S), JN=HN/log(S)) %>% ungroup()
evennesstoplot10copies=evenness %>% filter(nCopies>9)%>% group_by(genome) %>% mutate(total_bp=sum(bp), p_i=bp/total_bp, totalcopies=sum(nCopies), p_iN=nCopies/totalcopies)%>% summarize(nfam=dplyr::n(), H=-sum(p_i*log(p_i)), S=dplyr::n(), HN=-sum(p_iN*log(p_iN))) %>% mutate(J=H/log(S), JN=HN/log(S)) %>% ungroup()



evennesstoplot$ploidy=asize$ploidy[match(evennesstoplot$genome, asize$V2)]
evennesstoplot$haploidRepeatSize=asize$haploidRepeatSize[match(evennesstoplot$genome, asize$V2)]
evennesstoplot10copies$ploidy=asize$ploidy[match(evennesstoplot10copies$genome, asize$V2)]
evennesstoplot10copies$haploidRepeatSize=asize$haploidRepeatSize[match(evennesstoplot10copies$genome, asize$V2)]

evennesstoplot$mya=asize$mya[match(evennesstoplot$genome, asize$V2)]
evennesstoplot10copies$mya=asize$mya[match(evennesstoplot10copies$genome, asize$V2)]


pdf('~/transfer/evenness_tefam.pdf',8,8)
ggplot(evennesstoplot[evennesstoplot$genome!='zmB735',], aes(x=J, y=nfam, color=ploidy)) + geom_point() + scale_color_manual(values=ploidycolors)
ggplot(evennesstoplot[evennesstoplot$genome!='zmB735',], aes(x=JN, y=nfam, color=ploidy)) + geom_point() + scale_color_manual(values=ploidycolors)
ggplot(evennesstoplot[evennesstoplot$genome!='zmB735',], aes(x=J, y=nfam, color=ploidy)) + geom_point() + scale_color_manual(values=ploidycolors)+geom_text(aes(label=genome))
ggplot(evennesstoplot[evennesstoplot$genome!='zmB735',], aes(x=JN, y=nfam, color=ploidy)) + geom_point() + scale_color_manual(values=ploidycolors)+geom_text(aes(label=genome))
ggplot(evennesstoplot[evennesstoplot$genome!='zmB735',], aes(x=J, y=nfam, color=ploidy)) + geom_point(aes(size=haploidRepeatSize)) + scale_color_manual(values=ploidycolors)
ggplot(evennesstoplot[evennesstoplot$genome!='zmB735',], aes(x=JN, y=nfam, color=ploidy)) + geom_point(aes(size=haploidRepeatSize)) + scale_color_manual(values=ploidycolors)

ggplot(evennesstoplot10copies[evennesstoplot10copies$genome!='zmB735',], aes(x=J, y=nfam, color=ploidy)) + geom_point(aes(size=haploidRepeatSize)) + scale_color_manual(values=ploidycolors)
ggplot(evennesstoplot10copies[evennesstoplot10copies$genome!='zmB735',], aes(x=JN, y=nfam, color=ploidy)) + geom_point(aes(size=haploidRepeatSize)) + scale_color_manual(values=ploidycolors)

ggplot(evennesstoplot[evennesstoplot$genome!='zmB735',], aes(y=J, x=mya, color=ploidy)) + geom_point(aes(size=haploidRepeatSize)) + scale_color_manual(values=ploidycolors)
ggplot(evennesstoplot[evennesstoplot$genome!='zmB735',], aes(y=JN, x=mya, color=ploidy)) + geom_point(aes(size=haploidRepeatSize)) + scale_color_manual(values=ploidycolors)

ggplot(evennesstoplot10copies[evennesstoplot10copies$genome!='zmB735',], aes(y=J, x=mya, color=ploidy)) + geom_point(aes(size=haploidRepeatSize)) + scale_color_manual(values=ploidycolors)
ggplot(evennesstoplot10copies[evennesstoplot10copies$genome!='zmB735',], aes(y=JN, x=mya, color=ploidy)) + geom_point(aes(size=haploidRepeatSize)) + scale_color_manual(values=ploidycolors)

ggplot(ages[ages$genome!='zmB735',], aes(x=mya, y=meanage, color=ploidy)) + geom_point() + scale_color_manual(values=ploidycolors)
ggplot(ages[ages$genome!='zmB735',], aes(x=mya, y=medianage, color=ploidy)) + geom_point() + scale_color_manual(values=ploidycolors)
ggplot(ages[ages$genome!='zmB735',], aes(x=mya, y=propIdentical, color=ploidy)) + geom_point() + scale_color_manual(values=ploidycolors)
ggplot(ages[ages$genome!='zmB735',], aes(x=mya, y=propYoung, color=ploidy)) + geom_point() + scale_color_manual(values=ploidycolors)


ggplot(ages[ages$genome!='zmB735',], aes(x=mya, y=(1-medianage)/6.5e-9/2/1e6, color=ploidy, size=copies)) + geom_point() + scale_color_manual(values=ploidycolors) + ylab('Median Insertion Time, LTR retrotransposon (Mya)')+xlab('Divergence Between Parental Genomes, (Mya)')

dev.off()



### okay this will work!!! make a stacked bar plot of each of these genome-wide, using these colors
# Retrotransposon colors (pastel reds and pinks)
retro_colors <- c(RLG="#F7B7B7", RLC="#F28E8E", RLX="#E67373")
# DNA transposon colors (cohesive blue gradient)
dna_colors <- c(DHH="#377EB8", DTA="#4B8EBA", DTC="#5DA9BD", DTH="#70C2BF", 
                DTM="#83D6C1", DTT="#96E2C3", DTS="#A8EEC5")
# Tandem repeat color (pastel lavender)
tandem_color <- c(TandemRepeat="#D4B4F4", TR='#D4B4F4')
# Combine all colors
te_colors <- c(retro_colors, dna_colors, tandem_color)

#### then, make a stacked bar histogram of the largest scaffold/chromosome in each assembly, to show the relative proportions


#### WATCH oUYT
### changing this here so parv is shortened!!!!! for TILs
shorttaxonnames=c("Z. mays ssp. parv. TIL11",  "Z. mays ssp. parv. TIL01", "Z. mays ssp. mays B73v5", "Z. mays ssp. mex. TIL25", "Z. mays ssp. mex. TIL18", "Z. mays ssp. huehue.", 
                  "Z. luxurians", "Z. nicaraguensis", "Z. diploperennis Momo", "Z. diploperennis Gigi", "T. zoloptense", "T. dactyloides FL", "T. dactyloides Southern Hap2", 
                  "T. dactyloides Northern Hap2", "T. dactyloides KS", "T. dactyloides tetraploid", "U. digitatum", "V. cuspidata", "R. rottboellioides", "R. tuberculosa", 
                  "H. compressa", "E. tripsacoides", "S. scoparium", "S. microstachyum", "A. virginicum", "A. chinensis", "A. gerardi", 
                  "C. refractus", "C. citratus", "H. contortus", "T. triandra", "B. laguroides", "P. paniceum", "S. bicolor", 
                  "I. rugosum", "S. nutans", '"A." burmanicus', "T. elegans", "C. serrulatus", "P. vaginatum")
names(shorttaxonnames)=c("zTIL11",  "zTIL01", "zmB735", "zTIL25", "zTIL18", "zmhuet", 
                         "zluxur", "znicar", "zdmomo", "zdgigi", "tzopol", "tdacs1", "tdacs2", 
                         "tdacn2", "tdacn1", "tdactm", "udigit", "vcuspi", "rrottb", "rtuber", 
                         "hcompr", "etrips", "sscopa", "smicro", "avirgi", "achine", "agerar", 
                         "crefra", "ccitra", "hconto", "ttrian", "blagur", "ppanic", "sbicol", 
                         "irugos", "snutan", "atenui", "telega", "cserru", "pvagin")

# List of substrings to exclude from italics
exceptions <- c('subsp.', 'ssp.', 'TIL11', 'TIL01', 
                'TIL25', 'TIL18', 'Momo', 'Gigi', 
                'Southern Hap1', 'Northern Hap1', 
                'FL', 'KS', '\\*', '\\"', 'B73v5')

# Custom function to apply italics with exceptions
format_labels <- function(labels, exceptions) {
  sapply(labels, function(label) {
    # Apply italics to everything initially
    parts <- unlist(strsplit(label, " "))  # Split label into words
    formatted_parts <- sapply(parts, function(word) {
      if (word %in% exceptions) {
        return(word)  # Keep exceptions in regular text
      } else {
        return(paste0("<i>", word, "</i>"))  # Italicize other words
      }
    })
    # Combine parts back into a single string
    return(paste(formatted_parts, collapse = " "))
  })
} 

pdf(paste0('~/transfer/chromosomes_panand.mbscaled.genes.', Sys.Date(), '.pdf'), 20, 12)
for(genome in all$V2){
#for(i in unique(genomecount$seqnames[genomecount$end>50e6 & genomecount$genome==genome])){
longest=genomecount[genomecount$genome==genome,]
longest=longest[longest$seqnames==longest$seqnames[which.max(longest$end)],]

shortSpeciesLabel=format_labels(shorttaxonnames[genome], exceptions)
if(genome=='aburma'){shortSpeciesLabel='<i>&quot;A.&quot;</i> <i>burmanicus</i>'}


## these are just the gene entry
genes=read.table(paste0('../repeats/coverage_results/', all$V1[all$V2==genome], '_filtered_genes.gff'))
genes=genes[genes$V1==longest$seqnames[which.max(longest$end)],]


print(
plot_grid(
ggplot(longest, aes(x=start, fill=factor(sup), weight=width)) + scale_fill_manual(values=te_colors)+
geom_histogram(binwidth=1e6, position='stack') + ggtitle(paste0(shortSpeciesLabel,' Longest Sequence: ', longest$seqnames[1]))+ #, 
xlab('Position on Sequence (1 Mb bins)')+ylab('Count')+
  theme(
    legend.position = "top",                # Place legend at bottom
    legend.title = element_blank(),            # Remove legend title
    legend.text = element_text(size = 8),     # Adjust text size
    legend.key.height = unit(0.5, "cm"),       # Reduce height of keys
    legend.key.width = unit(1, "cm"),           # Adjust width of keys
    axis.text.y=element_markdown(size=9),
    legend.justification='center',
    plot.title=element_markdown(hjust=0.5)
  ) +
  guides(fill = guide_legend(
    title.position = "top",                    # Moves legend title to top (optional)
    label.position = "bottom",                    # Places labels above boxes
    nrow = 1,                                  # Arrange in one row
    byrow = TRUE,                               # Fill row by row
    reverse=T
  )),
  
               ggplot(genes, aes(x=V2)) + geom_histogram(binwidth=1e6)+xlab('Position on Sequence (1 Mb bins)')+ ylab('Count')+ ggtitle('Helixer Genes') + theme(plot.title=element_text(hjust=0.5), axis.text.y=element_markdown(size=9)),
               align='hv', ncol=1, rel_heights=c(1,0.3)))
               
          #    }
}
dev.off()

### okay exactly the same, but plot all the individual files
for(genome in all$V2){
pdf(paste0('~/transfer/', genome, '_chrdist.pdf'), 15, 8)

longest=genomecount[genomecount$genome==genome,]
longest=longest[longest$seqnames==longest$seqnames[which.max(longest$end)],]

shortSpeciesLabel=format_labels(shorttaxonnames[genome], exceptions)
if(genome=='aburma'){shortSpeciesLabel='<i>&quot;A.&quot;</i> <i>burmanicus</i>'}


## these are just the gene entry
genes=read.table(paste0('../repeats/coverage_results/', all$V1[all$V2==genome], '_filtered_genes.gff'))
genes=genes[genes$V1==longest$seqnames[which.max(longest$end)],]


print(
plot_grid(
ggplot(longest, aes(x=start, fill=factor(sup), weight=width)) + scale_fill_manual(values=te_colors)+
geom_histogram(binwidth=1e6, position='stack') + ggtitle(paste0(shortSpeciesLabel,' Longest Sequence: ', longest$seqnames[1]))+ #, 
xlab('Position on Sequence (1 Mb bins)')+ylab('Count')+
  theme(
    legend.position = "top",                # Place legend at bottom
    legend.title = element_blank(),            # Remove legend title
    legend.text = element_text(size = 8),     # Adjust text size
    legend.key.height = unit(0.5, "cm"),       # Reduce height of keys
    legend.key.width = unit(1, "cm"),           # Adjust width of keys
    axis.text.y=element_markdown(size=9),
    legend.justification='center',
    plot.title=element_markdown(hjust=0.5)
  ) +
  guides(fill = guide_legend(
    title.position = "top",                    # Moves legend title to top (optional)
    label.position = "bottom",                    # Places labels above boxes
    nrow = 1,                                  # Arrange in one row
    byrow = TRUE,                               # Fill row by row
    reverse=T
  )),
  
               ggplot(genes, aes(x=V2)) + geom_histogram(binwidth=1e6)+xlab('Position on Sequence (1 Mb bins)')+ ylab('Count')+ ggtitle('Helixer Genes') + theme(plot.title=element_text(hjust=0.5), axis.text.y=element_markdown(size=9)),
               align='hv', ncol=1, rel_heights=c(1,0.4)))
               
   dev.off()

}



pdf('~/transfer/barplots_tesup.pdf',8,8)

bardata=genomecount%>% group_by(sup, genome) %>% summarize(bp=sum(width),.groups='drop')

bardata$shortSpeciesLabel=shorttaxonnames[match(bardata$genome, names(shorttaxonnames))]
bardata$shortSpeciesLabel=factor(bardata$shortSpeciesLabel, levels=rev(shorttaxonnames))

# Apply the formatting to y-axis labels
formatted_labels <- format_labels(levels(bardata$shortSpeciesLabel), exceptions)
formatted_labels[4]='<i>&quot;A.&quot;</i> <i>burmanicus</i>'
formatted_labels['Z. luxurians']='<i>Z. luxurians</i>'

bardata$supLegend <- factor(bardata$sup, levels = c(
  "DHH", "DTA", "DTC", "DTH", "DTM", "DTT", "RLC", "RLG", "RLX", "TandemRepeat"),
  labels = c("DHH", "DTA", "DTC", "DTH", "DTM", "DTT", "RLC", "RLG", "RLX", "TR")  # Change label to TR
)
ggplot(bardata[!bardata$genome%in%c('zmB735', 'tdacn2', 'tdacs2'),], aes(y=genome, x=bp, fill=supLegend)) + geom_bar(stat='identity') + scale_fill_manual(values=te_colors)+ylab('Genome')+xlab('Total bp')
ggplot(bardata[!bardata$genome%in%c('zmB735', 'tdacn2', 'tdacs2'),], aes(y=shortSpeciesLabel, x=bp, fill=supLegend)) + geom_col(position='fill') + scale_fill_manual(values=te_colors)+ylab('Genome')+xlab('Proportion of Repeat Base Pairs')+
  geom_vline(xintercept=c(0.25,0.5,0.75), color='snow', alpha=0.1)+
  theme(
    legend.position = "top",                # Place legend at bottom
    legend.title = element_blank(),            # Remove legend title
    legend.text = element_text(size = 8),     # Adjust text size
    legend.key.height = unit(0.5, "cm"),       # Reduce height of keys
    legend.key.width = unit(1, "cm"),           # Adjust width of keys
    axis.text.y=element_markdown(size=9),
    legend.justification='center'
  ) +
  scale_y_discrete(labels=formatted_labels)+
  guides(fill = guide_legend(
    title.position = "top",                    # Moves legend title to top (optional)
    label.position = "bottom",                    # Places labels above boxes
    nrow = 1,                                  # Arrange in one row
    byrow = TRUE,                               # Fill row by row
    reverse=T
  ))+ylab('')

dev.off()


write.table(bardata, 'te_barplot_data.txt', row.names=F, col.names=T, sep='\t', quote=F)
write.table(evennesstoplot10copies, 'evenness10copies_plot_data.txt', row.names=F, col.names=T, sep='\t', quote=F)
write.table(ages, 'ages_plot_data.txt', row.names=F, col.names=T, sep='\t', quote=F)

