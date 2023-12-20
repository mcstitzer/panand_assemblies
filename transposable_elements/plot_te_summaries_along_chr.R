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

write.table(data.frame(genome=names(repeatlengthslist), repeatbp=unlist(repeatlengthslist)),'total_repeat_bp.txt', row.names=F, col.names=T, sep='\t')

## add a superfamily based on matchign up to the classification field - note most of these NAs are relics from the B73 annotation being included!!!!!
genomecount$sup=c(NA, 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', 'DHH', NA,NA,NA,NA,NA,NA,'RLC', 'RLG', 'RLG', 'RLX', 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', NA,NA,NA)[match(genomecount$Classification, c("Cent/CentC", "DNA/DTA", "DNA/DTC", "DNA/DTH", "DNA/DTM", "DNA/DTT", 
"DNA/Helitron", "knob/knob180", "knob/TR-1", "LINE/L1", "LINE/RTE", 
"LINE/unknown", "Low_complexity", "LTR/Copia", "LTR/CRM", "LTR/Gypsy", 
"LTR/unknown", "MITE/DTA", "MITE/DTC", "MITE/DTH", "MITE/DTM", 
"MITE/DTT", "rDNA/spacer", "Simple_repeat", "subtelomere/4-12-1"))]

## tandem repeat maybe want to do by length class???
#genomecount$sup[genomecount$source=='TRASH']=genomecount$Name[genomecount$source=='TRASH']
genomecount$sup[genomecount$source=='TRASH']='TandemRepeat'



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
cor.test(gs$haploidAssemblySize, gs$haploidRepeatSize)
gs$polyploid=gs$ploidy!='Diploid'
t.test(gs$haploidAssemblySize[gs$polyploid], gs$haploidAssemblySize[!gs$polyploid])
t.test(gs$haploidRepeatSize[gs$polyploid]/gs$haploidAssemblySize[gs$polyploid], gs$haploidRepeatSize[!gs$polyploid]/gs$haploidAssemblySize[!gs$polyploid])
t.test(gs$haploidRepeatSize[gs$polyploid & gs$ploidy!='Paleotetraploid']/gs$haploidAssemblySize[gs$polyploid & gs$ploidy!='Paleotetraploid'], gs$haploidRepeatSize[!gs$polyploid & gs$ploidy!='Paleotetraploid']/gs$haploidAssemblySize[!gs$polyploid & gs$ploidy!='Paleotetraploid'])

t.test(rowSums(nfam[,-1])[nfam$genome%in%gs$V2[gs$polyploid]], rowSums(nfam[,-1])[nfam$genome%in%gs$V2[!gs$polyploid]])
t.test(rowSums(ncopy[,-1])[nfam$genome%in%gs$V2[gs$polyploid]], rowSums(ncopy[,-1])[ncopy$genome%in%gs$V2[!gs$polyploid]])
t.test(rowSums(nbp[,-1])[nfam$genome%in%gs$V2[gs$polyploid]], rowSums(nbp[,-1])[nbp$genome%in%gs$V2[!gs$polyploid]])
t.test(rowSums(nfam[,-1])[nfam$genome%in%gs$V2[gs$polyploid& gs$ploidy!='Paleotetraploid']], rowSums(nfam[,-1])[nfam$genome%in%gs$V2[!gs$polyploid& gs$ploidy!='Paleotetraploid']])
t.test(rowSums(ncopy[,-1])[nfam$genome%in%gs$V2[gs$polyploid& gs$ploidy!='Paleotetraploid']], rowSums(ncopy[,-1])[ncopy$genome%in%gs$V2[!gs$polyploid& gs$ploidy!='Paleotetraploid']])
t.test(rowSums(nbp[,-1])[nfam$genome%in%gs$V2[gs$polyploid& gs$ploidy!='Paleotetraploid']], rowSums(nbp[,-1])[nbp$genome%in%gs$V2[!gs$polyploid& gs$ploidy!='Paleotetraploid']])

## make summaries of age
ages=genomecount %>% dplyr::group_by(genome, Method) %>% dplyr::filter(grepl('LTR', type)) %>% dplyr::summarize(copies=dplyr::n(), meanage=mean(as.numeric(Identity), na.rm=T), meanage90=mean(as.numeric(Identity)[as.numeric(Identity)>0.9], na.rm=T), nIdentical=sum(as.numeric(Identity)==1, na.rm=T), nYoung=sum(as.numeric(Identity)>0.99, na.rm=T)) %>% data.frame()                                  
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



pdf(paste0('~/transfer/te_panand_fig.', Sys.Date(), '.pdf'), 15,4)
                                      
bp1=ggplot(gs, aes(x=ploidy, y=haploidRepeatSize, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=4e9) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('TE base pairs') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
pr2=ggplot(gs, aes(x=ploidy, y=haploidRepeatSize/haploidAssemblySize*100, color=ploidy)) + geom_boxplot(outlier.shape=NA) +geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=100) + 
                                      scale_color_manual(values=ploidycolors) + xlab('Ploidy') + ylab('TE proportion')+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
fa3=ggplot(nfam, aes(x=ploidy, y=totFam, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=210000) + 
                                      scale_color_manual(values=ploidycolors) + xlab('Ploidy') + ylab('TE family count')+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
fa4=ggplot(nfam100, aes(x=ploidy, y=totFam, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=23000) + 
                                      scale_color_manual(values=ploidycolors) + xlab('Ploidy') + ylab('TE family count (>100 copies)')+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ag5=ggplot(gs, aes(x=ploidy, y=meanage, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid",label.y=0.89) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('LTR sequence identity') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ag590=ggplot(gs, aes(x=ploidy, y=meanage90, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid",label.y=0.89) + 
                                      scale_color_manual(values=ploidycolors, name='Ploidy') + xlab('Ploidy') + ylab('LTR sequence identity') +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
id6=ggplot(gs, aes(x=ploidy, y=propYoung, color=ploidy)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) +ggpubr::stat_compare_means(label = 'p.signif', show.legend = F,ref.group = "Diploid",label.y=0.011) +  
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

