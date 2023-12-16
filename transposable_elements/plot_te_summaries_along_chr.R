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



  agc=unlist(reduce(split(as_granges(genomecount, keep_mcols=T), ~c(Name,genome)))) ## merge bookended or overlapping TRs if they're the same repeat consensus


## simple count per genome
gcn=genomecount %>% group_by(genome, Classification) %>% summarize(nfam=length(unique(Name)), ncopy=n(), nbp=sum(width))
gcng=dcast(gcn, genome~Classification, value.var='nfam')
## remove na that are nam specific cols
gcng=gcng[,!is.na(gcng[1,])]
nfam=data.frame(genome=gcng$genome, DTA=gcng$`DNA/DTA`+gcng$`MITE/DTA`, DTC=gcng$`DNA/DTC`+gcng$`MITE/DTC`, DTH=gcng$`DNA/DTH`+gcng$`MITE/DTH`,
                  DTM=gcng$`DNA/DTM`+gcng$`MITE/DTM`, DTT=gcng$`DNA/DTT`+gcng$`MITE/DTT`, DHH=gcng$`DNA/Helitron`, 
                  RLC=gcng$`LTR/Copia`, RLG=gcng$`LTR/Gypsy`, RLX=gcng$`LTR/unknown`)
## ncopy
gcng=dcast(gcn, genome~Classification, value.var='ncopy')
## remove na that are nam specific cols
gcng=gcng[,!is.na(gcng[1,])]
ncopy=data.frame(genome=gcng$genome, DTA=gcng$`DNA/DTA`+gcng$`MITE/DTA`, DTC=gcng$`DNA/DTC`+gcng$`MITE/DTC`, DTH=gcng$`DNA/DTH`+gcng$`MITE/DTH`,
                  DTM=gcng$`DNA/DTM`+gcng$`MITE/DTM`, DTT=gcng$`DNA/DTT`+gcng$`MITE/DTT`, DHH=gcng$`DNA/Helitron`, 
                  RLC=gcng$`LTR/Copia`, RLG=gcng$`LTR/Gypsy`, RLX=gcng$`LTR/unknown`)
##nbp
gcng=dcast(gcn, genome~Classification, value.var='nbp')
## remove na that are nam specific cols
gcng=gcng[,!is.na(gcng[1,])]
nbp=data.frame(genome=gcng$genome, DTA=gcng$`DNA/DTA`+gcng$`MITE/DTA`, DTC=gcng$`DNA/DTC`+gcng$`MITE/DTC`, DTH=gcng$`DNA/DTH`+gcng$`MITE/DTH`,
                  DTM=gcng$`DNA/DTM`+gcng$`MITE/DTM`, DTT=gcng$`DNA/DTT`+gcng$`MITE/DTT`, DHH=gcng$`DNA/Helitron`, 
                  RLC=gcng$`LTR/Copia`, RLG=gcng$`LTR/Gypsy`, RLX=gcng$`LTR/unknown`)

## get statistics for the paper!
gs=read.table('../panand_')
polyploids=gs$V2[gs$ploidy!='Diploid']
          
t.test(rowSums(nfam[,-1])[nfam$genome%in%polyploids], rowSums(nfam[,-1])[!nfam$genome%in%polyploids])
t.test(rowSums(ncopy[,-1])[nfam$genome%in%polyploids], rowSums(ncopy[,-1])[!ncopy$genome%in%polyploids])
t.test(rowSums(nbp[,-1])[nfam$genome%in%polyploids], rowSums(nbp[,-1])[!nbp$genome%in%polyploids])

         
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
