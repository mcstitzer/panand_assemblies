library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(rtracklayer)
library(stringr)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(ape)
library(ggtree)


all=read.table('../panand_sp_ploidy.txt')
all=all[!all$V2 %in% c('pprate', 'tdactm', 'tzopol', 'osativ', 'bdista', 'agerjg', 'svirid', 'eophiu'),]

trash=lapply(all$V2, function(x) {
 a=read.csv(paste0(x, '/Summary.of.repetitive.regions.', all$V1[all$V2==x], '.fasta.csv'), header=T)
 a$genome=x
 return(a)
 })


at=do.call(rbind, trash)
at=at[at$consensus.primary!='none_identified' & at$consensus.count!=0,]

atl=at[at$most.freq.value.N>40,]

## by species, what are the repeats
ros=atl %>% group_by(consensus.primary, genome, most.freq.value.N) %>% summarize(n=n(), repeats=sum(repeats.identified)) %>% arrange(desc(repeats)) %>% filter(repeats>100)
ros$totbp=ros$most.freq.value.N*ros$repeats

ros%>% arrange(desc(totbp)) 

## my filters - remove repeats with sequence not found at least five positions in the genome
tanr=ros[ros$n>5,] ## n is positions per species


rout=list()

## write a fasta of unique sequences for each class - there will be some single snp diffs, but for now okay
for(genome in all$V2){
 to=atl[atl$genome==genome & atl$consensus.primary %in% tanr$consensus.primary[tanr$genome==genome],]
 if(nrow(to)>0){ ## filtered away everything for c serrulatus
 to=to %>% group_by(consensus.primary, most.freq.value.N) %>% summarize(bp=sum(width), repeats.identified=sum(repeats.identified))
# to=to[!duplicated(to$consensus.primary),]
 rownames(to)=1:nrow(to)
 # go=data.frame(seqid=to$name, source='TRASH', type='satellite_DNA', start=to$start+1, end=to$end, score=to$ave.score, strand='.', phase='.', ## start +1 for gff format
 #               attributes=paste0('ID=TandemRepeat', rownames(to), ';Name=', to$most.freq.value.N, 'bp_repeat;Note=consensus_',to$consensus.primary))
 # write.table(go, paste0(genome, '_filteredTandemRepeats.TRASH.gff3'), row.names=F, col.names=F, sep='\t', quote=F)
  to$genome=genome
  rout[[genome]]=to
 fa=data.frame(name=paste0('>', genome, '_TandemRepeat', rownames(to), '_', to$most.freq.value.N, 'bp_repeat'), seq=to$consensus.primary)
   write.table(fa, paste0(genome, '_filteredTandemRepeats.TRASH.fa'), row.names=F, col.names=F, sep='\n', quote=F)
  }
 }

trs=do.call(rbind, rout)
write.table(trs, 'tandem_repeats_by_species.txt', row.names=F, col.names=T, sep='\t', quote=F)


