library(ggplot2) ## cbsu or scinet - but need 96 cpu to run in parallel
library(cowplot)
theme_set(theme_cowplot())
library(stringr)
library(tidyr)
library(dplyr)
library(reshape2)

### make sure on scinet you've loaded modules first!!
## module load samtools
## module load bedtools
## module load r ## to get r to work on scinet

all=data.frame(V2=c('taesti'), V3=c(3))
reffa="../Phytozome/PhytozomeV13/Pvaginatum/v3.1/assembly/Pvaginatum_672_v3.0.fa"
hap1="Triticum_aestivum.IWGSC.dna.toplevel.fa"

all$fasta=c(hap1)


anchors=lapply(all$V2, function(x) {
 a=read.table(paste0(x, '-Pv-', all$V3[all$V2==x]), header=T)
 a$genome=x
 return(a)
 })

ab=Reduce(function(...) merge(..., all=T), anchors)

b=ab[ab$gene!='interanchor',] ## keep only genes
b$quickgene=substr(b$gene,1,14)


## just do them all, but note this includes genes with very little conservation
          
faf=b %>% group_by(quickgene) %>% summarize(n=n())
#write.table(faf, 'paspalum_gene_tree_tipcount.txt', sep='\t', quote=F, row.names=F, col.names=T)


library(doParallel)
detectCores()
registerDoParallel(9)
          
## start with 100 genes in gene tree??
foreach(gene=faf$quickgene[faf$n==3])  %dopar% { #[faf$n==100]){ #### for virginicus included!
#foreach(gene=faf$quickgene[faf$n==6])  %dopar% { #[faf$n==100]){

 #for(gene in faf$quickgene){ #[faf$n==100]){
for(sp in unique(b$genome)){

bed=b[b$quickgene==gene & b$genome==sp,c('queryChr', 'queryStart', 'queryEnd', 'quickgene', 'strand', 'genome')]
if(nrow(bed)!=0){
bed[,4]=paste0(bed$genome, '_', bed$quickgene,'_', bed$queryChr, '_', bed$queryStart, '-', bed$queryEnd)
bed$genome=NULL
}

outfile=paste0('gene_tree_beds/', gene, '_', sp, '.bed')

write.table(bed, outfile, row.names=F, col.names=F, sep='\t', quote=F)
system(paste0('bedtools getfasta -s -nameOnly -fi ', all$fasta[all$V2==sp], ' -bed ', outfile, ' >> gene_tree_unalignedfa/', gene, '.fa'))
#### bedtools getfasta -s -fi test.fa -bed test.bed
}
 ## also do paspalum genomic region!
 sp='pvagin'
 bed=b[b$quickgene==gene,c('refChr', 'referenceStart', 'referenceEnd', 'quickgene', 'strand', 'genome')]
if(nrow(bed)!=0){
 bed$referenceStart=min(bed$referenceStart)
 bed$referenceEnd=max(bed$referenceEnd)
bed[,4]=paste0('pvagin', '_', bed$quickgene,'_', bed$refChr, '_', bed$referenceStart, '-', bed$referenceEnd)
bed$genome=NULL
 bed=bed[1,]
}
outfile=paste0('gene_tree_beds/', gene, '_', sp, '.bed')
write.table(bed, outfile, row.names=F, col.names=F, sep='\t', quote=F)
system(paste0('bedtools getfasta -s -nameOnly -fi ../Phytozome/PhytozomeV13/Pvaginatum/v3.1/assembly/Pvaginatum_672_v3.0.fa -bed ', outfile, ' >> gene_tree_unalignedfa/', gene, '.fa'))

 print(paste0(gene, ' the ', which(faf$quickgene==gene), ' gene'))
}