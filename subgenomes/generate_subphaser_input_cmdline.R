
library(dplyr)
library(tidyr)
library(stringr)
library(tidyr)
## not quite right, but gets general picture of scaffolds

args <- commandArgs(trailingOnly = TRUE)
## anchorwaveanchorfile ploidy(numbercolumns) faifile
aw=read.table(args[1], header=T)
ploidy=as.numeric(args[2])
fai=read.table(args[3], header=F)

b=aw[aw$gene!='interanchor',]


## building on the assumption that most chromosomes are syntenic...
## jsut get the biggest n syntenic scaffolds to assign kmers! then add restt of scaffolds back in
## where n is subgenome count (ploidy, or 2*ploidy if allelic)

  ## for now, not worrying about allelic assemblies and hoping chromosomes will be random length by subgenome (not all longest will be alleles or something)

if(ploidy==2){
  
grouped=b %>% group_by(refChr, queryChr) %>% summarize(n=n()) %>% top_n(n=ploidy) %>% data.frame() %>% filter(grepl('Chr', refChr))
out=grouped %>% group_by(refChr) %>% summarize(first=first(queryChr), second=nth(queryChr, 2))
## get all contigs with synteny! 
#singletons=unique(b$queryChr)[! unique(b$queryChr) %in% grouped$queryChr]
#singletons=unique(fai$V1[!fai$V1%in%grouped$queryChr])

output=c(paste(out$first, out$second, sep='\t'))#, singletons)
}
  if(ploidy==3){
    grouped=b %>% group_by(refChr, queryChr) %>% summarize(n=n()) %>% top_n(n=ploidy) %>% data.frame() %>% filter(grepl('Chr', refChr))
out=grouped %>% group_by(refChr) %>% summarize(first=first(queryChr), second=nth(queryChr, 2), third=nth(queryChr,3))
## get all contigs with synteny! 
#singletons=unique(b$queryChr)[! unique(b$queryChr) %in% grouped$queryChr]
#singletons=unique(fai$V1[!fai$V1%in%grouped$queryChr])

output=c(paste(out$first, out$second, out$third, sep='\t'))#, singletons)
    }
  


if(ploidy==4){
  
grouped=b %>% group_by(refChr, queryChr) %>% summarize(n=n()) %>% top_n(n=ploidy) %>% data.frame() %>% filter(grepl('Chr', refChr))
out=grouped %>% group_by(refChr) %>% summarize(first=first(queryChr), second=nth(queryChr, 2), third=nth(queryChr,3), fourth=nth(queryChr,4))
## get all contigs with synteny! 
#singletons=unique(b$queryChr)[! unique(b$queryChr) %in% grouped$queryChr]

#singletons=unique(fai$V1[!fai$V1%in%grouped$queryChr])
output=c(paste(out$first, out$second, out$third, out$fourth, sep='\t'))#, singletons)
}

if(ploidy==5){
  
grouped=b %>% group_by(refChr, queryChr) %>% summarize(n=n()) %>% top_n(n=ploidy) %>% data.frame() %>% filter(grepl('Chr', refChr))
out=grouped %>% group_by(refChr) %>% summarize(first=first(queryChr), second=nth(queryChr, 2), third=nth(queryChr,3), fourth=nth(queryChr,4), fifth=nth(queryChr,5))
## get all contigs with synteny! 
#singletons=unique(b$queryChr)[! unique(b$queryChr) %in% grouped$queryChr]

#singletons=unique(fai$V1[!fai$V1%in%grouped$queryChr])
output=c(paste(out$first, out$second, out$third, out$fourth, out$fifth, sep='\t'))#, singletons)
}

  if(ploidy==6){
    grouped=b %>% group_by(refChr, queryChr) %>% summarize(n=n()) %>% top_n(n=ploidy) %>% data.frame() %>% filter(grepl('Chr', refChr))
out=grouped %>% group_by(refChr) %>% summarize(first=first(queryChr), second=nth(queryChr, 2), third=nth(queryChr,3), fourth=nth(queryChr,4), fifth=nth(queryChr,5), sixth=nth(queryChr,6))
## get all contigs with synteny! 
## add in step to get the unplaced scaffolds?!?
#bl=read.table('blagur/blagur/04.build/blagur_aggressivecorrection.FINAL.fa.fai', header=F)
#singletons=unique(b$queryChr)[! unique(b$queryChr) %in% grouped$queryChr]
#singletons=unique(fai$V1[!fai$V1%in%grouped$queryChr])

output=c(paste(out$first, out$second, out$third, out$fourth, out$fifth, out$sixth, sep='\t'))#, singletons)
    }
  


fileout=str_split_fixed(args[1], '-', 2)[,1]

write.table(output, paste0(fileout, '_subphaserinput.txt'), col.names=F, row.names=F, quote=F)

