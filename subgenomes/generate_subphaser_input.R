
library(dplyr)
library(tidyr)
library(stringr)
library(tidyr)
## not quite right, but gets general picture of scaffolds

all=read.table('../panand_sp_ploidy.txt', header=F)

for( i in 1:nrow(all)){
if(file.exists(paste0('../',all$V2[i],'-Pv-', all$V3[i]*2))){
a=read.table(paste0('../',all$V2[i],'-Pv-', all$V3[i]*2), header=T)

b=a[a$gene!='interanchor',]



## building on the assumption that most chromosomes are syntenic...
## jsut get the biggest n syntenic scaffolds to assign kmers! then add restt of scaffolds back in
## where n is subgenome count (ploidy, or 2*ploidy if allelic)

  ## for now, not worrying about allelic assemblies and hoping chromosomes will be random length by subgenome (not all longest will be alleles or something)

ploidy=all$V3[i]
if(ploidy==2){
  
grouped=b %>% group_by(refChr, queryChr) %>% summarize(n=n()) %>% top_n(n=ploidy) %>% data.frame() %>% filter(grepl('Chr', refChr))
out=grouped %>% group_by(refChr) %>% summarize(first=first(queryChr), second=nth(queryChr, 2))
## get all contigs with synteny! 
singletons=unique(b$queryChr)[! unique(b$queryChr) %in% grouped$queryChr]
output=c(paste(out$first, out$second, sep='\t'), singletons)
  write.table(output, paste0(all$V2[i],'_subphaserinput.aw.txt'), col.names=F, row.names=F, quote=F)
}
  if(ploidy==3){
    grouped=b %>% group_by(refChr, queryChr) %>% summarize(n=n()) %>% top_n(n=ploidy) %>% data.frame() %>% filter(grepl('Chr', refChr))
out=grouped %>% group_by(refChr) %>% summarize(first=first(queryChr), second=nth(queryChr, 2), third=nth(queryChr,3))
## get all contigs with synteny! 
singletons=unique(b$queryChr)[! unique(b$queryChr) %in% grouped$queryChr]
output=c(paste(out$first, out$second, out$third, sep='\t'), singletons)
    write.table(output, paste0(all$V2[i],'_subphaserinput.aw.txt'), col.names=F, row.names=F, quote=F)
    }
  

}
}
