
library(dplyr)
library(tidyr)
library(stringr)
library(tidyr)
## not quite right, but gets general picture of scaffolds

all=read.table('../panand_sp_ploidy.txt', header=F)

a=read.table('achine-Pv-4', header=T)

b=a[a$gene!='interanchor',]



## building on the assumption that most chromosomes are syntenic...
## jsut get the biggest n syntenic scaffolds to assign kmers! then add restt of scaffolds back in
## where n is subgenome count (ploidy, or 2*ploidy if allelic)

ploidy=2
grouped=b %>% group_by(refChr, queryChr) %>% summarize(n=n()) %>% top_n(n=ploidy) %>% data.frame() %>% filter(grepl('Chr', refChr))

out=grouped %>% group_by(refChr) %>% summarize(first=first(queryChr), second=nth(queryChr, 2))

## get all contigs with synteny! 
singletons=unique(b$queryChr)[! unique(b$queryChr) %in% grouped$queryChr]

output=c(paste(out$first, out$second, sep='\t'), singletons)

write.table(output, 'achine_subphaserinput.aw.txt', col.names=F, row.names=F, quote=F)
