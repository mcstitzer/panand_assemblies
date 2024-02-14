library(rtracklayer)
library(stringr)
library(plyranges)
library(dplyr)
library(tidyverse)

all=read.table('../panand_sp_ploidy.txt')
all=all[!all$V2 %in% c('pprate', 'tdactm', 'tzopol', 'osativ', 'bdista', 'agerjg', 'svirid', 'eophiu'),]

yw=lapply(all$V2, function(x) {
 a=import.gff(paste0(x, '.ywminiprot.gff'))
 a=a[a$type=='mRNA',]
  ## eventually deal with intersections (rank by score)
# a=
 a$genome=x
 return(a)
 })


at=do.call(c, yw)
at$name=str_split_fixed(at$Target, ' ', 3)[,1]


tegroup=read.table('yuan_wessler_2011.supnames.txt', header=T)
at$sup=tegroup$superfamily[match(at$name, tegroup$fasta_entry)]

data.frame(at) %>% dplyr::group_by(sup, genome) %>% dplyr::summarize(n=n()) %>% pivot_wider(names_from=sup, values_from=n, values_fill=0)
