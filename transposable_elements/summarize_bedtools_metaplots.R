library(data.table)
library(dplyr)
a=fread('Ac-Pasquet1232-DRAFT-PanAnd-1.0_coverage.txt')
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
a$id=sub("_[0-9]+$", "", a$V4)
g=read.table('Ac-Pasquet1232-DRAFT-PanAnd-1.0_genes_only.gff')
a$window=sub(".*_([0-9]+)$", "\\1", a$V4)
a$window=as.numeric(a$window)
b=a%>% group_by(id)%>% filter(n()==400, strand=='+') %>% group_by(window) %>% summarize(med=median(V8))
a$window[duplicated(a$window)]=a$window[duplicated(a$window)]+200
b=a%>% group_by(id)%>% filter(n()==400, strand=='+') %>% group_by(window) %>% summarize(med=median(V8))




## next need to filter the orthogroup ids by light orthogroup sizes
## and run with syntenic anchors

## also run with CDS start/end