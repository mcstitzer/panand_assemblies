library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

all=read.table('../panand_sp_ploidy.txt', header=F)

## read in each blast zdgigi.DUFs.tblastn.txt

blast=lapply(all$V2, function(x) {
 a=read.table(paste0(x, '.DUFs.tblastn.txt'), header=F)
 a$genome=x
 return(a)
 })


ab=Reduce(function(...) merge(..., all=T), blast)


abb=ab[ab$V11<1e-100,]

pdf('duf_tblastn_plot.pdf',14,10)
for(genome in unique(abb$genome)){
print(
  ggplot(abb[abb$genome==genome,], aes(x=V9, xend=V10, y=factor(V1), yend=factor(V1), color=factor(V1), linewidth=as.numeric(V3))) + geom_segment() + facet_wrap(~V2, ncol=1, scales='free') + ggtitle(genome)
  )
  }
 dev.off()
          
