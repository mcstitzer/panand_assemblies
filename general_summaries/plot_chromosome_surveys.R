library(stringr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(RColorBrewer)
library(dplyr)

a=read.table('tdacs1-Pv-2', header=T)
a=a[a$gene!='interanchor',]
a$genecount=as.numeric(table(a$gene)[a$gene])
d=read.table('genomes/Td-FL_9056069_6-REFERENCE-PanAnd-2.0a.fasta.fai', header=F)
a$maxChrLen=d$V2[match(a$queryChr, d$V1)]

tand=read.table('~/transfer/tandem_repeats_panand_rm_from_genomecount.txt', header=T)    
tand$rl=str_split_fixed(tand$Name,'_',4)[,3]
tand$class=NA
tand$class[tand$rl%in%c('180bp', '179bp')]='knob180'
tand$class[tand$rl%in%c('359.5bp', '359bp')]='knobTr1'
tand$class[tand$rl%in%c('156bp','155bp')]='centromere'       
tand=tand[!is.na(tand$class),]


a$queryBin100k=round(a$queryStart,digits=-5)
a$queryBin=round(a$queryStart,digits=-6)
aa=a %>% group_by(queryBin, queryChr, maxChrLen) %>% summarize(mean(genecount))

pdf('~/transfer/tripsacum_chrs.pdf', 12,8)
ggplot(a[substr(a$queryChr,1,1)=='c',], aes(x=queryStart, y=queryChr, color=factor(genecount))) + geom_point(pch='|') + geom_segment(aes(y=queryChr, yend=queryChr, x=0, xend=maxChrLen), color='black') 
ggplot(aa[substr(aa$queryChr,1,1)=='c',], aes(x=queryBin, y=`mean(genecount)`)) + facet_wrap(~queryChr, ncol=1) + geom_point() + geom_segment(aes(y=0, yend=0, x=0, xend=maxChrLen), color='black') + scale_color_brewer(palette='Dark2')
ggplot(a[substr(a$queryChr,1,1)=='c',], aes(x=queryStart, y=queryChr, color=factor(genecount))) + geom_point(pch='|') + geom_segment(aes(y=queryChr, yend=queryChr, x=0, xend=maxChrLen), color='black') + geom_point(data=tand[tand$genome=='tdacs1' & substr(tand$seqnames,1,1)=='c',], aes(x=start, y=seqnames, color=class), pch='|') + scale_color_brewer(palette='Dark2')
dev.off()


