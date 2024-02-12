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



## plot big contigs of each assembly
all=read.table('../panand_sp_ploidy.txt')
all=all[!all$V2 %in% c('pprate', 'tdactm', 'tzopol', 'osativ', 'bdista', 'agerjg', 'svirid', 'eophiu'),]

clen=lapply(all$V2, function(x) {
# a=read.table(paste0(x, '-Pv-', 2*all$V3[all$V2==x]), header=T)
  a=read.table(paste0('genomes/', all$V1[all$V2==x], '.fasta.fai'), header=F)
 a$genome=x
 return(a)
 })


at=do.call(rbind, clen)

at10=at[at$V2>10e6,]
at10=at10 %>% arrange(-V2)


pdf('~/transfer/panand_explore_chrs.pdf', 16,10)
for(i in 1:nrow(all)){
print(
  ggplot(at10[at10$genome==all$V2[i],], aes(x=V2, y=reorder(V1, V2))) + geom_point() + ggtitle(all$V2[i]) + geom_vline(xintercept=c(0, 416000000))
  )
  }

#ggplot(a[substr(a$queryChr,1,1)=='c',], aes(x=queryStart, y=queryChr, color=factor(genecount))) + geom_point(pch='|') + geom_segment(aes(y=queryChr, yend=queryChr, x=0, xend=maxChrLen), color='black') 
dev.off()
