 library(ggplot2)
 library(cowplot)
 theme_set(theme_cowplot())

## wait do i want two to two??
a=read.table('../b73_anchors/tdacs1-zmB735-2', header=T)
d=a[a$gene!='interanchor',]
d=d[!grepl('localAlignment', d$gene),]
d=d[d$queryChr %in% paste0('chr',1:18) & d$refChr %in% paste0('chr', 1:10),]


## 1:1
a=read.table('~/transfer/Td-FL-A_ZmB73V5_R1Q1_anchors', header=T)
b=a[a$gene!='interanchor',]
b=b[!grepl('localAlignment', b$gene),]
b=b[b$queryChr %in% paste0('chr',1:10) & b$refChr %in% paste0('chr', 1:18),]


## subgenomes
sg=read.table('~/transfer/panand_all_subgenomes_10052023.txt', header=F)
sgb73=sg[grepl('zB73v5', sg$V1),]
sgb73$queryChr=gsub('zB73v5\\.','', sgb73$V1)
## this was a bad idea: sgb73$referenceChr=gsub('zB73v5\\.','', sgb73$V1) ## jsut do both so i can use it for both :)
## poor mans chromosome lengths
chrlen=b %>% group_by(queryChr) %>% summarize(xmax=max(queryEnd))
sgb73$xmax= chrlen$xmax[match(sgb73$queryChr,chrlen$queryChr)]
sg$referenceStart=NA
sg$queryStart=NA

pdf('~/transfer/B73v5-tripsacumS.onepathSG.pdf', 16,9)
##  trip on x, maize on y
ggplot(b, aes(x=referenceStart/1e6, y=queryStart/1e6))  + facet_grid(factor(queryChr, levels=c(paste0('chr', 1:10)))~factor(refChr, levels=c(paste0('chr', 1:18))),space='free', scales='free', switch='y')  + xlab('Tripsacum dactyloides chromosome position (Mb)') + ylab('Zea mays chromosome position (Mb)')+theme(strip.text.y.left = element_text(angle = 0), axis.text.x=element_text(angle=45, size=8), axis.text.y=element_text(size=8)) + 
       geom_rect(data=sgb73, inherit.aes=F, aes(xmin=0, xmax=xmax/1e6, ymin=V2/1e6, ymax=V3/1e6, fill=V13), alpha=0.1) + geom_point()+ scale_fill_manual(values=c('red', 'blue'))

ggplot(d, aes(y=referenceStart/1e6, x=queryStart/1e6))  + facet_grid(factor(refChr, levels=c(paste0('chr', 1:10)))~factor(queryChr, levels=c(paste0('chr', 1:18))),space='free', scales='free', switch='y')  + xlab('Tripsacum dactyloides chromosome position (Mb)') + ylab('Zea mays chromosome position (Mb)')+theme(strip.text.y.left = element_text(angle = 0), axis.text.x=element_text(angle=45, size=8), axis.text.y=element_text(size=8)) + 
       geom_rect(data=sgb73, inherit.aes=F, aes(xmin=0, xmax=xmax/1e6, ymin=V2/1e6, ymax=V3/1e6, fill=V13), alpha=0.1) + geom_point()+ scale_fill_manual(values=c('red', 'blue'))


dev.off()


