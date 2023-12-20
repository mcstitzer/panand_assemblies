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
chrlen=b %>% group_by(queryChr) %>% summarize(xmax=max(referenceEnd))
sgb73$xmax= chrlen$xmax[match(sgb73$queryChr,chrlen$queryChr)]
sgb73$referenceStart=NA
sgb73$queryStart=NA


## try both v pasp?
tp=read.table("../tdacs1-Pv-2", header=T)
mp=read.table('../zmB735-Pv-4', header=T)
tp=tp[tp$gene!='interanchor',]
tp=tp[!grepl('localAlignment', tp$gene),]
tp=tp[tp$queryChr %in% paste0('chr',1:18) & tp$refChr %in% c(paste0('Chr0', 1:9), 'Chr10'),]
mp=mp[mp$gene!='interanchor',]
mp=mp[!grepl('localAlignment', mp$gene),]
mp=mp[mp$queryChr %in% paste0('chr',1:10) & mp$refChr %in% c(paste0('Chr0', 1:9), 'Chr10'),]



pdf('~/transfer/B73v5-tripsacumS.onepathSG.pdf', 16,9)
##  trip on x, maize on y
ggplot(b, aes(x=referenceStart/1e6, y=queryStart/1e6))  + facet_grid(factor(queryChr, levels=c(paste0('chr', 1:10)))~factor(refChr, levels=c(paste0('chr', 1:18))),space='free', scales='free', switch='y')  + xlab('Tripsacum dactyloides chromosome position (Mb)') + ylab('Zea mays chromosome position (Mb)')+theme(strip.text.y.left = element_text(angle = 0), axis.text.x=element_text(angle=45, size=8), axis.text.y=element_text(size=8)) + 
       geom_rect(data=sgb73, inherit.aes=F, aes(xmin=0, xmax=xmax/1e6, ymin=V2/1e6, ymax=V3/1e6, fill=V13), alpha=0.1) + geom_point()+ scale_fill_manual(values=c('red', 'blue'))

ggplot(b, aes(x=referenceStart/1e6, y=queryStart/1e6))  + facet_grid(factor(queryChr, levels=c(paste0('chr', 1:10)))~factor(refChr, levels=c(paste0('chr', c(1,5,8,2,9,14,17,7,10,12,13,3,6,16,4,15,11,18)))),space='free', scales='free', switch='y')  + xlab('Tripsacum dactyloides chromosome position (Mb)') + ylab('Zea mays chromosome position (Mb)')+theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0), axis.text.x=element_text(angle=45, size=8), axis.text.y=element_text(size=8)) + 
       geom_rect(data=sgb73, inherit.aes=F, xmax=Inf, aes(xmin=0,  ymin=V2/1e6, ymax=V3/1e6, fill=V13), alpha=0.1) + geom_point()+ scale_fill_manual(values=c('red', 'blue'))



ggplot(d, aes(y=referenceStart/1e6, x=queryStart/1e6))  + facet_grid(factor(refChr, levels=c(paste0('chr', 1:10)))~factor(queryChr, levels=c(paste0('chr', 1:18))),space='free', scales='free', switch='y')  + xlab('Tripsacum dactyloides chromosome position (Mb)') + ylab('Zea mays chromosome position (Mb)')+theme(strip.text.y.left = element_text(angle = 0), axis.text.x=element_text(angle=45, size=8), axis.text.y=element_text(size=8)) + 
       geom_rect(data=sgb73, inherit.aes=F, xmax=Inf, aes(xmin=0, ymin=V2/1e6, ymax=V3/1e6, fill=V13), alpha=0.1) + geom_point()+ scale_fill_manual(values=c('red', 'blue'))


## color these by genespace paspalum chromsoomes??
# gsb=paleotetraploid$plotData$sourceData$blocks
# gsb[gsb$genome2=='zB73v5',]
# gsb[gsb$genome2=='tdacs1',]
## plot as bars above facets??


mpp=ggplot(mp, aes(x=referenceStart/1e6, y=queryStart/1e6))+ facet_grid(factor(queryChr, levels=paste0('chr', 1:10))~factor(refChr, levels=c(paste0('Chr0', 1:9), 'Chr10')),space='free', scales='free', switch='y')  + xlab('Paspalum vaginatum chromosome position (Mb)') + ylab('Zea mays chromosome position (Mb)')+theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0), axis.text.x=element_text(angle=45, size=8), axis.text.y=element_text(size=8)) + 
       geom_rect(data=sgb73, inherit.aes=F, xmax=Inf, aes(xmin=0,  ymin=V2/1e6, ymax=V3/1e6, fill=V13), alpha=0.1) + geom_point()+ scale_fill_manual(values=c('red', 'blue'))
tpp=ggplot(tp, aes(x=referenceStart/1e6, y=queryStart/1e6))+ facet_grid(factor(queryChr, levels=c(paste0('chr', 1:18)))~factor(refChr, levels=c(paste0('Chr0', 1:9), 'Chr10')),space='free', scales='free', switch='y')  + xlab('Paspalum vaginatum chromosome position (Mb)') + ylab('Tripsacum dactyloides chromosome position (Mb)')+theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0), axis.text.x=element_text(angle=45, size=8), axis.text.y=element_text(size=8)) + 
    #   geom_rect(data=sgb73, inherit.aes=F, xmax=Inf, aes(xmin=0,  ymin=V2/1e6, ymax=V3/1e6, fill=V13), alpha=0.1) + 
           geom_point()+ scale_fill_manual(values=c('red', 'blue'))

plot_grid(mpp, tpp, ncol=1, align='hv', labels='auto')
dev.off()


