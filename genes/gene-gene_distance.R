library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(rtracklayer)
library(stringr)

all=read.table('../panand_sp_ploidy.txt')
all$V3[all$V2=='rtuber']=2
all$V3[all$V2=='telega']=2

anchors=lapply(all$V2, function(x) {
  a=read.table(paste0('../syntenic_anchors/anchors/',x, '-Pv-', 2*all$V3[all$V2==x]), header=T)
  a$genome=x
  return(a)
})

ab=Reduce(function(...) merge(..., all=T), anchors)

b=ab[ab$gene!='interanchor',] ## keep only genes

genedists=data.frame(genome=all$V2, meandist=NA, mediandist=NA)    
gdd=vector(mode = "list", length = length(genedists$genome))
names(gdd)=genedists$genome

gddna=gdd

synt=read.table('~/Downloads/sharedSyntenicAnchors.txt', header=F)



for(x in all$V2){
  bg=b[b$genome==x,]
#  a1=makeGRangesFromDataFrame(bg[!duplicated(bg$gene),1:3], seqnames.field='refChr')
  a2=makeGRangesFromDataFrame(bg[,4:6], seqnames.field='queryChr')
  a2$gene=bg$gene
  a2$blockIndex=bg$blockIndex
  a2synt=makeGRangesFromDataFrame(bg[bg$gene %in%paste0(synt$V1,'.1.v3.1'),4:6], seqnames.field='queryChr')
  a2synt$gene=bg$gene[bg$gene %in%paste0(synt$V1,'.1.v3.1')]
    #print(summary(mcols(distanceToNearest(a1))$distance))
  print(x)
  print(summary(mcols(distanceToNearest(a2))$distance))
  genedists$meandist[genedists$genome==x]=mean(mcols(distanceToNearest(a2))$distance)
  genedists$mediandist[genedists$genome==x]=median(mcols(distanceToNearest(a2))$distance)
  genedists$sddist[genedists$genome==x]=sd(mcols(distanceToNearest(a2))$distance)
  # genedists$meandistNoalt[genedists$genome==x]=mean(mcols(distanceToNearest(a2na))$distance)
  # genedists$mediandistNoalt[genedists$genome==x]=median(mcols(distanceToNearest(a2na))$distance)
  # genedists$sddistNoalt[genedists$genome==x]=sd(mcols(distanceToNearest(a2na))$distance)
  # 
  gddtemp=data.frame(genome=x, dists=mcols(distanceToNearest(a2))$distance, gene=a2$gene)
  gdd[[x]]=gddtemp
  
  gddtemp=data.frame(genome=x, dists=mcols(distanceToNearest(a2synt))$distance, gene=a2synt$gene[!is.na(nearest(a2synt))])
  gddna[[x]]=gddtemp
}

gdda=do.call(rbind, gdd)
gddasynt=do.call(rbind, gddna)
gddasynt=gddasynt[!gddasynt$genome %in% c('tdacn2', 'tdacs2', 'pprate'),]

# flow=read.table('panand_flow.txt', header=T)      
# all$flow=flow$flow[match(paste0(str_split_fixed(all$V1,'-',3)[,1], '-', str_split_fixed(all$V1,'-',3)[,2]), flow$genome)]
# 
# 
# genedists$flow=all$flow[match(genedists$genome, all$V2)]         

gdda$ploidy=asize$ploidy[match(gdda$genome, asize$V2)]
gdda$haploidSize=(asize$haploidAssemblySize-asize$haploidNCount)[match(gdda$genome, asize$V2)]
ggplot(gdda, aes(x=dists+1, y=reorder(genome, desc(haploidSize)), fill=ploidy)) + geom_density_ridges()  + scale_x_log10() + scale_fill_manual(values=ploidycolors)+ geom_vline(xintercept=c(101,1001,10001,100001, 1000001))

gddasynt$ploidy=asize$ploidy[match(gddasynt$genome, asize$V2)]
gddasynt$haploidSize=(asize$haploidAssemblySize-asize$haploidNCount)[match(gddasynt$genome, asize$V2)]
gddasynt$medianKs=asize$medianKs[match(gddasynt$genome, asize$V2)]

tempDist=gddasynt %>% group_by(genome) %>% summarize(medGGdist=median(dists))
asize$medGGdist=tempDist$medGGdist[match(asize$V2, tempDist$genome)]

## plot relationship between genome size and gene-gene distance!!
ggplot(asize, aes(x=haploidAssemblySize-haploidNCount, y=medGGdist)) + geom_point(aes( color=ploidy)) + scale_color_manual(values=ploidycolors) + geom_text(aes(label=V2, color=ploidy)) + stat_smooth(data=asize[asize$ploidy!='Paleotetraploid',], method='lm', lty='longdash', color='gray', se=F)
summary(lm(medGGdist~haploidAssemblySize-haploidNCount, data=asize[asize$ploidy!='Paleotetraploid',]))

## they're all 5 kb!!!
asize %>% group_by(ploidy) %>% summarize(medgg=median(medGGdist))

### got to make a asizezt for this one!
ggplot(asize, aes(x=haploidRepeatSize-haploidNCount, y=medGGdist)) + geom_point(aes( color=ploidy)) + scale_color_manual(values=ploidycolors) + geom_text(aes(label=V2, color=ploidy)) + stat_smooth(data=asize[asize$ploidy!='Paleotetraploid',], method='lm', lty='longdash', color='gray', se=F)+ stat_smooth(data=asizezt, method='lm',color='gray', se=F)
summary(lm(medGGdist~c(haploidRepeatSize/(haploidAssemblySize-haploidNCount)), data=asize[asize$ploidy!='Paleotetraploid',]))


#### finalish
rpgg=ggplot(asize, aes(x=haploidRepeatSize/(haploidAssemblySize-haploidNCount), y=medGGdist)) + geom_point(aes( color=ploidy), size=3) + scale_color_manual(values=ploidycolors) + #geom_text(aes(label=V2, color=ploidy)) + 
  stat_smooth(data=asize[asize$ploidy!='Paleotetraploid',], method='lm', lty='longdash', color='gray', se=F)+ 
#  stat_smooth(data=asize, method='lm',color='gray', se=F) + 
  ylab('Median Syntenic\nGene-Gene Distance') + xlab('Repeat Proportion') + theme(legend.position='NULL')
rpgg #repeat proportion gg

ggplot(asize, aes(x=haploidRepeatSize, y=medGGdist)) + geom_point(aes( color=ploidy), size=3) + scale_color_manual(values=ploidycolors) + #geom_text(aes(label=V2, color=ploidy)) + 
  stat_smooth(data=asize[asize$ploidy!='Paleotetraploid',], method='lm', lty='longdash', color='gray', se=F)+ 
  #  stat_smooth(data=asize, method='lm',color='gray', se=F) + 
  ylab('Median Syntenic\nGene-Gene Distance') + xlab('Repeat Bases') + theme(legend.position='NULL') + geom_
##10k genes, 

gsgg=ggplot(asize, aes(x=(haploidAssemblySize-haploidNCount), y=medGGdist)) + geom_point(aes( color=ploidy), size=3) + scale_color_manual(values=ploidycolors) + #geom_text(aes(label=V2, color=ploidy)) + 
  stat_smooth(data=asize[asize$ploidy!='Paleotetraploid',], method='lm', lty='longdash', color='gray', se=F)+ 
  #  stat_smooth(data=asize, method='lm',color='gray', se=F) + 
  ylab('Median Syntenic\nGene-Gene Distance') + xlab('Genome Size') + theme(legend.position='NULL') 
gsgg ## genome size gene dist
##

ggplot(asize[!asize$V2%in%lowQualAssemblies,], aes(x=medFrac, y=medGGdist)) + geom_point(aes( color=ploidy), size=3) + scale_color_manual(values=ploidycolors) + #geom_text(aes(label=V2, color=ploidy)) + 
  stat_smooth(data=asize[asize$ploidy!='Paleotetraploid' & !asize$V2%in%lowQualAssemblies,], method='lm', lty='longdash', color='gray', se=F)+ 
  stat_smooth(data=asize[!asize$V2%in%lowQualAssemblies,], method='lm',color='gray', se=F) + 
  ylab('Median Syntenic\nGene-Gene Distance') + xlab('Retained Gene Proportion') + theme(legend.position='NULL')


## maybe suplement?? ridgeline of diff gene-gene dists
ggplot(gddasynt, aes(x=dists+1, y=genome, fill=ploidy)) + geom_density_ridges()  + scale_x_log10() + scale_fill_manual(values=ploidycolors)+ geom_vline(xintercept=c(101,1001,10001,100001, 1000001)) 



ggplot(gddasynt, aes(x=dists+1, y=reorder(genome, desc(haploidSize)), fill=ploidy)) + geom_density_ridges()  + scale_x_log10() + scale_fill_manual(values=ploidycolors)+ geom_vline(xintercept=c(101,1001,10001,100001, 1000001))

ggplot(gddasynt, aes(x=dists+1, y=reorder(genome, desc(haploidSize)), fill=ploidy, color='black'))   + scale_x_log10() + scale_fill_manual(values=ploidycolors)+ geom_vline(xintercept=c(101,1001,10001,100001, 1000001)) +stat_density_ridges(quantile_lines = TRUE, linetype='solid', color='lightgray')

ggplot(gddasynt, aes(x=dists+1, y=reorder(genome, desc(medianKs)), fill=ploidy))   + scale_x_log10(limits=c(11,1e7+1)) + scale_fill_manual(values=ploidycolors)+ geom_vline(xintercept=c(11,101,1001,10001,100001, 1000001,10000001), lty='dotted', color='lightgray')+ geom_density_ridges()+ ylab('Genome, sorted by increasing ks between duplicate genes') + xlab('Distance Between Syntenic Genes') 


ggplot(gddasynt, aes(x=dists+1, y=reorder(genome, desc(haploidSize)), fill=ploidy))   + scale_x_log10(limits=c(11,1e7+1)) + scale_fill_manual(values=ploidycolors)+ geom_vline(xintercept=c(11,101,1001,10001,100001, 1000001,10000001), lty='dotted', color='lightgray')+ geom_density_ridges() + ylab('Genome, sorted by increasing haploid genome size') + xlab('Distance Between Syntenic Genes') 


taxonnames=c("Z. mays ssp. parviglumis TIL11", "Z. mays ssp. mays B73v5", "Z. mays ssp. parviglumis TIL01", "Z. mays ssp. mexicana TIL25", "Z. mays ssp. mexicana TIL18", "Z. mays ssp. huehuetenangensis", 
             "Z. luxurians", "Z. nicaraguensis", "Z. diploperennis Momo", "Z. diploperennis Gigi", "T. zoloptense", "T. dactyloides FL", "T. dactyloides Southern Hap2", 
             "T. dactyloides Northern Hap2", "T. dactyloides KS", "T. dactyloides tetraploid", "U. digitatum", "V. cuspidata", "R. rottboellioides", "R. tuberculosa", 
             "H. compressa", "E. tripsacoides", "S. scoparium", "S. microstachyum", "A. virginicum", "A. chinensis", "A. gerardi", 
             "C. refractus", "C. citratus", "H. contortus", "T. triandra", "B. laguroides", "P. paniceum", "S. bicolor", 
             "I. rugosum", "S. nutans", '"A." burmanicus', "T. elegans", "C. serrulatus", "P. vaginatum")
names(taxonnames)=c("zTIL11", "zmB735", "zTIL01", "zTIL25", "zTIL18", "zmhuet", 
                    "zluxur", "znicar", "zdmomo", "zdgigi", "tzopol", "tdacs1", "tdacs2", 
                    "tdacn2", "tdacn1", "tdactm", "udigit", "vcuspi", "rrottb", "rtuber", 
                    "hcompr", "etrips", "sscopa", "smicro", "avirgi", "achine", "agerar", 
                    "crefra", "ccitra", "hconto", "ttrian", "blagur", "ppanic", "sbicol", 
                    "irugos", "snutan", "atenui", "telega", "cserru", "pvagin")
gddasynt$taxon=taxonnames[match(gddasynt$genome, names(taxonnames))]
ggplot(gddasynt, aes(x=dists+1, y=reorder(taxon, desc(medianKs)), fill=ploidy))   + scale_x_log10(limits=c(11,1e7+1)) + scale_fill_manual(values=ploidycolors)+ geom_vline(xintercept=c(11,101,1001,10001,100001, 1000001,10000001), lty='dotted', color='lightgray')+ geom_density_ridges()+ ylab('Genome, sorted by increasing ks between duplicate genes') + xlab('Distance Between Syntenic Genes') 
## phyloorder
ggplot(gddasynt, aes(x=dists+1, y=factor(taxon, levels=rev(taxonnames)), fill=ploidy))   + scale_x_log10(limits=c(11,1e7+1)) + scale_fill_manual(values=ploidycolors)+ geom_vline(xintercept=c(11,101,1001,10001,100001, 1000001,10000001), lty='dotted', color='lightgray')+ geom_density_ridges()+ ylab('Genome, sorted by phylogeny') + xlab('Distance Between Syntenic Genes') 



closeby=gddasynt%>%group_by(gene)%>%summarize(prop=sum(dists<5e2)/n(), count=n())
cb=run_GO_analysis(gsub('.1.v3.1', '.1', closeby[closeby$prop>0.95,]$gene), genomeid='closeby', background=syntbackground, GODB=GODB[names(GODB)%in%syntbackground])

faraway=gddasynt%>%group_by(gene)%>%summarize(prop=sum(dists>5e4)/n(), count=n())
fa=run_GO_analysis(gsub('.1.v3.1', '.1', faraway[faraway$prop>0.9,]$gene), genomeid='faraway', background=syntbackground, GODB=GODB[names(GODB)%in%syntbackground])


closeby=gddasynt%>%filter(dists<5e2)
cb=run_GO_analysis(gsub('.1.v3.1', '.1', closeby$gene), genomeid='closeby', background=syntbackground, GODB=GODB[names(GODB)%in%syntbackground])

## so 500 bp away is stress and there's a terpenoid
tmp=factor(as.integer(background%in%gsub('.1.v3.1', '.1', closeby$gene)))
names(tmp)=background
tgd1=new( "topGOdata", ontology="BP", allGenes = tmp, nodeSize=5,annot=annFUN.gene2GO, gene2GO = GODB)
allGO=genesInTerm(tgd1)
lapply(allGO[cb$GO.ID[1:4]] , function(x) x[x%in%gsub('.1.v3.1', '.1', closeby$gene)])




faraway=gddasynt%>%filter(dists>1e5)
fa=run_GO_analysis(gsub('.1.v3.1', '.1', faraway$gene), genomeid='faraway', background=syntbackground, GODB=GODB[names(GODB)%in%syntbackground])


#genedists=rbind(genedists, data.frame(genome='pvagin', meandist=mean(mcols(distanceToNearest(a1))$distance), mediandist=median(mcols(distanceToNearest(a1))$distance), sddist=sd(mcols(distanceToNearest(a1))$distance),
#                                      meandistNoalt=mean(mcols(distanceToNearest(a1))$distance), mediandistNoalt=median(mcols(distanceToNearest(a1))$distance), sddistNoalt=sd(mcols(distanceToNearest(a1))$distance),flow=593)) ## they say 593 Mb in the genome paper




medianRepeat=median(asize$haploidRepeatSize[asize$ploidy=='Diploid']-asize$haploidNCount[asize$ploidy=='Diploid'])

summary(lm((asize$haploidRepeatSize-asize$haploidNCount)~asize$medianKs))
ggplot(asize, aes(x=medianKs, y=(haploidRepeatSize-haploidNCount)/1e6, color=ploidy)) + geom_point(size=3) + scale_color_manual(values=ploidycolors) + theme(legend.position = 'NA')+ stat_smooth(aes(group=1), method = 'lm', se=F) + xlab('Median ks between gene duplicates') + ylab('Repeat sequence (Mb)') + geom_hline(yintercept=c(1,2,3)*medianRepeat/1e6, color=ploidycolors[1:3])
ggplot(asize, aes(x=medianKs, y=(haploidRepeatSize-haploidNCount)/1e6, color=ploidy)) + geom_point(size=3) + scale_color_manual(values=ploidycolors) + theme(legend.position = 'NA')+ stat_smooth(aes(group=1), method = 'lm', se=F) + xlab('Median ks between gene duplicates') + ylab('Repeat sequence (Mb)') + geom_hline(yintercept=c(1,2,3)*medianRepeat/1e6, color=ploidycolors[1:3]) + geom_text(aes(label=V2))

ggplot(asize, aes(x=medianKs, y=(haploidRepeatSize-haploidNCount)/(haploidAssemblySize-haploidNCount), color=ploidy)) + geom_point(size=3) + scale_color_manual(values=ploidycolors) + theme(legend.position = 'NA') + stat_smooth(aes(group=1), method = 'lm', se=F)



## ask whether fractionation distance bigger or smaller in homeologs!!

## groups of three in b, to see which are 2,1,2
### b needs to be the presence/absence table :(
a2$copynumber=table(a2$gene)[a2$gene]

posn=gregexpr('212', paste0(b$etrips,collapse=''))[[1]]
indices <- unlist(lapply(posn, function(x) x:(x+2)))
subset_b <- b[indices,]
subset_a2 <- a2[a2$gene %in% subset_b$gene,]
subset_by_copy_number <- subset_a2[subset_a2$copynumber == c(2,2, 1,2, 2),]

results_list <- list()
for (i in seq(1, n, by = 5)) {
if ((i + 4) <= n) {
group <- subset_by_copy_number[i:(i + 4)]
chrs=unique(seqnames(group))
# Step 5: Calculate min and max start positions for the group
result <- data.frame(
  seqnames = chrs,
  range = sapply(chrs, function(x) max(start(group[seqnames(group)==x,]))-min(start(group[seqnames(group)==x,]))),
  copynumber = sapply(chrs, function(x) paste(group[seqnames(group)==x,]$copynumber, collapse=',')),
  genes = paste(unique(group$gene), collapse = ", "),
  stringsAsFactors = FALSE
)
# Append to the results list
if(length(unique(group$blockIndex))==2){
results_list[[length(results_list) + 1]] <- result
}
}}
results_df <- do.call(rbind, results_list)




a2synt$pvChr=substr(a2synt$gene,1,7)
a2synt$geneIndex=as.numeric(factor(a2synt$gene))


fractBias=data.frame(a2synt) %>% group_by(window=geneIndex%/%50, pvChr, seqnames) %>% summarize(fractBias=n()/50)
fractTotal=data.frame(a2synt) %>% group_by(geneIndex, pvChr) %>% summarize(geneCount=n()) %>% group_by(window=geneIndex%/%50, pvChr) %>% summarize(fractTotal=sum(geneCount)/50)

## maize
##fractBias=fractBias[fractBias$seqnames%in%paste0('chr',1:10),]
ggplot(fractBias, aes(x=window, y=fractBias, group=seqnames, color=seqnames)) + geom_line() + facet_wrap(~pvChr, ncol=1)
ggplot(fractBias, aes(x=window, y=fractBias, group=seqnames, color=seqnames)) + geom_line() 
ggplot(fractBias[fractBias$pvChr=='Pavag01',], aes(x=window, y=fractBias, group=seqnames, color=seqnames)) + geom_line() 

## total
ggplot(fractTotal[fractTotal$pvChr=='Pavag01',], aes(x=window, y=fractTotal)) + geom_line() 


## try for every genome!!!
pdf('~/Downloads/fractionationBias.pdf',8,4)
for(chr in c(paste0('Pavag0', 1:9), 'Pavag10')){
for(x in all$V2){
  bg=b[b$genome==x,]
  #  a1=makeGRangesFromDataFrame(bg[!duplicated(bg$gene),1:3], seqnames.field='refChr')
  a2=makeGRangesFromDataFrame(bg[,4:6], seqnames.field='queryChr')
  a2$gene=bg$gene
  a2$blockIndex=bg$blockIndex
  a2synt=makeGRangesFromDataFrame(bg[bg$gene %in%paste0(synt$V1,'.1.v3.1'),4:6], seqnames.field='queryChr')
  a2synt$gene=bg$gene[bg$gene %in%paste0(synt$V1,'.1.v3.1')]
  
  a2$pvChr=substr(a2$gene,1,7)
  a2$geneIndex=as.numeric(factor(a2$gene))
  
  a2synt$pvChr=substr(a2synt$gene,1,7)
  a2synt$geneIndex=as.numeric(factor(a2synt$gene))
  fractBias=data.frame(a2synt) %>% group_by(window=geneIndex%/%100, pvChr, seqnames) %>% summarize(fractBias=n()/100)
  print(ggplot(fractBias[fractBias$pvChr==chr,], aes(x=window, y=fractBias, group=seqnames, color=seqnames)) + geom_line() + ggtitle(paste(chr, x)) + theme(legend.position='NA'))
}}
dev.off()  

library(zoo)









pdf('~/Downloads/fractionationBias.rolling.pdf',8,4)
for(chr in c(paste0('Pavag0', 1:9), 'Pavag10')){
  for(x in all$V2){
    bg=b[b$genome==x,]
    #  a1=makeGRangesFromDataFrame(bg[!duplicated(bg$gene),1:3], seqnames.field='refChr')
    a2=makeGRangesFromDataFrame(bg[,4:6], seqnames.field='queryChr')
    a2$gene=bg$gene
    a2$blockIndex=bg$blockIndex
    a2synt=makeGRangesFromDataFrame(bg[bg$gene %in%paste0(synt$V1,'.1.v3.1'),4:6], seqnames.field='queryChr')
    a2synt$gene=bg$gene[bg$gene %in%paste0(synt$V1,'.1.v3.1')]
    
    a2$pvChr=substr(a2$gene,1,7)
    a2$geneIndex=as.numeric(factor(a2$gene))
    
    a2synt$pvChr=substr(a2synt$gene,1,7)
    a2synt$geneIndex=as.numeric(factor(a2synt$gene))
    fractBias=data.frame(a2synt) %>% group_by(window=geneIndex%/%100, pvChr, seqnames) %>% summarize(fractBias=n()/100)
    # Convert synt to a vector
    synt_genes <- synt$V1
    
    # Extract the relevant gene identifier from bg
    bg <- bg %>%
      mutate(gene_id = sub("\\..*", "", gene))  # Removing version numbers from gene IDs
    
    # Define a custom function to count genes in synt
    count_synt_genes <- function(x, synt_genes) {
      sum(x %in% synt_genes)
    }
    
    # Apply the rolling window function grouped by queryChr
    bg <- bg %>%
      group_by(queryChr) %>%
      arrange(referenceStart) %>%
      mutate(synt_gene_count = rollapply(gene_id, width = 100, FUN = function(x) count_synt_genes(x, synt_genes), fill = NA, align = "right")) %>%
      ungroup()
    bg$pvChr=substr(bg$gene,1,7)
    
    bg$fractBias=bg$synt_gene_count/100
    bg <- bg %>%
      mutate(gene_id = sub("\\..*", "", gene)) %>%  # Ensure gene_id is stripped of version numbers
      mutate(geneIndex = match(gene_id, synt_genes))
    
    # View the result
    head(bg)
    print(
      ggplot(bg[bg$pvChr==chr,], aes(x=geneIndex, y=fractBias, group=queryChr, color=queryChr)) + geom_line() + ggtitle(paste(chr, x)) + theme(legend.position='NA')
      )
  }}
dev.off()  

#### that was all genes, now switch to syntenic!!!!!!!!!!

results_list=vector(mode='list', length=length(all$V2)) #empty_list <- vector(mode = "list", length = desired_length)
names(results_list)=all$V2

pdf('~/Downloads/fractionationBias.rollingSyntenic.pdf',8,4)
  for(x in all$V2){
    bg=b[b$genome==x,]
    #  a1=makeGRangesFromDataFrame(bg[!duplicated(bg$gene),1:3], seqnames.field='refChr')
    a2=makeGRangesFromDataFrame(bg[,4:6], seqnames.field='queryChr')
    a2$gene=bg$gene
    a2$blockIndex=bg$blockIndex
    a2synt=makeGRangesFromDataFrame(bg[bg$gene %in%paste0(synt$V1,'.1.v3.1'),4:6], seqnames.field='queryChr')
    a2synt$gene=bg$gene[bg$gene %in%paste0(synt$V1,'.1.v3.1')]
    
    a2$pvChr=substr(a2$gene,1,7)
    a2$geneIndex=as.numeric(factor(a2$gene))
    
    a2synt$pvChr=substr(a2synt$gene,1,7)
    a2synt$geneIndex=as.numeric(factor(a2synt$gene))
    fractBias=data.frame(a2synt) %>% group_by(window=geneIndex%/%100, pvChr, seqnames) %>% summarize(fractBias=n()/100)
    # Convert synt to a vector
    synt_genes <- synt$V1
    
    # Extract the relevant gene identifier from bg
    bg <- bg %>%
      mutate(gene_id = sub("\\..*", "", gene)) %>% # Removing version numbers from gene IDs
      filter(gene_id%in%synt_genes)
    # Define a custom function to count bg genes in each window of synt_genes and track additional information
    count_bg_genes_with_index <- function(synt_window, synt_index, bg_df, query_chr) {
      bg_subset <- bg_df %>% filter(queryChr == query_chr)
      count <- sum(bg_subset$gene_id %in% synt_window)
      first_gene <- synt_window[1]
      return(c(count = count, first_gene_index = synt_index, first_gene = first_gene))
    }
    
    # Apply the rolling window function on synt_genes and store results
    result <- data.frame()
    
    for (chr in unique(bg$queryChr)) {
      rolling_results <- rollapply(1:length(synt_genes), width = 100, FUN = function(idx) {
        synt_window <- synt_genes[idx]
        count_bg_genes_with_index(synt_window, idx[1], bg, chr)
      }, by.column = FALSE, fill = NA, align = "right")
      
      rolling_results_df <- as.data.frame(rolling_results)
      rolling_results_df$queryChr <- chr
      rolling_results_df$synt_window_start_index <- as.numeric(rolling_results_df$first_gene_index)
      rolling_results_df$first_gene <- as.character(rolling_results_df$first_gene)
      rolling_results_df$count <- as.numeric(rolling_results_df$count)
      
      result <- rbind(result, rolling_results_df)
    }
    result$pvChr=substr(result$first_gene,1,7)
    
    result$fractBias=result$count/100
    results_list[[x]]=result
    # View the result
 #   head(bg)
    for(chr in c(paste0('Pavag0', 1:9), 'Pavag10')){
      
    print(
      ggplot(result[result$pvChr==chr,], aes(x=synt_window_start_index, y=fractBias, group=queryChr, color=queryChr)) + geom_line() + ggtitle(paste(chr, x)) + theme(legend.position='NA')
    )
  }}
dev.off()  

outlist=results_list
for(i in 1:length(results_list)){
  temp=results_list[[i]]
  temp$genome=all$V2[i]
  outlist[[i]]=temp
}

out=do.call(rbind, outlist)

pdf('~/Downloads/fractionationBias.rollingSyntenic.pdf',8,4)
for(x in 1:length(all$V2)){
  for(chr in c(paste0('Pavag0', 1:9), 'Pavag10')){
    
    print(
      ggplot(results_list[[x]][results_list[[x]]$pvChr==chr,], aes(x=synt_window_start_index, y=fractBias, group=queryChr, color=queryChr)) + geom_line() + ggtitle(paste(chr, all$V2[x])) + theme(legend.position='NA')
    )
  }}

chr='Pavag01'
print(
  ggplot(out[out$pvChr==chr,], aes(x=synt_window_start_index, y=fractBias, group=queryChr, color=queryChr)) + geom_line() + ggtitle(chr) + theme(legend.position='NA') + facet_wrap(~genome)
)

dev.off()


pdf('~/Downloads/fractionationBias.rollingSyntenicByChr.pdf',20,20)

outFiltered=out %>% filter(count!=0) %>% group_by(queryChr, pvChr, genome) %>% mutate(min=min(as.numeric(first_gene_index), na.rm=T)+100, max=max(as.numeric(first_gene_index), na.rm=T)-100)
outFiltered= outFiltered %>% filter(as.numeric(first_gene_index)<=max & as.numeric(first_gene_index)>=min)
outFiltered= outFiltered %>% filter(genome!='pprate')

for(chr in c(paste0('Pavag0', c(1,2,3,4,6,7,8,9)), 'Pavag10')){
  print(
    ggplot(outFiltered[outFiltered$pvChr==chr,], aes(x=synt_window_start_index, y=fractBias, group=queryChr, color=queryChr)) + xlim(min(max(out[out$pvChr==chr,]$synt_window_start_index)+100),max(out[out$pvChr==chr,]$synt_window_start_index)-100)+ geom_line() + ggtitle(chr) + theme(legend.position='NA') + facet_wrap(~genome, ncol=10)
  )
}
  dev.off()


  
  outFiltered$ploidy=asize$ploidy[match(outFiltered$genome, asize$V2)]
  
  outFiltered %>% filter(pvChr=='Pavag01') %>% ggplot(aes(x=fractBias, y=factor(genome, levels=rev(names(taxonnames))), fill=ploidy)) + geom_density_ridges() + scale_fill_manual(values=ploidycolors)
  outFiltered %>% ggplot(aes(x=fractBias, y=factor(genome, levels=rev(names(taxonnames))), fill=ploidy)) + geom_density_ridges() + scale_fill_manual(values=ploidycolors)
  

  off=outFiltered %>% group_by(ploidy, genome) %>% summarize(m=median(fractBias)) 
  off$medianKs=asize$medianKs[match(off$genome, asize$V2)]
  ## older polyploidies are more fractionated
  ggplot(off, aes(x=medianKs, y=m, color=ploidy)) + geom_point() + scale_color_manual(values=ploidycolors) + geom_text(aes(label=genome))
  
########

  ## test for biased fractionation - rank alleles by dNdS and calculate average dnds in windows along each subgenome??
  
  ks=fread('~/Downloads/ksp_to_pasp.txt.gz')
  outFiltered$dnds=ks$V19[match(paste(outFiltered$genome, outFiltered$first_gene), paste(ks$genome, ks$gene))]
  
  ggplot(outFiltered[outFiltered$pvChr==chr,], aes(y=dnds, x=fractBias, color=queryChr)) + ylim(0,1.5) + geom_point() + ggtitle(chr) + theme(legend.position='NA') + facet_wrap(~genome, ncol=10)
  
  
  ggplot(outFiltered[outFiltered$pvChr==chr,], aes(y=dnds, x=fractBias, color=queryChr)) + ylim(0,1.5) + geom_point() + ggtitle(chr) + theme(legend.position='NA') + facet_wrap(~genome, ncol=10)
  
  ## by chrom, look at median fractionation and median dnds - are these correlated
  mof=outFiltered %>% group_by(queryChr, genome) %>% summarize(fractBias=median(fractBias), dnds=median(dnds, na.rm=T), ngenes=n())
  ggplot(mof[mof$ngenes>100,], aes(y=dnds, x=fractBias))  + geom_point()  + theme(legend.position='NA') + facet_wrap(~genome, ncol=10)
  ggplot(mof[mof$ngenes>100,], aes(y=dnds, x=fractBias))  + geom_point()  + theme(legend.position='NA') + facet_wrap(~genome, ncol=10, scales='free')
## weak negative relationship - more genes retained, lower dNdS....
    cor.test(mof$dnds[mof$ngenes>100], mof$fractBias[mof$ngenes>100])
    
    ## get chatgpt to make me a table... sigs are mostly negative?? bothriochloa, hconto, udigit, crefra is weird - pretty positive
    ## ccitra hcompr, irugos neg
    ## could this be a thing of chromsome length????
    for(i in unique(mof$genome)){
      print(i)
      print(cor.test(mof$dnds[mof$ngenes>100 & mof$genome==i], mof$fractBias[mof$ngenes>100 & mof$genome==i]))
    }
    
  
  ####
##### ended up not using the below...

library(parallel)


pdf('~/Downloads/fractionationBias.rollingSyntenicParallel.pdf',8,4)
# Parallel processing setup

# Parallel processing setup
num_cores <- detectCores() - 1  # Use all but one core
cl <- makeCluster(num_cores)
clusterEvalQ(cl, {
  library(dplyr)
  library(zoo)
  library(GenomicRanges)
  library(data.table)
})

process_genome <- function(x, synt, b) {
  bg <- b[b$genome == x,]
  
  if (nrow(bg) == 0) return(NULL)  # Skip if no data for the genome
  
  # Create GRanges objects
  a2 <- makeGRangesFromDataFrame(bg[, 4:6], seqnames.field = 'queryChr')
  a2$gene <- bg$gene
  a2$blockIndex <- bg$blockIndex
  
  a2synt <- makeGRangesFromDataFrame(bg[bg$gene %in% paste0(synt$V1, '.1.v3.1'), 4:6], seqnames.field = 'queryChr')
  a2synt$gene <- bg$gene[bg$gene %in% paste0(synt$V1, '.1.v3.1')]
  
  a2$pvChr <- substr(a2$gene, 1, 7)
  a2$geneIndex <- as.numeric(factor(a2$gene))
  
  a2synt$pvChr <- substr(a2synt$gene, 1, 7)
  a2synt$geneIndex <- as.numeric(factor(a2synt$gene))
  
  fractBias <- data.frame(a2synt) %>%
    group_by(window = geneIndex %/% 100, pvChr, seqnames) %>%
    summarize(fractBias = n() / 100)
  
  # Filter bg to include only synt genes
  synt_genes <- synt$V1
  bg <- bg %>%
    mutate(gene_id = sub("\\..*", "", gene)) %>%
    filter(gene_id %in% synt_genes)
  
  if (length(synt_genes) < 100) return(NULL)  # Skip if not enough genes for a window
  
  # Function to count bg genes in each window of synt_genes
  count_bg_genes_with_index <- function(synt_window, synt_index, bg_df, query_chr) {
    bg_subset <- bg_df %>% filter(queryChr == query_chr)
    if (nrow(bg_subset) == 0) {
      first_gene <- if (length(synt_window) > 0) synt_window[1] else NA
      return(list(count = 0, first_gene_index = synt_index, first_gene = first_gene))
    }
    
    count <- sum(bg_subset$gene_id %in% synt_window)
    first_gene <- if (length(synt_window) > 0) synt_window[1] else NA
    return(list(count = count, first_gene_index = synt_index, first_gene = first_gene))
  }
  
  # Apply the rolling window function on synt_genes and store results
  result <- data.table()
  
  for (chr in unique(bg$queryChr)) {
    if (length(unique(bg$gene_id[bg$queryChr == chr])) < 100) next  # Skip if not enough unique genes
    
    rolling_results <- tryCatch({
      rollapply(1:length(synt_genes), width = 100, FUN = function(idx) {
        synt_window <- synt_genes[idx]
        count_bg_genes_with_index(synt_window, idx[1], bg, chr)
      }, by.column = FALSE, fill = NA, align = "right")
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(rolling_results)) next
    
    rolling_results_df <- rbindlist(rolling_results, use.names = TRUE, fill = TRUE)
    rolling_results_df[, queryChr := chr]
    
    result <- rbind(result, rolling_results_df)
  }
  
  if (nrow(result) > 0) {
    result[, pvChr := ifelse(!is.na(first_gene), substr(first_gene, 1, 7), NA)]
    result[, fractBias := count / 100]
    result[, genome := x]  # Add genome identifier to the result
  } else {
    return(NULL)
  }
  
  return(result)
}

# Export necessary objects to cluster
clusterExport(cl, list("synt", "b", "process_genome"))

# Apply the function across all genomes in parallel
results_list <- parLapply(cl, all$V2, function(x) {
  process_genome(x, synt, b)
})

stopCluster(cl)

# Combine results and plot
results <- rbindlist(results_list, use.names = TRUE, fill = TRUE)

# Check if 'pvChr' exists before plotting
if (!"pvChr" %in% names(results)) {
  stop("Column 'pvChr' not found in results")
}

for (genome in unique(results$genome)) {
  genome_results <- results[results$genome == genome, ]
  
  for (chr in c(paste0('Pavag0', 1:9), 'Pavag10')) {
    chr_results <- genome_results[genome_results$pvChr == chr, ]
    
    if (nrow(chr_results) > 0) {
      print(
        ggplot(chr_results, aes(x = first_gene_index, y = fractBias, group = queryChr, color = queryChr)) +
          geom_line() +
          ggtitle(paste("Genome:", genome, "| Chromosome:", chr)) +
          theme(legend.position = 'NA')
      )
    }
  }
}

dev.off()






