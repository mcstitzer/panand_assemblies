library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(rtracklayer)
library(stringr)
library(tidyr)
library(dplyr)

all=read.table('../panand_sp_ploidy.txt')
all$V3[all$V2=='rtuber']=2
all$V3[all$V2=='telega']=2

anchors=lapply(all$V2, function(x) {
  a=import.gff3(paste0('../helixer_annotations/',all$V1[all$V2==x], '_helixer.gff'))
  a$genome=x
  a=a[a$type=='gene',]
  return(a)
})

ab=do.call('c', anchors)


genedists=data.frame(genome=all$V2, meandist=NA, mediandist=NA)    
gdd=vector(mode = "list", length = length(genedists$genome))
names(gdd)=genedists$genome


### got to get rid of alll the helixer mess - filter by orthogroup size?
og=read.table('../genespace/orthofinder/Results_Nov14/Orthogroups/Orthogroups.tsv', header=T, sep='\t')
# Assuming your data frame is called 'og'
# Step 1: Convert the data frame to long format
og_long <- og %>%
  pivot_longer(
    cols = -Orthogroup, # Keep Orthogroup as the identifier
    names_to = "species", # Column name for species
    values_to = "genes" # Column name for genes
  )

# Step 2: Separate the comma-separated genes into individual rows
og_melted <- og_long %>%
  separate_rows(genes, sep = ", ") # Split by comma and space

# View the result
head(og_melted)

a=read.table('../genespace/orthofinder/Results_Nov14/Orthogroups/Orthogroups.GeneCount.tsv', header=T)

b=a[rowSums(a[,2:39]==0)==0 & rowSums(a[,2:39]>12)==0,]
## simpler
b=a[a$Total>40 & a$Total<200,]

filtog=og_melted[og_melted$Orthogroup%in%b$Orthogroup,]



for(x in all$V2){
  a2=ab[ab$genome==x,]
  a2$gene=paste0(a2$genome, '_', a2$ID)
  a2=a2[a2$gene%in%filtog$genes,]
  print(x)
  print(summary(mcols(distanceToNearest(a2, ignore.strand=T))$distance))
  genedists$meandist[genedists$genome==x]=mean(mcols(distanceToNearest(a2, ignore.strand=T))$distance)
  genedists$mediandist[genedists$genome==x]=median(mcols(distanceToNearest(a2, ignore.strand=T))$distance)
  genedists$sddist[genedists$genome==x]=sd(mcols(distanceToNearest(a2, ignore.strand=T))$distance)
  # genedists$meandistNoalt[genedists$genome==x]=mean(mcols(distanceToNearest(a2na))$distance)
  # genedists$mediandistNoalt[genedists$genome==x]=median(mcols(distanceToNearest(a2na))$distance)
  # genedists$sddistNoalt[genedists$genome==x]=sd(mcols(distanceToNearest(a2na))$distance)
  # 
  gddtemp=data.frame(genome=x, dists=mcols(distanceToNearest(a2, ignore.strand=T))$distance, gene=a2$gene[queryHits(distanceToNearest(a2, ignore.strand=T))])
  gdd[[x]]=gddtemp

}

gdda=do.call(rbind, gdd)

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
summary(lm(medGGdist~haploidRepeatSize-haploidNCount, data=asize[asize$ploidy!='Paleotetraploid',]))


#### finalish
ggplot(asize, aes(x=haploidRepeatSize/(haploidAssemblySize-haploidNCount), y=medGGdist)) + geom_point(aes( color=ploidy), size=3) + scale_color_manual(values=ploidycolors) + #geom_text(aes(label=V2, color=ploidy)) + 
  stat_smooth(data=asize[asize$ploidy!='Paleotetraploid',], method='lm', lty='longdash', color='gray', se=F)+ 
  stat_smooth(data=asize, method='lm',color='gray', se=F) + ylab('Median Gene-Gene Distance') + xlab('Genomic Repeat Proportion')

