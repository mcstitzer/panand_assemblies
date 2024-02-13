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
## is specific contig a chromsoomes?!?!??!
at10$chromLevel=grepl('chr', at10$V1)
at$chromLevel=grepl('chr', at$V1)

gs=read.table('~/transfer/panand_assembly_sizes.txt', header=T, sep='\t')

## plot nicely :)
taxonnames=c("Z. mays ssp. parviglumis TIL11", "Z. mays ssp. mays B73v5", "Z. mays ssp. parviglumis TIL01", "Z. mays ssp. mexicana TIL25", "Z. mays ssp. mexicana TIL18", "Z. mays ssp. huehuetengensis", 
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

gs$species=taxonnames[match(gs$V2, names(taxonnames))]
gs$species=factor(gs$species, levels=taxonnames)
gs$speciesLabel=ifelse(gs$haploid, paste0(gs$species, '*'), as.character(gs$species))
gs$speciesLabel=gs$species
levels(gs$speciesLabel)[levels(gs$speciesLabel) %in% gs$species[gs$haploid]]=paste0(levels(gs$speciesLabel)[levels(gs$speciesLabel) %in% gs$species[gs$haploid]], '*')

ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')



at10$ploidy=gs$ploidy[match(at10$genome, gs$V2)]
at$ploidy=gs$ploidy[match(at$genome, gs$V2)]

at10$speciesLabel=gs$speciesLabel[match(at10$genome, gs$V2)]
at$speciesLabel=gs$speciesLabel[match(at$genome, gs$V2)]


## switch zeas back to "diploid" for chromsoome numbers becasue they effectively are
all$V3[grepl('z', all$V2)]=1

pdf('~/transfer/panand_explore_chrs.pdf', 16,10)
for(i in 1:nrow(all)){
print(
  ggplot(at10[at10$genome==all$V2[i],], aes(x=V2, y=reorder(V1, V2), color=ploidy)) + geom_point() + ggtitle(paste0(all$V2[i], '-', ifelse(gs$haploid[gs$V2==all$V2[i]], 'haploid', 'allelic'))) + geom_vline(xintercept=c(0, 416000000)) + geom_vline(data=gs[gs$V2==all$V2[i],], aes(xintercept=ifelse(haploid, haploidAssemblySize/(10*all$V3[i]), rawAssemblySize/(20*all$V3[i]))), lty='dashed')+ scale_color_manual(values=ploidycolors)
  )
  }

ggplot(at10[!at10$genome %in% c('tdacn2', 'tdacs2'),], aes(x=V2, y=speciesLabel, color=ploidy)) + geom_point(aes(shape=chromLevel)) + geom_vline(xintercept=c(0, 416000000)) + scale_color_manual(values=ploidycolors)
ggplot(at[!at$genome %in% c('tdacn2', 'tdacs2'),], aes(x=V2, y=speciesLabel, color=ploidy)) + geom_point(aes(shape=chromLevel)) + geom_vline(xintercept=c(0, 416000000))+ scale_color_manual(values=ploidycolors)

ggplot(at[!at$genome %in% c('tdacn2', 'tdacs2'),], aes(x=V2, y=speciesLabel, color=ploidy)) + geom_point(aes(shape=chromLevel)) + geom_vline(xintercept=c(0, 416000000))+ scale_color_manual(values=ploidycolors) + geom_point(data=gs, aes(x=ifelse(haploid, haploidAssemblySize/(10*all$V3[i]), rawAssemblySize/(20*all$V3[i]))), pch='|', color='black', size=3)


#ggplot(a[substr(a$queryChr,1,1)=='c',], aes(x=queryStart, y=queryChr, color=factor(genecount))) + geom_point(pch='|') + geom_segment(aes(y=queryChr, yend=queryChr, x=0, xend=maxChrLen), color='black') 
dev.off()

pdf('~/transfer/supplemental_assembly_contiguity.pdf', 16,10)
ggplot(at[!at$genome %in% c('tdacn2', 'tdacs2'),], aes(x=V2, y=speciesLabel, color=ploidy)) + geom_point() + geom_vline(xintercept=c(0, 416000000))+ scale_color_manual(values=ploidycolors) + geom_point(data=gs, aes(x=ifelse(haploid, haploidAssemblySize/(10*all$V3[i]), rawAssemblySize/(20*all$V3[i]))), pch='|', color='black', size=3) + xlab('Contig Length (bp)')
dev.off()

