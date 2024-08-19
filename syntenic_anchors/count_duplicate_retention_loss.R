library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(rtracklayer)
library(topGO)
library(dplyr)
#setwd('~/Documents/GitHub/panand_assemblies/syntenic_anchors/')

## count genes as boxplot, show relaxed selection


ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')



b=read.table('../syntenic_anchors/gene_anchorTable_AnchorWave_Paspalum.txt', header=T)
b=b[,-which(colnames(b)%in%c('svirid',"pprate", "osativ", "eophiu", "bdista", "agerjg", 'tdactm', 'tzopol'))]
b=b[,-which(colnames(b)%in%c('tdacs2', 'tdacn2'))]

## themeda might loose abn arm like contortus?
## serrulatus is nanopore
#b=b[,-which(colnames(b)%in%c('ttrian', 'cserru'))]


table(rowSums(b[,-1]>0))
bb=b[rowSums(b[,-1]>0)>=32,]



asize=fread('../general_summaries/panand_assembly_sizes.txt', header=T, quote="", fill=T)
asize$syntAnchors=colSums(bb[,-1])[match(asize$V2, names(colSums(bb[,-1])))]
asize$syntAnchorsCount=colSums(bb[,-1]>0)[match(asize$V2, names(colSums(bb[,-1]>0)))]
asize$syntAnchors=ifelse(asize$haploid, asize$syntAnchors, asize$syntAnchors/2)
asize$ploidy=factor(asize$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))
asize$diploidEquivalent=ifelse(asize$ploidy=='Diploid', asize$syntAnchors, ifelse(asize$ploidy%in%c('Tetraploid', 'Paleotetraploid'), asize$syntAnchors/2, ifelse(asize$ploidy=='Hexaploid', asize$syntAnchors/3, NA)))

## use this to look at ks and dnds for each!
ks=fread('../general_summaries/ks_to_look_for_mixtures.txt', header=T, quote='')
ks=data.frame(ks)[ks$ks>0.001,]
mks=ks[ks$ploidy %in%c('Tetraploid', 'Paleotetraploid', 'Hexaploid'),] %>% group_by(genome, ploidy, species, haploid) %>% dplyr::summarize(median=median(ks, na.rm=T), nonallelic=median(ks[ks>0.005], na.rm=T))
mks$medianNonAllelicCorr=ifelse(mks$nonallelic-mks$median>0.01, mks$nonallelic, mks$median)


mu=6.5e-9
mks$mya=mks$medianNonAllelicCorr/2/mu/1e6
ks$mya=ks$ks/2/mu/1e6
ggplot(data.frame(ks)[ks$ploidy %in% c('Tetraploid', 'Paleotetraploid', 'Hexaploid') & !ks$genome%in%c('zdmomo', 'zdgigi', 'zdnica'),], aes(x=mya, color=ploidy, group=genome)) + geom_density() + scale_color_manual(values=ploidycolors) + facet_wrap(~ploidy, ncol=1, scales='free_y') + geom_vline(data=mks[!mks$genome %in% c('zdmomo', 'zdgigi', 'zdnica'),], aes(xintercept=mya))


asize$medianKs=mks$medianNonAllelicCorr[match(asize$V2, mks$genome)]
asize$mya=mks$mya[match(asize$V2, mks$genome)]
asize$mya[asize$ploidy=='Diploid']=0
asize$medianKs[asize$ploidy=='Diploid']=0


## make a version where haploid individuals are doubled, so I can check by ploidy quickly
bbd=bb
bbd[,asize$V2[asize$haploid]]=bbd[,asize$V2[asize$haploid]]*2

asize$ploidyNumber=ifelse(asize$ploidy=='Diploid', 2, ifelse(asize$ploidy%in%c('Tetraploid', 'Paleotetraploid'), 4, ifelse(asize$ploidy=='Hexaploid', 6, NA)))

### for each species, which genes are kept in duplicate, which are reduced to single copy?
## what to do about hexaploids?? are 1, 2 haploid copies equivalent?

## this is whether the gene is retained in duplicate in that individual or not!
bbdd=bbd
#bbdd[,-1]=lapply(colnames(bbdd[,-1]), function(x) bbdd[,x]>= (asize$ploidyNumber[asize$V2==x]-1)) ## -1 because alleles?
bbdd[,-1]=lapply(colnames(bbdd[,-1]), function(x) bbdd[,x]>= 3) ## so just say there's more than a diploid's copies??

asize$retainedInDup=sapply(asize$V2, function(x) sum(bbdd[,x]))


ggplot(asize, aes(x=medianKs, y=retainedInDup, group=ploidy, color=ploidy)) + 
  #  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(seed = 1), size=3)+ 
  scale_color_manual(values=ploidycolors) + #ylim(0,110000) + 
  # ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=110000) + 
  #  geom_hline(yintercept=c(median(asize$diploidEquivalent[asize$ploidy=='Diploid'], na.rm=T)*c(1)), lty='dotted', color='darkgray')+
  #  annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  # annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  # annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  xlab('Ks') + ylab('Retained In Duplicate') + theme(legend.position='NULL') + geom_text(aes(label=V2))

### so now get out gene names and do go~!
## save to file for each gene set
GO_table_path='~/Downloads/Pv_GOTable.txt'
GODB <- readMappings(GO_table_path, IDsep = ",")
background <- names(GODB)

## new background from syntenic genes
syntbackground=gsub('.1.v3.1', '.1',bbdd$gene)

run_GO_analysis <- function(goi, genomeid, GO_table_path='~/Downloads/Pv_GOTable.txt', background=background, GODB=GODB) {
  # Read the GO mapping file
  # Create factor indicating genes of interest
  tmp <- factor(as.integer(background %in% goi))
  names(tmp) <- background
  # Create a topGOdata object
  tgd1 <- new("topGOdata", ontology = "BP", allGenes = tmp, nodeSize = 5, annot = annFUN.gene2GO, gene2GO = GODB)
  # Run GO enrichment tests
  resTopGO.classic <- runTest(tgd1, algorithm = "classic", statistic = "Fisher")
  resTopGO.weight01 <- runTest(tgd1, algorithm = "weight01", statistic = "Fisher")
  # Generate results table
  GO_res_table <- GenTable(tgd1, 
                           Fisher.classic = resTopGO.classic, 
                           Fisher.weight01 = resTopGO.weight01, 
                           orderBy = "Fisher.weight01", 
                           ranksOf = "Fisher.classic", 
                           topNodes = length(resTopGO.classic@score), 
                           numChar = 100)
  GO_res_table$genome=genomeid
  GO_res_table$Fisher.weight01=as.numeric(GO_res_table$Fisher.weight01)
  return(GO_res_table)
}

# Example usage:
# goi <- c("gene1", "gene2", "gene3")  # replace with your list of genes
# GO_table_path <- '~/Downloads/Pv_GOTable.txt'
# results <- run_GO_analysis(goi, GO_table_path)
# print(results)

snutdup=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[bbdd$snutan]), genomeid='snutan', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])
snutsc=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[!bbdd$snutan]), genomeid='snutan', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])


stayduplist=lapply(asize$V2[asize$ploidy!='Diploid'], function(x) run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[bbdd[,x]]), genomeid=x, background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground]))
staydup=do.call(rbind, stayduplist)
sort(table(staydup[as.numeric(staydup$Fisher.weight01)<0.05,]$Term))
## so transcription factors!!!

singlecopylist=lapply(asize$V2[asize$ploidy!='Diploid'], function(x) run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[!bbdd[,x]]), genomeid=x, background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground]))
singlecopy=do.call(rbind, singlecopylist)
sort(table(singlecopy[singlecopy$Fisher.weight01<0.05,]$Term))


## then across all (most?) species
##### oh wait we need to do just polyploids!!!
table(rowSums(bbdd[,asize$V2[asize$ploidy!='Diploid']]))
#bbdd$gene[rowSums(bbdd[,asize$V2[asize$ploidy!='Diploid']])>=20]
dupinmost=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(bbdd[,asize$V2[asize$ploidy!='Diploid']])>=24]), genomeid='dupinmost', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])
singlecopyinmost=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(!bbdd[,asize$V2[asize$ploidy!='Diploid']])>=19]), genomeid='singlecopyinmost', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])

## try just tetraploid, just hexaploid
table(rowSums(bbdd[,asize$V2[asize$ploidy=='Tetraploid']]))
dupintetra=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(bbdd[,asize$V2[asize$ploidy=='Tetraploid']])>=9]), genomeid='dupintetra', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])
singlecopyintetra=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(!bbdd[,asize$V2[asize$ploidy=='Tetraploid']])>=6]), genomeid='singlecopyintetra', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])
table(rowSums(bbdd[,asize$V2[asize$ploidy=='Hexaploid']]))
dupinhexa=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(bbdd[,asize$V2[asize$ploidy=='Hexaploid']])>=4]), genomeid='dupinhexa', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])
singlecopyinhexa=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(!bbdd[,asize$V2[asize$ploidy=='Hexaploid']])>=3]), genomeid='singlecopyinhexa', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])
table(rowSums(bbdd[,asize$V2[asize$ploidy=='Paleotetraploid']]))
dupinpaleo=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(bbdd[,asize$V2[asize$ploidy=='Paleotetraploid']])>=12]), genomeid='dupinpaleo', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])
singlecopyinpaleo=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(!bbdd[,asize$V2[asize$ploidy=='Paleotetraploid']])>=12]), genomeid='singlecopyinpaleo', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])

## now try youngest three polyploids, oldest four, middle6
table(rowSums(bbdd[,asize$V2[asize$ploidy!='Diploid' & asize$medianKs<0.03]]))
dupinyoung=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(bbdd[,asize$V2[asize$ploidy!='Diploid' & asize$medianKs<0.03]])>=3]), genomeid='dupinyoung', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])
singlecopyinyoung=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(!bbdd[,asize$V2[asize$ploidy!='Diploid' & asize$medianKs<0.03]])>=3]), genomeid='singlecopyinyoung', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])

table(rowSums(bbdd[,asize$V2[asize$ploidy!='Diploid' & asize$medianKs>0.05 & asize$medianKs<0.1]]))
dupinold=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(bbdd[,asize$V2[asize$ploidy!='Diploid' & asize$medianKs>0.05 & asize$medianKs<0.1]])>=4]), genomeid='dupinold', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])
singlecopyinold=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(!bbdd[,asize$V2[asize$ploidy!='Diploid' & asize$medianKs>0.05 & asize$medianKs<0.1]])>=4]), genomeid='singlecopyinold', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])

## from below, where i loook for ALL zea/trip to have duplicates!!
##dupinall=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(bbddz)==14]), genomeid='dupinall')


## start plotting!
dupcombo=rbind(staydup, dupinmost)
dupcombo=rbind(dupcombo, dupintetra)
dupcombo=rbind(dupcombo,dupinhexa)
dupcombo=rbind(dupcombo,dupinpaleo)
dupcombo=rbind(dupcombo, dupinold)
dupcombo=rbind(dupcombo,dupinyoung)
dupcombo=rbind(dupcombo,dupinalls)
dupcombo=rbind(dupcombo,dupintrips)
dupcombo$ploidy=asize$ploidy[match(dupcombo$genome, asize$V2)]
dupcombo$medianKs=asize$medianKs[match(dupcombo$genome, asize$V2)]
dim(dupcombo)
dupcombo$Fisher.weight01=as.numeric(dupcombo$Fisher.weight01)


data_summary <- dupcombo %>%
  group_by(Term) %>%
  summarize(Significant_Count_05 = sum(substr(genome,1,3)!='dup' & Fisher.weight01<5e-2),
            Significant_Count_B = sum(substr(genome,1,3)!='dup' & Fisher.weight01<(5e-2/26))) %>%
  arrange(desc(Significant_Count_05))

# Step 2: Merge the counts back to the original data
data <- merge(dupcombo, data_summary, by = "Term")

# Step 3: Sort the data by the significant count
data <- data %>%
  arrange(desc(Significant_Count_05))

ggplot(data[data$Fisher.weight01<(5e-2/26),], aes(y = reorder(paste(Term,GO.ID), Significant_Count_B), x=reorder(genome, medianKs), color=ploidy, size=-log10(as.numeric(Fisher.weight01)))) + geom_point() + scale_color_manual(values=ploidycolors) + theme(legend.position = 'None') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(data[data$Fisher.weight01<(5e-2) ,], aes(y = reorder(Term, Significant_Count_05), x=reorder(genome, medianKs), color=ploidy, size=-log10(as.numeric(Fisher.weight01)))) + geom_point() + scale_color_manual(values=ploidycolors) + theme(legend.position = 'None') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(data[data$Fisher.weight01<(5e-2) ,], aes(y = reorder(Term, Significant_Count_05), x=reorder(genome, medianKs), color=ploidy, size=-log10(as.numeric(Fisher.weight01)))) + geom_point() + scale_color_manual(values=ploidycolors) + theme(legend.position = 'None') + theme(axis.text.x = element_text(angle = 45, hjust = 1))



## try with single copy - less convergence!
sccombo=rbind(singlecopy, singlecopyinmost)
sccombo=rbind(sccombo, singlecopyintetra)
sccombo=rbind(sccombo,singlecopyinhexa)
sccombo=rbind(sccombo,singlecopyinpaleo)
sccombo=rbind(sccombo, singlecopyinyoung)
sccombo=rbind(sccombo, singlecopyinold)
sccombo$ploidy=asize$ploidy[match(sccombo$genome, asize$V2)]
sccombo$medianKs=asize$medianKs[match(sccombo$genome, asize$V2)]

dim(sccombo)
sccombo$Fisher.weight01=as.numeric(sccombo$Fisher.weight01)

data_summary <- sccombo %>%
  group_by(Term) %>%
  summarize(Significant_Count_05 = sum(substr(genome,1,3)!='sin' & Fisher.weight01<5e-2),
            Significant_Count_B = sum(substr(genome,1,3)!='sin' & Fisher.weight01<(5e-2/26))) %>%
  arrange(desc(Significant_Count_05))

# Step 2: Merge the counts back to the original data
data <- merge(sccombo, data_summary, by = "Term")

# Step 3: Sort the data by the significant count
data <- data %>%
  arrange(desc(Significant_Count_05))
ggplot(data[data$Fisher.weight01<(5e-2/26),], aes(y=reorder(paste(Term,GO.ID), Significant_Count_B), x=reorder(genome, medianKs), color=ploidy, size=-log10(as.numeric(Fisher.weight01)))) + geom_point() + scale_color_manual(values=ploidycolors) + theme(legend.position = 'None') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(data[data$Fisher.weight01<(5e-2),], aes(y=reorder(paste(Term,GO.ID), Significant_Count_05), x=reorder(genome, medianKs), color=ploidy, size=-log10(as.numeric(Fisher.weight01)))) + geom_point() + scale_color_manual(values=ploidycolors) + theme(legend.position = 'None') + theme(axis.text.x = element_text(angle = 45, hjust = 1))



### plot retained in dup from this analysis
ggplot(asize, aes(x=medianKs, y=retainedInDup, group=ploidy, color=ploidy)) + 
  #  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(seed = 1), size=4)+ 
 scale_color_manual(values=ploidycolors) + #ylim(0,110000) + 
  # ggpubr::stat_compare_means(aes(group=ploidy, x=ploidy), label = 'p.signif', show.legend = F,ref.group = "Diploid", label.y=110000) + 
  #  geom_hline(yintercept=c(median(asize$diploidEquivalent[asize$ploidy=='Diploid'], na.rm=T)*c(1)), lty='dotted', color='darkgray')+
  #  annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T), label = "\u00D71", vjust = -0.5) + 
  # annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*2, label = "\u00D72", vjust = -0.5)+ 
  # annotate("text", x = 4.45, y = median(asize$haploidGeneCount[asize$ploidy=='Diploid'], na.rm=T)*3, label = "\u00D73", vjust = -0.5) + 
  xlab('Ks') + ylab('Genes Retained in Duplicate') + theme(legend.position='NULL') #+ geom_text(aes(label=V2))






### simple descriptor plot of shared gene duplicates across taxa
### is there a way to do it with 'average' zea/tripsacum??
### whatever just start with ALL zea/trip having duplicates
bbddz=bbdd[,asize$V2[asize$ploidy %in% c('Hexaploid', 'Tetraploid')]]
bbddz$tripsacineae=sapply(1:nrow(bbdd), function(x) all(bbdd[x,asize$V2[asize$ploidy=='Paleotetraploid']]))

bbddt=bbdd[,asize$V2[asize$ploidy %in% c('Hexaploid', 'Tetraploid') ]]
bbddt$tripsacum=sapply(1:nrow(bbdd), function(x) all(bbdd[x,asize$V2[asize$V2%in%c('tdacn1', 'tdacs1')]]))


dupgenecount=data.frame(table(rowSums(bbddz)))
colnames(dupgenecount)=c('numbergenomes', 'numbergenes')
dupgenecount$numbergenomes=as.numeric(as.character(dupgenecount$numbergenomes))
dupgenecount=rbind(dupgenecount, data.frame(numbergenomes=0:3, numbergenes=0))
ggplot(dupgenecount, aes(x=as.numeric(as.character(numbergenomes)), y=numbergenes)) + geom_bar(stat='identity') + xlab('Number Polyploidy Events with Gene Duplicated') + ylab('Number Genes')  + scale_x_continuous(breaks = 0:14, labels = 0:14) + geom_text(aes(label = numbergenes), size=3,vjust = -0.5, position = position_dodge(0.9))

## using just tripsacum
dupgenecountt=data.frame(table(rowSums(bbddt)))
colnames(dupgenecountt)=c('numbergenomes', 'numbergenes')
dupgenecountt$numbergenomes=as.numeric(as.character(dupgenecountt$numbergenomes))
dupgenecountt=rbind(dupgenecountt, data.frame(numbergenomes=0:3, numbergenes=0))
ggplot(dupgenecountt, aes(x=as.numeric(as.character(numbergenomes)), y=numbergenes)) + geom_bar(stat='identity') + xlab('Number Polyploidy Events with Gene Duplicated') + ylab('Number Genes')  + scale_x_continuous(breaks = 0:14, labels = 0:14) + geom_text(aes(label = numbergenes), size=3,vjust = -0.5, position = position_dodge(0.9))





## redo zea count for sc
bbddzs=bbdd[,asize$V2[asize$ploidy %in% c('Hexaploid', 'Tetraploid')]]
bbddzs=data.frame(lapply(bbddzs, function(x) !x))
bbddzs$tripsacineae=sapply(1:nrow(bbdd), function(x) all(!bbdd[x,asize$V2[asize$ploidy=='Paleotetraploid']]))
scgenecount=data.frame(table(rowSums(bbddzs)))
colnames(scgenecount)=c('numbergenomes', 'numbergenes')
scgenecount$numbergenomes=as.numeric(as.character(scgenecount$numbergenomes))
scgenecount=rbind(scgenecount, data.frame(numbergenomes=c(10:14), numbergenes=0))
ggplot(scgenecount, aes(x=as.numeric(as.character(numbergenomes)), y=numbergenes)) + geom_bar(stat='identity') + xlab('Number Polyploidy Events with Gene Single Copy') + ylab('Number Genes')  + scale_x_continuous(breaks = 0:14, labels = 0:14) + geom_text(aes(label = numbergenes), size=3,vjust = -0.5, position = position_dodge(0.9))



## oh cool, check these 156 genes always retained in duplicate!!!
dupinall=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(bbddz)>=14]), genomeid='dupinall', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])
dupin13=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(bbddz)>=13]), genomeid='dupin13', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])
dupin12=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(bbddz)>=12]), genomeid='dupin12', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])
dupin11=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(bbddz)>=11]), genomeid='dupin11', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])
dupin10=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(bbddz)>=10]), genomeid='dupin10', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])

dupintrip=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(bbddt)>=14]), genomeid='dupinalltrip')

nullset=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene), genomeid='nullset')



scin4=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(bbddzs)>=4]), genomeid='scin4')

##table(rowSums(bbdd[,asize$V2[asize$ploidy!='Diploid']]))

##########################################
## now be smarter and use syntenic genes as background
dupinalls=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(bbddz)>=14]), genomeid='dupinalls', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])
dupintrips=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(bbddt)>=14]), genomeid='dupinalltrips', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])

dupin13trip=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(bbddt)>=13]), genomeid='dupin13trips', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])


table(rowSums(bbdd[,asize$V2[asize$ploidy!='Diploid' & asize$medianKs<0.03]]))
dupinyoungs=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(bbdd[,asize$V2[asize$ploidy!='Diploid' & asize$medianKs<0.03]])>=3]), genomeid='dupinyoungs', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])
singlecopyinyoungs=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(!bbdd[,asize$V2[asize$ploidy!='Diploid' & asize$medianKs<0.03]])>=3]), genomeid='singlecopyinyoungs', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])

table(rowSums(bbdd[,asize$V2[asize$ploidy!='Diploid' & asize$medianKs>0.05 & asize$medianKs<0.1]]))
dupinolds=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(bbdd[,asize$V2[asize$ploidy!='Diploid' & asize$medianKs>0.05 & asize$medianKs<0.1]])>=4]), genomeid='dupinolds', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])
singlecopyinolds=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(!bbdd[,asize$V2[asize$ploidy!='Diploid' & asize$medianKs>0.05 & asize$medianKs<0.1]])>=3]), genomeid='singlecopyinolds', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])

table(rowSums(bbdd[,asize$V2[asize$ploidy!='Diploid' & asize$medianKs>0.05 & asize$medianKs<0.1]]&!bbdd[,asize$V2[asize$ploidy!='Diploid' & asize$medianKs<0.03]]))
dupinoldnotyoung=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(bbdd[,asize$V2[asize$ploidy!='Diploid' & asize$medianKs>0.05 & asize$medianKs<0.1]]&!bbdd[,asize$V2[asize$ploidy!='Diploid' & asize$medianKs<0.03]])>=4]), genomeid='dupinolds', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])

dupinyoungnotold=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(!bbdd[,asize$V2[asize$ploidy!='Diploid' & asize$medianKs>0.05 & asize$medianKs<0.1]]&bbdd[,asize$V2[asize$ploidy!='Diploid' & asize$medianKs<0.03]])>=2]), genomeid='dupinolds', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])



table(rowSums(bbdd[,asize$V2[asize$ploidy=='Paleotetraploid']]))
dupintrips=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(bbdd[,asize$V2[asize$ploidy=='Paleotetraploid' ]])>=12]), genomeid='dupintrips', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])
singlecopyintrips=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(!bbdd[,asize$V2[asize$ploidy=='Paleotetraploid']])>=12]), genomeid='singlecopyintrips', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])
dupinzeas=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(bbdd[,asize$V2[asize$ploidy=='Paleotetraploid' & !asize$V2 %in% c('tdacn1', 'tdacs1') ]])>=10]), genomeid='dupintrips', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])
singlecopyinzeas=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(!bbdd[,asize$V2[asize$ploidy=='Paleotetraploid' & !asize$V2 %in% c('tdacn1', 'tdacs1')]])>=10]), genomeid='singlecopyintrips', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])


## two oldest blagur, sscopa

table(rowSums(bbddz[,c('blagur', 'sscopa', 'tripsacineae')]))
dupintripsolds=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(bbddz[,c('blagur', 'sscopa', 'tripsacineae')])>=3]), genomeid='dupintripolds', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])
singlecopyintripsolds=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(!bbddz[,c('blagur', 'sscopa', 'tripsacineae')])>=3]), genomeid='singlecopyintripolds', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])

dupintripsoldsscinyoungs=run_GO_analysis(gsub('.1.v3.1', '.1', bbdd$gene[rowSums(!bbddz[,c('blagur', 'sscopa', 'tripsacineae')] & bbddz[,c('vcuspi', 'hcompr', 'atenui')])>=3]), genomeid='dupintripolds', background = syntbackground, GODB=GODB[names(GODB)%in%syntbackground])



## read in all dnds!!!!1

dnds=fread('~/Downloads/ksp_to_pasp.txt.gz', header=T)

### regular vs staydup
ggplot(dnds[dnds$genome%in%asize$V2[!asize$ploidy%in%c('Diploid','Paleotetraploid')],], aes(x=V19, fill=paste0(gene, '.1.v3.1')%in%bbdd$gene[rowSums(bbddz)>=14])) + geom_histogram(binwidth=0.01) + xlim(0,1.2) + facet_grid(genome~paste0(gene, '.1.v3.1')%in%bbdd$gene[rowSums(bbddz)>=14], scales='free_y') + theme(legend.position='NULL')

## within staydup, difference between gene copies
dc=dnds %>% group_by(genome, gene) %>% summarize(dndsdiff=max(V19)-min(V19), mindnds=min(V19), copynumber=n())
dc$staydup=paste0(dc$gene,'.1.v3.1') %in% bbdd$gene[rowSums(bbddz)>=14]
  
ggplot(dc, aes(x=dndsdiff, fill=staydup)) + geom_histogram(binwidth=0.01)  + xlim(0,1) + ylim(0,6e4)
  
  
  
tail(sort(table(dc[ dc$copynumber>1 & dc$dndsdiff>0.3 & dc$mindnds<0.5 & dc$genome%in%asize$V2[asize$ploidy%in%c('Tetraploid','Hexaploid')],]$gene))  )
 



annuals=c('zmB735', 'sbicol', 'zTIL11', 'zTIL01', 'zTIL25', 'zTIL18', 'amhuet', 'telega')

ggplot(dnds, aes(x=V19, fill=genome%in%annuals)) + geom_histogram(binwidth=0.01) + xlim(0,1.2) 
median(dnds$V19[dnds$genome%in%annuals], na.rm=T)
median(dnds$V19[!dnds$genome%in%annuals], na.rm=T)
median(dnds$V19[dnds$genome%in%asize$V2[asize$ploidy=='Diploid']], na.rm=T)
median(dnds$V19[dnds$genome%in%asize$V2[!asize$ploidy=='Diploid']], na.rm=T)
