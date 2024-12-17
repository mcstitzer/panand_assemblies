library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

og=read.table('~/Downloads/genespace_dotplots/Orthogroups.GeneCount.tsv', header=T)

## duhhhh scale by haploid!!!
asize=read.table('../general_summaries/panand_assembly_sizes.txt', header=T)

dog=og
dog[,c(asize$V2[asize$haploid], 'tdacn2', 'tdacs2')]=dog[,asize$V2[asize$haploid]]*2

bog=dog[rowSums(dog[,2:39]==0)==0 & rowSums(dog[,2:39]>12)==0,]
bog$oneTrip=rowSums(bog[,c(2:22,26:29)]) ## only one tripsacinae representative, tdacs1


## set up with median tripsacinae value for og

tripsacinae=c("tdacn2", "tdacs2", "tdacn1", "tdacs1", "zdgigi", "zdmomo", 
              "zluxur", "zmhuet", "zTIL18", "zTIL25", "zTIL01", "zTIL11", "znicar", 
              "zmB735")

bog$medianTZ=sapply(1:nrow(bog), function(x) median(as.numeric(bog[x, tripsacinae])))




annual=c('sbicol', 'zmhuet', 'zTIL01', 'zTIL11', 'zTIL18', 'zTIL25','telega', 'zmB735')
perennial=setdiff(colnames(bog)[-c(1,40,41:ncol(bog))], annual)


results=apply(bog[,-c(1,40,41)],1,function(row){
  pheno=as.numeric(row[annual])
  other=as.numeric(row[perennial])
  t.test(pheno, other)$p.value
})
bog$pvalue=results

meds=apply(bog[,-c(1,40,41,42)],2, median)
results=apply(bog[,-c(1,40:42)], 1, function(row){
  deviations=row-meds
  phen=deviations[annual]
  other=deviations[perennial]
  wilcox.test(phen, other)$p.value
})
#warnings()
bog$med_pvalue=results


head(bog[order(bog$med_pvalue),c('Orthogroup', annual, 'med_pvalue', perennial)])
dput(head(bog[order(bog$med_pvalue),c('Orthogroup', annual, 'med_pvalue', perennial)]$Orthogroup, 100))

## while waiting for blast2go,...
## OG0016237 snf1 kinase which is response to engergy starvation in arabidopsis??
## OG0009289 methyl cpg binding protein 4 - DNA glycosylase that excises deaminated cytosines
## OG0013534 metacapsase - cystine protease involved in programmed cell death during development - look at this grant "fantasy" of transforming annuals into perennials using metacapsases?? https://kaw.wallenberg.org/en/research/mapping-hidden-world-ancient-enzymes
## OG0002582 ORR9 type a response regulator cytokinin regulator, light response?
## OG0010854 RALF - rapid alkanization factor quickly changes ph and affects development
## OG0012253 SH2 domain containing protein - signal transduction from binding to phosphorylated tyrosines
## OG0015142 protein argonaute 5-like or Mel1-like - meiosis and phasirnas !! (lots of variation in protein length)
## OG0009333 nuclear transcription factor y subunit A 10-like - center of many de3velopmental stress-responsive processes, incl. drought
## OG0008347 aaa-type atpase family protein- represses flowering in perennial Arabis alpina
##### i can make a story here i think!!!

meds=apply(bog[,-c(1,40,41,42)],2, median)
results=apply(bog[,-c(1,40:42)], 1, function(row){
  deviations=row/meds
  phen=deviations[annual]
  other=deviations[perennial]
  wilcox.test(phen, other, mu=1)$p.value
})
#warnings()
bog$div_pvalue=results

head(bog[order(bog$div_pvalue),c('Orthogroup', annual, 'div_pvalue', perennial)])
head(bog[order(bog$div_pvalue),c(annual,perennial)]-meds[c(annual,perennial)])

dput(head(bog[order(bog$div_pvalue),c('Orthogroup', annual, 'div_pvalue', perennial)]$Orthogroup, 100))



### broader annuals
annualish=c('irugos', 'zluxur', 'znicar','sbicol', 'zmhuet', 'zTIL01', 'zTIL11', 'zTIL18', 'zTIL25','telega', 'zmB735')
perennialish=setdiff(colnames(bog)[-c(1,40:ncol(bog))], annualish)


results=apply(bog[,-c(1,40,41:ncol(bog))],1,function(row){
  pheno=as.numeric(row[annualish])
  other=as.numeric(row[perennialish])
  t.test(pheno, other)$p.value
})
bog$annualish_pvalue=results

meds=apply(bog[,-c(1,40,41,42:ncol(bog))],2, median)
results=apply(bog[,-c(1,40:ncol(bog))], 1, function(row){
  deviations=row-meds
  phen=deviations[annualish]
  other=deviations[perennialish]
  wilcox.test(phen, other)$p.value
})
bog$annualish_med_pvalue=results

head(bog[order(bog$annualish_med_pvalue),c('Orthogroup', annualish, 'annualish_med_pvalue', perennialish)])

dput(head(bog[order(bog$annualish_med_pvalue),c('Orthogroup', annual, 'annualish_med_pvalue', perennial)]$Orthogroup, 100))


## try for tannins
tannin=c( "agerar", "atenui", "avirgi", "blagur", 
          "ccitra", "crefra", "cserru", "etrips", "hcompr", "hconto", 
          "ppanic",  "smicro", "snutan", 
          "sscopa", "ttrian", 
          "udigit", "vcuspi")
notannin=setdiff(colnames(bog)[-c(1,40,41,42,43:ncol(bog))], tannin)

results=apply(bog[,-c(1,40,41:ncol(bog))],1,function(row){
  pheno=as.numeric(row[tannin])
  other=as.numeric(row[notannin])
  t.test(pheno, other)$p.value
})
bog$tannin_pvalue=results

meds=apply(bog[,-c(1,40:ncol(bog))],2, median)

results=apply(bog[,-c(1,40:ncol(bog))], 1, function(row){
  deviations=row-meds
  phen=deviations[tannin]
  other=deviations[notannin]
  wilcox.test(phen, other)$p.value
})
bog$tannin_med_pvalue=results

head(bog[order(bog$tannin_med_pvalue),c('Orthogroup', tannin, 'tannin_med_pvalue', notannin)])
### OG0006062 flowering-promoting factor1 like - small peptide that promotes flowering but who knows what it is
### OG0004244 apo protein 4 mitochondrial (don't know what this is - google confuses alzheimers)
### OG0008122 serine carboxypeptidase-like 34, diversification of metabolites (appears to biosynthesize tannins???!?!!!)
### OG0009316 10 kDa chaperonin 1 chloroplastic
### OG0012258 E3 ubiquitin ligase RHY1a
### OG0010087 ethylene responsive transcription factor

## got to account for all the tripsacinae sampling - 
dput(head(bog[order(bog$tannin_med_pvalue),c('Orthogroup', annual, 'tannin_med_pvalue', perennial)]$Orthogroup, 100))


### self compatibility
sc=c( "avirgi", "blagur",  "irugos", 
      "ppanic", "pvagin", "rtuber", "sbicol",  
      "sscopa", "tdacn1", "tdacn2", "tdacs1", "tdacs2",  
      "udigit", "zTIL01", "zTIL11", "zTIL18", "zTIL25", "zdgigi", 
      "zdmomo", "zluxur", "zmB735", "zmhuet", "znicar")
si=c("agerar", 'crefra', 'ccitra', 'snutan', 'telega')

results=apply(bog[,-c(1,40,41:ncol(bog))],1,function(row){
  pheno=as.numeric(row[sc])
  other=as.numeric(row[si])
  t.test(pheno, other)$p.value
})
bog$self_pvalue=results

meds=apply(bog[,-c(1,40:ncol(bog))],2, median)

results=apply(bog[,-c(1,40:ncol(bog))], 1, function(row){
  deviations=row-meds
  phen=deviations[sc]
  other=deviations[si]
  wilcox.test(phen, other)$p.value
})
bog$self_med_pvalue=results

head(bog[order(bog$self_med_pvalue),c('Orthogroup', sc, 'self_med_pvalue', si)])
##
## OG0021131

dput(head(bog[order(bog$self_med_pvalue),c('Orthogroup', annual, 'self_med_pvalue', perennial)]$Orthogroup, 100))





####### polyploids!!!

diploids=c("cserru", "irugos", "sbicol", "ppanic", "ttrian", "crefra", 
           "avirgi", "smicro", "rtuber", 'telega')
tetraploids=c("snutan", "hconto", "ccitra", "achine", "sscopa", "etrips", 
              "vcuspi", 'atenui', 'rrottb')
hexaploids=c("udigit", "agerar", "hcompr", "blagur")
paleotetraploids=asize$V2[!asize$V2 %in% c(diploids, tetraploids, hexaploids)]

results=apply(bog[,-c(1,40,41:ncol(bog))],1,function(row){
  pheno=as.numeric(row[diploids])
  other=as.numeric(row[c(tetraploids, hexaploids, paleotetraploids)])
  t.test(pheno, other)$p.value
})
bog$polyploid_pvalue=results

meds=apply(bog[,-c(1,40:ncol(bog))],2, median)

results=apply(bog[,-c(1,40:ncol(bog))], 1, function(row){
  deviations=row-meds
  phen=deviations[diploids]
  other=deviations[c(tetraploids, hexaploids, paleotetraploids)]
  wilcox.test(phen, other)$p.value
})
bog$polyploid_med_pvalue=results


head(bog[order(bog$polyploid_med_pvalue),c('Orthogroup', diploids, 'polyploid_med_pvalue', c(tetraploids, hexaploids, paleotetraploids))])
### omfg these are so interpretable!?!>!>!>!!??!?!?
# OG0001629 bc10 glycosyltransferase, for brittle culm 10 of rice, so cell wall, affects cellulose and mechanical properties
# OG0001506 eukaryotic initiation factor 4a, RNA helicase needed for cell size homeostasis in arabidopsis, ovule size - variants between A. virgi and A. gerardi are at DEAD box helicase motif?!?!?!?
# OG0001702 bHLH130 transcription factor, lignin biosynthesis under nitrogen stress in apple
# OG0001695 snrk7 salt and drought stress
# OG0001848 H/ACA ribonucleoprotein complex subunit 4-like, pseudouridylation of ribosomal dna, telomeric maintenance, giant lysine homopolymer stretches in protein
# OG0001563 60S ribosomal protein L5-1
# OG0001403 serine/threonine phosphatase
# OG0001653 ELMO domain engulfment and cell motility, pollen aperturues???
# OG0001589 E3 ubiquitin-protein ligase SINAT5 - ubiquitin degradation of NAC1 to modulate auxin response, root morphogenesis
# OG0001931 calcium dependent protein kinase


## top 100 orthogroups -- send to go enrichment
dput(head(bog[order(bog$polyploid_med_pvalue),c('Orthogroup', diploids, 'polyploid_med_pvalue', c(tetraploids, hexaploids, paleotetraploids))]$Orthogroup,100))

## so for root hairs:
bog[bog$Orthogroup%in%c("OG0001626", "OG0001829", "OG0002131"),c('Orthogroup', diploids, 'polyploid_med_pvalue', c(tetraploids, hexaploids, paleotetraploids))]


##just tet hex because what if paleotets have doen it all
results=apply(bog[,-c(1,40,41:ncol(bog))],1,function(row){
  pheno=as.numeric(row[diploids])
  other=as.numeric(row[c(tetraploids, hexaploids)])
  t.test(pheno, other)$p.value
})
bog$tethexaploid_pvalue=results

meds=apply(bog[,-c(1,40:ncol(bog))],2, median)

results=apply(bog[,-c(1,40:ncol(bog))], 1, function(row){
  deviations=row-meds
  phen=deviations[diploids]
  other=deviations[c(tetraploids, hexaploids)]
  wilcox.test(phen, other)$p.value
})
bog$tethexaploid_med_pvalue=results


head(bog[order(bog$tethexaploid_med_pvalue),c('Orthogroup', diploids, 'tethexaploid_med_pvalue', c(tetraploids, hexaploids, paleotetraploids))])
## get top 100
dput(head(bog[order(bog$tethexaploid_med_pvalue),c('Orthogroup', diploids, 'tethexaploid_med_pvalue', c(tetraploids, hexaploids, paleotetraploids))]$Orthogroup,100))


### top hit is the same as with paleos!! third is also in top10

#OG0001629
#OG0002099
#OG0001589
#OG0001754 glutamate receptor 3,4 like receptor, lateral root initiation in arabidopsis
#OG0001911 nicotianamine synthase 3, long distance iron transport
#OG0002143 mitogen activated protein kinase, innate immunity


### wait a sec should i be dividing??
##just tet hex because what if paleotets have doen it all

meds=apply(bog[,-c(1,40:ncol(bog))],2, median)

results=apply(bog[,-c(1,40:ncol(bog))], 1, function(row){
  deviations=row/meds
  phen=deviations[diploids]
  other=deviations[c(tetraploids)]
  wilcox.test(phen, other, mu=1, alternative='less')$p.value
})
bog$tethexaploid_div_pvalue=results


head(bog[order(bog$tethexaploid_div_pvalue),c('Orthogroup', diploids, 'tethexaploid_div_pvalue', c(tetraploids, hexaploids, paleotetraploids))])


