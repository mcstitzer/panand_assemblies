
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())



a=read.table('Orthogroups.GeneCount.tsv', header=T)


b=a[rowSums(a[,2:39]==0)==0 & rowSums(a[,2:39]>12)==0,]
## simpler
b=a[a$Total>40 & a$Total<200,]

b$oneTrip=rowSums(b[,c(2:22,26:29)])



ggplot(b, aes(x=oneTrip)) + geom_histogram(binwidth=1)
ggplot(b, aes(x=Total)) + geom_histogram(binwidth=1)


annual=c('sbicol', 'zmhuet', 'zTIL01', 'zTIL11', 'zTIL18', 'zTIL25','telega', 'zmB735')
perennial=setdiff(colnames(b)[-c(1,40,41)], annual)


results=apply(b[,-c(1,40,41)],1,function(row){
pheno=as.numeric(row[annual])
other=as.numeric(row[perennial])
t.test(pheno, other)$p.value
})
b$pvalue=results

meds=apply(b[,-c(1,40,41,42)],2, median)
results=apply(b[,-c(1,40:42)], 1, function(row){
deviations=row-meds
phen=deviations[annual]
other=deviations[perennial]
wilcox.test(phen, other)$p.value
})
warnings()
b$med_pvalue=results


head(b[order(b$med_pvalue),])



### check the amino acid sequence of some of these orthogroups!!!
## i can at least support a lot of kinesins!!!

## try annual with only one zea anjual rep
annual=c('sbicol', 'zmhuet', 'telega')
results=apply(b[,-c(1,40:ncol(b))],1,function(row){
pheno=as.numeric(row[annual])
other=as.numeric(row[perennial])
t.test(pheno, other)$p.value
})
b$stictannual_pvalue=results

results=apply(b[,-c(1,40:ncol(b))], 1, function(row){
deviations=row-meds
phen=deviations[annual]
other=deviations[perennial]
wilcox.test(phen, other)$p.value
})
warnings()
b$strictannual_med_pvalue=results


head(b[order(b$strictannual_med_pvalue),])



### broader annuals
annualish=c('irugos', 'zluxur', 'znicar','sbicol', 'zmhuet', 'zTIL01', 'zTIL11', 'zTIL18', 'zTIL25','telega', 'zmB735')
perennialish=setdiff(colnames(b)[-c(1,40,41)], annualish)


results=apply(b[,-c(1,40,41)],1,function(row){
pheno=as.numeric(row[annualish])
other=as.numeric(row[perennialish])
t.test(pheno, other)$p.value
})
b$annualish_pvalue=results

meds=apply(b[,-c(1,40,41,42)],2, median)
results=apply(b[,-c(1,40:42)], 1, function(row){
deviations=row-meds
phen=deviations[annualish]
other=deviations[perennialish]
wilcox.test(phen, other)$p.value
})
warnings()
b$annualish_med_pvalue=results

head(b[order(b$annualish_med_pvalue),])



## try for tannins
tannin=c( "agerar", "atenui", "avirgi", "blagur", 
"ccitra", "crefra", "cserru", "etrips", "hcompr", "hconto", 
"ppanic",  "smicro", "snutan", 
"sscopa", "ttrian", 
"udigit", "vcuspi")
notannin=setdiff(colnames(b)[-c(1,40,41,42,43)], tannin)

results=apply(b[,-c(1,40,41)],1,function(row){
pheno=as.numeric(row[tannin])
other=as.numeric(row[notannin])
t.test(pheno, other)$p.value
})
b$tannin_pvalue=results

results=apply(b[,-c(1,40:42)], 1, function(row){
deviations=row-meds
phen=deviations[tannin]
other=deviations[notannin]
wilcox.test(phen, other)$p.value
})
warnings()
b$tannin_med_pvalue=results

head(b[order(b$tannin_med_pvalue),])


## got to account for all the tripsacinae sampling - 


### self compatibility
sc=c( "avirgi", "blagur",  "irugos", 
"ppanic", "pvagin", "rtuber", "sbicol",  
"sscopa", "tdacn1", "tdacn2", "tdacs1", "tdacs2",  
"udigit", "zTIL01", "zTIL11", "zTIL18", "zTIL25", "zdgigi", 
"zdmomo", "zluxur", "zmB735", "zmhuet", "znicar")
si=c("agerar", 'crefra', 'ccitra', 'snutan', 'telega')

results=apply(b[,-c(1,40,41)],1,function(row){
pheno=as.numeric(row[sc])
other=as.numeric(row[si])
t.test(pheno, other)$p.value
})
b$self_pvalue=results

results=apply(b[,-c(1,40:42)], 1, function(row){
deviations=row-meds
phen=deviations[sc]
other=deviations[si]
wilcox.test(phen, other)$p.value
})
warnings()
b$self_med_pvalue=results

head(b[order(b$self_med_pvalue),])

