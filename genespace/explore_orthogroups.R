
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