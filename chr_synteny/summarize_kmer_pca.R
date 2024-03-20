library(stringr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(ggfortify)

a=Sys.glob('*/*.kmer.mat')

b=data.frame(genome=substr(a, 1,6), kmer=gsub('k','', unlist(str_extract_all(a, "k(\\d+)"))), q=gsub('q','',str_split_fixed(a, '_',3)[,2]))

pdf('biplot_all_subphaser.pdf',10,10)

for(i in 1:length(a)){
mat=read.table(a[i], header=T)
mat2=t(mat[,-1])
rownames(mat2)=colnames(mat)[-1]
colnames(mat2)=mat$kmer
if(ncol(mat2)>1){
pca=prcomp(mat2)
print(autoplot(pca, data = mat2, label = TRUE, label.size = 3, loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3) + ggtitle(a[i]))
pcprops=pca$sdev^2/sum(pca$sdev^2)
b$PC1prop[i]=pcprops[1]
b$PC2prop[i]=pcprops[2]
}
}
dev.off()
