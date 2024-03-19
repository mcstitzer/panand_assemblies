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
if(nrow(mat)>1){
pca=prcomp(mat[,-1])
print(autoplot(pca, data = mat, label = TRUE, label.size = 3, loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3) + ggtitle(a[i]))
pcprops=pca$sdev^2/sum(pca$sdev^2)
b$PC1prop[i]=pcprops[1]
b$PC2prop[i]=pcprops[2]
}
}
dev.off()
