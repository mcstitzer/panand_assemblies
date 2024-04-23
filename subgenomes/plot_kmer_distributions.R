library(stringr)
library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())


a=fread('udigit.5.tsv')
rownames(a)=a$region_id
k=kmeans(data.frame(a[,-1]), centers=3)
aa=cbind(a, cluster=k$cluster)

pca=prcomp(a[,-1])

df=data.frame(cbind(pca$x[,1:2], aa$cluster))

df$PC1=as.numeric(df$PC1)/(pca$sdev[1]*sqrt(nrow(a)))
df$PC2=as.numeric(df$PC2)/(pca$sdev[2]*sqrt(nrow(a)))
df$pos=a$region_id



pdf('udigit.5.pdf',10,10)
ggplot(df, aes(x=PC1, y=PC2, color=as.factor(V3))) + geom_point()
ggplot(df, aes(x=PC1, y=PC2, color=as.factor(V3))) + geom_text(aes(label=pos))
ggplot(df, aes(x=PC1, y=PC2, color=str_split_fixed(pos, '_',3)[,3]=='0')) + geom_text(aes(label=pos))
dev.off()


