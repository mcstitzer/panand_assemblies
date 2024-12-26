library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(stringr)
library(tidyr)
library(dplyr)


#a=fread('Ac-Pasquet1232-DRAFT-PanAnd-1.0_coverage.txt')
a=fread('Ac-Pasquet1232-DRAFT-PanAnd-1.0_syntcoverage.txt')
a$id=sub("_[0-9]+$", "", a$V4)
#g=read.table('Ac-Pasquet1232-DRAFT-PanAnd-1.0_genes_only.gff')
g=read.table('../../syntenic/Ac-Pasquet1232-DRAFT-PanAnd-1.0.syntenicAnchors.gff3')
g$id=gsub('ID=', '', g$V9)
a$strand=g$V7[match(a$id, g$id)]
a$window=sub(".*_([0-9]+)$", "\\1", a$V4)
a$window=as.numeric(a$window)
#a$window[duplicated(a$window)]=a$window[duplicated(a$window)]+200
a=a%>% group_by(id, window) %>% mutate(windowadj=ifelse(row_number()==2, window+200, window)) %>% ungroup()
b=a%>% group_by(id)%>% filter(dplyr::n()==400)%>%filter(strand=='+') 
brev=a%>% group_by(id)%>% filter(dplyr::n()==400)%>%filter(strand=='-') %>% group_by(id) %>% mutate(windowadj=max(windowadj)+1-windowadj) %>% ungroup() 

bb=rbind(b, brev)

bbmed=bb %>% group_by(windowadj) %>% dplyr::summarize(med=median(V8), mean=mean(V8))

pdf('~/transfer/achine_metaplot_syntenic.pdf',14,8)
ggplot(bbmed, aes(x=windowadj, y=med)) + geom_line()
ggplot(bbmed, aes(x=windowadj, y=mean)) + geom_line()
ggplot(bb, aes(x=windowadj, y=V8, group=windowadj)) + geom_violin()

dev.off()

## next need to filter the orthogroup ids by light orthogroup sizes
## and run with syntenic anchors

## also run with CDS start/end



## chatgpt loop
combined_data <- data.frame()
for(i in seq_len(nrow(all))){
  genome_id <- all$V2[i]
  base_name <- all$V1[i]

  # Read coverage and gene data
  a <- fread(paste0(base_name, '_syntcoverage.txt'))
  a$id <- sub("_[0-9]+$", "", a$V4)
  g <- read.table(paste0('../../syntenic/', base_name, '.syntenicAnchors.bed'))
  g$id <- g$V4

  # Process coverage data
  a$strand <- g$V6[match(a$id, g$id)]
  a$window <- sub(".*_([0-9]+)$", "\\1", a$V4)
  a$window <- as.numeric(a$window)

  a <- a %>% group_by(id, window) %>% mutate(windowadj = ifelse(row_number() == 2, window + 200, window)) %>% ungroup()

  b <- a %>% group_by(id) %>% filter(dplyr::n() == 400) %>% filter(strand == '+')
  brev <- a %>% group_by(id) %>% filter(dplyr::n() == 400) %>% filter(strand == '-') %>% group_by(id) %>% mutate(windowadj = max(windowadj) + 1 - windowadj) %>% ungroup()

  bb <- rbind(b, brev)

  bbmed <- bb %>% group_by(windowadj) %>% dplyr::summarize(med = median(V8), mean = mean(V8))
  # Add genome_id column
  bbmed <- bbmed %>% mutate(genome_id = genome_id)

  # Combine results
  combined_data <- rbind(combined_data, bbmed)

#   # Write summarized data
#   write.table(bbmed, file = paste0(base_name, '_summary.txt'), sep = "\t", quote = FALSE, row.names = FALSE)
#   cat("Written:", paste0(base_name, '_summary.txt'), "\n")
}


ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')
asize=fread('../..//panand_assembly_sizes.txt', header=T, quote="", fill=T)
## factor so correct order
asize$ploidy=factor(asize$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))

combined_data$ploidy=asize$ploidy[match(combined_data$genome_id, asize$V2)]
pdf('~/transfer/combo_metaplot_syntenic.pdf',14,8)
ggplot(combined_data, aes(x=windowadj, y=med, color=ploidy, group=genome_id)) + geom_line() + scale_color_manual(values=ploidycolors)
ggplot(combined_data, aes(x=windowadj, y=mean, color=ploidy, group=genome_id)) + geom_line() + scale_color_manual(values=ploidycolors)

dev.off()



## nonsynt - all helixer
## chatgpt loop
combined_data_helixer <- data.frame()
for(i in seq_len(nrow(all))){
  genome_id <- all$V2[i]
  base_name <- all$V1[i]

  # Read coverage and gene data
  a <- fread(paste0(base_name, '_coverage.txt'))
  a$id <- sub("_[0-9]+$", "", a$V4)
  g <- read.table(paste0(base_name, '_filtered_genes.gff'))
  g$id <- g$V4

  # Process coverage data
  a$strand <- g$V6[match(a$id, g$id)]
  a$window <- sub(".*_([0-9]+)$", "\\1", a$V4)
  a$window <- as.numeric(a$window)

  a <- a %>% group_by(id, window) %>% mutate(windowadj = ifelse(row_number() == 2, window + 200, window)) %>% ungroup()

  b <- a %>% group_by(id) %>% filter(dplyr::n() == 400) %>% filter(strand == '+')
  brev <- a %>% group_by(id) %>% filter(dplyr::n() == 400) %>% filter(strand == '-') %>% group_by(id) %>% mutate(windowadj = max(windowadj) + 1 - windowadj) %>% ungroup()

  bb <- rbind(b, brev)

  bbmed <- bb %>% group_by(windowadj) %>% dplyr::summarize(med = median(V8), mean = mean(V8))
  # Add genome_id column
  bbmed <- bbmed %>% mutate(genome_id = genome_id)

  # Combine results
  combined_data_helixer <- rbind(combined_data_helixer, bbmed)

#   # Write summarized data
#   write.table(bbmed, file = paste0(base_name, '_summary.txt'), sep = "\t", quote = FALSE, row.names = FALSE)
#   cat("Written:", paste0(base_name, '_summary.txt'), "\n")
}


ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')
asize=fread('../..//panand_assembly_sizes.txt', header=T, quote="", fill=T)
## factor so correct order
asize$ploidy=factor(asize$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))

combined_data_helixer$ploidy=asize$ploidy[match(combined_data_helixer$genome_id, asize$V2)]

write.table(combined_data_helixer, 'combined_tedist_helixer_genes.txt', row.names=F, col.names=T, sep='\t', quote=F)

pdf('~/transfer/combo_metaplot_helixer.pdf',14,8)
ggplot(combined_data_helixer, aes(x=windowadj, y=med, color=ploidy, group=genome_id)) + geom_line() + scale_color_manual(values=ploidycolors)
ggplot(combined_data_helixer, aes(x=windowadj, y=mean, color=ploidy, group=genome_id)) + geom_line() + scale_color_manual(values=ploidycolors)


### upstream
uph=ggplot(combined_data_helixer[combined_data_helixer$windowadj%in%1:200,], aes(group=genome_id, x=(windowadj-200)*100, y=mean, color=ploidy))+ geom_vline(xintercept=c(0:20*10*-100), color='whitesmoke', lty='dotted')+ geom_hline(yintercept=c(0.2,0.4,0.6,0.8), color='seashell2', lty='longdash') + geom_line() + scale_color_manual(values=ploidycolors)  + xlab('Base pairs away from TSS') + ylab('Mean TE base pairs in 100bp window') + theme(legend.position='NULL')
### downstream
downh=ggplot(combined_data_helixer[combined_data_helixer$windowadj%in%201:400,], aes(group=genome_id, x=(windowadj-200)*100, y=mean, color=ploidy))+ geom_vline(xintercept=c(0:20*10*100), color='whitesmoke', lty='dotted')+ geom_hline(yintercept=c(0.2,0.4,0.6,0.8), color='seashell2', lty='longdash') + geom_line() + scale_color_manual(values=ploidycolors)  + xlab('Base pairs away from TTS') + ylab('Mean TE base pairs in 100bp window') + theme(legend.position='NULL')

plot_grid(uph, downh, align='tb', axis='hv', ncol=2, labels=c('B', 'C'))

dev.off()


#### filter to orthogroups
## nonsynt - orthogroup helixer
## chatgpt loop


og=read.table('../../../genespace/orthofinder/Results_Nov14/Orthogroups/Orthogroups.tsv', header=T, sep='\t')
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

a=read.table('../../../genespace/orthofinder/Results_Nov14/Orthogroups/Orthogroups.GeneCount.tsv', header=T)

b=a[rowSums(a[,2:39]==0)==0 & rowSums(a[,2:39]>12)==0,]
## simpler
b=a[a$Total>40 & a$Total<200,]

filtog=og_melted[og_melted$Orthogroup%in%b$Orthogroup,]




combined_data_oghelixer <- data.frame()
for(i in seq_len(nrow(all))){
  genome_id <- all$V2[i]
  base_name <- all$V1[i]

  # Read coverage and gene data
  a <- fread(paste0(base_name, '_coverage.txt'))
  a$id <- sub("_[0-9]+$", "", a$V4)
  a=a[paste0( genome_id, '_', gsub('ID=', '', a$id)) %in% filtog$genes,]
  
  g <- read.table(paste0(base_name, '_filtered_genes.gff'))
  g$id <- g$V4

  # Process coverage data
  a$strand <- g$V6[match(a$id, g$id)]
  a$window <- sub(".*_([0-9]+)$", "\\1", a$V4)
  a$window <- as.numeric(a$window)

  a <- a %>% group_by(id, window) %>% mutate(windowadj = ifelse(row_number() == 2, window + 200, window)) %>% ungroup()

  b <- a %>% group_by(id) %>% filter(dplyr::n() == 400) %>% filter(strand == '+')
  brev <- a %>% group_by(id) %>% filter(dplyr::n() == 400) %>% filter(strand == '-') %>% group_by(id) %>% mutate(windowadj = max(windowadj) + 1 - windowadj) %>% ungroup()

  bb <- rbind(b, brev)

  bbmed <- bb %>% group_by(windowadj) %>% dplyr::summarize(med = median(V8), mean = mean(V8))
  # Add genome_id column
  bbmed <- bbmed %>% mutate(genome_id = genome_id)

  # Combine results
  combined_data_oghelixer <- rbind(combined_data_oghelixer, bbmed)

#   # Write summarized data
#   write.table(bbmed, file = paste0(base_name, '_summary.txt'), sep = "\t", quote = FALSE, row.names = FALSE)
#   cat("Written:", paste0(base_name, '_summary.txt'), "\n")
}


ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')
asize=fread('../..//panand_assembly_sizes.txt', header=T, quote="", fill=T)
## factor so correct order
asize$ploidy=factor(asize$ploidy, levels=c('Diploid', 'Tetraploid', 'Paleotetraploid', 'Hexaploid'))

combined_data_oghelixer$ploidy=asize$ploidy[match(combined_data_oghelixer$genome_id, asize$V2)]

write.table(combined_data_helixer, 'combined_tedist_helixer_genes.txt', row.names=F, col.names=T, sep='\t', quote=F)

pdf('~/transfer/combo_metaplot_orthogroup_helixer.pdf',14,8)
ggplot(combined_data_oghelixer, aes(x=windowadj, y=med, color=ploidy, group=genome_id)) + geom_line() + scale_color_manual(values=ploidycolors)
ggplot(combined_data_oghelixer, aes(x=windowadj, y=mean, color=ploidy, group=genome_id)) + geom_line() + scale_color_manual(values=ploidycolors)


### upstream
uph=ggplot(combined_data_helixer[combined_data_helixer$windowadj%in%1:200,], aes(group=genome_id, x=(windowadj-200)*100, y=mean, color=ploidy))+ geom_vline(xintercept=c(0:20*10*-100), color='whitesmoke', lty='dotted')+ geom_hline(yintercept=c(0.2,0.4,0.6,0.8), color='seashell2', lty='longdash') + geom_line() + scale_color_manual(values=ploidycolors)  + xlab('Base pairs away from TSS') + ylab('Mean TE base pairs in 100bp window') + theme(legend.position='NULL')
### downstream
downh=ggplot(combined_data_helixer[combined_data_helixer$windowadj%in%201:400,], aes(group=genome_id, x=(windowadj-200)*100, y=mean, color=ploidy))+ geom_vline(xintercept=c(0:20*10*100), color='whitesmoke', lty='dotted')+ geom_hline(yintercept=c(0.2,0.4,0.6,0.8), color='seashell2', lty='longdash') + geom_line() + scale_color_manual(values=ploidycolors)  + xlab('Base pairs away from TTS') + ylab('Mean TE base pairs in 100bp window') + theme(legend.position='NULL')

plot_grid(uph, downh, align='tb', axis='hv', ncol=2, labels=c('B', 'C'))

dev.off()


### measure distance to 0.5 coverage

combined_data_helixer %>%
+     group_by(genome_id, ploidy) %>%
+     filter(mean < 0.5) %>%          # Filter rows where mean < 0.5
+     slice_min(windowadj) %>%       # Select the first occurrence based on windowadj
+     select(genome_id, windowadj)%>% group_by(ploidy) %>% summarize(median(windowadj))  # Keep only relevant columns


## supplemental figure: repeat proportion doesn't explain rearrangments
pdf('../figures/supp_fig_repeat_prop_vs_rearrangements.pdf',8,8)
ggplot(asize, aes(x=repeatProp, y=scaledTransloc, color=ploidy))+geom_point(size=3)+xlab('Proportion Repeats')+ylab('Rearrangements per Diploid Equivalent')+scale_color_manual(values=ploidycolors)+theme(legend.position='NULL')
dev.off()


