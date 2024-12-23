library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())


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
pdf('~/transfer/combo_metaplot_helixer.pdf',14,8)
ggplot(combined_data_helixer, aes(x=windowadj, y=med, color=ploidy, group=genome_id)) + geom_line() + scale_color_manual(values=ploidycolors)
ggplot(combined_data_helixer, aes(x=windowadj, y=mean, color=ploidy, group=genome_id)) + geom_line() + scale_color_manual(values=ploidycolors)

dev.off()


