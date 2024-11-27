library(seqinr)
library(ape)
library(rphast)
library(reshape2)
library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(stringr)

# Function to check if a sequence ends with "CAT"
check_sequence_end <- function(msa, pattern = "Pavag", suffix = "CAT") {
  # Find sequences that start with the pattern
  seqs <- msa$seq[grep(paste0("^", pattern), msa$names)]
  # Check if any of these sequences end with the suffix
  any(sapply(seqs, function(seq) endsWith(seq, suffix)))
}

rotate_clades <- function(tree) {
  # Helper function to check if a node is a leaf
  is.leaf <- function(node) {
    return(!(node %in% tree$edge[, 1]))
  }
  
  # Helper function to get all tips of a clade
  get_tips <- function(tree, node) {
    if (is.leaf(node)) {
      return(tree$tip.label[node])
    } else {
      children <- tree$edge[tree$edge[, 1] == node, 2]
      return(unlist(sapply(children, function(x) get_tips(tree, x))))
    }
  }
  
  # Function to recursively sort the labels in a clade
  sort_clade <- function(node) {
    if (is.leaf(node)) {
      return(node)
    } else {
      children <- which(tree$edge[, 1] == node)
      child_nodes <- tree$edge[children, 2]
      
      # Recursively sort the child nodes
      sorted_children <- sapply(child_nodes, sort_clade)
      
      # Sort the child nodes alphabetically based on their tips
      tips <- sapply(sorted_children, function(x) min(get_tips(tree, x)))
      sorted_order <- order(tips)
      
      tree$edge[children, 2] <<- sorted_children[sorted_order]
      return(node)
    }
  }
  
  # Start sorting from the root
  root <- Ntip(tree) + 1
  sort_clade(root)
  return(tree)
}

# Define a function to process each gene
process_gene <- function(gene_id) {
  # Read MSA and tree files with error handling
  msa_file <- paste0(gene_id, ".withPvCDS.aln.fa")
  tree_file <- paste0("RAxML_bipartitions.", gene_id)
  
  x <- tryCatch({
    read.msa(msa_file, format = "FASTA")
  }, error = function(e) {
    # Print the error message and return NULL
    message("Error reading MSA for ", gene_id, ": ", e$message)
    return(NULL)
  })
  
  tree <- tryCatch({
    read.tree(tree_file)
  }, error = function(e) {
    # Print the error message and return NULL
    message("Error reading tree for ", gene_id, ": ", e$message)
    return(NULL)
  })
  
  # Skip processing if reading files failed
  if (is.null(x) || is.null(tree)) {
    return(NULL)
  }
  

  
  # Process the tree
  treer <- root(tree, tree$tip.label[grepl("pvagin", tree$tip.label)], resolve.root = TRUE)
  treer$tip.label[grepl("taesti", treer$tip.label)] <- paste0(substr(treer$tip.label[grepl("taesti", treer$tip.label)], 1, 6), substr(treer$tip.label[grepl("taesti", treer$tip.label)], 24, 24))
 treer$tip.label[grepl('pvagin', treer$tip.label)]='pvagin'
#treea=keep.tip(treer, treer$tip.label[treer$tip.label %in% c("taestiA", "taestiB", "taestiD", "pvagin")])
treea=drop.tip(treer, treer$tip.label[substr(treer$tip.label,1,5)=='Pavag'])

    
  stripped_tree <- treea
  stripped_tree$edge.length <- NULL
  stripped_tree$node.label <- NULL
  
  sorted_tree <- rotate_clades(stripped_tree)
  topology <- write.tree(sorted_tree)
  
  # Extract genomic positions
  # gerarditips <- tree$tip.label[grepl("agerjg", tree$tip.label)]
  # matches <- regexec("_Chr(\\d+[A-Z]?)_(\\d+)-(\\d+)", gerarditips)
  # matches_list <- regmatches(gerarditips, matches)
    androtips=tree$tip.label[grepl('^t', tree$tip.label)]            
    matches <- regexec("hr(\\d+[A-Z]?)_(\\d+)-(\\d+)", androtips)
    matches_list <- regmatches(androtips, matches)
  
  chrom <- sapply(matches_list, function(x) paste0("Chr", x[2]))
  start <- sapply(matches_list, function(x) as.numeric(x[3]))
  end <- sapply(matches_list, function(x) as.numeric(x[4]))
  
  chrom=str_split_fixed(androtips, '_',4)[,3]
  start=as.numeric(str_split_fixed(str_split_fixed(androtips, '_', 4)[,4], '-', 2)[,1])
  end=as.numeric(str_split_fixed(str_split_fixed(androtips, '_', 4)[,4], '-', 2)[,2])
  
  tip_index=sapply(androtips, function(x) which(tree$tip.label==x))
  
  tbl=sapply(tip_index, function(x) tree$edge.length[tree$edge[,2]==x])
  
  treelength=sum(drop.tip(tree, tree$tip.label[!tree$tip.label%in%androtips])$edge.length)
  
  ultra=chronos(tree, lambda=0)  
  
  ultratbl=sapply(tip_index, function(x) ultra$edge.length[ultra$edge[,2]==x])
  
  outpos <- data.frame(
    gene = paste0(gene_id, ".1.v3.1"),
    copy = paste0(substr(androtips, 1, 7), substr(androtips, 29, 29)),
    chrom = chrom,
    start = start,
    end = end,
    tbl=tbl,
    ultratbl=ultratbl,
    treelength=treelength
  )
  
  out <- data.frame(
    gene = paste0(gene_id, ".1.v3.1"),
#    t(ks_values),
    topology = topology
  )
  
  return(list(out = out, outpos = outpos))
}



# List of gene IDs to process
#gene_ids <- c("Pavag01G000700", "Pavag01G000800", "Pavag01G000900")
genes=read.table('Pv_falist.txt')
gene_ids=gsub('.fa','',genes$V1)
                
# Process each gene and combine results, skipping NULLs
results <- lapply(gene_ids, process_gene)
results <- Filter(Negate(is.null), results)
out_list <- lapply(results, function(x) x$out)
outpos_list <- lapply(results, function(x) x$outpos)

# Combine into final data frames
out_combined <- do.call(rbind, out_list)
outpos_combined <- do.call(rbind, outpos_list)

# Print the resulting data frames
#print(out_combined)
#print(outpos_combined)

outpos_combined$topology=out_combined$topology[match(outpos_combined$gene, out_combined$gene)]
outpos_combined$common=outpos_combined$topology%in%names(tail(sort(table(outpos_combined$topology)),8))
outpos_combined$common15=outpos_combined$topology%in%names(tail(sort(table(outpos_combined$topology)),15))

outpos_combined =  outpos_combined %>% group_by(gene) %>% mutate(samechr=length(unique(substr(chrom,1,1)))==1)
                      
color_palette <- c(
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # red
  "#CC79A7", # pink
  "#999999"  # grey
)

color_palette15 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#8E44AD", "#2ECC71", "#3498DB", "#E74C3C", "#F1C40F", "#34495E", "#1ABC9C")

outmelt=melt(out_combined[,-11])

## add virginicus


# Function to determine sister lineage based on the topology
### stupid chat gpt let me just hard code this dumb thing

oo=outpos_combined[outpos_combined$common & outpos_combined$samechr,]
oo$sister_lineage=NA
oo$sister_lineage[oo$topology=='(pvagin,((taestiA,taestiB),taestiD));' & substr(oo$chrom,2,2)=='A']='B'
oo$sister_lineage[oo$topology=='(pvagin,((taestiA,taestiB),taestiD));' & substr(oo$chrom,2,2)=='B']='A'
oo$sister_lineage[oo$topology=='(pvagin,((taestiA,taestiB),taestiD));' & substr(oo$chrom,2,2)=='D']='(A,B)'
oo$sister_lineage[oo$topology=='(pvagin,((taestiA,taestiD),taestiB));' & substr(oo$chrom,2,2)=='A']='D'
oo$sister_lineage[oo$topology=='(pvagin,((taestiA,taestiD),taestiB));' & substr(oo$chrom,2,2)=='B']='(A,D)'
oo$sister_lineage[oo$topology=='(pvagin,((taestiA,taestiD),taestiB));' & substr(oo$chrom,2,2)=='D']='A'
oo$sister_lineage[oo$topology=='(pvagin,(taestiA,(taestiB,taestiD)));' & substr(oo$chrom,2,2)=='A']='(B,D)'
oo$sister_lineage[oo$topology=='(pvagin,(taestiA,(taestiB,taestiD)));' & substr(oo$chrom,2,2)=='B']='D'
oo$sister_lineage[oo$topology=='(pvagin,(taestiA,(taestiB,taestiD)));' & substr(oo$chrom,2,2)=='D']='B'


## okay, may be easier to think about "who" is sister





pdf('taesti_subgenome_proportions.pdf',14,8)
ggplot(outpos_combined[outpos_combined$common15 & !duplicated(outpos_combined$gene),], aes(x = topology, fill = topology)) +
  geom_histogram(stat='count') +
  scale_fill_manual(values = color_palette15) +
  labs(x = "Topology", y = "Count", title = "Count of Topologies") 
ggplot(outmelt, aes(x=value, group=variable, color=variable, fill=variable)) + geom_density(alpha=0.2) + scale_color_manual(values=c(color_palette, 'black'))+ scale_fill_manual(values=c(color_palette, 'black')) + xlim(0,0.2)
ggplot(outmelt, aes(x=value, group=variable, color=variable, fill=variable)) + geom_density(alpha=0.2) + scale_color_manual(values=c(color_palette, 'black'))+ scale_fill_manual(values=c(color_palette, 'black')) + xlim(0,0.2) + facet_wrap(~variable, ncol=1)
                     
dev.off()

pdf('taesti_trees_chr.pdf',14,24)
ggplot(outpos_combined[outpos_combined$common15 & grepl('taesti',outpos_combined$copy),], aes(x=start,color=topology,y=topology)) + geom_point()+facet_wrap(~chrom, ncol=1) + scale_color_manual(values=color_palette15)
ggplot(outpos_combined[outpos_combined$common15 & grepl('taesti',outpos_combined$copy),], aes(x = start, fill = topology)) +
  geom_histogram(binwidth = 1e6, position = "stack") +
  facet_wrap(~ chrom, ncol = 1) +
  scale_fill_manual(values = color_palette15) +
  labs(x = "Genomic Position (Mb)", y = "Count", title = "Stacked Histogram of Topologies in 1 Mb Windows") 
ggplot(outpos_combined[outpos_combined$samechr & outpos_combined$common15 & grepl('taesti',outpos_combined$copy),], aes(x = start, fill = topology)) +
  geom_histogram(binwidth = 1e6, position = "fill") +
  facet_wrap(~ chrom, ncol = 1) +
  scale_fill_manual(values = color_palette15) +
  labs(x = "Genomic Position (Mb)", y = "Proportion", title = "Proportional Stacked Histogram of Topologies in 1 Mb Windows") +
  theme_minimal()
ggplot(outpos_combined[outpos_combined$samechr & outpos_combined$common15 & grepl('taesti',outpos_combined$copy),], aes(x = start, fill = topology)) +
  geom_histogram(binwidth = 10e6, position = "fill") +
  facet_wrap(~ chrom, ncol = 1) +
  scale_fill_manual(values = color_palette15) +
  labs(x = "Genomic Position (Mb)", y = "Proportion", title = "Proportional Stacked Histogram of Topologies in 10 Mb Windows") +
  theme_minimal()
ggplot(outpos_combined[outpos_combined$samechr & outpos_combined$common15 & grepl('taesti',outpos_combined$copy),], aes(x = start, fill = topology)) +
  geom_histogram(binwidth = 50e6, position = "fill") +
  facet_wrap(~ chrom, ncol = 1) +
  scale_fill_manual(values = color_palette15) +
  labs(x = "Genomic Position (Mb)", y = "Proportion", title = "Proportional Stacked Histogram of Topologies in 50 Mb Windows") +
  theme_minimal()
  
  ggplot(outpos_combined[outpos_combined$samechr & 
                       outpos_combined$common15 & 
                       grepl('taesti', outpos_combined$copy),], 
       aes(x = substr(chrom,1,1), fill = topology)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = color_palette15) +
  labs(x = "", y = "Proportion", title = "Overall Proportion of Topologies Across All Chromosomes") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 40))+
  geom_text(stat = "count", aes(label = ..count../3), 
            position = position_fill(vjust = 0.5), size = 30)


ggplot(outpos_combined[outpos_combined$samechr & 
                       outpos_combined$common15 & 
                       grepl('taesti', outpos_combined$copy),], 
       aes(x = "All Chromosomes", fill = topology)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = color_palette15) +
  labs(x = "", y = "Proportion", title = "Overall Proportion of Topologies Across All Chromosomes") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 40))+
  geom_text(stat = "count", aes(label = ..count../3), 
            position = position_fill(vjust = 0.5), size = 30)


ggplot(oo, aes(x = start, fill = sister_lineage)) +
  geom_histogram(binwidth = 10e6, position = "fill") +
  facet_wrap(~ chrom, ncol = 1) +
  scale_fill_manual(values = color_palette15) +
  labs(x = "Genomic Position (Mb)", y = "Proportion", title = "Proportional Stacked Histogram of Topologies in 10 Mb Windows") +
  theme_minimal()

ggplot(oo, aes(x = chrom, fill = sister_lineage)) +
  geom_bar(position='fill') + scale_fill_manual(values=color_palette15) +
  labs(x = "Chromosome", y = "Proportion", title = "Proportional Stacked Histogram of Topologies by Chromosome") +
  theme_minimal()

ggplot(oo, aes(x = chrom, fill = sister_lineage)) +
  geom_bar(position='dodge') + scale_fill_manual(values=color_palette15) +
  labs(x = "Chromosome", y = "Proportion", title = "Topologies by Chromosome") +
  theme_minimal()

ggplot(oo%>% group_by(gene)%>% mutate(maxtbl=max(tbl))%>% filter(maxtbl<0.1), aes(x = chrom, fill = sister_lineage, y=tbl)) +
  geom_violin() + scale_fill_manual(values=color_palette15) +
  labs(x = "Chromosome", y = "Proportion", title = "Topologies by Chromosome") +
  theme_minimal()

ggplot(oo%>% group_by(gene)%>% mutate(maxtbl=max(tbl))%>% filter(maxtbl<0.1), aes(x = chrom, fill = sister_lineage, y=tbl)) +
  geom_boxplot(position='dodge') + scale_fill_manual(values=color_palette15) +
  labs(x = "Chromosome", y = "Proportion", title = "Topologies by Chromosome") +
  theme_minimal()

dev.off()
         
pdf('taesti_trees_tbl.pdf',6,4)         
ggplot(oo%>% group_by(gene)%>% mutate(maxtbl=max(tbl))%>% filter(maxtbl<0.1), aes(x = chrom, fill = sister_lineage, y=tbl)) +
  geom_violin() + scale_fill_manual(values=color_palette15) +
  labs(x = "Chromosome", y = "Proportion", title = "Topologies by Chromosome") +
  theme_minimal()

ggplot(oo%>% group_by(gene)%>% mutate(maxtbl=max(tbl))%>% filter(maxtbl<0.1), aes(x = chrom, fill = sister_lineage, y=tbl)) +
  geom_boxplot(position='dodge') + scale_fill_manual(values=color_palette15) +
  labs(x = "Chromosome", y = "Proportion", title = "Topologies by Chromosome") +
  theme_minimal()     
ggplot(oo%>% group_by(gene)%>% mutate(maxtbl=max(tbl))%>% filter(maxtbl<0.1), aes(x = chrom, fill = sister_lineage, y=tbl)) +
  geom_boxplot(position='dodge') + scale_fill_manual(values=color_palette15) +
  labs(x = "Chromosome", y = "Terminal Branch Length of Focal Chromosome") +
  theme_minimal()   + facet_wrap(substr(chrom,2,2)~substr(chrom,1,1), ncol=7, scale='free_x')
  
ggplot(oo%>% group_by(gene)%>% mutate(maxtbl=max(tbl))%>% filter(maxtbl<0.1), aes(x = substr(chrom,2,2), fill = sister_lineage, y=tbl)) +
  geom_boxplot(position='dodge') + scale_fill_manual(values=color_palette15) +
  labs(x = "Subgenome", y = "Terminal Branch Length of Focal Chromosome") +
  theme_minimal()                    

ggplot(oo, aes(x = chrom, fill = sister_lineage, y=ultratbl)) +
  geom_boxplot(position='dodge') + scale_fill_manual(values=color_palette15) +
  labs(x = "Chromosome", y = "Ultrametric Terminal Branch Length of Focal Chromosome") +
  theme_minimal()   + facet_wrap(substr(chrom,2,2)~substr(chrom,1,1), ncol=7, scale='free_x')
  
ggplot(oo, aes(x = substr(chrom,2,2), fill = sister_lineage, y=ultratbl)) +
  geom_boxplot(position='dodge') + scale_fill_manual(values=color_palette15) +
  labs(x = "Subgenome", y = "Ultrametric Terminal Branch Length of Focal Chromosome") +
  theme_minimal()    

ggplot(oo, aes(x = chrom, fill = sister_lineage, y=treelength)) +
  geom_boxplot(position='dodge') + scale_fill_manual(values=color_palette15) +
  labs(x = "Chromosome", y = "Tree Length of Focal Chromosome") + coord_cartesian(y=c(0,0.15))+
  theme_minimal()   + facet_wrap(substr(chrom,2,2)~substr(chrom,1,1), ncol=7, scale='free_x')
  
ggplot(oo, aes(x = substr(chrom,2,2), fill = sister_lineage, y=treelength)) +
  geom_boxplot(position='dodge') + scale_fill_manual(values=color_palette15) +
  labs(x = "Subgenome", y = "Tree Length of Focal Chromosome") + coord_cartesian(y=c(0,0.15))+
  theme_minimal()    


dev.off()


anova_result <- aov(treelength ~ paste(sister_lineage), data = oo[substr(oo$chrom,2,2)=='A',])
pairwise.t.test(oo$treelength[substr(oo$chrom,2,2)=='D'], oo$sister_lineage[substr(oo$chrom,2,2)=='D'], p.adjust.method = "bonferroni")

####
# Define a function to process each gene
process_trees_thomas <- function(tree_file='wheat_all_trees.txt') {
  # Read MSA and tree files with error handling
#  tree_file <- "wheat_all_trees.txt"
  
  
  tree <- tryCatch({
    read.tree(tree_file)
  }, error = function(e) {
    # Print the error message and return NULL
    message("Error reading tree  ", ": ", e$message)
    return(NULL)
  })
  
  
  out <- data.frame(
    gene = 1:length(tree),
   chrom = NA,
   topology=NA
  )

  for(i in 1:length(tree)){
  if(sum(grepl('TAES', tree[[i]]$tip.label))>2){
#  print('tree processed')
  # Process the tree
  tre=tree[[i]]
  treer=tre
  treer$tip.label[grepl("TAES", treer$tip.label)] <- paste0("TAES_", substr(treer$tip.label[grepl("TAES", treer$tip.label)], 8, 8))
# treer$tip.label[grepl('pvagin', treer$tip.label)]='pvagin'
#treea=keep.tip(treer, treer$tip.label[treer$tip.label %in% c("taestiA", "taestiB", "taestiD", "pvagin")])
treea=keep.tip(treer, treer$tip.label[grepl("TAES", treer$tip.label)])

    
  stripped_tree <- treea
  stripped_tree$edge.length <- NULL
  stripped_tree$node.label <- NULL
  
  sorted_tree <- rotate_clades(stripped_tree)
  topology <- write.tree(sorted_tree)
  
  # Extract genomic positions
  # gerarditips <- tree$tip.label[grepl("agerjg", tree$tip.label)]
  # matches <- regexec("_Chr(\\d+[A-Z]?)_(\\d+)-(\\d+)", gerarditips)
  # matches_list <- regmatches(gerarditips, matches)
    androtips=treer$tip.label[grepl('TAES', treer$tip.label)]     
 #   print(androtips)       
    fullnames=tre$tip.label[grepl('TAES', treer$tip.label)]
  

  if(all(androtips %in% c('TAES_A', 'TAES_B', 'TAES_D')) & length(androtips)==3){
  chrom=substr(str_split_fixed(fullnames, '_',4)[,2],1,1)


  out$chrom[i]=chrom
  out$topology[i]=topology
}
  
  }}
  return(out)
}
            
            
             
# Process each gene and combine results, skipping NULLs
out=process_trees_thomas()
# Print the resulting data frames
#print(out_combined)
#print(outpos_combined)

out$common=out$topology%in%names(tail(sort(table(out$topology)),8))
out$common15=out$topology%in%names(tail(sort(table(out$topology)),3))

library(tidyr)
oo <- out[out$common15, ] %>%
  uncount(3) %>%
  mutate(new_column = rep(c("A", "B", "D"), length.out = n()),   # Create the new column
         chrSG = paste0(chrom, new_column))  # Concatenate 'chrom' with new_column

oo$sister_lineage=NA
oo$sister_lineage[oo$topology=='((TAES_A,TAES_B),TAES_D);' & substr(oo$chrSG,2,2)=='A']='B'
oo$sister_lineage[oo$topology=='((TAES_A,TAES_B),TAES_D);' & substr(oo$chrSG,2,2)=='B']='A'
oo$sister_lineage[oo$topology=='((TAES_A,TAES_B),TAES_D);' & substr(oo$chrSG,2,2)=='D']='(A,B)'
oo$sister_lineage[oo$topology=='((TAES_A,TAES_D),TAES_B);' & substr(oo$chrSG,2,2)=='A']='D'
oo$sister_lineage[oo$topology=='((TAES_A,TAES_D),TAES_B);' & substr(oo$chrSG,2,2)=='B']='(A,D)'
oo$sister_lineage[oo$topology=='((TAES_A,TAES_D),TAES_B);' & substr(oo$chrSG,2,2)=='D']='A'
oo$sister_lineage[oo$topology=='(TAES_A,(TAES_B,TAES_D));' & substr(oo$chrSG,2,2)=='A']='(B,D)'
oo$sister_lineage[oo$topology=='(TAES_A,(TAES_B,TAES_D));' & substr(oo$chrSG,2,2)=='B']='D'
oo$sister_lineage[oo$topology=='(TAES_A,(TAES_B,TAES_D));' & substr(oo$chrSG,2,2)=='D']='B'

#outpos_combined =  outpos_combined %>% group_by(gene) %>% mutate(samechr=length(unique(substr(chrom,1,1)))==1)
 outpos_combined=out         
          
pdf('taesti_thomas_trees_chr.pdf',14,10)
ggplot(out[outpos_combined$common15,], aes(x=chrom,fill=topology)) + 
geom_bar(position='fill') + scale_fill_manual(values=color_palette15) + geom_hline(yintercept=c(0.33,0.66), lty='dashed', color='gray')

ggplot(oo, aes(x=chrSG,fill=sister_lineage)) + 
geom_bar(position='fill') + scale_fill_manual(values=color_palette15) 

ggplot(outpos_combined[ 
                       outpos_combined$common15,], 
       aes(x = "All Chromosomes", fill = topology)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = color_palette15) +
  labs(x = "", y = "Proportion", title = "Overall Proportion of Topologies Across All Chromosomes") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 40))+
  geom_text(stat = "count", aes(label = ..count..), 
            position = position_fill(vjust = 0.5), size = 30)

dev.off()
                        
          
          