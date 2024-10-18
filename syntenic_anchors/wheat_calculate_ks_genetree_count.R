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
  
  outpos <- data.frame(
    gene = paste0(gene_id, ".1.v3.1"),
    copy = paste0(substr(androtips, 1, 7), substr(androtips, 29, 29)),
    chrom = chrom,
    start = start,
    end = end
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
       aes(x = "All Chromosomes", fill = topology)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = color_palette15) +
  labs(x = "", y = "Proportion", title = "Overall Proportion of Topologies Across All Chromosomes") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 40))+
  geom_text(stat = "count", aes(label = ..count../3), 
            position = position_fill(vjust = 0.5), size = 30)

dev.off()
                      

############# OLD - GERARDI ONLY!!                      
pdf('gerardi_subgenome_proportions.pdf',14,8)
ggplot(outpos_combined[outpos_combined$common & !duplicated(outpos_combined$gene),], aes(x = topology, fill = topology)) +
  geom_histogram(stat='count') +
  scale_fill_manual(values = color_palette) +
  labs(x = "Haplotype", y = "Proportion", title = "Proportional Stacked Histogram of Topologies") 
ggplot(outmelt, aes(x=value, group=variable, color=variable, fill=variable)) + geom_density(alpha=0.2) + scale_color_manual(values=c(color_palette, 'black'))+ scale_fill_manual(values=c(color_palette, 'black')) + xlim(0,0.2)
ggplot(outmelt, aes(x=value, group=variable, color=variable, fill=variable)) + geom_density(alpha=0.2) + scale_color_manual(values=c(color_palette, 'black'))+ scale_fill_manual(values=c(color_palette, 'black')) + xlim(0,0.2) + facet_wrap(~variable, ncol=1)
                     
dev.off()

## make a nice figure
## plot each topology in the color it's plotted
top=read.tree(text='(((A.1,A.2),(B.1,B.2)),(C.1,C.2));')
mid=read.tree(text='(((A.1,A.2),(C.1,C.2)),(B.1,B.2));')               
low=read.tree(text='((A.1,A.2),((B.1,B.2),(C.1,C.2)));')
genecount=outpos_combined[!duplicated(outpos_combined$gene),] %>% group_by(topology) %>% summarize(genes=n())
gco=genecount %>% filter(genes>1000) %>% add_row(topology='Other', genes=sum(genecount$genes[genecount$genes<=1000]))
genecount=genecount %>% arrange(desc(genes))
                      
pdf('gerardi_subgenome_fig.pdf',8,3)
tp=ggtree(top, layout='slanted')+ geom_tree(color = "#009E73", layout='slanted')+ scale_color_manual("#009E73") + geom_tiplab(size=1.5) + xlim(-3,9)
mp=ggtree(mid, layout='slanted')+ geom_tree(color = "#F0E442", layout='slanted') + scale_color_manual("#F0E442") + geom_tiplab(size=1.5)+ xlim(-3,9)
bp=ggtree(low, layout='slanted')+ geom_tree(color = "#CC79A7", layout='slanted') + scale_color_manual("#CC79A7") + geom_tiplab(size=1.5)+ xlim(-3,9)
empty_plot <- ggplot() + 
  theme_void() + 
  theme(
    panel.background = element_blank(), 
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
trees=plot_grid(tp,mp,bp, ncol=1)
treesr=plot_grid(tp,mp,bp, empty_plot, nrow=1)
genecount$topology_highlight='gray'
genecount$topology_highlight[1]="#009E73"
genecount$topology_highlight[2]="#F0E442"
genecount$topology_highlight[3]="#CC79A7"
tops=ggplot(genecount, aes(x=reorder(topology, -genes), y=genes, fill=topology_highlight, color=topology_highlight)) + geom_bar(stat='identity',position='identity') + scale_color_identity()  +scale_fill_identity()                   
plot_grid(tops, trees, ncol=2, rel_widths=c(1,0.1))             

p2_grob <- ggplotGrob(treesr)

# Overlay the grob on the first plot
tops + annotation_custom(grob = p2_grob, xmin = 10, xmax = 130, ymin = 1000, ymax = 2000)

tops2=ggplot(genecount, aes(x=reorder(topology_highlight, -genes), group=topology, y=genes, fill=topology_highlight)) + geom_bar(stat='identity',position='stack') + scale_color_identity()  +scale_fill_identity()                   
tops2
plot_grid(tops2 + theme(axis.title.x=element_blank()) + scale_x_discrete(labels=c('((A,B),C)', '((A,C),B)', '((B,C),A)', 'Other')), treesr, nrow=2, rel_heights=c(1,0.6), align='hv')             

plot_grid(tops2 + theme(axis.title.x=element_blank()) + scale_x_discrete(labels=c('((A,B),C)', '((A,C),B)', '((B,C),A)', 'Other')) + theme(plot.margin = margin(6, 0, 0, 0)), treesr+ theme(plot.margin = margin(0, 0, 0, 0)), nrow=2, rel_heights=c(1,0.6), align='hv')             
                     
                      
dev.off()

                      
                      
pdf('gerardi_trees_chr.pdf',14,24)
ggplot(outpos_combined[outpos_combined$common & grepl('1',outpos_combined$copy),], aes(x=start,color=topology,y=topology)) + geom_point()+facet_wrap(~chrom, ncol=1) + scale_color_manual(values=color_palette)
ggplot(outpos_combined[outpos_combined$common & grepl('1',outpos_combined$copy),], aes(x = start, fill = topology)) +
  geom_histogram(binwidth = 1e6, position = "stack") +
  facet_wrap(~ chrom, ncol = 1) +
  scale_fill_manual(values = color_palette) +
  labs(x = "Genomic Position (Mb)", y = "Count", title = "Stacked Histogram of Topologies in 1 Mb Windows") 
ggplot(outpos_combined[outpos_combined$common & grepl('1',outpos_combined$copy),], aes(x = start, fill = topology)) +
  geom_histogram(binwidth = 1e6, position = "fill") +
  facet_wrap(~ chrom, ncol = 1) +
  scale_fill_manual(values = color_palette) +
  labs(x = "Genomic Position (Mb)", y = "Proportion", title = "Proportional Stacked Histogram of Topologies in 1 Mb Windows") +
  theme_minimal()
ggplot(outpos_combined[outpos_combined$common & grepl('2',outpos_combined$copy),], aes(x=start,color=topology,y=topology)) + geom_point()+facet_wrap(~chrom, ncol=1) + scale_color_manual(values=color_palette)
ggplot(outpos_combined[outpos_combined$common & grepl('2',outpos_combined$copy),], aes(x = start, fill = topology)) +
  geom_histogram(binwidth = 1e6, position = "stack") +
  facet_wrap(~ chrom, ncol = 1) +
  scale_fill_manual(values = color_palette) +
  labs(x = "Genomic Position (Mb)", y = "Count", title = "Stacked Histogram of Topologies in 1 Mb Windows") 
ggplot(outpos_combined[outpos_combined$common & grepl('2',outpos_combined$copy),], aes(x = start, fill = topology)) +
  geom_histogram(binwidth = 1e6, position = "fill") +
  facet_wrap(~ chrom, ncol = 1) +
  scale_fill_manual(values = color_palette) +
  labs(x = "Genomic Position (Mb)", y = "Proportion", title = "Proportional Stacked Histogram of Topologies in 1 Mb Windows") 

dev.off()
                      