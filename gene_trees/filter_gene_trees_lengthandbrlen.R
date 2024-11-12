# Load necessary packages
library(ape)
library(phytools)
library(Biostrings)

## this is all chatgpt!!!

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check for required arguments
if (length(args) < 2) {
  stop("Usage: Rscript script_name.R <tree_file> <alignment_file>")
}

# Assign arguments to variables
tree_file <- args[1]
alignment_file <- args[2]

# Load input files
tree <- read.tree(tree_file)
fasta <- readDNAStringSet(alignment_file)

# Step 1: Reroot the tree to tips starting with "pvagin"
pvagin_tips <- grep("^pvagin", tree$tip.label, value = TRUE)

if (length(pvagin_tips) > 0) {
  tree <- reroot(tree, which(tree$tip.label == pvagin_tips[1]))
} else {
  stop("No tips found starting with 'pvagin'")
}

# Step 2: Calculate root-to-tip distances
node_heights <- nodeHeights(tree)
root_to_tip_distances <- node_heights[tree$edge[, 2] <= length(tree$tip.label), 2]
names(root_to_tip_distances) <- tree$tip.label

# Step 3: Identify tips with root-to-tip distances > 3 standard deviations
distance_threshold <- mean(root_to_tip_distances) + 3 * sd(root_to_tip_distances)
long_distance_tips <- names(root_to_tip_distances)[root_to_tip_distances > distance_threshold]

# Step 4: Filter sequences based on nucleotide length
nucleotide_counts <- rowSums(letterFrequency(fasta, letters = c("A", "T", "G", "C")))
average_nucleotide_length <- mean(nucleotide_counts)
length_threshold <- 0.25 * average_nucleotide_length
short_sequences <- names(fasta)[nucleotide_counts < length_threshold]

# Step 5: Combine filters and prune
tips_to_remove <- union(long_distance_tips, short_sequences)
pruned_tree <- drop.tip(tree, tips_to_remove)
filtered_fasta <- fasta[!(names(fasta) %in% tips_to_remove)]

# Generate output filenames based on input
tree_output_file <- sub("$", "_filtered.tre", tree_file)
alignment_output_file <- sub("\\.fa$", "_filtered.fasta", alignment_file)

# Save outputs
write.tree(pruned_tree, tree_output_file)
cat("Filtered and rooted tree saved as:", tree_output_file, "\n")

if(length(tips_to_remove)>0){
writeXStringSet(filtered_fasta, alignment_output_file)
cat("Filtered alignment saved as:", alignment_output_file, "\n")

}

print(tips_to_remove)


