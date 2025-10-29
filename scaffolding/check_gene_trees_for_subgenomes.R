library(ape)
library(stringr)

## on blfs1

input_dir='.'
output_file='~/transfer/all_hconto_trees.tre'
tree_files <- list.files(input_dir, pattern = "00$", full.names = TRUE)

process_tree_file <- function(file) {
  # Read the tree
  trees <- read.tree(file)

  # Handle cases where the file contains a single tree (not a list)
  if (class(trees) == "phylo") {
    trees <- list(trees)
  }

  # Subset trees: keep tips with "hconto" or "pvagin" in the name
  subset_trees <- lapply(trees, function(tree) {
    tips_to_keep <- tree$tip.label[str_detect(tree$tip.label, "hconto|pvagin")]
    drop.tip(tree, setdiff(tree$tip.label, tips_to_keep))
  })

  # Combine subset trees into a single list
  return(subset_trees)
}



all_subset_trees <- lapply(tree_files, process_tree_file)

# Flatten the list of lists into a single list of trees
all_subset_trees <- do.call(c, all_subset_trees)

# Save all subset trees into a single output file
write.tree(all_subset_trees, file = output_file)

cat("Subset trees written to:", output_file, "\n")

## on cbsuxm01

subset_trees=read.tree('~/transfer/all_hconto_trees.tre')
# Function to check if a tree has any `hconto` tips
has_hconto_tips <- function(tree) {
  any(str_detect(tree$tip.label, "hconto"))
}

# Filter trees that contain at least one `hconto` tip
filtered_trees <- subset_trees[sapply(subset_trees, has_hconto_tips)]

# If no matching trees remain after filtering
if (length(filtered_trees) == 0) {
  stop("No trees contain any `hconto` tips after filtering.")
}

cat("Number of trees with `hconto` tips:", length(filtered_trees), "\n")





agp_mapping <- read.table("hconto.agp", header = F, stringsAsFactors = FALSE)

# Example of expected column names: agp_mapping$Scaffold, agp_mapping$Chromosome
scaffold_to_chr <- setNames(agp_mapping$V1, agp_mapping$V6)

convert_scaffold_to_chr <- function(tree) {
  # Extract all unique scaffolds in the tree's tip labels
  scaffolds_in_tree <- unique(str_extract(tree$tip.label, "scaf_[0-9]+"))
  
#   # Check if all scaffolds are in the AGP mapping
#   if (any(is.na(scaffolds_in_tree)) || !all(scaffolds_in_tree %in% names(scaffold_to_chr))) {
#     return(NULL)  # Drop the tree if any scaffold is missing
#   }

  # Simplify and standardize tip labels
  tree$tip.label <- sapply(tree$tip.label, function(label) {
    # Extract scaffold name
    scaffold_match <- str_extract(label, "scaf_[0-9]+")
    
    # Check if the label is a pvagin outgroup OR an hconto ingroup
    if (str_detect(label, "^pvagin")) {
      # Extract chromosome if directly part of the `pvagin` tip label
      pvagin_chromosome <- str_extract(label, "Chr[0-9A-Za-z]+")
      if (!is.na(pvagin_chromosome)) {
        # If a chromosome is present, keep it in the label
        return(paste0("pvagin_", pvagin_chromosome))
      } else {
        # If no `ChrXX` is found, fallback to NA (shouldn't happen often here)
        return("pvagin_NA")
      }
    } else {
      # Otherwise, handle hconto or other species that rely on scaffold-to-chromosome mapping
      chromosome <- scaffold_to_chr[scaffold_match]
      species_name <- str_extract(label, "^[^_]+")  # Extract species name
      return(paste0(species_name, "_", chromosome))
    }
  })

  return(tree)
}


# Apply the conversion to all subset trees
trees_with_chr <- lapply(subset_trees, convert_scaffold_to_chr)
trees_with_chr <- trees_with_chr[!sapply(trees_with_chr, is.null)]


# Ensure pvagin is used as the root for all trees
root_trees_on_pvagin <- function(tree) {
  # Find the outgroup (tip containing pvagin)
  outgroup <- tree$tip.label[grepl("pvagin", tree$tip.label)]
  
  if (length(outgroup) != 1) {
    stop("Each tree must contain exactly one `pvagin` tip for rooting.")
  }
  
  # Re-root the tree on the pvagin outgroup
  tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)
  
  return(tree)
}

rooted_trees <- lapply(trees_with_chr, root_trees_on_pvagin)

cat("All trees successfully rooted on pvagin.\n")



# Convert rooted trees to Newick strings
topology_strings <- sapply(rooted_trees, write.tree)

clean_newick <- function(newick_string) {
  # Remove all branch lengths (anything after `:`) using regex
  no_branch_lengths <- gsub(":.*?(\\)|,)", "\\1", newick_string)
  # Remove node labels (e.g., `Root` or `100`) using regex
  cleaned_string <- gsub("\\)\\w+", ")", no_branch_lengths)
  return(cleaned_string)
}
# Clean all topology strings (remove branch lengths and node labels)
cleaned_topology_strings <- sapply(topology_strings, clean_newick)


cleaned_topology_counts <- as.data.frame(table(cleaned_topology_strings))
colnames(cleaned_topology_counts) <- c("topology", "count")  # Rename columns

cleaned_topology_counts <- cleaned_topology_counts %>%
  mutate(
    pvagin_chr = str_extract(topology, "Chr[0-9A-Za-z]+")
  )

# View the updated table
cleaned_topology_counts %>% filter(pvagin_chr=='Chr01')
  arrange(-count) %>%
  head(20)




# Step 1: Normalize pvagin naming in topologies
filtered_topologies <- cleaned_topology_counts %>%
  mutate(topology = as.character(topology)) %>%  # Ensure topology is character
  mutate(topology = str_replace_all(topology, "(?i)pvagin", "pvagin")) %>%  # Standardize pvagin
  filter(str_detect(topology, "pvagin"))  # Keep only topologies with pvagin

# Step 2: Filter for topologies with exactly four hconto tips
filter_four_hconto <- function(topology) {
  hconto_tips <- str_extract_all(topology, "hconto_chr[0-9]+_[0-9]+")[[1]]
  return(length(hconto_tips) == 4)
}

filtered_topologies <- filtered_topologies %>%
  filter(sapply(topology, filter_four_hconto))

# Step 3: Canonicalize topologies while keeping the root
canonicalize_rooted_topology <- function(topology) {
  # Parse the tree
  tree <- tryCatch({
    read.tree(text = topology)
  }, error = function(e) {
    message("Invalid topology detected and skipped:\n", topology)
    return(NULL)  # Skip invalid topology
  })
  
  # Check if any tip label starts with "pvagin"
  if (is.null(tree) || !any(grepl("^pvagin", tree$tip.label))) {
    message("Skipping topology: pvagin not found as a tip\n", topology)
    return(NULL)
  }
  
  # Canonicalize the tree: retain root and reorder branches
  tree <- ape::ladderize(tree, right = FALSE)  # Order branches consistently
  
  # Rewrite into canonical Newick string
  return(write.tree(tree))
}

filtered_topologies <- filtered_topologies %>%
  mutate(
    canonical_topology = sapply(topology, canonicalize_rooted_topology)
  ) %>%
  filter(!is.na(canonical_topology))  # Remove rows with NULL canonical_topology

# Step 4: Group by canonicalized topologies and count
canonical_topology_counts <- filtered_topologies %>%
  group_by(canonical_topology, pvagin_chr) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  arrange(pvagin_chr,desc(count))





### now convert to a mrc

canonicalize_tree <- function(tree) {
  # Internal function to reorder the edges of a phylo object recursively
  reorder_subtree <- function(node, tree) {
    if (node <= length(tree$tip.label)) {
      # Return terminal node (tip)
      return(tree$tip.label[node])
    }
    
    # Identify the indices of children (subtrees)
    children <- which(tree$edge[, 1] == node)
    child_nodes <- tree$edge[children, 2]
    
    # Recursively order child nodes
    ordered_children <- lapply(child_nodes, reorder_subtree, tree = tree)
    
    # Sort children lexicographically and rebuild subtree
    ordered_children <- ordered_children[order(sapply(ordered_children, function(x) {
      if (is.character(x)) return(x)  # Sort by tip label lexicographically
      return(paste(sort(unlist(x)), collapse = ""))  # Internal node
    }))]
    
    return(ordered_children)
  }
  
  # Start from the root and reorder the entire tree recursively
  reordered_tree <- reorder_subtree(length(tree$tip.label) + 1, tree)
  
  # Rewrite the reordered tree into a newick-compatible object
  return(write.tree(tree))
}
canonicalize_topologies <- function(topologies) {
  # Convert Newick strings into ape trees
  trees <- lapply(topologies, function(topology) read.tree(text = topology))
  
  # Canonicalize each tree to remove rotational ambiguities
  canonicalized_trees <- lapply(trees, function(tree) {
    tryCatch({
      canonicalize_tree(tree)
    }, error = function(e) {
      message("Error canonicalizing tree: ", e$message)
      return(NULL)
    })
  })
  
  # Return canonicalized Newick strings, removing any NULL results
  return(canonicalized_trees)
}

expand_topologies <- canonical_topology_counts %>%
  group_by(pvagin_chr) %>%
  summarise(
    expanded_topologies = list(rep(canonical_topology, count)),  # Replicate trees by count
    .groups = "drop"
  )

standardize_tip_labels <- function(trees) {
  # Get all unique tip labels across all trees
  all_tips <- unique(unlist(lapply(trees, function(tree) tree$tip.label)))
  
  # Add missing tips as "zero-length" polytomies to each tree
  standardized_trees <- lapply(trees, function(tree) {
    # Identify missing tips
    missing_tips <- setdiff(all_tips, tree$tip.label)
    
    # If no tips are missing, return the tree as-is
    if (length(missing_tips) == 0) {
      return(tree)
    }
    
    # Create dummy trees for missing tips
    dummy_trees <- lapply(missing_tips, function(tip) {
      read.tree(text = paste0("(", tip, ");"))  # Create single-tip trees
    })
    
    # Add the missing tips to the original tree
    tree <- Reduce(function(x, y) bind.tree(x, y, where = "root"), dummy_trees, init = tree)
    return(tree)
  })
  
  return(standardized_trees)
}

# Canonicalize topologies for each chromosome
canonical_species_trees <- expand_topologies %>%
  mutate(expanded_topologies = lapply(expanded_topologies, canonicalize_topologies))
  
  
  build_species_tree <- function(canonicalized_topologies) {
  # Convert canonicalized topologies back to trees
  trees <- lapply(canonicalized_topologies, function(topology) read.tree(text = topology))
  
  # Standardize tips if necessary (optional, ensure consistency across trees)
  standardized_trees <- standardize_tip_labels(trees)
  
  # Compute consensus tree
  return(consensus(standardized_trees, p = 0.5))  # Majority-rule consensus tree
}

# Apply consensus tree building per chromosome
species_trees <- canonical_species_trees %>%
  mutate(species_tree = lapply(expanded_topologies, build_species_tree))
  
  
  
  
  
# Step 2: Build a weighted species tree for each chromosome
build_weighted_species_tree <- function(expanded_topologies) {
  # Convert Newick strings to ape trees
  trees <- lapply(expanded_topologies, function(topology) read.tree(text = topology))
  
  # Compute the majority-rule consensus tree
  consensus <- consensus(trees, p = 0.5, check.labels=T)  # Majority-rule consensus (adjust p as needed)
  return(consensus)
}

weighted_species_trees <- expand_topologies %>%
  mutate(species_tree = lapply(expanded_topologies, build_weighted_species_tree))



#### OMG DO THIS MYSELF OR MAKE CHATGPT THINK ABOUT IT AGAIN, IT"S OFF ON SOME WEIRD TANGENT
