library(ape)
library(stringr)
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(purrr)


## on blfs1 /users/mcs368/panand_gene_trees/gene_trees_with_paspalumCDS_Dec2023
species='blagur'
ploidy=6
input_dir='.'
output_file=paste0('~/transfer/all_', species,'_trees.tre')
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
    tips_to_keep <- tree$tip.label[str_detect(tree$tip.label, paste0(species,"|pvagin"))]
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

### chatgpt got out of control on this. i'll try again with thinkikng a bit more about individual topologies...

## essentially, I need to canonnicalize rooted topology, and then count
## maybe not even care about majority rule consensus tree?
## but it would be nice to have on a chromosome by chromsome basis...



subset_trees=read.tree(paste0('~/transfer/all_',species,'_trees.tre'))
# Function to check if a tree has any `hconto` tips
has_hconto_tips <- function(tree) {
  any(str_detect(tree$tip.label, species))
}

# Filter trees that contain at least one `hconto` tip
filtered_trees <- subset_trees[sapply(subset_trees, has_hconto_tips)]

# If no matching trees remain after filtering
if (length(filtered_trees) == 0) {
  stop("No trees contain any species tips after filtering.")
}

cat("Number of trees with species tips:", length(filtered_trees), "\n")





agp_mapping <- read.table(paste0(species,".agp"), header = F, stringsAsFactors = FALSE)

# Example of expected column names: agp_mapping$Scaffold, agp_mapping$Chromosome
scaffold_to_chr <- setNames(agp_mapping$V1, agp_mapping$V6)

convert_scaffold_to_chr <- function(tree) {
  # Extract all unique scaffolds in the tree's tip labels
  scaffolds_in_tree <- unique(str_extract(tree$tip.label, "(alt-)?scaf_[0-9]+"))
  
#   # Check if all scaffolds are in the AGP mapping
#   if (any(is.na(scaffolds_in_tree)) || !all(scaffolds_in_tree %in% names(scaffold_to_chr))) {
#     return(NULL)  # Drop the tree if any scaffold is missing
#   }

  # Simplify and standardize tip labels
  tree$tip.label <- sapply(tree$tip.label, function(label) {
    # Extract scaffold name
     scaffold_match <- str_extract(label, "(alt-)?scaf_[0-9]+")
    
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
           # Handle cases where the scaffold is missing from mapping
      if (is.na(chromosome)) {
        return(paste0(species_name, "_NA"))  # Assign NA if no scaffold-chromosome mapping is found
      } else {
        return(paste0(species_name, "_", chromosome))  # Combine species name with chromosome
      }

    }
  })

  return(tree)
}


# Apply the conversion to all subset trees
trees_with_chr <- lapply(subset_trees, convert_scaffold_to_chr)


convert_scaffold_to_chr_metadata_optimized <- function(tree, agp_mapping) {
  # Ensure scaffolds in AGP mapping are formatted consistently
#   agp_mapping <- agp_mapping %>%
#     mutate(V6 = str_replace(V6, "scaf_0?([0-9]+)", "scaf_\\1"))
#   
  # Check if tree is invalid and handle gracefully
  if (is.null(tree) || !("tip.label" %in% names(tree)) || length(tree$tip.label) == 0) {
    return(NULL)
  }

  # Extract metadata for all tip labels
  df <- data.frame(
    tip_label = tree$tip.label,
    scaffold_match = str_extract(tree$tip.label, "(alt-)?scaf_[0-9]+"),  # Extract scaffold
    scaffold_start = as.numeric(str_match(tree$tip.label, "_([0-9]+)-([0-9]+)")[, 2]),  # Extract start
    scaffold_end = as.numeric(str_match(tree$tip.label, "_([0-9]+)-([0-9]+)")[, 3])     # Extract end
  )
  
  # Add extracted coordinates for scaffold-based or chromosome-based positions
  df <- df %>%
    mutate(
      pvagin_chromosome = str_extract(tip_label, "Chr[0-9A-Za-z]+")
    )
  
  # Perform a left join with the AGP mapping file to merge scaffold-based metadata
  df <- df %>%
    left_join(agp_mapping, by = c("scaffold_match" = "V6")) %>%  # Join by scaffold column in AGP Mapping ("V6")
    mutate(
      species = ifelse(str_detect(tip_label, "^pvagin"), "pvagin", str_extract(tip_label, "^[^_]+")),
      chromosome = ifelse(str_detect(tip_label, "^pvagin"), pvagin_chromosome, V1),   # Use pvagin chromosome or AGP mapping chromosome
      chrom_start_adjusted = case_when(
        str_detect(tip_label, "^pvagin") ~ scaffold_start,  # For pvagin, use direct chromosome start position
        !is.na(scaffold_start) & !is.na(V2) ~ V2 + scaffold_start - 1,  # Adjust scaffold start position to chromosome start
        TRUE ~ NA_real_
      ),
      chrom_end_adjusted = case_when(
        str_detect(tip_label, "^pvagin") ~ scaffold_end,    # For pvagin, use direct chromosome end position
        !is.na(scaffold_end) & !is.na(V2) ~ V2 + scaffold_end - 1,      # Adjust scaffold end position to chromosome start
        TRUE ~ NA_real_
      )
    ) %>%
    select(
      tip_label, species, chromosome, chrom_start_adjusted, chrom_end_adjusted, scaffold_start, scaffold_end
    ) # Keep only relevant columns
  
  # Filter out rows with invalid data if necessary
  if (nrow(df) == 0) {
    return(NULL)
  }
  
  return(df)
}

chr_with_chr <- lapply(subset_trees, function(x) convert_scaffold_to_chr_metadata_optimized(x, agp_mapping))


## remove from both!
chr_with_chr <- chr_with_chr[!sapply(trees_with_chr, is.null)]
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

# Canonicalize tree ensuring sorted clades
canonicalize_tree <- function(tree) {
  # Recursive function to process clades
  process_clades <- function(node) {
    # Get descendant nodes for the current node
    descendants <- tree$edge[tree$edge[, 1] == node, 2]
    
    # Check if the node is a tip (leaf)
    if (length(descendants) == 0) {
      return(tree$tip.label[node])  # Return the tip label directly
    }
    
    # Process each descendant recursively
    child_clades <- sapply(descendants, process_clades)
    
    # Convert to atomic structure: combine child clades and sort them
    sorted_clades <- sort(child_clades)
    
    # Combine sorted clades into Newick representation
    return(paste0("(", paste(sorted_clades, collapse = ","), ")"))
  }
  
  # Identify the root node (unique parent that isn't listed as a child)
  root_node <- setdiff(tree$edge[, 1], tree$edge[, 2])[1]
  
  # Generate canonical Newick string starting from the root
  canonical_tree <- process_clades(root_node)
  
  # Return the complete Newick string
  return(paste0(canonical_tree, ";"))
}





# Convert rooted trees to Newick strings
### tip count depends on ploidy, so ploidy+pvagin
fivetips=sapply(rooted_trees, function(x){ length(x$tip.label)==ploidy+1})
topology_strings <- sapply(rooted_trees[fivetips], canonicalize_tree)
chr_with_chr_final=chr_with_chr[fivetips]

# Extract only the `pvagin` rows from each data frame in the list
pvagin_rows <- lapply(chr_with_chr_final, function(df) {
  df[df$species == "pvagin", ]
})

# Combine all extracted `pvagin` rows into a single data frame
pvagin_df <- do.call(rbind, pvagin_rows)

#### AAHAAHAHAH THIS IS AMAZING!!!! IT'S WORKIHNG!!!!

## now i just have to plot the topologies along each chromosome!!

























### now convert to a mrc
#topologies=data.frame(topology=topology_strings, pvaginChr=sub(".*pvagin_([A-Za-z0-9]{5}).*", "\\1", topology_strings))

# chr1=lapply(topologies[topologies$pvaginChr=='Chr01',]$topology, function(x) read.tree(text=x))
# 
# chr1tips=chr1[unlist(lapply(chr1, function(x) all(x$tip.label%in%c("hconto_chr3_1", "hconto_chr3_2", "hconto_chr3_3", "hconto_chr3_4", 
# "pvagin_Chr01"))))]
# 
# # Combine all tree objects into a single `multiPhylo` object
# chr1mrc=consensus(chr1tips, p=0.5, check.labels=T)

topologies=data.frame(topology=topology_strings, pvaginChr=sub(".*pvagin_([A-Za-z0-9]{5}).*", "\\1", topology_strings), pvaginStart=pvagin_df$chrom_start_adjusted)

chromosomes <- paste0("Chr", sprintf("%02d", 1:10)) # Generates Chr01, Chr02, ..., Chr10


# Iterate over all chromosomes
chromosome_trees <- lapply(chromosomes, function(chrom) {
  # Filter rows for the current chromosome
  topologies_chr <- topologies[topologies$pvaginChr == chrom, ]
  newick_strings <- topologies_chr$topology
  
  # Parse all Newick strings into tree objects
  trees <- lapply(newick_strings, function(x) read.tree(text = x))
  
  # Convert trees to Newick format for topology comparison
  topology_strings <- sapply(trees, write.tree)  # Get the Newick string for each tree
  
  # Identify the most common topology
  most_common_topology <- names(sort(table(topology_strings), decreasing = TRUE))[1]
  
  # Parse the most common topology back into a tree object
  most_common_tree <- read.tree(text = most_common_topology)
  
  # Extract the tip labels from the most common topology
  most_common_tips <- most_common_tree$tip.label
  
  # Filter trees containing exactly the tips in the most common topology
filter_vector <- unlist(lapply(trees, function(x) 
  all(x$tip.label %in% most_common_tips) && all(most_common_tips %in% x$tip.label)
))

# Use `filter_vector` to filter trees
filtered_trees <- trees[filter_vector]

  # If filtered_trees is empty, skip
  if (length(filtered_trees) == 0) {
    return(NULL)
  }
  
  filtered_topology_strings <- sapply(filtered_trees, write.tree)

  # Create the filtered data frame
  filtered_topologies <- data.frame(
    topology = filtered_topology_strings,
    pvaginChr = sub(".*pvagin_([A-Za-z0-9]{5}).*", "\\1", filtered_topology_strings),
    pvaginStart=topologies_chr$pvaginStart[filter_vector]
  )


  return(filtered_topologies) # Return the filtered data frame
})

# Combine all filtered data frames into a single unified data frame
too <- rbindlist(chromosome_trees, fill = TRUE)

#### aaah somethign went wrong with chromsome_trees there are way too many here!?!?!


 too %>% group_by(topology, pvaginChr)%>%summarize(n=n())%>%arrange(-n)


## yikes this is a messy string, but gets at proportion of topologies with paired clades...
 too %>%group_by(pvaginChr)%>%mutate(totaltrees=n())%>%ungroup()%>% group_by(topology, pvaginChr, totaltrees)%>%filter(grepl('\\)\\,\\(', topology))%>%summarize(n=n())%>%ungroup()%>%group_by(pvaginChr)%>%mutate(total=sum(n))%>%ungroup()%>%group_by(topology,pvaginChr)%>%mutate(prop=n/total, proptotal=n/totaltrees)%>%arrange(pvaginChr,-n)%>%data.frame


## okay, just plot!

stoo= too %>%group_by(pvaginChr)%>%mutate(totaltrees=n())%>% group_by(topology, pvaginChr, totaltrees)%>%summarize(n=n())

pdf(paste0('~/transfer/',species,'_topologies.pdf'),16,10)
for(i in chromosomes){
alongchr <- ggplot(too[too$pvaginChr == i, ], aes(x = pvaginStart, y = topology, color = topology)) +
  geom_point(alpha = 0.1) +
  theme(legend.position = 'NULL') +
  labs(x = "pvaginStart", y = "Topology", title=i)

# Bar plot (fix axes and alignment for combining with scatter plot)
bp <- ggplot(stoo[stoo$pvaginChr == i, ], aes(x = n/totaltrees, y = factor(topology), fill = factor(topology))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Proportion", y = "Topology", title = i) +
  theme_minimal() +
  guides(fill = guide_legend(title = "Topology")) +
  theme(legend.position = 'NULL', 
        axis.text.y = element_blank(), # To remove redundant y-labels
        axis.title.y = element_blank())

# Combine plots side by side ensuring alignment
combined_plot <- plot_grid(alongchr, bp, nrow = 1, rel_widths = c(2, 1), align = 'h', axis = 'b')

# Display final combined plot
print(combined_plot)
}
dev.off()

## deconvolute to pairs in trees?

# Step 1: Create a function to extract pairs
extract_pairs <- function(topology) {
  # Clean the topology string by removing text outside the relevant tree structure
  topology_clean <- gsub("pvagin_Chr01.*", "", topology)  # Remove ending text
  topology_clean <- gsub(";|\\s+", "", topology)          # Remove semicolons and spaces
  
  # Match valid leaf pairs anywhere in the string
  matches <- str_extract_all(topology_clean, paste0(species,paste0("_chr\\d+_\\d+,",species,"_chr\\d+_\\d+")))[[1]]

  # Return the matches or NULL if none are found
  if (length(matches) > 0) {
    return(matches)
  } else {
    return(NULL)
  }
}


# Step 2: Apply the function and expand rows
too_expanded <- too %>%                  # Filter for Chr01
  mutate(pairs = map(topology, extract_pairs)) %>%  # Apply the function to extract pairs
  unnest(pairs) %>%                                 # Expand rows for multiple pairs
  mutate(pairs = str_remove_all(pairs, "[()]")) %>%    # Clean out parentheses in pairs
  mutate(
    simplified_pair = str_replace_all(pairs, ".*_(\\d+),.*_(\\d+)", "\\1,\\2")  # Extract final digits
  )

## this may not work for hexaploids?? or need to think ahrder
too_expanded$sortedpairs=NA
too_expanded$sortedpairs[too_expanded$simplified_pair%in%c('1,2','3,4')]='1,2-3,4'
too_expanded$sortedpairs[too_expanded$simplified_pair%in%c('1,3','2,4')]='1,3-2,4'
too_expanded$sortedpairs[too_expanded$simplified_pair%in%c('1,4','2,3')]='1,4-2,3'


## switch to simplified_pair for now
#stoo_expanded= too_expanded %>%group_by(pvaginChr)%>%mutate(totaltrees=n())%>% group_by(sortedpairs, pvaginChr, totaltrees)%>%summarize(n=n())
stoo_expanded= too_expanded %>%group_by(pvaginChr)%>%mutate(totaltrees=n())%>% group_by(simplified_pair, pvaginChr, totaltrees)%>%summarize(n=n())

pdf(paste0('~/transfer/',species,'_pairs.pdf'),16,10)
for(i in chromosomes){
alongchr <- ggplot(too_expanded[too_expanded$pvaginChr == i, ], aes(x = pvaginStart, y = simplified_pair, color = simplified_pair)) +
  geom_point(alpha = 0.1) +
  theme(legend.position = 'NULL') +
  labs(x = "pvaginStart", y = "Topology", title=i)

# Bar plot (fix axes and alignment for combining with scatter plot)
bp <- ggplot(stoo_expanded[stoo_expanded$pvaginChr == i, ], aes(x = n/totaltrees, y = factor(simplified_pair), fill = factor(simplified_pair))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Proportion", y = "Topology", title = i) +
  theme_minimal() +
  guides(fill = guide_legend(title = "Topology")) +
  theme(legend.position = 'NULL', 
        axis.text.y = element_blank(), # To remove redundant y-labels
        axis.title.y = element_blank())

# Combine plots side by side ensuring alignment
combined_plot <- plot_grid(alongchr, bp, nrow = 1, rel_widths = c(2, 1), align = 'h', axis = 'b')

# Display final combined plot
print(combined_plot)
}
dev.off()




########## all below doesn't work!





topologies%>%group_by(pvaginChr)%>%summarize(sptree=consensus(topology, p=0.5, check.labels=T))

  
  
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
