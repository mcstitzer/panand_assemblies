library(ggtree)  ## /programs/R-4.2.1-r9/bin/R
library(ggplot2) ## on cbsu
library(cowplot)
theme_set(theme_cowplot())
library(ape)
library(treeio)
library(phytools)
library(reshape2)
library(dplyr)
library(ggridges)
library(tidypaleo) ## facet species names in italics!!!!!
library(ggh4x) ## facet strips spanning groups (subtribe)

## need to be in `scinet_trees/` ###AGGagfag
#filenames=list.files('.', pattern='RAxML_bipartitionsBranchLabels.')
dd=read.tree('../sp_tree/paspalum_anchors_aster.dated.2024-11-25.tre')                                              
 outdf=data.frame(pvgene=sapply(dd, function(tree){str_split_fixed(tree$tip.label, '_', 3)[1,2]}) )                                     
#outdf=data.frame(filenames=filenames)
ddd=lapply(dd, function(tree){
 nt=tree
 nt$tip.label=substr(tree$tip.label,1,6)
 return(nt)
})
ddd=do.call('c', ddd)




all=read.table('../panand_sp_ploidy.txt')                                
#for(sp in c('vcuspi', 'rtuber', 'blagur', 'hcompr', 'udigit', 'telega')){

## do southerns because they're in the gene tree
tripsacinae=c('tdacn2', 'tdacs2',"tdacn1", "tdacs1", "zdgigi", "zdmomo", "zluxur", "zmhuet", "zTIL18", "zTIL25", "zTIL01", "zTIL11", "znicar", "zmB735")
  
for(sp in all$V2){

  outdf[,paste0(sp, 'Count')]=NA
  outdf[,paste0(sp, 'Monophyletic')]=NA
#brlenlist=list()
# 
# for(i in (1:length(filenames))){ ## ep2 is not there
# 
# awt=read.raxml(paste0('',filenames[i]))
# #tryCatch(("reroot" (
# #awt=root(awt, awt$tip.label[substr(awt$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate')]))), error = function(err) print("Outgroups already at the root"))
# awt=as.phylo(awt)
# awt$tip.label=gsub('_R_', '', awt$tip.label)
#   ## add to make sure outgroup is there, and is monophyletic - otherwise skip this tree
#   if(any(substr(as.phylo(awt)$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate'))){
#     if(is.monophyletic(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate')])){
# awt=root(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate')])
# 
# vc=sum(substr(as.phylo(awt)$tip.label,1,6) %in% sp)
# outdf[i, paste0(sp, 'Count')]=vc
# if(vc>1){
# mrcanode=getMRCA(as.phylo(awt), as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c(sp)])
# getDescendants(as.phylo(awt), mrcanode)
# tips=as.phylo(awt)$tip.label[getDescendants(as.phylo(awt), mrcanode)]
# tips=tips[!is.na(tips)]
# outdf[i, paste0(sp, 'Monophyletic')]=all(substr(tips,1,6)==sp)
#   ## get terminal branch length for each copy
#   ## first get the node numbers of the tips
# nodes<-sapply(as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c(sp)],function(x,y) which(y==x),y=as.phylo(awt)$tip.label)
# ## then get the edge lengths for those nodes
# edge.lengths<-setNames(awt$edge.length[sapply(nodes,function(x,y) which(y==x),y=awt$edge[,2])],names(nodes))
# #brlenlist[i]=edge.lengths
# }}
# }
# }
#                                               }
               
   for(i in seq_along(ddd)){  ## loop over trees
  awt = dd[[i]]
  awt = as.phylo(awt)
  awt$tip.label = gsub('_R_', '', awt$tip.label)
  
  ## if sp is in tripsacinae, drop other tripsacinae tips
  if(sp %in% tripsacinae){
    awt = drop.tip(awt, awt$tip.label[substr(awt$tip.label,1,6) %in% tripsacinae[-which(tripsacinae == sp)]])
  }
  
  ## count number of focal tips
  vc = sum(substr(awt$tip.label, 1, 6) == sp)
  outdf[i, paste0(sp, 'Count')] = vc
  
  if(vc > 0){
    ## Define the focal node: if only one tip, use that tip index;
    ## if more than one, get the MRCA of all focal tips.
    if(vc == 1){
      focal_node = which(substr(awt$tip.label, 1, 6) == sp)[1]
    } else {
      focal_node = getMRCA(awt, awt$tip.label[substr(awt$tip.label, 1, 6) == sp])
    }
    
    ## Find the parent of the focal node
    parent_ind = which(awt$edge[,2] == focal_node)
    if(length(parent_ind) == 0){
      ## Focal node is the root; no sister exists.
      sister_tips = NA
    } else {
      parent_node = awt$edge[parent_ind, 1]
      ## Get the children of the parent
      children = awt$edge[awt$edge[,1] == parent_node, 2]
      ## The sister node is the one that is not the focal node
      sister_node = children[children != focal_node]
      
      if(length(sister_node) == 0){
        sister_tips = NA
      } else {
        ## If the sister node is a tip:
        if(sister_node <= length(awt$tip.label)){
          sister_tips = awt$tip.label[sister_node]
        } else {
          ## Otherwise, extract the clade defined by the sister node
          sister_clade = extract.clade(awt, node = sister_node)
          sister_tips = sister_clade$tip.label
        }
      }
    }
    
    ## Record the sister tips (as a comma-separated string)
    outdf[i, paste0(sp, 'SisterTips')] = paste(sister_tips, collapse = ",")
  }
}

}

#outdf now has columns with sisters - but this is for the mrca of all species tips! 
## i want to get for each tip to do stuff with



# Assuming dd is your list of trees read via read.tree()
results_list <- list()
tree_counter <- 1

for(tree in dd) {
  # Convert to phylo object if needed.
  current_tree <- as.phylo(tree)
  
  # Optionally, you might want to remove specific patterns,
  # but we're not truncating the tip names.
  # current_tree$tip.label <- gsub('_R_', '', current_tree$tip.label)
  
  ntips <- length(current_tree$tip.label)
  
  # Loop over every tip in the current tree.
  for(tip_idx in 1:ntips) {
    tip_label <- current_tree$tip.label[tip_idx]  # full tip name
    
    # Find the edge(s) where the child equals this tip.
    parent_edge_idx <- which(current_tree$edge[,2] == tip_idx)
    
    # If the tip has no parent (should not happen for a tip), assign NA.
    if(length(parent_edge_idx) == 0) {
      sister_str <- NA
    } else {
      parent_node <- current_tree$edge[parent_edge_idx, 1]
      
      # Get all children of the parent node.
      children <- current_tree$edge[current_tree$edge[,1] == parent_node, 2]
      
      # Exclude the focal tip to get the sister node(s).
      sister_nodes <- children[children != tip_idx]
      
      sister_tips <- c()
      # For each sister node, check if it is a tip or an internal node.
      for(node in sister_nodes) {
        if(node <= ntips) {
          # The node is a tip.
          sister_tips <- c(sister_tips, current_tree$tip.label[node])
        } else {
          # The node is internal: extract its clade to get all descendant tip labels.
          clade <- extract.clade(current_tree, node)
          sister_tips <- c(sister_tips, clade$tip.label)
        }
      }
      
      # Create a comma-separated string of unique sister tip labels.
      sister_str <- paste(unique(sister_tips), collapse = ",")
    }
    
    # Record the result for the current tip.
    results_list[[length(results_list) + 1]] <- data.frame(
      tree = tree_counter,
      tip = tip_label,
      sisters = sister_str,
      stringsAsFactors = FALSE
    )
  }
  
  tree_counter <- tree_counter + 1
}


# Combine all results into one data frame.
tree_tips <- do.call(rbind, results_list)

head(tree_tips)



# For the focal tip, extract the first six characters.
tree_tips$tip_species <- substr(tree_tips$tip, 1, 6)

# For the sisters, split the comma-separated list, extract the first six characters,
# remove duplicates, and then collapse back to a comma-separated string.
tree_tips$sister_species <- sapply(tree_tips$sisters, function(s) {
  if (is.na(s) || s == "") return(NA)
  # Split into individual tip names, trim whitespace, and extract first six characters.
  sp_names <- trimws(unlist(strsplit(s, ",")))
  sp_names <- substr(sp_names, 1, 6)
  # Keep only unique species names.
  sp_names <- sort(unique(sp_names))
  # Collapse back into a comma-separated string.
  paste(sp_names, collapse = ",")
})


# Define a helper function to parse a single tip string.
parse_tip <- function(tip) {
  # Split by underscore.
  parts <- unlist(strsplit(tip, "_"))
  n <- length(parts)
  
  # The last part should be the start-end string.
  pos <- parts[n]
  pos_parts <- unlist(strsplit(pos, "-"))
  
  # Get start and end positions.
  tip_start <- pos_parts[1]
  tip_end <- pos_parts[2]
  
  # Anchor is always the second token.
  tip_anchor <- parts[2]
  
  # The contig is made from tokens 3 to (n-1).
  # This will work whether there's one token or more.
  tip_contig <- paste(parts[3:(n-1)], collapse = "_")
  
  # Return a data frame row.
  return(data.frame(tip_anchor = tip_anchor,
                    tip_contig = tip_contig,
                    tip_start = tip_start,
                    tip_end = tip_end,
                    stringsAsFactors = FALSE))
}

# Example: suppose output_df is your current data frame with a column "tip"
# (which contains the full tip names).

# Apply the function to every tip and combine the results.
parsed_tips <- do.call(rbind, lapply(tree_tips$tip, parse_tip))

# Combine the new columns with your original data frame.
tree_tips <- cbind(tree_tips, parsed_tips)

# Check the output:
head(tree_tips)



# Now define a function to parse a comma-separated list of sister tips.
parse_sisters <- function(s) {
  # If s is NA or empty, return NA for all components.
  if (is.na(s) || s == "") {
    return(c(anchor = NA, contig = NA, start = NA, end = NA))
  }
  
  # Split the string by comma and trim any whitespace.
  tip_list <- trimws(unlist(strsplit(s, ",")))
  
  # Parse each sister tip.
  parsed_list <- lapply(tip_list, parse_tip)
  
  # Extract each component and take unique values.
  anchors <- unique(sapply(parsed_list, function(x) x$tip_anchor))
  contigs <- unique(sapply(parsed_list, function(x) x$tip_contig))
  starts  <- unique(sapply(parsed_list, function(x) x$tip_start))
  ends    <- unique(sapply(parsed_list, function(x) x$tip_end))
  
  # Collapse unique values back to comma separated strings.
  out <- c(anchor = paste(anchors, collapse = ","),
           contig = paste(contigs, collapse = ","),
           start  = paste(starts, collapse = ","),
           end    = paste(ends, collapse = ","))
  return(out)
}

# Now apply parse_sisters() to the sisters column in your output data frame.
# (Assuming your data frame is named output_df and has a column named 'sisters'.)
sister_parsed <- t(sapply(tree_tips$sisters, parse_sisters))

# Bind the new columns (sister_anchor, sister_contig, sister_start, sister_end) to your output_df.
tree_tips$sister_anchor <- sister_parsed[, "anchor"]
tree_tips$sister_contig <- sister_parsed[, "contig"]
tree_tips$sister_start  <- sister_parsed[, "start"]
tree_tips$sister_end    <- sister_parsed[, "end"]

# View the updated data frame.
head(tree_tips)

tree_tips$pvchr=substr(tree_tips$tip_anchor,1,7)

ks2p=fread('../syntenic_anchors/ksp_to_pasp.txt.gz')
tree_tips$ksp=ks2p$V17[match(tree_tips$tip, ks2p$V3)]
tree_tips$dndsp=ks2p$V19[match(tree_tips$tip, ks2p$V3)]





tree_tips %>% group_by(tip_species, sister_species, tip_contig, sister_contig, pvchr=substr(tip_anchor,1,7)) %>% summarize(n=n()) %>% arrange(-n)%>% filter(tip_species=='blagur' & pvchr=='Pavag01') %>% print(n=30)

## try with blagur
filtered_pairs <- tree_tips %>%
    group_by(tip_species, sister_species, tip_contig, sister_contig, 
             pvchr = substr(tip_anchor, 1, 7)) %>%
    summarize(n = n(), .groups = "drop") %>%
    arrange(-n) %>%
    filter(tip_species == 'blagur') %>%   # or remove this filter to do all pvchr
    filter(n > 30)


clusters_by_pvchr <- filtered_pairs %>%
    group_by(pvchr, tip_species) %>%  # cluster separately for each pvchr and species
    group_modify(~ {
        group_data <- .x
        
        # Build an edge list from the tip_contig and sister_contig columns.
        edge_list <- group_data %>%
            select(tip_contig, sister_contig) %>%
            distinct()
        
        # Build an undirected graph from the edge list.
        g <- graph_from_data_frame(edge_list, directed = FALSE)
        comps <- components(g)
        
        # Create a lookup table mapping each contig (node) to its cluster.
        contig_clusters <- tibble(
            contig  = names(comps$membership),
            cluster = comps$membership
        )
        
        # Compute the support for each cluster.
        # Here we join the original group_data with contig_clusters by tip_contig.
        # (Since tip_species == sister_species, both nodes belong to the same cluster.)
        support_df <- group_data %>%
            left_join(contig_clusters, by = c("tip_contig" = "contig")) %>%
            group_by(cluster) %>%
            summarize(support = sum(n), .groups = "drop")
        
        # Build a lookup of clusters: for each cluster, get a comma-separated list of contigs.
        cluster_lookup <- contig_clusters %>%
            group_by(cluster) %>%
            summarize(contigs = paste(unique(contig), collapse = ","), .groups = "drop")
        
        # Merge the cluster membership with the support counts.
        clusters <- cluster_lookup %>%
            left_join(support_df, by = "cluster") %>%
            mutate(pvchr = unique(group_data$pvchr),
                   tip_species = unique(group_data$tip_species))
        
        clusters
    }) %>%
    ungroup()

# View the final result.
print(clusters_by_pvchr)



#### this gets me close to "alleels" for each!!!!!!!



