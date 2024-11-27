library(ggtree)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ape)




at=read.tree('../sp_tree/paspalum_anchors_aster.2024-11-22.ASTRALPRO3OUT.support.tre')


node_data <- data.frame(
  node = c(length(at$tip.label)+1):c(length(at$tip.label)+length(at$node.label)),
  label = at$node.label
) %>%
  mutate(label = gsub("'", "", label)) %>%  # Remove single quotes
  mutate(label = gsub("\\]$", "", label)) %>%  # Remove trailing ] from q3
  separate_rows(label, sep = ";") %>%  # Split key-value pairs
  separate(label, into = c("key", "value"), sep = "=", fill = "right") %>%  # Split into key and value
  mutate(value = as.numeric(value)) %>%  # Convert to numeric (will drop non-numeric)
  filter(key %in% c("q1", "q2", "q3")) %>%  # Keep only q1, q2, q3
  pivot_wider(names_from = key, values_from = value)  # Reshape into wide format




# Prepare data for nodepie
pie_data <- node_data %>%
  select(node, q1, q2, q3) %>%
  pivot_longer(cols = c(q1, q2, q3), names_to = "group", values_to = "value")

# Create pie chart plots for each node
 pies=nodepie(node_data, cols=c('q1', 'q2', 'q3'), color=c('#2c7fb8', '#7fcdbb', '#edf8b1'),  alpha=0.95)
 
 
# Plot the tree with pies
p <- ggtree(at) +
geom_tiplab()+
#  geom_text2(aes(label = node), subset = .(isTip == FALSE)) +
  geom_inset(pies, x = "node", height=0.05, width=0.05)

print(p)


#chronogram <- chronos(at, lambda = 0.1, calibration = list(node = 39, age.min = 14.6, age.max = 27.6))
noog=drop.tip(at, c('pvagin', 'tdacn2', 'tdacs2'))

chronogram <- chronos(noog, lambda = 0.005, calibration = list(node = length(noog$tip.label)+1, age.min = 14.6, age.max = 27.6))


node_data <- data.frame(
  node = c(length(noog$tip.label)+1):c(length(noog$tip.label)+length(noog$node.label)),
  label = noog$node.label
) %>%
  mutate(label = gsub("'", "", label)) %>%  # Remove single quotes
  mutate(label = gsub("\\]$", "", label)) %>%  # Remove trailing ] from q3
  separate_rows(label, sep = ";") %>%  # Split key-value pairs
  separate(label, into = c("key", "value"), sep = "=", fill = "right") %>%  # Split into key and value
  mutate(value = as.numeric(value)) %>%  # Convert to numeric (will drop non-numeric)
  filter(key %in% c("q1", "q2", "q3")) %>%  # Keep only q1, q2, q3
  pivot_wider(names_from = key, values_from = value)  # Reshape into wide format




# Prepare data for nodepie
pie_data <- node_data %>%
  select(node, q1, q2, q3) %>%
  pivot_longer(cols = c(q1, q2, q3), names_to = "group", values_to = "value")

# Create pie chart plots for each node
 pies=nodepie(node_data, cols=c('q1', 'q2', 'q3'), color=c('#2c7fb8', '#7fcdbb', '#edf8b1'),  alpha=0.95)
 

chronopie <- ggtree(chronogram) +
geom_tiplab()+
#  geom_text2(aes(label = node), subset = .(isTip == FALSE)) +
  geom_inset(pies, x = "node", height=0.1, width=0.1)+
  theme_tree2()+ scale_x_continuous(expand = expansion(mult = 1.5), labels = abs)
chronopie


chronopienotip <- ggtree(chronogram) +
#geom_tiplab()+
#  geom_text2(aes(label = node), subset = .(isTip == FALSE)) +
  geom_inset(pies, x = "node", height=0.15, width=0.15)+
  theme_tree2()+ scale_x_continuous( labels = abs)


## maek a legend
 pie_labels=c('Main Topology', 'First Alternative', 'Second Alternative')
 pie_colors=c('#2c7fb8', '#7fcdbb', '#edf8b1')
dummy_data <- data.frame(
  group = pie_labels,
  value = 1:3
)
chronopienotip <- ggtree(chronogram) +
 # geom_tiplab() +
 geom_point(data = dummy_data, aes(x = -Inf, y = -Inf, fill = group), 
             shape = 21, size = 0.000001, color = "black") +  # Use shape 21 for fill support
  scale_fill_manual(
    values = pie_colors,
    labels = pie_labels,
    name = ""
  ) +
  theme_tree2()+ scale_x_continuous( labels = abs)+
  theme(
    legend.position = c(0.05, 1),  # Top-left corner (x, y)
    legend.justification = c(0, 1),  # Align legend by top-left corner
    legend.title = element_text(size = 1),
    legend.text = element_text(size = 7)
  )+
    geom_inset(pies, x = "node", height=0.15, width=0.15)+
    guides(
    fill = guide_legend(override.aes = list(size = 5))  # Make legend points larger
  )



ddd=lapply(dd, function(tree){
 nt=tree
 nt$tip.label=substr(tree$tip.label,1,6)
 return(nt)
})
ddd=do.call('c', ddd)


ggdensitree(trees.fort, layout="rectangular",tip.order=names(rev(taxonnames)), align.tips=T, aes(color=tree, alpha=alpha)) + geom_tiplab(cex=1) + scale_color_manual(values=c('black', 'snow4'))

ggdensitree(dd, alpha=.3, colour='steelblue') + 
    geom_tiplab(size=3) + hexpand(.35)
