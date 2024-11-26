library(ggtree)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ape)
library(stringr)



at=read.tree('../sp_tree/paspalum_anchors_aster.2024-11-22.ASTRALPRO3OUT.u3.tre')


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



alt=read.table('../sp_tree/freqQuad.csv', comment.char='') ## tab not csv :(

alt$left=str_split_fixed(alt$V3, '#',2)[,1]
alt$right=str_split_fixed(alt$V3, '#',2)[,2]

alt$A=str_split_fixed(alt$left, '\\|', 2)[,1]
alt$B=str_split_fixed(alt$left, '\\|', 2)[,2]
alt$C=str_split_fixed(alt$right, '\\|', 2)[,1]
alt$D=str_split_fixed(alt$right, '\\|', 2)[,2]


## alt$V1 node number matches node_data node number






quar=data.frame(xstart=c(0.25,0.25,0.4,0.6,0.6), xend=c(0.4,0.4,0.6,0.75,0.75), ystart=c(0.4,0.6,0.5,0.5,0.5), yend=c(0.5,0.5,0.5,0.6,0.4))


first=ggplot(alt[1,]) + geom_text(aes(label=gsub(',','\n',A), x=0, y=1))+geom_text(aes(label=gsub(',','\n',B), x=0,y=0))+geom_text(aes(label=gsub(',','\n',C), x=1,y=1))+geom_text(aes(label=gsub(',','\n',D), x=1,y=0)) + geom_segment(data=quar, inherit.aes=F,aes(x=xstart, xend=xend, y=ystart,yend=yend)) + geom_text(aes(label=V5), x=0.5,y=0.55)+ geom_text(aes(label=round(V5/V6*100, digits=2)), x=0.5,y=0.45)
second=ggplot(alt[2,]) + geom_text(aes(label=gsub(',','\n',A), x=0, y=1))+geom_text(aes(label=gsub(',','\n',B), x=0,y=0))+geom_text(aes(label=gsub(',','\n',C), x=1,y=1))+geom_text(aes(label=gsub(',','\n',D), x=1,y=0)) + geom_segment(data=quar, inherit.aes=F,aes(x=xstart, xend=xend, y=ystart,yend=yend)) + geom_text(aes(label=V5), x=0.5,y=0.55)+ geom_text(aes(label=round(V5/V6*100, digits=2)), x=0.5,y=0.45)
third=ggplot(alt[3,]) + geom_text(aes(label=gsub(',','\n',A), x=0, y=1))+geom_text(aes(label=gsub(',','\n',B), x=0,y=0))+geom_text(aes(label=gsub(',','\n',C), x=1,y=1))+geom_text(aes(label=gsub(',','\n',D), x=1,y=0)) + geom_segment(data=quar, inherit.aes=F,aes(x=xstart, xend=xend, y=ystart,yend=yend)) + geom_text(aes(label=V5), x=0.5,y=0.55)+ geom_text(aes(label=round(V5/V6*100, digits=2)), x=0.5,y=0.45)
plot_grid(first,second,third, ncol=3)



generate_plot <- function(row) {
  ggplot() +
    geom_text(aes(label = gsub(",", "\n", row$A), x = 0, y = 1)) +
    geom_text(aes(label = gsub(",", "\n", row$B), x = 0, y = 0)) +
    geom_text(aes(label = gsub(",", "\n", row$C), x = 1, y = 1)) +
    geom_text(aes(label = gsub(",", "\n", row$D), x = 1, y = 0)) +
    geom_segment(data = quar, inherit.aes = FALSE,
                 aes(x = xstart, xend = xend, y = ystart, yend = yend)) +
    geom_text(aes(label = row$V5), x = 0.5, y = 0.55) +
    geom_text(aes(label = round(row$V5 / row$V6 * 100, digits = 2)), x = 0.5, y = 0.45)+
    xlim(-0.5,1.5)+ylim(-0.5,1.5)
}

# Generate plots for every 3 rows
plot_list <- lapply(seq(1, nrow(alt), by = 3), function(i) {
  rows <- alt[i:(i+2), , drop = FALSE] # Get the group of three rows
  plot_list <- lapply(seq_len(nrow(rows)), function(j) generate_plot(rows[j, ]))
  do.call(plot_grid, c(plot_list, ncol = 3))
})

pdf('../sp_tree/alternative_quartets.pdf',14,9)
for(i in plot_list){
print(i)}
dev.off()

# Combine all plots into a single plot
final_plot <- plot_grid(plotlist = plot_list, ncol = 1)
print(final_plot)


