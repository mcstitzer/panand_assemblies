library(ape)
library(ggtree)
library(phytools)
library(tidypaleo) ## facet species names in italics!!!!!
library(ggh4x) ## facet strips spanning groups (subtribe)
library(treeio)


taxonnames=c("Zea mays ssp. parviglumis TIL11", "Zea mays ssp. mays B73v5", "Zea mays ssp. parviglumis TIL01", "Zea mays ssp. mexicana TIL25", "Zea mays ssp. mexicana TIL18", "Zea mays ssp. huehuetengensis", 
             "Zea luxurians", "Zea nicaraguensis", "Zea diploperennis Momo", "Zea diploperennis Gigi", "Tripsacum zoloptense", "Tripsacum dactyloides FL", "Tripsacum dactyloides Southern Hap2", 
             "Tripsacum dactyloides Northern Hap2", "Tripsacum dactyloides KS", "Tripsacum dactyloides tetraploid", "Urelytrum digitatum", "Vossia cuspidata", "Rhytachne rottboellioides", "Rottboellia tuberculosa", 
             "Hemarthria compressa", "Elionurus tripsacoides", "Schizachyrium scoparium", "Schizachyrium microstachyum", "Anatherum virginicum", "Andropogon chinensis", "Andropogon gerardi", 
             "Cymbopogon refractus", "Cymbopogon citratus", "Heteropogon contortus", "Themeda triandra", "Bothriochloa laguroides", "Pogonatherum paniceum", "Sorghum bicolor", 
             "Ischaemum rugosum", "Sorghastrum nutans", '"Andropogon" burmanicus', "Thelepogon elegans", "Chrysopogon serrulatus", "Paspalum vaginatum")
names(taxonnames)=c("zTIL11", "zmB735", "zTIL01", "zTIL25", "zTIL18", "zmhuet", 
                    "zluxur", "znicar", "zdmomo", "zdgigi", "tzopol", "tdacs1", "tdacs2", 
                    "tdacn2", "tdacn1", "tdactm", "udigit", "vcuspi", "rrottb", "rtuber", 
                    "hcompr", "etrips", "sscopa", "smicro", "avirgi", "achine", "agerar", 
                    "crefra", "ccitra", "hconto", "ttrian", "blagur", "ppanic", "sbicol", 
                    "irugos", "snutan", "atenui", "telega", "cserru", "pvagin")


#rsp=read.tree(text='((((((((((((zTIL11:1.0,zmB735:1.0):0.00399,zTIL01:1.0):0.123008,(zTIL25:1.0,zTIL18:1.0):0.051481):0.078959,zmhuet:1.0):0.315942,(zluxur:1.0,znicar:1.0):1.343):0.035313,(zdmomo:1.0,zdgigi:1.0):1.085725):4.966642,(((((tdacs1:1.0,tdacs2:1.0):0.104529,tdacn2:1.0):0.344325,tdacn1:1.0):0.085028,tdactm:1.0):3.804823,tzopol:1.0):0.163575):3.340689,((udigit:1.0,vcuspi:1.0):0.130867,rrottb:1.0):0.2067):0.532529,((rtuber:1.0,hcompr:1.0):1.047055,etrips:1.0):0.311621):0.382165,(((((((((((sscopa:1.0,smicro:1.0):0.375104,avirgi:1.0):0.464688,(achine:1.0,agerar:1.0):0.197343):2.589164,(crefra:1.0,ccitra:1.0):4.376438):0.274712,((hconto:1.0,ttrian:1.0):0.601483,blagur:1.0):0.896994):0.456651,ppanic:1.0):0.093559,sbicol:1.0):0.131326,irugos:1.0):0.006381,snutan:1.0):0.144602,atenui:1.0):0.524102,(telega:1.0,cserru:1.0):0.446963):0.337256):2.735692,((bdista:1.0,osativ:1.0):8.038372,svirid:1.0):0.322258):1.0,paspal:1.0);')
rsp=read.tree('paspalum_anchors_aster.2024-11-22.ASTRALOUT.tre')
rsp=root(rsp, outgroup='pvagin', resolve.root=T)
rsp=drop.tip(rsp, c('svirid', 'bdista', 'osativ', 'tdacs2', 'tdacn2', 'tdactm', 'tzopol', 'paspal'))
rsp=rotateConstr(rsp, rev(rsp$tip.label)) 
rsp$tip.label=taxonnames[rsp$tip.label]



## asize from other plot :(
#hgs=ggplot(asize, aes(x=(haploidAssemblySize-haploidNCount)/1e9, y=1, color=ploidy)) + geom_segment(aes(y=1,yend=1, x=0, xend=doubledAssembly/1e9/2)) + geom_vline(xintercept=c(2,4), color='snow2', linetype='dotted') + geom_vline(xintercept=c(1,3,5), color='snow3', linetype='dotted') + scale_color_manual(values=ploidycolors) + scale_fill_manual(values=ploidycolors) + geom_point(aes(x=haploidAssemblySize/1e9), color='snow3', size=2)+ geom_point(size=4)+ geom_point(aes(x=haploidRepeatSize/1e9, bg=ploidy, y=1.6), shape=25, size=3)+ facet_wrap(~speciesLabel, ncol=1, strip.position='left', labeller=purrr::partial(label_species, dont_italicize=c('subsp.', 'ssp.', 'TIL11', 'TIL01', 'TIL25', 'TIL18', 'Momo', 'Gigi', 'Southern Hap1', 'Northern Hap1', 'FL', 'KS',  '\\*', '\\"', 'B73v5'))) + theme(strip.placement = "outside",   strip.text.y.left = element_text(angle=0), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9), background_x = elem_list_rect(fill = ploidycolors[asize$ploidy]))+ theme(legend.position = "none") + ylab('')+ xlab('Haploid Size (Gb)') + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_text(color=ploidycolors[asize$ploidy])) + ylim(0,2)
#plot_grid(densit, NULL, hgs, align='hv',axis='tb', ncol=3, rel_widths=c(0.2,-0.05,0.7))

densit=ggtree(force.ultrametric(rsp, method='extend'))

data.frame(label=rsp$tip.label, color=asize$ploidy[match(names(rsp$tip.label), asize$V2)])

densit + geom_tiplab(aes(color=ploidycolors[asize$ploidy[match(rsp$tip.label, asize$V2)]])) + coord_cartesian(clip="off") +  scale_x_continuous(expand = expansion(mult = 14.9))
densit + geom_tiplab(fontface='italic') + coord_cartesian(clip="off") +  scale_x_continuous(expand = expansion(mult = 3))

  ## from bianconi, 21.1 (14.6–27.6) for macrofossils
  ## 34.0 (23.5–44.4) for microfossils
rsp$edge.length[1]=1e-8
## node nubmer is tips (36)+1
chronogram <- chronos(force.ultrametric(rsp, method = "extend"), lambda = 0.1, calibration = list(node = 37, age.min = 14.6, age.max = 27.6))

## it said i needed force ultrametric - i don't seem to need to on the gene trees???
### is this messing anything up???
chronogram <- chronos(force.ultrametric(rsp, method = "extend"), lambda = 0.01, calibration = list(node = 37, age.min = 12.2, age.max = 23.7))

p=ggtree(chronogram) + geom_tiplab(fontface='italic') + coord_cartesian(clip="off") +  scale_x_continuous(expand = expansion(mult = 3), labels = abs) + theme_tree2()

revts(p)



### time calibrate each gene tree
trees=read.tree('paspalum_anchors_aster.2024-11-22.tre')
dated <- vector("list", length(trees))  # Initialize list for dated trees

filenames=list.files('../gene_trees/trees/', pattern='RAxML_bipartitionsBranchLabels.*0')

#all=read.table('../panand_sp_ploidy.txt')                                


b=lapply(c(1:length(filenames))[-40], function(i){ ## ep2 is not there

print(i)
awto=read.raxml(paste0('../gene_trees/trees/',filenames[i]))
awt=as.phylo(awto)
awt$node.label <- awto@data$bootstrap
awt$tip.label=gsub('_R_', '', awt$tip.label)
  ## add to make sure outgroup is there, and is monophyletic - otherwise skip this tree
  if(any(substr(as.phylo(awt)$tip.label,1,6) %in% c('pvagin'))){
#    if(is.monophyletic(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('pvagin')])){
  awt=root(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('pvagin')], resolve.root=T)

## now keep only six digit code
####awt$tip.label=substr(awt$tip.label,1,6)

return(awt)
}
#}
                                              })
trees=do.call("c",b)



# Loop through each tree
chronograms=lapply(1:length(trees), function(tree) {
  # Root the tree
 # print(tree)
  #rsp <- root(trees[[tree]], outgroup = 'pvagin', resolve.root = TRUE)
  rsp=trees[[tree]]
  # Try different lambda values, fallback to NA if both fail
  chronogram <- tryCatch({
    # First attempt with lambda = 0.01
    chronos(rsp, lambda = 0.01, calibration = list(node = length(rsp$tip.label) + 1, age.min = 12.2, age.max = 23.7))
  }, error = function(e) {
    tryCatch({
      # Second attempt with lambda = 0.1
      chronos(rsp, lambda = 0.1, calibration = list(node = length(rsp$tip.label) + 1, age.min = 12.2, age.max = 23.7))
    }, error = function(e2) {
      # If it fails again, return NA
      NA
    })
  })
  
  # Store the result
  return(chronogram)
})

valid_trees <- chronograms[!vapply(chronograms, function(x) {
  # Ensure x is scalar and test its properties
  if (length(x) != 1) {
    TRUE  # Exclude non-scalar objects
  } else {
    is.na(x) || !inherits(x, "phylo")
  }
}, logical(1))]

# Check if any valid trees exist
if (length(valid_trees) > 0) {
  dd <- do.call(c, valid_trees)
  class(dd) <- "multiPhylo"
} else {
  stop("No valid phylogenetic trees available for export.")
}

# Check the results
dd=do.call("c",c(chronograms, na.rm=T))
class(dd)='multiPhylo'

### not working - need to fix missing trees, and i stopped it early

##rsp$tip.label=taxonnames[rsp$tip.label]

## some trees are wrong???
dd <- dd[!sapply(dd, function(tree) {
  !inherits(tree, "phylo") || is.null(tree$edge.length)
})]
class(dd) <- "multiPhylo"

write.tree(dd, paste0('paspalum_anchors_aster.dated.', Sys.Date(), '.tre'))

ggtree(dated[[2]]) + geom_tiplab(fontface='italic') + coord_cartesian(clip="off") +  scale_x_continuous(expand = expansion(mult = 1.5), labels = abs) + theme_tree2()

#revts(p)




# Initialize a vector to store MRCA-to-tip values
mrca_to_tip_values <- numeric(length(dd))

# Loop through all trees in dd
for (i in seq_along(dd)) {
  tree <- dd[[i]]  # Get the tree
  
  # Ensure the tree is valid
  if (inherits(tree, "phylo")) {
    # Calculate tree height
    tree_height <- max(node.depth.edgelength(tree))
    
    # Find the MRCA of tips matching a specific condition
    group_tips <- tree$tip.label[substr(tree$tip.label, 1, 6) == "ttrian"]
    
    # Check if there are enough tips to calculate an MRCA
    if (length(group_tips) > 1) {
      mrca_node <- getMRCA(tree, group_tips)
      
      # Calculate the MRCA depth
      mrca_depth <- node.depth.edgelength(tree)[mrca_node]
      
      # Calculate MRCA-to-tip length
      mrca_to_tip_values[i] <- tree_height - mrca_depth
    } else {
      mrca_to_tip_values[i] <- NA  # Not enough tips
    }
  } else {
    mrca_to_tip_values[i] <- NA  # Invalid tree
  }
}

summary(mrca_to_tip_values)

## oh dang i should be doing this with the whole tip name, so i can go back to subgenomes!!!!!


