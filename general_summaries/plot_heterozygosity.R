library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(stringr)




synt=read.table('~/Downloads/sharedSyntenicAnchors.txt', header=F)


b=read.table('../syntenic_anchors/gene_anchorTable_AnchorWave_Paspalum.txt', header=T)
b=b[,-which(colnames(b)%in%c('svirid',"pprate", "osativ", "eophiu", "bdista", "agerjg", 'tdactm', 'tzopol'))]
b=b[,-which(colnames(b)%in%c('tdacs2', 'tdacn2'))]



table(rowSums(b[,-1]>0))
bb=b[rowSums(b[,-1]>0)>=32,]

sum(synt$V1==str_split_fixed(bb$gene, '\\.',2)[,1])


## now get positions in the genome for thsese
all=read.table('../panand_sp_ploidy.txt')
all$V3[all$V2=='rtuber']=2
all$V3[all$V2=='telega']=2



#### now do for each 
merged_results <- list()

# Loop through each file, process it, and append it to the merged_results list
for (file in all$V2) {
  # Read in the ks file
  ks <- fread(paste0('../syntenic_anchors/ks_self/', file, '_matches.txt'))
  
  # Extract gene information
  ks$gene <- str_split_fixed(ks$V3, '_', 3)[, 2]
  
  # Filter based on 'synt' gene list
  kss <- ks[ks$gene %in% synt$V1, ]
  
  # Process queryChr
  kss$queryChr <- str_split_fixed(kss$V3, '_', 4)[, 3]
  kss$queryChr[substr(kss$queryChr, 1, 1) %in% c('a', 's')] <- str_extract(kss$V3, "(scaf_[0-9]+|alt-scaf_[0-9]+)")[substr(kss$queryChr, 1, 1) %in% c('a', 's')]
  kss$queryChr[substr(kss$queryChr, 1, 3) %in% c('ctg')] <- str_extract(kss$V3, "(ctg_[0-9]+)")[substr(kss$queryChr, 1, 3) %in% c('ctg')]
  
  # Extract queryStart and queryEnd
  kss$queryStart <- as.integer(str_match(kss$V3, "_(\\d+)-(\\d+)$")[, 2])
  kss$queryEnd <- as.integer(str_match(kss$V3, "_(\\d+)-(\\d+)$")[, 3])
  
  # Process reference information
  kss$refGenome <- substr(kss$V4, 1, 6)
  kss$refChr <- str_split_fixed(kss$V4, '_', 4)[, 3]
  kss$refChr[substr(kss$refChr, 1, 1) %in% c('a', 's')] <- str_extract(kss$V4, "(scaf_[0-9]+|alt-scaf_[0-9]+)")[substr(kss$refChr, 1, 1) %in% c('a', 's')]
  kss$refChr[substr(kss$refChr, 1, 3) %in% c('ctg')] <- str_extract(kss$V4, "(ctg_[0-9]+)")[substr(kss$refChr, 1, 3) %in% c('ctg')]
  
  # Extract refStart and refEnd
  kss$refStart <- as.integer(str_match(kss$V4, "_(\\d+)-(\\d+)$")[, 2])
  kss$refEnd <- as.integer(str_match(kss$V4, "_(\\d+)-(\\d+)$")[, 3])
  
  # Add paspChr column
  kss$paspChr <- substr(kss$gene, 6, 7)
  
  # Add a column for the genome from the file name (assuming 'blagur' is part of the file name)
  kss$genome <- file
  
  # Append the result to the merged_results list
  merged_results[[file]] <- kss
}

# Combine all processed data frames into one
final_merged_data <- rbindlist(merged_results)
final_merged_data$ploidy=asize$ploidy[match(final_merged_data$genome, asize$V2)]
final_merged_data$speciesLabel=asize$speciesLabel[match(final_merged_data$genome, asize$V2)]

final_merged_data %>% group_by(queryChr, refChr, paspChr, genome) %>%summarize(ks=median(V17), n=n()) %>% filter(n>10) %>% ggplot() + geom_histogram(aes(x=ks), binwidth=0.001) + facet_wrap(~paspChr, ncol=1) 

shorttaxonnames=c("Z. mays ssp. parviglumis TIL11", "Z. mays ssp. mays B73v5", "Z. mays ssp. parviglumis TIL01", "Z. mays ssp. mexicana TIL25", "Z. mays ssp. mexicana TIL18", "Z. mays ssp. huehuetengensis", 
                  "Z. luxurians", "Z. nicaraguensis", "Z. diploperennis Momo", "Z. diploperennis Gigi", "T. zoloptense", "T. dactyloides FL", "T. dactyloides Southern Hap2", 
                  "T. dactyloides Northern Hap2", "T. dactyloides KS", "T. dactyloides tetraploid", "U. digitatum", "V. cuspidata", "R. rottboellioides", "R. tuberculosa", 
                  "H. compressa", "E. tripsacoides", "S. scoparium", "S. microstachyum", "A. virginicum", "A. chinensis", "A. gerardi", 
                  "C. refractus", "C. citratus", "H. contortus", "T. triandra", "B. laguroides", "P. paniceum", "S. bicolor", 
                  "I. rugosum", "S. nutans", '"A." burmanicus', "T. elegans", "C. serrulatus", "P. vaginatum")
names(shorttaxonnames)=c("zTIL11", "zmB735", "zTIL01", "zTIL25", "zTIL18", "zmhuet", 
                         "zluxur", "znicar", "zdmomo", "zdgigi", "tzopol", "tdacs1", "tdacs2", 
                         "tdacn2", "tdacn1", "tdactm", "udigit", "vcuspi", "rrottb", "rtuber", 
                         "hcompr", "etrips", "sscopa", "smicro", "avirgi", "achine", "agerar", 
                         "crefra", "ccitra", "hconto", "ttrian", "blagur", "ppanic", "sbicol", 
                         "irugos", "snutan", "atenui", "telega", "cserru", "pvagin")


## from ks_plotting
final_merged_data$shortspeciesLabel=ifelse(final_merged_data$genome%in%asize$V2[asize$haploid], paste0(shorttaxonnames[match(final_merged_data$genome, names(shorttaxonnames))], '*'), shorttaxonnames[match(final_merged_data$genome, names(shorttaxonnames))])


### put in all the het/homeol results!!!! 
homeol=final_merged_data %>% filter(!is.na(shortspeciesLabel)) %>% group_by(queryChr, refChr, paspChr, genome,ploidy, speciesLabel, shortspeciesLabel) %>%mutate(ks=mean(V17), n=length(unique(gene))) %>% filter(n>10)   %>% group_by(genome, shortspeciesLabel, ploidy) %>% summarize(n=n(), ks=mean(V17, na.rm=T))

het=data.frame(genome=unique(final_merged_data$refGenome), shortspeciesLabel=unique(final_merged_data$shortspeciesLabel))
het$ks=homeol$ks[match(het$genome, homeol$genome)]
het$ploidy=asize$ploidy[match(het$genome, asize$V2)]




het$ks[is.na(het$ks)]=0
het$het=NA
het$het[het$genome=='avirgi']=2.2e-6
het$het[het$genome=='ppanic']=0.006
het$het[het$genome=='rrottb']=0.002
het$het[het$genome=='smicro']=0.004
het$het[het$genome=='udigit']=2.81e-6
het$het[het$genome=='achine']=0.004859949
het$het[het$genome=='tdacs1']=0.004768089
het$het[het$genome=='tdacn1']=0.01200872
het$het[het$genome=='rtuber']=1.202468e-06
het$het[het$genome=='irugos']=0.003000871
het$het[het$genome=='etrips']=0.007441668
het$het[het$genome=='zmhuet']=0.001369575
het$het[het$genome=='zTIL25']=2.602244e-05 ### with just chromosomes, without is 1.629226e-06
het$het[het$genome=='zTIL18']=9.86714e-05  ### with just chromosomes, witout is 9.86714e-05
het$het[het$genome=='zTIL01']=0.0001058948 ## just chromosmes, didn't calculate without..
het$het[het$genome=='ccitra']=het$ks[het$genome=='ccitra']
het$het[het$genome=='vcuspi']=het$ks[het$genome=='vcuspi']
het$het[het$genome=='atenui']=het$ks[het$genome=='atenui']
het$het[het$genome=='hconto']=het$ks[het$genome=='hconto']
het$het[het$genome=='hcompr']=het$ks[het$genome=='hcompr']
het$het[het$genome=='zdgigi']=final_merged_data %>% filter(ploidy!='Diploid', !is.na(shortspeciesLabel), genome=='zdgigi') %>% group_by(queryChr, refChr, paspChr, genome,ploidy, speciesLabel, shortspeciesLabel) %>%summarize(ks=mean(V17), n=length(unique(gene))) %>% filter(n>10) %>% filter(ks<0.05) %>% group_by(genome, shortspeciesLabel, ploidy) %>% summarize(lowks=mean(ks)) %>% pull(lowks)
het$het[het$genome=='blagur']=final_merged_data %>% filter(ploidy!='Diploid', !is.na(shortspeciesLabel), genome=='blagur') %>% group_by(queryChr, refChr, paspChr, genome,ploidy, speciesLabel, shortspeciesLabel) %>%summarize(ks=mean(V17), n=length(unique(gene))) %>% filter(n>10) %>% filter(ks<0.05) %>% group_by(genome, shortspeciesLabel, ploidy) %>% summarize(lowks=mean(ks)) %>% pull(lowks)
het$het[het$genome=='agerar']=final_merged_data %>% filter(ploidy!='Diploid', !is.na(shortspeciesLabel), genome=='agerar') %>% group_by(queryChr, refChr, paspChr, genome,ploidy, speciesLabel, shortspeciesLabel) %>%summarize(ks=mean(V17), n=length(unique(gene))) %>% filter(n>10) %>% filter(ks<0.05) %>% group_by(genome, shortspeciesLabel, ploidy) %>% summarize(lowks=mean(ks)) %>% pull(lowks)
het$het[het$genome=='sscopa']=final_merged_data %>% filter(ploidy!='Diploid', !is.na(shortspeciesLabel), genome=='sscopa') %>% group_by(queryChr, refChr, paspChr, genome,ploidy, speciesLabel, shortspeciesLabel) %>%summarize(ks=mean(V17), n=length(unique(gene))) %>% filter(n>10) %>% filter(ks<0.02) %>% group_by(genome, shortspeciesLabel, ploidy) %>% summarize(lowks=mean(ks)) %>% pull(lowks)
het$het[het$genome=='snutan']=final_merged_data %>% filter(ploidy!='Diploid', !is.na(shortspeciesLabel), genome=='snutan') %>% group_by(queryChr, refChr, paspChr, genome,ploidy, speciesLabel, shortspeciesLabel) %>%summarize(ks=mean(V17), n=length(unique(gene))) %>% filter(n>10) %>% filter(ks<0.02) %>% group_by(genome, shortspeciesLabel, ploidy) %>% summarize(lowks=mean(ks)) %>% pull(lowks)
## add back diploids
het$het[het$genome=='cserru']=NA ## ignore nanopore  het$ks[het$genome=='cserru']
het$het[het$genome=='telega']=het$ks[het$genome=='telega']
het$het[het$genome=='ttrian']=het$ks[het$genome=='ttrian']

het$ks[het$ploidy=='Diploid']=0
#het$het[het$genome%in% c('zTIL01', 'zmB735', 'zTIL11', 'zTIL25', 'zTIL18', 'zmhuet', 'zluxur')]=0

het$type=ifelse(het$het==het$ks, 'autopolyploid', 'allopolyploid')
het$type[het$genome%in%c('achine', 'ccitra', 'blagur', 'agerar', 'sscopa')| het$ploidy=='Paleotetraploid']='non-monophyletic allopolyploid'


library(ggrepel)
heterozygosity=ggplot(het, aes(x = het, y = 0, color=ploidy)) +  # Keep y=0 to represent the number line , pch=type
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +  # Horizontal number line
  geom_point(size = 3, alpha=0.8) +  # Plot points on the number line
  geom_text_repel(aes(label = shortspeciesLabel), 
                  nudge_y = 0.2,       # Adjust the vertical position of the labels
                  size = 2,            # Label size
                  direction = "y",     # Labels are repelled vertically
                  arrow = arrow(length = unit(0.02, "npc"), type = "closed"),  # Add arrows from labels to points
                  box.padding = 0.5,   # Space between label and the point
                  point.padding = 0.5, max.overlaps = 15) +  # Space between points and labels
  scale_color_manual(values=ploidycolors)+theme(legend.position='none') +xlab('Heterozygosity')+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #     axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),axis.line.y=element_blank(),
        panel.grid = element_blank())  # Remove y-axis and other clutter

hethomeo=ggplot(het[!is.na(het$type),], aes(x=het, y=ks, color=ploidy, pch=ifelse(type=='autopolyploid', 'autopolyploid', ifelse(ploidy=='Diploid', 'diploid', 'allopolyploid'))))+geom_abline(slope=1, intercept=0, lwd=5, color='gray90') +
  geom_point(size=4, alpha=0.8)+
  scale_color_manual(values=ploidycolors)+theme(legend.position=c(0.55,0.9), legend.text = element_text(size=10)) +
  ylab('Ks between\nduplicates') + xlab('Heterozygosity') + labs(color='Ploidy', pch='') +guides(color='none')+
  theme() # Remove y-axis and other clutter


hethomeomya=ggplot(het[!is.na(het$type),], aes(x=het, y=ks/2/6.5e-9/1e6, color=ploidy, pch=ifelse(type=='autopolyploid', 'autopolyploid', ifelse(ploidy=='Diploid', 'diploid', 'allopolyploid'))))+geom_abline(slope=76.92, intercept=0, lwd=5, color='gray90') +
  geom_point(size=3, alpha=0.8)+
  scale_color_manual(values=ploidycolors)+theme(legend.position=c(0.61,0.9), legend.text = element_text(size=8)) +
  ylab('Parental Divergence (Mya)') + xlab('Heterozygosity') + labs(color='Ploidy', pch='') +guides(color='none')+
  theme() # Remove y-axis and other clutter




