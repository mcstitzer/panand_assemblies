

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(rtracklayer)
library(stringr)
library(zoo)
library(patchwork)

all=read.table('../panand_sp_ploidy.txt')
all$V3[all$V2=='rtuber']=2
all$V3[all$V2=='telega']=2

anchors=lapply(all$V2, function(x) {
  a=read.table(paste0('../syntenic_anchors/anchors/',x, '-Pv-', 2*all$V3[all$V2==x]), header=T)
  a$genome=x
  return(a)
})

ab=Reduce(function(...) merge(..., all=T), anchors)

b=ab[ab$gene!='interanchor',] ## keep only genes

genedists=data.frame(genome=all$V2, meandist=NA, mediandist=NA)    
gdd=vector(mode = "list", length = length(genedists$genome))
names(gdd)=genedists$genome

gddna=gdd

synt=read.table('~/Downloads/sharedSyntenicAnchors.txt', header=F)


taxonnames=c("Z. mays ssp. parviglumis TIL11", "Z. mays ssp. mays B73v5", "Z. mays ssp. parviglumis TIL01", "Z. mays ssp. mexicana TIL25", "Z. mays ssp. mexicana TIL18", "Z. mays ssp. huehuetenangensis", 
             "Z. luxurians", "Z. nicaraguensis", "Z. diploperennis Momo", "Z. diploperennis Gigi", "T. zoloptense", "T. dactyloides FL", "T. dactyloides Southern Hap2", 
             "T. dactyloides Northern Hap2", "T. dactyloides KS", "T. dactyloides tetraploid", "U. digitatum", "V. cuspidata", "R. rottboellioides", "R. tuberculosa", 
             "H. compressa", "E. tripsacoides", "S. scoparium", "S. microstachyum", "A. virginicum", "A. chinensis", "A. gerardi", 
             "C. refractus", "C. citratus", "H. contortus", "T. triandra", "B. laguroides", "P. paniceum", "S. bicolor", 
             "I. rugosum", "S. nutans", '"A." burmanicus', "T. elegans", "C. serrulatus", "P. vaginatum")
names(taxonnames)=c("zTIL11", "zmB735", "zTIL01", "zTIL25", "zTIL18", "zmhuet", 
                    "zluxur", "znicar", "zdmomo", "zdgigi", "tzopol", "tdacs1", "tdacs2", 
                    "tdacn2", "tdacn1", "tdactm", "udigit", "vcuspi", "rrottb", "rtuber", 
                    "hcompr", "etrips", "sscopa", "smicro", "avirgi", "achine", "agerar", 
                    "crefra", "ccitra", "hconto", "ttrian", "blagur", "ppanic", "sbicol", 
                    "irugos", "snutan", "atenui", "telega", "cserru", "pvagin")



results_list=vector(mode='list', length=length(all$V2)) #empty_list <- vector(mode = "list", length = desired_length)
names(results_list)=all$V2



#### set up output, make fractionation bias and fractionation totals for each genome

##pdf('~/Downloads/fractionationBias.rollingSyntenic.pdf',8,4)
for(x in all$V2){
  bg=b[b$genome==x,]
  bg <- bg %>% filter(substr(gene, 1, 7) == "Pavag01")
  #  a1=makeGRangesFromDataFrame(bg[!duplicated(bg$gene),1:3], seqnames.field='refChr')
  a2=makeGRangesFromDataFrame(bg[,4:6], seqnames.field='queryChr')
  a2$gene=bg$gene
  a2$blockIndex=bg$blockIndex
  a2synt=makeGRangesFromDataFrame(bg[bg$gene %in%paste0(synt$V1,'.1.v3.1'),4:6], seqnames.field='queryChr')
  a2synt$gene=bg$gene[bg$gene %in%paste0(synt$V1,'.1.v3.1')]
  
  a2$pvChr=substr(a2$gene,1,7)
  a2$geneIndex=as.numeric(factor(a2$gene))
  
  a2synt$pvChr=substr(a2synt$gene,1,7)
  a2synt$geneIndex=as.numeric(factor(a2synt$gene))
  #fractBias=data.frame(a2synt) %>% group_by(window=geneIndex%/%100, pvChr, seqnames) %>% summarize(fractBias=n()/100)
  fractBias <- data.frame(a2synt) %>%
    group_by(window = geneIndex %/% 100, pvChr, seqnames) %>%
    summarize(fractBias = n()/100, .groups = 'drop') %>%
    left_join(
      data.frame(a2synt) %>%
        group_by(geneIndex, pvChr) %>%
        summarize(geneCount = n(), .groups = 'drop') %>%
        group_by(window = geneIndex %/% 100, pvChr) %>%
        summarize(fractTotal = sum(geneCount)/100, fractMax=max(geneCount), .groups = 'drop'),
      by = c("window", "pvChr")
    )%>% filter(pvChr=='Pavag01')
  # Convert synt to a vector
  synt_genes <- synt$V1
  
  # Extract the relevant gene identifier from bg
  bg <- bg %>%
    mutate(gene_id = sub("\\..*", "", gene)) %>% # Removing version numbers from gene IDs
    filter(gene_id%in%synt_genes)
  # Define a custom function to count bg genes in each window of synt_genes and track additional information
  count_bg_genes_with_index <- function(synt_window, synt_index, bg_df, query_chr) {
    bg_subset <- bg_df %>% filter(queryChr == query_chr)
    count <- sum(bg_subset$gene_id %in% synt_window)
    first_gene <- synt_window[1]
    total = sum(bg_df$gene_id %in% synt_window) ## sort of double counting this because doing for each query chr????? let's see how long it takes
    max = max(table(bg_df$gene_id))
    return(c(count = count, first_gene_index = synt_index, first_gene = first_gene, total=total, max=max))
  }
  
  # Apply the rolling window function on synt_genes and store results
  result <- data.frame()
  
for (chr in unique(bg$queryChr)) {
  chr_length <- length(synt_genes)
rolling_results <- rollapply(1:chr_length, width = 100, FUN = function(idx) {
  synt_window <- synt_genes[idx]
  # Determine the actual size of the window
  actual_window_size <- min(100, chr_length - idx[1] + 1)
  
  result <- count_bg_genes_with_index(synt_window, idx[1], bg, chr)
  
  # Extract numeric values for calculations
  count <- as.numeric(result["count"])
  total <- as.numeric(result["total"])
  
  # Adjust counts for the actual window size
  result["adjusted_fractBias"] <- count / actual_window_size
  result["adjusted_fractTotal"] <- total / actual_window_size
  
  return(result)
}, by.column = FALSE, fill = NA, align = "right")

  rolling_results_df <- as.data.frame(rolling_results)
  rolling_results_df$queryChr <- chr
  rolling_results_df$synt_window_start_index <- as.numeric(rolling_results_df$first_gene_index)
  rolling_results_df$first_gene <- as.character(rolling_results_df$first_gene)
  rolling_results_df$count <- as.numeric(rolling_results_df$count)
  rolling_results_df$total <- as.numeric(rolling_results_df$total)
  rolling_results_df$adjusted_fractBias <- as.numeric(rolling_results_df$adjusted_fractBias)
  rolling_results_df$adjusted_fractTotal <- as.numeric(rolling_results_df$adjusted_fractTotal)
  result <- rbind(result, rolling_results_df)
}

  
  # Add pvChr and fractBias columns
  result$pvChr <- substr(result$first_gene, 1, 7)
  result$fractBias <- result$count / 100
  result$fractTotal <- result$total / 100
  result$copyCount=table(bg$gene_id)[result$first_gene]

  # Display the final result

  results_list[[x]]=result
  # View the result
  #   head(bg)
  # for(chr in c(paste0('Pavag0', 1:9), 'Pavag10')){
  #   
  #   print(
  #     ggplot(result[result$pvChr==chr,], aes(x=synt_window_start_index, y=fractBias, group=queryChr, color=queryChr)) + geom_line() + ggtitle(paste(chr, x)) + theme(legend.position='NA')
  #   )
  # }
  }
#dev.off()  

outlist=results_list
for(i in 1:length(results_list)){
  temp=results_list[[i]]
  temp$genome=all$V2[i]
  outlist[[i]]=temp
}

out=do.call(rbind, outlist)

pdf('~/Downloads/fractionationBias.rollingSyntenic.pdf',8,4)
for(x in 1:length(all$V2)){
  for(chr in c(paste0('Pavag0', 1:9), 'Pavag10')){
    
    print(
      ggplot(results_list[[x]][results_list[[x]]$pvChr==chr,], aes(x=synt_window_start_index, y=fractBias, group=queryChr, color=queryChr)) + geom_line() + ggtitle(paste(chr, all$V2[x])) + theme(legend.position='NA')
    )
  }}

chr='Pavag01'
print(
  ggplot(out[out$pvChr==chr,], aes(x=synt_window_start_index, y=fractBias, group=queryChr, color=queryChr)) + geom_line() + ggtitle(chr) + theme(legend.position='NA') + facet_wrap(~genome)
)

dev.off()


ogout=out
## got to filter the end of the chr, since they're walking off by one gene intervals
## do it by flagging subsequent rows that decrease by 1??
chr='Pavag01'
out=out[!is.na(out$count),]
out$ploidy=asize$ploidy[match(out$genome, asize$V2)]
out=out%>% mutate(delta1=abs(count-lag(count))%in%c(1,2)) %>% filter(!is.na(delta1) & !is.na(copyCount))  ## one or two off??? will this be too permissive?
#out$stretch=rle(out$delta1)$values==T&rle(out$delta1)$lengths>50
out$stretch=data.table::rleid(out$delta1) ### oh dang this is just an id for each
out$stretchLen=as.numeric(table(out$stretch)[out$stretch]) ## this is dumb but i am lazy - just want to see how long each stretch is!!!
out$filterSwitch=out$stretchLen>5 & out$delta1
  


#out=out %>% group_by(genome, queryChr) %>% filter(slice(100:c(n()-100)))

ggplot(out[out$pvChr==chr & out$genome=='crefra'& !out$genome%in%lowQualAssemblies & !out$filterSwitch & out$fractBias>0.2 & out$synt_window_start_index<1981 & out$stretchLen>1,], aes(x=synt_window_start_index, y=fractBias, group=queryChr, color=queryChr)) + geom_line() + ggtitle(chr) + facet_wrap(ploidy~genome) + xlim(0,2081-100)+ylim(0,1)


fig3sp=c('smicro', 'avirgi', 'crefra', 'ppanic', 'etrips', 'achine', 'hconto', 'snutan', 'udigit', 'agerar', 'blagur', 'hcompr', 'zTIL01', 'zluxur', 'tdacs1', 'tdacn1')

ggplot(out[out$genome%in% fig3sp & out$pvChr==chr& !out$genome%in%lowQualAssemblies & !out$filterSwitch & out$fractBias>0.2 & out$synt_window_start_index<1981 & out$stretchLen>1,], 
       aes(x=synt_window_start_index, y=fractBias, group=queryChr, color=ploidy)) + scale_color_manual(values=ploidycolors) + geom_line(alpha=0.8) + theme(legend.position='NA')+ facet_wrap(ploidy~genome) + xlim(0,2081-100)+ylim(0,1)+
       geom_hline(yintercept=c(0,0.25,0.5,0.75,1), lty='dotted', color='gray90') + ylab('Proportion Syntenic Genes Present (per 100 gene window)') + xlab('Paspalum gene index (Chr 1)') + theme(strip.background=element_blank())


#### WATCH oUYT
### changing this here so parv is shortened!!!!! for TIL01
shorttaxonnames=c("Z. mays ssp. parviglumis TIL11",  "Z. mays ssp. parv. TIL01", "Z. mays ssp. mays B73v5", "Z. mays ssp. mexicana TIL25", "Z. mays ssp. mexicana TIL18", "Z. mays ssp. huehuetenangensis", 
                  "Z. luxurians", "Z. nicaraguensis", "Z. diploperennis Momo", "Z. diploperennis Gigi", "T. zoloptense", "T. dactyloides FL", "T. dactyloides Southern Hap2", 
                  "T. dactyloides Northern Hap2", "T. dactyloides KS", "T. dactyloides tetraploid", "U. digitatum", "V. cuspidata", "R. rottboellioides", "R. tuberculosa", 
                  "H. compressa", "E. tripsacoides", "S. scoparium", "S. microstachyum", "A. virginicum", "A. chinensis", "A. gerardi", 
                  "C. refractus", "C. citratus", "H. contortus", "T. triandra", "B. laguroides", "P. paniceum", "S. bicolor", 
                  "I. rugosum", "S. nutans", '"A." burmanicus', "T. elegans", "C. serrulatus", "P. vaginatum")
names(shorttaxonnames)=c("zTIL11",  "zTIL01", "zmB735", "zTIL25", "zTIL18", "zmhuet", 
                         "zluxur", "znicar", "zdmomo", "zdgigi", "tzopol", "tdacs1", "tdacs2", 
                         "tdacn2", "tdacn1", "tdactm", "udigit", "vcuspi", "rrottb", "rtuber", 
                         "hcompr", "etrips", "sscopa", "smicro", "avirgi", "achine", "agerar", 
                         "crefra", "ccitra", "hconto", "ttrian", "blagur", "ppanic", "sbicol", 
                         "irugos", "snutan", "atenui", "telega", "cserru", "pvagin")


## from ks_plotting
out$shortSpeciesLabel=shorttaxonnames[match(out$genome, names(shorttaxonnames))]
out$shortSpeciesLabel=factor(out$shortSpeciesLabel, levels=shorttaxonnames)


fractchr=ggplot(out[out$genome%in% fig3sp & out$pvChr==chr& !out$genome%in%lowQualAssemblies & !out$filterSwitch & out$fractBias>0.2 & out$synt_window_start_index<1981 & out$stretchLen>1,], 
       aes(x=synt_window_start_index, y=fractBias, group=queryChr, color=ploidy)) + scale_color_manual(values=ploidycolors) + geom_line(alpha=0.8) + 
  theme(strip.background = element_blank(), 
        strip.text.y.right = element_text(angle = 0, hjust=0, size=10, vjust=0.5),
        panel.spacing = unit(3, "pt"), 
        axis.text.y = element_blank(), 
        axis.text.x = element_text(size = 9),
        legend.position = 'NULL',
 #       axis.ticks.y=element_blank(),
#        axis.line.y=element_blank()
 ) +
  facet_wrap(ploidy~shortSpeciesLabel, nrow=4, strip.position = 'top', labeller=purrr::partial(label_species, dont_italicize=c('subsp.', 'ssp.', 'TIL11', 'TIL01', 'TIL25', 'TIL18', 'Momo', 'Gigi', 'Southern Hap1', 'Northern Hap1', 'FL', 'KS',  '\\*', '\\"', 'B73v5'))) + 
  xlim(0,2081-100)+ylim(0,1)+
  geom_hline(yintercept=c(0,0.25,0.5,0.75,1), lty='dotted', color='gray90') + ylab('Proportion Syntenic Genes Present (per 100 gene window)') + xlab('Paspalum gene index (Chr 1)') + theme(strip.background=element_blank())


dipFracChrPlot=ggplot(out[out$ploidy=='Diploid' & out$genome%in% fig3sp & out$pvChr==chr& !out$genome%in%lowQualAssemblies & !out$filterSwitch & out$fractBias>0.2 & out$synt_window_start_index<1981 & out$stretchLen>1,], 
                      aes(x=synt_window_start_index, y=fractBias, group=queryChr, color=ploidy)) + scale_color_manual(values=ploidycolors) + geom_line(alpha=0.8) + 
  theme(strip.background = element_blank(), 
        strip.text.y.right = element_text(angle = 0, hjust=0, size=10, vjust=0.5),
        panel.spacing = unit(3, "pt"), 
        axis.text.y = element_blank(), 
        axis.text.x = element_text(size = 9),
        legend.position = 'NULL',
        #       axis.ticks.y=element_blank(),
        #        axis.line.y=element_blank()
  ) +
  facet_wrap(~shortSpeciesLabel, nrow=1, strip.position = 'top', labeller=purrr::partial(label_species, dont_italicize=c('subsp.', 'ssp.', 'TIL11', 'TIL01', 'TIL25', 'TIL18', 'Momo', 'Gigi', 'Southern Hap1', 'Northern Hap1', 'FL', 'KS',  '\\*', '\\"', 'B73v5'))) + 
  xlim(0,2081-100)+ylim(0,1)+
  geom_hline(yintercept=c(0,0.25,0.5,0.75,1), lty='dotted', color='gray90') + ylab('Proportion Syntenic Genes Present (per 100 gene window)') + xlab('Paspalum gene index (Chr 1)') + theme(strip.background=element_blank())



## get medians to plot against
medFrac=out[ out$pvChr==chr& !out$filterSwitch & out$fractBias>0.2 & out$synt_window_start_index<1981 & out$stretchLen>1,]%>% group_by(genome, ploidy) %>% summarize(median=median(fractBias), mean=mean(fractBias))
asize$medFrac=medFrac$median[match(asize$V2, medFrac$genome)]
asize$meanFrac=medFrac$mean[match(asize$V2, medFrac$genome)]

ggplot(asize, aes(x=mya, y=medFrac, color=ploidy)) + geom_point() + scale_color_manual(values=ploidycolors)

asizezt=asize
asizezt$V2zt=ifelse(asize$V2 %in% c('tdacn1', 'tdacs1'), 'trips', asize$V2)
asizezt$V2zt=ifelse(asize$V2 %in% c("zdgigi", "zdmomo", "zluxur", "zmhuet", "zTIL18", "zTIL25", "zTIL01", "zTIL11", "znicar", "zmB735"), 'zea', asize$V2)
asizezt=asizezt%>% group_by(V2zt, ploidy) %>% summarize(medFrac=median(medFrac), meanFrac=median(meanFrac), syntAnchorsCount=median(syntAnchorsCount), syntAnchors=median(syntAnchors), diploidEquivalentsyntAnchors=median(diploidEquivalentsyntAnchors), mya=median(mya))

cor.test(asizezt$medFrac[!asizezt$medFrac%in%lowQualAssemblies], asizezt$mya)
## negative -0.5081658 
cor.test(asizezt$medFrac[!asizezt$ploidy%in%c('Diploid', 'Paleotetraploid')], asizezt$mya[!asizezt$ploidy%in%c('Diploid', 'Paleotetraploid')])
## positive 0.5994824 

medFracPlot=ggplot(asize[!asize$V2%in%lowQualAssemblies,], aes(x=mya, y=medFrac, color=ploidy)) + geom_vline(xintercept=c(2,4,6,8,10,12,14), color='gray90', lty='dotted', alpha=0.3) + geom_vline(xintercept=c(1,3,5,7,9,11,13), color='gray80', lty='dashed', alpha=0.3)+
  geom_point(size=3) + scale_color_manual(values=ploidycolors)+ ylim(0,1)+ 
  theme(legend.position='none') + 
  xlab('Divergence between\nParental Subgenomes (Mya)') + 
  ylab('Median Proportion Genes Retained') +
  stat_smooth(data=asizezt[!asizezt$V2zt%in%lowQualAssemblies,], method='lm', aes(group=NA), se=F, color='gray80')#+
#  stat_smooth(data=asizezt[!asizezt$ploidy%in%c('Diploid', 'Paleotetraploid'),], method='lm', aes(group=NA), se=F, color='gray80', lty='longdash')






ggplot(asize, aes(x=mya, y=meanFrac, color=ploidy)) + geom_point() + scale_color_manual(values=ploidycolors) + geom_text(aes(label=V2))
## maybe filter by contig length?????? onlly count things >10 mb?? 
#### MUST DO THIS FILTER!!! gross genomes are ruining it all...

cor.test(asize$mya, asize$medFrac)
## so negative correlation....
## these high fractionation diploids worry me... are these actually tetraploids??? i guess they are the nanorpore and the permanent translocation heteroygote....
## i think filtering by contig lenght shoudl clean up the teosintes that are getting blumped by the duplicate genes from extra scaffolds

### oh alos i should calculate these values for the genome as a whole!!
## can i put these histograms alongside the y axis of the chromosome plots??
## so maybe i do want a vertical column
asize$genome=asize$V2
ggplot(out[ out$pvChr==chr& !out$genome%in%lowQualAssemblies & !out$filterSwitch & out$fractBias>0.2 & out$synt_window_start_index<1981 & out$stretchLen>1,],
       aes(x=fractBias, color=ploidy, fill=ploidy))+ geom_density() + facet_wrap(~genome, scales='free_x') + xlim(0,1) + scale_color_manual(values=ploidycolors) + scale_fill_manual(values=ploidycolors) +coord_flip() + geom_vline(data=asize, aes(xintercept=medFrac), lty='dashed', alpha=0.5)


ggplot(out[out$pvChr==chr & !out$genome%in%lowQualAssemblies & !out$filterSwitch,], aes(x=synt_window_start_index, y=fractBias, group=queryChr, color=queryChr)) + geom_line() + ggtitle(chr) + theme(legend.position='NA') + facet_wrap(ploidy~genome) + xlim(0,2081-100)

## try grouping by ploicy??
ggplot(out[out$pvChr==chr & !out$genome%in%lowQualAssemblies & out$first_gene_index<2081-100,], aes(x=synt_window_start_index, y=fractBias, group=queryChr, color=queryChr)) + geom_line() + ggtitle(chr) + theme(legend.position='NA') + facet_wrap(~ploidy, ncol=1)


################# combine linea nd density


label_species <- function(species_name, dont_italicize = c()) {
 species_name=as.character(species_name)
  # Split the species name into words
  words <- unlist(strsplit(species_name, " "))
  
  # Process each word: italicize unless it matches dont_italicize
  formatted_words <- sapply(words, function(word) {
    if (word %in% dont_italicize) {
      return(paste0("'", word, "'"))  # Leave plain with quotes
    } else {
      return(paste0("italic('", word, "')"))  # Wrap in italic() with quotes
    }
  })
  
  # Combine words into a single expression with spaces explicitly added
  combined_label <- paste(formatted_words, collapse = " * ' ' * ")
  return(parse(text = combined_label))  # Convert to R expression
}




# Filter data
filtered_data <- out %>%
  filter(ploidy == 'Diploid',
         genome %in% fig3sp,
         pvChr == chr,
         !genome %in% lowQualAssemblies,
         !filterSwitch,
         fractBias > 0.2,
         synt_window_start_index < 1981,
         stretchLen > 1)


# Unique species
species_list <- as.character(unique(filtered_data$shortSpeciesLabel))

# Calculate medians for each species
medians <- filtered_data %>%
  group_by(shortSpeciesLabel) %>%
  summarize(medFrac = median(fractBias, na.rm = TRUE))

# Create alternating line and density plots
plot_list <- list()
for (i in seq_along(species_list)) {
  species <- species_list[i]
  # Subset data for the current species
  species_data <- filtered_data %>% filter(shortSpeciesLabel == species)
  species_median <- medians %>% filter(shortSpeciesLabel == species)
  species_label <- label_species(species, dont_italicize = c('subsp.', 'ssp.', 'TIL11', 'TIL01', 
                                                             'TIL25', 'TIL18', 'Momo', 'Gigi', 
                                                             'Southern Hap1', 'Northern Hap1', 
                                                             'FL', 'KS', '\\*', '\\"', 'B73v5'))
  
  ## filter edges:
  species_data <- species_data %>%
  group_by(queryChr) %>% # Group by queryChr
  filter(row_number() > 90 & row_number() <= (n() - 90)) %>% # Keep entries outside the first and last 100
  ungroup() # Ungroup after filtering

  
  # Line plot
  line_plot <- ggplot(species_data, aes(x = synt_window_start_index, y = fractBias, group = queryChr, color = ploidy)) +
    geom_line(alpha = 0.8) +
    scale_color_manual(values = ploidycolors) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 8),
          axis.title.y = if (i == 1) element_text(size = 10) else element_blank(), # Y-axis title for first line plot only
          axis.text.y = if (i == 1) element_text(size = 8) else element_blank(),  # Y-axis labels only for first plot
          axis.ticks.y = if (i == 1) element_line() else element_blank(),         # Y-axis ticks only for first plot
          plot.title = element_text(size = 10, hjust = 0.5)) +
    geom_hline(yintercept=c(0,0.25,0.5,0.75,1), lty='dotted', color='gray90') +
    ggtitle(species_label) +
    ylim(0, 1) +
    labs(x = "Gene index", y = if (i == 1) "Proportion Syntenic\nGenes Present" else NULL)
  
  # Density plot with median line
  density_plot <- ggplot(species_data, aes(x = fractBias, color = ploidy, fill = ploidy)) +
    geom_density(alpha = 0.6) +
    geom_vline(xintercept = species_median$medFrac, linetype = "dashed", color = "gray40", alpha = 0.8) +
    scale_color_manual(values = ploidycolors) +
    scale_fill_manual(values = ploidycolors) +
    coord_flip() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),  # Remove x-axis tick labels
          axis.ticks.x = element_blank(), # Remove x-axis ticks
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line=element_blank()) +
    xlim(0, 1) +
    labs(y = NULL, x = NULL)
  
  # Combine line plot and density plot horizontally with width ratio
  combined <- line_plot + density_plot + plot_layout(widths = c(6, 2), guides = "collect") & theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
  
  # Add combined plot to the list
  plot_list[[species]] <- combined
}

# Combine all plots in one row
final_plot_dip <- wrap_plots(plot_list, nrow = 1)
final_plot_dip



#### tetraploid

# Filter data
filtered_data <- out %>%
  filter(ploidy == 'Tetraploid',
         genome %in% fig3sp,
         pvChr == chr,
         !genome %in% lowQualAssemblies,
         !filterSwitch,
         fractBias > 0.2,
         synt_window_start_index < 1981,
         stretchLen > 1)




# Unique species
species_list <- as.character(unique(filtered_data$shortSpeciesLabel))

# Calculate medians for each species
medians <- filtered_data %>%
  group_by(shortSpeciesLabel) %>%
  summarize(medFrac = median(fractBias, na.rm = TRUE))

# Create alternating line and density plots
plot_list <- list()
for (i in seq_along(species_list)) {
  species <- species_list[i]
  # Subset data for the current species
  species_data <- filtered_data %>% filter(shortSpeciesLabel == species)
  species_median <- medians %>% filter(shortSpeciesLabel == species)
  species_label <- label_species(species, dont_italicize = c('subsp.', 'ssp.', 'TIL11', 'TIL01', 
                                                             'TIL25', 'TIL18', 'Momo', 'Gigi', 
                                                             'Southern Hap1', 'Northern Hap1', 
                                                             'FL', 'KS', '\\*', '\\"', 'B73v5'))
    ## filter edges:
  species_data <- species_data %>%
  group_by(queryChr) %>% # Group by queryChr
  filter(row_number() > 90 & row_number() <= (n() - 90)) %>% # Keep entries outside the first and last 100
  ungroup() # Ungroup after filtering

  # Line plot
  line_plot <- ggplot(species_data, aes(x = synt_window_start_index, y = fractBias, group = queryChr, color = ploidy)) +
    geom_line(alpha = 0.8) +
    scale_color_manual(values = ploidycolors) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 8),
          axis.title.y = if (i == 1) element_text(size = 10) else element_blank(), # Y-axis title for first line plot only
          axis.text.y = if (i == 1) element_text(size = 8) else element_blank(),  # Y-axis labels only for first plot
          axis.ticks.y = if (i == 1) element_line() else element_blank(),         # Y-axis ticks only for first plot
          plot.title = element_text(size = 10, hjust = 0.5)) +
    geom_hline(yintercept=c(0,0.25,0.5,0.75,1), lty='dotted', color='gray90') +
    ggtitle(species_label) +
    ylim(0, 1) +
    labs(x = "Gene index", y = if (i == 1) "Proportion Syntenic\nGenes Present" else NULL)
  
  # Density plot with median line
  density_plot <- ggplot(species_data, aes(x = fractBias, color = ploidy, fill = ploidy)) +
    geom_density(alpha = 0.6) +
    geom_vline(xintercept = species_median$medFrac, linetype = "dashed", color = "gray40", alpha = 0.8) +
    scale_color_manual(values = ploidycolors) +
    scale_fill_manual(values = ploidycolors) +
    coord_flip() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),  # Remove x-axis tick labels
          axis.ticks.x = element_blank(), # Remove x-axis ticks
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line=element_blank()) +
    xlim(0, 1) +
    labs(y = NULL, x = NULL)
  
  # Combine line plot and density plot horizontally with width ratio
  combined <- line_plot + density_plot + plot_layout(widths = c(6, 2), guides = "collect") & theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
  
  # Add combined plot to the list
  plot_list[[species]] <- combined
}

# Combine all plots in one row
final_plot_tetra <- wrap_plots(plot_list, nrow = 1)
final_plot_tetra


####### hexaploid

# Filter data
filtered_data <- out %>%
  filter(ploidy == 'Hexaploid',
         genome %in% fig3sp,
         pvChr == chr,
         !genome %in% lowQualAssemblies,
         !filterSwitch,
         fractBias > 0.2,
         synt_window_start_index < 1981,
         stretchLen > 1)



# Unique species
species_list <- as.character(unique(filtered_data$shortSpeciesLabel))

# Calculate medians for each species
medians <- filtered_data %>%
  group_by(shortSpeciesLabel) %>%
  summarize(medFrac = median(fractBias, na.rm = TRUE))

# Create alternating line and density plots
plot_list <- list()
for (i in seq_along(species_list)) {
  species <- species_list[i]
  # Subset data for the current species
  species_data <- filtered_data %>% filter(shortSpeciesLabel == species)
  species_median <- medians %>% filter(shortSpeciesLabel == species)
  species_label <- label_species(species, dont_italicize = c('subsp.', 'ssp.', 'TIL11', 'TIL01', 
                                                             'TIL25', 'TIL18', 'Momo', 'Gigi', 
                                                             'Southern Hap1', 'Northern Hap1', 
                                                             'FL', 'KS', '\\*', '\\"', 'B73v5'))
    ## filter edges:
  species_data <- species_data %>%
  group_by(queryChr) %>% # Group by queryChr
  filter(row_number() > 90 & row_number() <= (n() - 90)) %>% # Keep entries outside the first and last 100
  ungroup() # Ungroup after filtering

  # Line plot
  line_plot <- ggplot(species_data, aes(x = synt_window_start_index, y = fractBias, group = queryChr, color = ploidy)) +
    geom_line(alpha = 0.8) +
    scale_color_manual(values = ploidycolors) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 8),
          axis.title.y = if (i == 1) element_text(size = 10) else element_blank(), # Y-axis title for first line plot only
          axis.text.y = if (i == 1) element_text(size = 8) else element_blank(),  # Y-axis labels only for first plot
          axis.ticks.y = if (i == 1) element_line() else element_blank(),         # Y-axis ticks only for first plot
          plot.title = element_text(size = 10, hjust = 0.5)) +
    geom_hline(yintercept=c(0,0.25,0.5,0.75,1), lty='dotted', color='gray90') +
    ggtitle(species_label) +
    ylim(0, 1) +
    labs(x = "Gene index", y = if (i == 1) "Proportion Syntenic\nGenes Present" else NULL)
  
  # Density plot with median line
  density_plot <- ggplot(species_data, aes(x = fractBias, color = ploidy, fill = ploidy)) +
    geom_density(alpha = 0.6) +
    geom_vline(xintercept = species_median$medFrac, linetype = "dashed", color = "gray40", alpha = 0.8) +
    scale_color_manual(values = ploidycolors) +
    scale_fill_manual(values = ploidycolors) +
    coord_flip() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),  # Remove x-axis tick labels
          axis.ticks.x = element_blank(), # Remove x-axis ticks
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line=element_blank()) +
    xlim(0, 1) +
    labs(y = NULL, x = NULL)
  
  # Combine line plot and density plot horizontally with width ratio
  combined <- line_plot + density_plot + plot_layout(widths = c(6, 2), guides = "collect") & theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
  
  # Add combined plot to the list
  plot_list[[species]] <- combined
}

# Combine all plots in one row
final_plot_hex <- wrap_plots(plot_list, nrow = 1)
final_plot_hex

### paleotetrploid

# Filter data
filtered_data <- out %>%
  filter(ploidy == 'Paleotetraploid',
         genome %in% fig3sp,
         pvChr == chr,
         !genome %in% lowQualAssemblies,
         !filterSwitch,
         fractBias > 0.2,
         synt_window_start_index < 1981,
         stretchLen > 1)




# Unique species
species_list <- as.character(unique(filtered_data$shortSpeciesLabel))

# Calculate medians for each species
medians <- filtered_data %>%
  group_by(shortSpeciesLabel) %>%
  summarize(medFrac = median(fractBias, na.rm = TRUE))

# Create alternating line and density plots
plot_list <- list()
for (i in seq_along(species_list)) {
  species <- species_list[i]
  # Subset data for the current species
  species_data <- filtered_data %>% filter(shortSpeciesLabel == species)
  species_median <- medians %>% filter(shortSpeciesLabel == species)
  species_label <- label_species(species, dont_italicize = c('subsp.', 'ssp.', 'TIL11', 'TIL01', 
                                                             'TIL25', 'TIL18', 'Momo', 'Gigi', 
                                                             'Southern Hap1', 'Northern Hap1', 
                                                             'FL', 'KS', '\\*', '\\"', 'B73v5'))
    ## filter edges:
  species_data <- species_data %>%
  group_by(queryChr) %>% # Group by queryChr
  filter(row_number() > 90 & row_number() <= (n() - 90)) %>% # Keep entries outside the first and last 100
  ungroup() # Ungroup after filtering

  # Line plot
  line_plot <- ggplot(species_data, aes(x = synt_window_start_index, y = fractBias, group = queryChr, color = ploidy)) +
    geom_line(alpha = 0.8) +
    scale_color_manual(values = ploidycolors) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 8),
          axis.title.y = if (i == 1) element_text(size = 10) else element_blank(), # Y-axis title for first line plot only
          axis.text.y = if (i == 1) element_text(size = 8) else element_blank(),  # Y-axis labels only for first plot
          axis.ticks.y = if (i == 1) element_line() else element_blank(),         # Y-axis ticks only for first plot
          plot.title = element_text(size = 10, hjust = 0.5)) +
    geom_hline(yintercept=c(0,0.25,0.5,0.75,1), lty='dotted', color='gray90') +
    ggtitle(species_label) +
    ylim(0, 1) +
    labs(x = "Gene index", y = if (i == 1) "Proportion Syntenic\nGenes Present" else NULL)
  
  # Density plot with median line
  density_plot <- ggplot(species_data, aes(x = fractBias, color = ploidy, fill = ploidy)) +
    geom_density(alpha = 0.6) +
    geom_vline(xintercept = species_median$medFrac, linetype = "dashed", color = "gray40", alpha = 0.8) +
    scale_color_manual(values = ploidycolors) +
    scale_fill_manual(values = ploidycolors) +
    coord_flip() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),  # Remove x-axis tick labels
          axis.ticks.x = element_blank(), # Remove x-axis ticks
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line=element_blank()) +
    xlim(0, 1) +
    labs(y = NULL, x = NULL)
  
  # Combine line plot and density plot horizontally with width ratio
  combined <- line_plot + density_plot + plot_layout(widths = c(6, 2), guides = "collect") & theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
  
  # Add combined plot to the list
  plot_list[[species]] <- combined
}

# Combine all plots in one row
# Combine all species plots into one row
final_plot_paleo <- wrap_plots(plot_list, nrow = 1)
final_plot_paleo


plot_grid(final_plot_dip, final_plot_tetra, final_plot_hex, final_plot_paleo + xlab('Paspalum Chr 1 Gene Position'), ncol=1, labels = 'AUTO', align='v', axis='lr' )

# Combine all final plots with plot_grid
combined_plots <- plot_grid(
  final_plot_dip, 
  final_plot_tetra, 
  final_plot_hex, 
  final_plot_paleo, 
  ncol = 1, 
  labels = 'AUTO', 
  align = 'v', 
  axis = 'lr'  # Align axes vertically
)

# Add shared x-axis label below the combined plots
final_with_xlab <- plot_grid(
  combined_plots,
  ggdraw() + draw_label("Paspalum Chr 1 Gene Position", size = 12, hjust = 0.5),
  ncol = 1,
  rel_heights = c(1, 0.05)  # Adjust relative heights: 95% plots, 5% x-axis label
)

# Display the final plot
final_with_xlab



## now try to plot points for each subgenome/chr
out$max_n=as.numeric(out$max)
## filter low qual
sgfrac=out[out$genome%in% fig3sp & out$pvChr==chr& !out$genome%in%lowQualAssemblies & !out$filterSwitch & out$fractBias>0.2 & out$synt_window_start_index<1981 & out$stretchLen>1,]%>% 
           group_by(shortSpeciesLabel, queryChr, genome, max_n) %>% 
summarize(medFrac=median(fractBias, na.rm=T), count=dplyr::n()) %>% ungroup()
## keep all
# sgfrac=out[ out$pvChr==chr& !out$genome%in%lowQualAssemblies & !out$filterSwitch & out$fractBias>0.2 & out$synt_window_start_index<1981 & out$stretchLen>1,]%>% 
#            group_by(shortSpeciesLabel, queryChr, genome, max_n) %>% 
# summarize(medFrac=median(fractBias, na.rm=T), count=dplyr::n()) %>% ungroup()


sgfrac <- sgfrac %>%
  group_by(genome) %>%  # Split by genome
  slice_max(order_by=count, n=6) %>% ## whatever this is not working dynamicaly
  filter(count > 50)  # Filter rows with count > 100
sgfrac$ploidy=factor(asize$ploidy[match(sgfrac$genome, asize$V2)], levels=c('Diploid', 'Tetraploid', 'Hexaploid', 'Paleotetraploid'))
sgfrac$mya=asize$mya[match(sgfrac$genome, asize$V2)]

sgfrac=sgfrac[sgfrac$queryChr!='alt-scaf_37',] ## remove luxurians alt scaffold

## cjhatgpt solution
## mnevermid sometimes i hate chatgpt
sgfrac_plot=ggplot(sgfrac, 
  aes(x=genome, y=medFrac, color=ploidy, group=genome)) +
  geom_point(size=3, alpha=0.8)+geom_line()+
scale_color_manual(values=ploidycolors)+scale_fill_manual(values=ploidycolors) +
  ylab('Median Proportion Retained Genes') + 
  labs(color='Ploidy', pch='') +guides(color='none', fill='none')+
geom_hline(yintercept=c(0,0.25,0.5,0.75,1), lty='dotted', color='gray90') +
ylim(0,1)+ facet_wrap(~ploidy, ncol=4, drop=T, scales='free_x')+
  theme() # Remove y-axis and other clutter




### plot individually for each species


for(genome in asize$V2){
  #ksplot
  if(genome %in% out$genome){
  pdf(paste0('../figures/fractionation_by_sp/', genome, '_fract.PvChr01.pdf'), 8,4)
  
  



# Filter data
filtered_data <- out[out$genome==genome,] %>%
  filter(
         pvChr == chr,
         !filterSwitch,
         fractBias > 0.2,
         synt_window_start_index < 1981,
         stretchLen > 1)
# Unique species
species_list <- unique(filtered_data$shortSpeciesLabel[filtered_data$genome==genome])

# Calculate medians for each species
medians <- filtered_data %>%
  group_by(shortSpeciesLabel) %>%
  summarize(medFrac = median(fractBias, na.rm = TRUE))

# Create alternating line and density plots
plot_list <- list()
for (i in seq_along(species_list)) {
  species <- species_list[i]
  # Subset data for the current species
  species_data <- filtered_data %>% filter(shortSpeciesLabel == species)
  species_median <- medians %>% filter(shortSpeciesLabel == species)
  species_label=species
#  species_label <- label_species(species, dont_italicize = c('subsp.', 'ssp.', 'TIL11', 'TIL01', 
#                                                             'TIL25', 'TIL18', 'Momo', 'Gigi', 
#                                                             'Southern Hap1', 'Northern Hap1', 
#                                                             'FL', 'KS', '\\*', '\\"', 'B73v5'))
  
  # Line plot
  line_plot <- ggplot(species_data, aes(x = synt_window_start_index, y = fractBias, group = queryChr, color = ploidy)) +
    geom_line(alpha = 0.8) +
    scale_color_manual(values = ploidycolors) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 8),
          axis.title.y = if (i == 1) element_text(size = 10) else element_blank(), # Y-axis title for first line plot only
          axis.text.y = if (i == 1) element_text(size = 8) else element_blank(),  # Y-axis labels only for first plot
          axis.ticks.y = if (i == 1) element_line() else element_blank(),         # Y-axis ticks only for first plot
          plot.title = element_text(size = 10, hjust = 0.5)) +
    geom_hline(yintercept=c(0,0.25,0.5,0.75,1), lty='dotted', color='gray90') +
    ggtitle(species_label) +
    ylim(0, 1) +
    labs(x = "Gene index", y = if (i == 1) "Proportion Syntenic\nGenes Present" else NULL)
  
  # Density plot with median line
  density_plot <- ggplot(species_data, aes(x = fractBias, color = ploidy, fill = ploidy)) +
    geom_density(alpha = 0.6) +
    geom_vline(xintercept = species_median$medFrac, linetype = "dashed", color = "gray40", alpha = 0.8) +
    scale_color_manual(values = ploidycolors) +
    scale_fill_manual(values = ploidycolors) +
    coord_flip() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),  # Remove x-axis tick labels
          axis.ticks.x = element_blank(), # Remove x-axis ticks
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line=element_blank()) +
    xlim(0, 1) +
    labs(y = NULL, x = NULL)
  
  # Combine line plot and density plot horizontally with width ratio
  combined <- line_plot + density_plot + plot_layout(widths = c(6, 1), guides = "collect") & theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))
  
  # Add combined plot to the list
  plot_list[[species]] <- combined
}

# Combine all plots in one row
# Combine all species plots into one row
final_plot_sp <- wrap_plots(plot_list, nrow = 1)


# Add shared x-axis label below the combined plots
final_with_xlab <- plot_grid(
  final_plot_sp,
  ggdraw() + draw_label("Paspalum Chr 1 Gene Position", size = 12, hjust = 0.5),
  ncol = 1,
  rel_heights = c(1, 0.05)  # Adjust relative heights: 95% plots, 5% x-axis label
)

print(final_with_xlab)
dev.off()
}}
