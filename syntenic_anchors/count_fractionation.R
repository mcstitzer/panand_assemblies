

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(rtracklayer)
library(stringr)

all=read.table('../panand_sp_ploidy.txt')
all$V3[all$V2=='rtuber']=2
all$V3[all$V2=='telega']=2

anchors=lapply(all$V2, function(x) {
  a=read.table(paste0('anchors/',x, '-Pv-', 2*all$V3[all$V2==x]), header=T)
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
    )
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
    rolling_results <- rollapply(1:length(synt_genes), width = 100, FUN = function(idx) {
      synt_window <- synt_genes[idx]
      count_bg_genes_with_index(synt_window, idx[1], bg, chr)
    }, by.column = FALSE, fill = NA, align = "right")
    
    rolling_results_df <- as.data.frame(rolling_results)
    rolling_results_df$queryChr <- chr
    rolling_results_df$synt_window_start_index <- as.numeric(rolling_results_df$first_gene_index)
    rolling_results_df$first_gene <- as.character(rolling_results_df$first_gene)
    rolling_results_df$count <- as.numeric(rolling_results_df$count)
    rolling_results_df$total=as.numeric(rolling_results_df$total)
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
