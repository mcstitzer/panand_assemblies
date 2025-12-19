

## plot tandem repeats found on each query scaffold for each reference chromosome - 
all=read.table('../panand_sp_ploidy.txt')
all=all[!all$V2 %in% c('pprate', 'tdactm', 'tzopol', 'osativ', 'bdista', 'agerjg', 'svirid', 'eophiu'),]


species='hconto'
#filepath='../syntenic_anchors/anchors/crefra-Pv-2'
filepath='../syntenic_anchors/anchors/hconto-Pv-4'
color_palette=muted_colors
minBlock=10
title=''
refChrs=c(paste0('Chr0', 1:9), 'Chr10')
queryChrs=''
queryChrtoFlip=''
ylabelspecies=''
ploidy=''

trfile=paste0('grep Tandem ~/Downloads/repeatmask_tandems/', all$V1[all$V2==species], '_EDTAandTandemRepeat.gff3')

tr=fread(trfile)
#crscaf=c('scaf_1', 'scaf_21', 'scaf_8', 'scaf_4', 'scaf_7', 'scaf_16', 'scaf_5', 'scaf_19', 'scaf_2', 'scaf_9', 'scaf_20', 'scaf_3', 'scaf_11', 'scaf_14', 'scaf_18', 'scaf_15', 'scaf_17', 'scaf_10', 'scaf_12', 'scaf_13')
tr$rp=paste0(str_split_fixed(tr$V9, '_', 4)[,2], '_', str_split_fixed(tr$V9, '_', 4)[,3])
#tr=tr[tr$rp!='TandemRepeat96_59bp',]
tr$queryChr=tr$V1


#process_anchors_to_dotplot_SUPP <- function(filepath, color_palette=muted_colors, minBlock=10, title='', refChrs=c(paste0('Chr0', 1:9), 'Chr10'), queryChrs='', 
#                                              queryChrtoFlip='', ylabelspecies='',ploidy='') {
  # Load data
  data <- read.table(filepath, header = TRUE)
  data <- data[data$gene != 'interanchor', ]
  
  ## get queryChrs if they aren't supplied
  if(queryChrs[1]==''){
    queryChrs=unique(data$queryChr)
  }
  
  
  # Reduce to blocks and calculate stats
  data <- data %>%
    group_by(blockIndex) %>%
    mutate(blockLength = dplyr::n()) %>%
    group_by(queryChr) %>%
    mutate(freqStrand = names(which.max(table(strand))),
           maxChr = max(queryStart),
           freqRef = names(which.max(table(refChr))))
  
  # Filter data based on block length
  data <- data[data$blockLength > minBlock, ]
  data$refChr <- factor(data$refChr, levels = c(paste0('Chr0', 1:9), 'Chr10'))
  
  # Reverse strand calculations
  data <- data %>%
    arrange(freqRef, referenceStart, queryStart)
  data$queryChr <- factor(data$queryChr, levels = rev(data$queryChr[!duplicated(data$queryChr)]))
  data$revQueryStart <- data$queryStart
  data$revQueryStart[data$freqStrand == '-'] <- abs(data$queryStart - data$maxChr)[data$freqStrand == '-']
  
  ## flip querychr that look funny
  if(queryChrtoFlip!=''){
    data$revQueryStart[data$queryChr%in%queryChrtoFlip]=abs(data$queryStart - data$maxChr)[data$queryChr%in%queryChrtoFlip]
  }
  
  ## clean up future facet labels - this is very specific to chr level s
  #  data$refLabel=gsub('Chr0', '', data$refChr)
  #  data$refLabel=gsub('Chr', '', data$refLabel)
  #  data$queryLabel=gsub('chr', '', data$queryLabel)
  
  ## set up the image
  #  im=data.frame(refChr='Chr01', queryChr=NA, referenceStart=1, revQueryStart=1, path=pathtokaryotype)
  #  im$queryChr=names(which.max(table(data$queryChr[data$refChr=='Chr01'])))
  
  # Combine the image and the border into a single grob
  #  imgA <- grobTree(border, imgA)
  
  ystriptextsize=ifelse(length(unique(queryChrs))>50, 6,9)
  
  
  
telo=read.csv('../repeats/AllAssembliesTelomeres_formatted_10TR_withLengths.csv', header=T)

telo=telo[telo$Code==species,]
telo$queryChr=telo$Contig
  
  telo$maxChr=data$maxChr[match(telo$queryChr, data$queryChr)]
  telo$teloStart=ifelse(telo$LeftTelomere, 1, NA)
  telo$teloEnd=ifelse(telo$RightTelomere,telo$maxChr,NA)
  
  
  # Create the plot
  p=ggplot(data[ data$refChr %in% names(color_palette) & data$refChr%in%refChrs & data$queryChr%in%queryChrs, ],
           aes(x = referenceStart / 1e6, y = queryStart / 1e6, color = refChr)) +
    geom_point(size=0.5) +
    facet_grid(queryChr ~ refChr, scales = 'free', space='free') +
    scale_color_manual(values = color_palette) +
    theme(legend.position = 'none') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    geom_hline(aes(yintercept=maxChr/1e6), lty='dashed', color='gray') +
    geom_vline(xintercept=0, color='gray95')+
    ggtitle(paste0(ylabelspecies,', ', ploidy)) + xlab('P. vaginatum position (Mb)') +
    ylab(paste0(ylabelspecies, ' position (Mb)')) + 
    #  ylab('Position (Mb)') +
    theme(#strip.text.y = element_text(angle = 0, hjust = 0), 
      strip.placement.y = "outside" , 
      strip.text = element_text(size = 8, color = "darkblue", face = "bold"),
      strip.background = element_rect(fill = "lightblue", color = "darkblue", linewidth = 1),
      strip.text.y = element_text(angle=0, size=ystriptextsize),
      axis.text.x=element_text(size=9),
      axis.text.y=element_text(size=5),
      panel.spacing = unit(0.1, 'lines'),
      plot.title=element_text(color=ploidycolors[ploidy], size=10))+
    scale_x_continuous(breaks = scales::breaks_pretty(n = 2), expand=c(0,0)) +  # Automatically choose 3 breaks for x-axis
    scale_y_continuous(breaks = scales::breaks_pretty(n = 2), expand=c(0,0))    # Automatically choose 3 breaks for y-axis
  
p  +new_scale_color()+geom_hline(data=tr[tr$queryChr%in%crscaf,], aes(yintercept=V4/1e6, color=rp), lty='dotted')+scale_color_brewer(palette='Set1')+theme(legend.position='bottom')

## with telomeres
if(nrow(tr)!=1){
p  +new_scale_color()+geom_hline(data=tr[tr$queryChr%in%queryChrs,], aes(yintercept=V4/1e6, color=rp), lty='dotted')+scale_color_brewer(palette='Set1')+theme(legend.position='bottom')+
    geom_hline(data=telo[telo$queryChr%in%queryChrs,], aes(yintercept=teloStart/1e6), color='purple')+
    geom_hline(data=telo[telo$queryChr%in%queryChrs,], aes(yintercept=teloEnd/1e6), color='purple')
}else{
p+    geom_hline(data=telo[telo$queryChr%in%queryChrs,], aes(yintercept=teloStart/1e6), color='purple')+
    geom_hline(data=telo[telo$queryChr%in%queryChrs,], aes(yintercept=teloEnd/1e6), color='purple')

}

#}

## hmm, these all feel subtelomeric... at least this 185/186 bp repeat






