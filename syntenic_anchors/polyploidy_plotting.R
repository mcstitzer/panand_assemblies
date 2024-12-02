


## /data1/users/mcs368/panand_gene_trees/gene_trees_with_paspalumCDS_Dec2023

## work done on xm01
##/local/workdir/mcs368/gene_trees_with_paspalumCDS_Dec2023

synt=read.table('~/transfer/sharedSyntenicAnchors.txt', header=F)
## synt=read.table('~/Downloads/sharedSyntenicAnchors.txt', header=F)

#RAxML_bipartitionsBranchLabels.Pavag04G184900
#RAxML_bipartitionsBranchLabels


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
library(tidytext)

## need to be in `scinet_trees/` ###AGGagfag
#filenames=list.files('.', pattern='RAxML_bipartitionsBranchLabels.')
filenames=paste0('RAxML_bipartitionsBranchLabels.', synt$V1)
filenames=filenames[file.exists(filenames)]

outdf=data.frame(filenames=filenames)
all=read.table('../panand_sp_ploidy.txt')                                
#for(sp in c('vcuspi', 'rtuber', 'blagur', 'hcompr', 'udigit', 'telega')){
for(sp in all$V2){
  
  outdf[,paste0(sp, 'Count')]=NA
  outdf[,paste0(sp, 'Monophyletic')]=NA
  #brlenlist=list()
  
  for(i in (1:length(filenames))){ ## ep2 is not there
    
    awt=read.raxml(paste0('',filenames[i]))
    #tryCatch(("reroot" (
    #awt=root(awt, awt$tip.label[substr(awt$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate')]))), error = function(err) print("Outgroups already at the root"))
    awt=as.phylo(awt)
    awt$tip.label=gsub('_R_', '', awt$tip.label)
    ## add to make sure outgroup is there, and is monophyletic - otherwise skip this tree
    if(any(substr(as.phylo(awt)$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate'))){
      if(is.monophyletic(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate')])){
        awt=root(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate')])
        awt=drop.tip(awt, tip = as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('agerjg', 'eophiu')])
        vc=sum(substr(as.phylo(awt)$tip.label,1,6) %in% sp)
        outdf[i, paste0(sp, 'Count')]=vc
        if(vc>1){
          mrcanode=getMRCA(as.phylo(awt), as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c(sp)])
          getDescendants(as.phylo(awt), mrcanode)
          tips=as.phylo(awt)$tip.label[getDescendants(as.phylo(awt), mrcanode)]
          tips=tips[!is.na(tips)]
          outdf[i, paste0(sp, 'Monophyletic')]=all(substr(tips,1,6)==sp)
          ## get terminal branch length for each copy
          ## first get the node numbers of the tips
          nodes<-sapply(as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c(sp)],function(x,y) which(y==x),y=as.phylo(awt)$tip.label)
          ## then get the edge lengths for those nodes
          edge.lengths<-setNames(awt$edge.length[sapply(nodes,function(x,y) which(y==x),y=awt$edge[,2])],names(nodes))
          #brlenlist[i]=edge.lengths
        }}
    }
  }
}




melt(outdf[,-1])                                              
## 80 columns, take and combine
outdflist=lapply((1:((ncol(outdf)-1)/2))*2, function(x){ 
  temp=outdf[,c(x, x+1)]
  sp=substr(colnames(temp),1,6)[1]
  colnames(temp)=gsub(sp, '', colnames(temp))
  temp$genome=sp
  return(temp)
})
tp=do.call(rbind, outdflist)
##dplyr is being stupid, or i am, so just make columsn
tp$Polyphyletic=!tp$Monophyletic
tp$NotApplicable=is.na(tp$Monophyletic)
tpp=tp[!is.na(tp$Count),] %>% group_by(Count, genome) %>% dplyr::summarize(Monophyletic=sum(Monophyletic, na.rm=T), Polyphyletic=sum(Polyphyletic, na.rm=T), NotApplicable=sum(NotApplicable))
#write.table(tpp, '../tpp_summarized_syntenic_counts.txt', row.names=F, col.names=T, quote=F, sep='\t')

####### above was on cbsuxm01
###### now on this computer

tpp=read.table('tpp_summarized_syntenic_counts.txt', header = T)


tppp=reshape2::melt(tpp[!is.na(tpp$Count) & !tpp$genome %in% c('tdactm', 'tzopol', 'osativ', 'pprate', 'tdacs2', 'tdacn2', 'pvagin'),], id.vars=c('genome', 'Count'))
tppp$genome=factor(tppp$genome, levels=names(taxonnames))
tppp$ploidy=asize$ploidy[match(tppp$genome, asize$V2)]
tppp$phyleticcol=paste0(tppp$ploidy, tppp$variable)
## make singletons "monophyletic"
tppp$phyleticcol[tppp$Count==1]=paste0(tppp$ploidy[tppp$Count==1], 'Monophyletic')
tppp$species=taxonnames[match(tppp$genome, names(taxonnames))]
tppp$species=factor(tppp$species, levels=taxonnames)
tppp=tppp[!is.na(tppp$genome),]




tppp$haploid=asize$haploid[match(tppp$genome, asize$V2)]
tppp$doubledCount=ifelse(tppp$haploid, tppp$Count*2, tppp$Count)

tppp$speciesLabel=ifelse(tppp$haploid, paste0(tppp$species, '*'), as.character(tppp$species))
tppp$speciesLabel=tppp$species
levels(tppp$speciesLabel)[levels(tppp$speciesLabel) %in% tppp$species[tppp$haploid]]=paste0(levels(tppp$speciesLabel)[levels(tppp$speciesLabel) %in% tppp$species[tppp$haploid]], '*')

tppp$linetype=NA
tppp$linetype[tppp$doubledCount%in%1:6]=rep(c('dotted', 'dashed'),3)[tppp$doubledCount[tppp$doubledCount%in%1:6]]

## switch it back, not thinking this through!
tppp$variable[tppp$Count==1]='NotApplicable'                                



aap=tppp %>% group_by(speciesLabel, variable) %>% summarize(n=n(), count=sum(value)) %>% filter(variable!='NotApplicable') %>% mutate(pct = count/sum(count)*100, width=sum(count))%>%
  ggplot(aes(x=width/2, y=pct, fill=variable, width=width)) + geom_bar(stat='identity', position='fill') + coord_polar(theta='y') + facet_wrap(~speciesLabel, ncol=1, strip.position='left')+   theme( strip.background = element_blank(), strip.text.y.left = element_blank(), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + theme(axis.ticks=element_blank(), axis.text=element_blank(), panel.grid=element_blank(), panel.border=element_blank()) + scale_fill_manual(values=c('#5F4B8BFF', '#E69A8DFF')) + ylab('Proportion\nDuplicates\nMonophyletic') + xlab('')
aapp=tppp[tppp$ploidy!='Diploid',] %>% group_by(speciesLabel, variable) %>% summarize(n=n(), count=sum(value)) %>% filter(variable!='NotApplicable') %>% mutate(pct = count/sum(count)*100, width=sum(count))%>%
  ggplot(aes(x=width/2, y=pct, fill=variable, width=width)) + geom_bar(stat='identity', position='fill') + coord_polar(theta='y') + facet_wrap(~speciesLabel, ncol=1, strip.position='left')+   theme( strip.background = element_blank(), strip.text.y.left = element_blank(), panel.spacing = unit(3, "pt"), axis.text=element_text(size=9))+ theme(legend.position = "none") + theme(axis.ticks=element_blank(), axis.text=element_blank(), panel.grid=element_blank(), panel.border=element_blank()) + scale_fill_manual(values=c('#5F4B8BFF', '#E69A8DFF')) + ylab('Proportion\nDuplicates\nMonophyletic') + xlab('')


#cpb=ggplot(tppp[tppp$doubledCount%in% 1:6 & !is.na(tppp$species) & tppp$ploidy!='Diploid',], aes(x=doubledCount, y=value, group=ploidy, fill=ploidy)) + geom_hline(yintercept=0, color='snow2', linetype='dotted') + geom_vline(xintercept=c(1,3,5), color='snow2', linetype='dotted') + geom_vline(xintercept=c(2,4,6), color='snow3', linetype='dotted') + geom_histogram(stat='identity', position='stack') + 
#  facet_wrap(~speciesLabel, ncol=1, strip.position='left', drop=T) + scale_fill_manual(values=ploidycolors) +   theme( strip.background = element_blank(), strip.text.y.left = element_blank(), panel.spacing = unit(3, "pt"), axis.text.y=element_blank(), axis.text.x=element_text(size=9))+ theme(legend.position = "none", plot.margin = margin(l = -15)) + ylab('')+ xlab('Syntenic Gene\nCopy Number')+ scale_y_continuous(n.breaks = 3)


tppp$shade <- ifelse(tppp$doubledCount %in% c(1, 3, 5), "light", "dark")
tppp$haploidCount=ifelse(tppp$haploid, tppp$Count, ifelse(tppp$doubledCount%in%c(1,3,5), (tppp$doubledCount+1)/2, tppp$doubledCount/2))
#tppp$haploidCount[tppp$ploidy=='Paleotetraploid' & tppp$haploidCount==3]=
### if it is haploid, we need odd numbers to be counted in the bin smaller for zea
#### if it is allelic, we need odd numbers to be counted in the next bin up
### haploid tetraploids and paleotetraploids with 3
tppp$haploidCount[tppp$ploidy%in%c('Paleotetraploid', 'Tetraploid')&tppp$Count==3]=2
tppp$heterozygous=(tppp$doubledCount%in%c(1,3,5))
tppp$heterozygous[tppp$ploidy%in%c('Paleotetraploid', 'Tetraploid')&tppp$Count==3]=T


ploidycolors_lighter <- c(
  DiploidLight = "#FFE4A1",       # Lighter version of "#FFC857"
  TetraploidLight = "#D6CFF0",    # Lighter version of "#A997DF"
  HexaploidLight = "#F58A8F",     # Lighter version of "#E5323B"
  OctaploidLight = "#637180",     # Lighter version of "#2E4052"
  PaleotetraploidLight = "#CFE5EF" # Lighter version of "#97cddf"
)

tppp$colorbar=ifelse(tppp$heterozygous,  paste0(as.character(tppp$ploidy), 'Light'), as.character(tppp$ploidy))

# Modify the plot code
cpb <- ggplot(tppp[tppp$doubledCount %in% 1:6 & !is.na(tppp$species) & tppp$ploidy != 'Diploid', ], 
              aes(x = factor(haploidCount*2), y = value, fill = colorbar)) +
  geom_hline(yintercept = 0, color = 'snow2', linetype = 'dotted') +
  geom_vline(xintercept = c(1, 2, 3), color = 'snow2', linetype = 'dotted') +
  geom_histogram(stat = 'identity', position = 'stack') +
  facet_wrap(~speciesLabel, ncol = 1, strip.position = 'left', drop = TRUE) +
  scale_fill_manual(values = c(ploidycolors, ploidycolors_lighter)) + # Adjust the color scheme here
  theme(strip.background = element_blank(), 
        strip.text.y.left = element_blank(), #element_text(angle = 180), 
        panel.spacing = unit(3, "pt"), 
        axis.text.y = element_blank(), 
        axis.text.x = element_text(size = 9)) +
  theme(legend.position = "none")+#, plot.margin = margin(l = -15)) +
  ylab('') +
  xlab('Syntenic Gene\nCopy Number') +
  scale_y_continuous(n.breaks = 3)


shorttaxonnames=c("Z. mays ssp. parviglumis TIL11",  "Z. mays ssp. parviglumis TIL01", "Z. mays ssp. mays B73v5", "Z. mays ssp. mexicana TIL25", "Z. mays ssp. mexicana TIL18", "Z. mays ssp. huehuetenangensis", 
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
#final_merged_data$shortspeciesLabel=ifelse(final_merged_data$genome%in%asize$V2[asize$haploid], paste0(shorttaxonnames[match(final_merged_data$genome, names(shorttaxonnames))], '*'), shorttaxonnames[match(final_merged_data$genome, names(shorttaxonnames))])
final_merged_data$shortspeciesLabel=shorttaxonnames[match(final_merged_data$genome, names(shorttaxonnames))]
final_merged_data$shortspeciesLabel=factor(final_merged_data$shortspeciesLabel, levels=shorttaxonnames)
#final_merged_data$speciesLabel=factor(final_merged_data$speciesLabel, levels=unique(tppp$speciesLabel))
#levels(final_merged_data$speciesLabel)=levels(tppp$speciesLabel)
####### what i really need ot do is consider this in 100 gene windows or something so that entire chromososomes don't skew edians!!!!!!
ksplot=final_merged_data %>% filter(ploidy!='Diploid', !is.na(shortspeciesLabel)) %>% group_by(queryChr, refChr, paspChr, genome,ploidy, speciesLabel, shortspeciesLabel) %>%summarize(ks=mean(V17), n=length(unique(gene))) %>% filter(n>10) %>% ggplot()+ 
  geom_vline(xintercept=seq(from=0,to=0.25,by=0.01), color='gray', lty='dashed', alpha=0.2) + geom_histogram(aes(x=ks, color=ploidy, fill=ploidy, weight=n), binwidth=0.001) + 
  facet_wrap(~shortspeciesLabel, ncol=1, scales='free_y', strip.position = 'left', drop = TRUE, labeller=purrr::partial(label_species, dont_italicize=c('subsp.', 'ssp.', 'TIL11', 'TIL01', 'TIL25', 'TIL18', 'Momo', 'Gigi', 'Southern Hap1', 'Northern Hap1', 'FL', 'KS',  '\\*', '\\"', 'B73v5'))) + scale_color_manual(values=ploidycolors)+ scale_fill_manual(values=ploidycolors)+
  theme(strip.background = element_blank(), 
        strip.text.y.left = element_text(angle = 0, hjust=1), 
        panel.spacing = unit(3, "pt"), 
        axis.text.y = element_blank(), 
        axis.text.x = element_text(size = 9),
        legend.position = 'NULL',
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank()) +
  ylab(label = '')+
#  xlim(0,0.25)
  ## from gerardi plot
  scale_x_continuous( "Mean Ks of Syntenic Homologs per Block", limits=c(0,0.25),  sec.axis = sec_axis(~ . /6.5e-9/2/1e6, name = "Syntenic Homolog Divergence\n(million years)"))+
  theme(axis.text.x.top=element_text(color='blue'), axis.title.x.top=element_text(color='blue'))

## labels on right
ksplotright=final_merged_data %>% filter(ploidy!='Diploid', !is.na(shortspeciesLabel)) %>% group_by(queryChr, refChr, paspChr, genome,ploidy, speciesLabel, shortspeciesLabel) %>%summarize(ks=mean(V17), n=length(unique(gene))) %>% filter(n>10) %>% ggplot()+ 
  geom_vline(xintercept=seq(from=0,to=0.25,by=0.01), color='gray', lty='dashed', alpha=0.2) + geom_histogram(aes(x=ks, color=ploidy, fill=ploidy, weight=n), binwidth=0.001) + 
  facet_wrap(~shortspeciesLabel, ncol=1, scales='free_y', strip.position = 'right', drop = TRUE, labeller=purrr::partial(label_species, dont_italicize=c('subsp.', 'ssp.', 'TIL11', 'TIL01', 'TIL25', 'TIL18', 'Momo', 'Gigi', 'Southern Hap1', 'Northern Hap1', 'FL', 'KS',  '\\*', '\\"', 'B73v5'))) + scale_color_manual(values=ploidycolors)+ scale_fill_manual(values=ploidycolors)+
  theme(strip.background = element_blank(), 
        strip.text.y.right = element_text(angle = 0, hjust=0, size=9, vjust=0), 
        panel.spacing = unit(3, "pt"), 
        axis.text.y = element_blank(), 
        axis.text.x = element_text(size = 9),
        legend.position = 'NULL',
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank()) +
  ylab(label = '')+
#  xlim(0,0.25)
  ## from gerardi plot
  scale_x_continuous( "Mean Ks of Syntenic Homologs per Block", limits=c(0,0.25),  sec.axis = sec_axis(~ . /6.5e-9/2/1e6, name = "Syntenic Homolog Divergence\n(million years)"))+
  theme(axis.text.x.top=element_text(color='blue'), axis.title.x.top=element_text(color='blue'))


#purrr::partial(label_species, dont_italicize=c('subsp.', 'ssp.', 'TIL11', 'TIL01', 'TIL25', 'TIL18', 'Momo', 'Gigi', 'Southern Hap1', 'Northern Hap1', 'FL', 'KS',  '\\*', '\\"', 'B73v5'))

plot_grid(cpb + theme(axis.title.x=element_text(angle = 45, hjust = 0.8, vjust = 1)), NULL, aapp+ theme(axis.title.x=element_text(angle = 45, hjust = 1, vjust=1)), ncol=4, rel_widths = c(0.1,-0.5,1), align='hv')


plot_grid(ksplot, cpb ,  aapp + theme(axis.line.y=element_blank()) , ncol=3,  rel_widths=c(1,0.2,0.2), align='h')





##### sister clades

polysisdf=data.frame(tree=1:nrow(outdf))
for(sp in unique(tp$genome)){  
  polysisdf[,paste0(sp, 'Clade')]=NA
  polysisdf[,paste0(sp, 'CountinClade')]=NA
  polysisdf[,paste0(sp, 'Bootstrap')]=NA
}

for(sp in unique(tp$genome)){                                            
  #ccpolysis=list() ## don't know how to organize this - jsut need vector of what species are there for each tree? but what about when two tips, so unique for each tree?
  ## 
  #sp='ccitra'
  for(i in which(!is.na(outdf$achineCount))){ ## for each tree
    awtog=read.raxml(paste0('',filenames[i]))
    awt=as.phylo(awtog)
    awtog=awt
    rename=data.frame(og=awt$tip.label, new=gsub('_R_', '', awt$tip.label))
    awtog=rename_taxa(awtog, rename, og, new)
    
    
    if(length(as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate')])==1){
      awtog=root(awtog, outgroup=as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate')])
    }else{
      awtog=root(awtog, node=getMRCA(as.phylo(awt), as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate')]))
    }
    awt$tip.label=gsub('_R_', '', awt$tip.label)
    awt=root(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate')], resolve.root=T)
    
    awt=drop.tip(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('agerjg', 'eophiu', 'tzopol')])
    ## now they're rooted, do things
    
    if(sum(substr(as.phylo(awt)$tip.label,1,6) %in% sp)>1){
      mrcanode=getMRCA(as.phylo(awt), as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c(sp)])
      getDescendants(as.phylo(awt), mrcanode)
      tips=as.phylo(awt)$tip.label[getDescendants(as.phylo(awt), mrcanode)]
      tips=tips[!is.na(tips)]
      nonselftips=tips[substr(tips,1,6)!=sp]
      funnames=paste(length(nonselftips), paste(sort(unique(substr(nonselftips,1,6))), collapse = ' '))
      funnamesnl=paste(paste(sort(unique(substr(nonselftips,1,6))), collapse = ' '))
      #    mrcaBootstrap=awtog@data$bootstrap[mrcanode - Ntip(awt)]
      
      ##  if(length(unique(nonselftips))==1){ ## keep if in a clade with one other species!
      #ccpolysis[[i]]=unique(substr(nonselftips,1,6))
      if(length(nonselftips)>0){
        #    ccpolysis[[i]]=funnamesnl ## see what it looks like without counts - this is just tree representation
        polysisdf[i,paste0(sp, 'Clade')]=funnamesnl
        polysisdf[i,paste0(sp, 'CountinClade')]=length(nonselftips)
  #      polysisdf[i,paste0(sp, 'Bootstrap')]=mrcaBootstrap
      }
      if(length(nonselftips)==0){
        polysisdf[i,paste0(sp, 'Clade')]='Monophyletic'
        polysisdf[i,paste0(sp, 'CountinClade')]=0
        polysisdf[i,paste0(sp, 'Bootstrap')]=NA
      }
      ##}
    }## if more than one copy of the species
    if(sum(substr(as.phylo(awt)$tip.label,1,6) %in% sp)==1){ ## single tip
      polysisdf[i,paste0(sp, 'Clade')]='Single Copy'
      polysisdf[i,paste0(sp, 'CountinClade')]=0
      polysisdf[i,paste0(sp, 'Bootstrap')]=NA
    }
    if(sum(substr(as.phylo(awt)$tip.label,1,6) %in% sp)==0){
      polysisdf[i,paste0(sp, 'Clade')]=NA
      polysisdf[i,paste0(sp, 'CountinClade')]=NA
      polysisdf[i,paste0(sp, 'Bootstrap')]=NA
    }
    
  }
  
  #print(sp)
  #  funsp=table(unlist(ccpolysis))
  #print(sort(funsp[funsp>=10]              )      )       
  
  
}    



for(x in unique(tp$genome)){ 
  funsp=table(polysisdf[,paste0(x, 'Clade')])
  print(x)
  print(sort(funsp[funsp>=10]))
}

library(stringr)

polysisdf$paspalumgene=str_split_fixed(outdf$filenames, '\\.', 4)[,2]
polysisdf$paspalumchr=substr(polysisdf$paspalumgene,6,7)
head(polysisdf[polysisdf$blagurClade=='ttrian' & !is.na(polysisdf$blagurClade),]) ## eventually want to look for where on chr these are - is there a spatial pattern?



polysisdf %>% group_by(sister=gsub('tdacn1 | tdacn1|tdacn1|zdgigi | zdgigi|zdgigi|zdmomo | zdmomo|zdmomo|zluxur | zluxur|zluxur|zmhuet | zmhuet|zmhuet|zTIL18 | zTIL18|zTIL18|zTIL25 | zTIL25|zTIL25|zTIL01 | zTIL01|zTIL01|zTIL11 | zTIL11|zTIL11|znicar | znicar|znicar|tdacs2 | tdacs2|tdacs2|tdacn2 | tdacn2|tdacn2|zmB735 | zmB735|zmB735|tdactm | tdactm|tdactm|tzopol | tzopol|tzopol', '', tdacs1Clade)) %>% summarize(n=n()) %>% filter(!is.na(sister), sister!='') %>% top_n(5) %>% arrange(n)
polysisdf %>% group_by(sister=gsub('', '', sscopaClade)) %>% summarize(n=n()) %>% filter(!is.na(sister)) %>% top_n(5) %>% arrange(n) 
polysisdf %>% group_by(sister=gsub('', '', achineClade)) %>% summarize(n=n()) %>% filter(!is.na(sister)) %>% top_n(5) %>% arrange(n)
polysisdf %>% group_by(sister=gsub('', '', agerarClade)) %>% summarize(n=n()) %>% filter(!is.na(sister)) %>% top_n(5) %>% arrange(n)
polysisdf %>% group_by(sister=gsub('', '', ccitraClade)) %>% summarize(n=n()) %>% filter(!is.na(sister)) %>% top_n(5) %>% arrange(n)
polysisdf %>% group_by(sister=gsub('', '', blagurClade)) %>% summarize(n=n()) %>% filter(!is.na(sister)) %>% top_n(5) %>% arrange(n)

# List of clades to process
clades <- c("tdacs1Clade", "sscopaClade", "achineClade", "agerarClade", "ccitraClade", "blagurClade")

# Initialize an empty data frame to store results
results_df <- data.frame()

# Loop over the clades and process each one
for (clade in clades) {
  # Remove specific substrings from the clade (modify gsub pattern as needed)
  # summary_df <- polysisdf %>%
  #   group_by(sister = gsub('tdacn1 | tdacn1|tdacn1|zdgigi | zdgigi|zdgigi|zdmomo | zdmomo|zdmomo|zluxur | zluxur|zluxur|zmhuet | zmhuet|zmhuet|zTIL18 | zTIL18|zTIL18|zTIL25 | zTIL25|zTIL25|zTIL01 | zTIL01|zTIL01|zTIL11 | zTIL11|zTIL11|znicar | znicar|znicar|tdacs2 | tdacs2|tdacs2|tdacn2 | tdacn2|tdacn2|zmB735 | zmB735|zmB735|tdactm | tdactm|tdactm|tzopol | tzopol|tzopol', '', get(clade))) %>%
  #   summarize(n = n()) %>%
  #   filter(!is.na(sister), sister != '', sister!='Single Copy') %>% ## second parts are to clean up tdacs1 because of zea/trip sharing
  #   top_n(5) %>%
  #   arrange(desc(n))
  summary_df <- polysisdf %>%
    group_by(sister = gsub('tdacn1 | tdacn1|tdacn1|zdgigi | zdgigi|zdgigi|zdmomo | zdmomo|zdmomo|zluxur | zluxur|zluxur|zmhuet | zmhuet|zmhuet|zTIL18 | zTIL18|zTIL18|zTIL25 | zTIL25|zTIL25|zTIL01 | zTIL01|zTIL01|zTIL11 | zTIL11|zTIL11|znicar | znicar|znicar|tdacs2 | tdacs2|tdacs2|tdacn2 | tdacn2|tdacn2|zmB735 | zmB735|zmB735|tdactm | tdactm|tdactm|tzopol | tzopol|tzopol', '', get(clade))) %>%
    mutate(sister = ifelse(sister == '', 'Monophyletic', sister)) %>%  # Rename empty sister to 'Monophyletic'
    summarize(n = n()) %>%
    filter(!is.na(sister), sister!='Single Copy') %>%
    top_n(5) %>%
    arrange(desc(n))
  
  # Add a column with the clade name
  summary_df <- summary_df %>%
    mutate(clade = clade)
  
  # Combine the results for this clade with the previous results
  results_df <- bind_rows(results_df, summary_df)
}

# Display or save the results
print(results_df)


# Define a color palette for non-Monophyletic species
color_palette <- colorRampPalette(c("#dd7a65", "#efc2b9"))  # Light to dark shades of #E69A8DFF

# Modify the dataframe to assign a fixed color for 'Monophyletic' and generate shades for others
results_df <- results_df %>%
  mutate(
    sister = ifelse(sister == '', 'Monophyletic', sister),  # Rename '' to 'Monophyletic'
    color = ifelse(sister == 'Monophyletic', '#5F4B8BFF', NA)  # Assign color to 'Monophyletic'
  ) %>%
  group_by(clade) %>%
  mutate(
    rank = rank(-n),  # Rank by count within each clade, with highest n getting rank 1
    shade = rank / max(rank),  # Normalize rank to get values between 0 and 1
    color = ifelse(is.na(color), color_palette(max(rank))[rank], color)  # Assign shades to non-Monophyletic
  )


# Create a function to replace taxon codes with their corresponding names from shorttaxonnames
replace_taxon_names <- function(sister_string, shorttaxonnames) {
  # Split the string by spaces into individual taxon codes
  taxon_codes <- unlist(strsplit(sister_string, " "))
  
  # Replace each taxon code with its corresponding name, if it exists in shorttaxonnames
  replaced_names <- sapply(taxon_codes, function(code) {
    if (code %in% names(shorttaxonnames)) {
      return(shorttaxonnames[[code]])
    } else {
      return(code)  # If the code is not found, keep the original
    }
  })
  if(length(replaced_names)==23){
    replaced_names="All taxa"
  }
  if(length(replaced_names)==22){ ## not generalizable!!!!! htis is just this stupid tripsacinae for now
    replaced_names='All taxa except H. compressa'
  }
  # Join the replaced names with '\n'
  return(paste(replaced_names, collapse = "\n"))
}

# Apply this function to each row in results_df$sister
results_df$sisterSpecies <- sapply(results_df$sister, replace_taxon_names, shorttaxonnames)


# Function to italicize genus.species parts of the labels
italicize_taxon <- function(label) {

  # Apply italic only to genus and species parts, leave other parts unchanged
  if (substr(label,1,3)!='All' & label!='Monophyletic') {
    return(bquote(italic(.(label))))
  } else {
    return(bquote(.(paste0('"', label, '"'))))  # If no genus/species parts, return as is
  }
}
## quote because of the color strings with #
#write.table(results_df, '~/transfer/polyphyletic_clades_counts.txt', row.names=F, col.names = T, sep='\t', quote=T)

results_df=fread('~/Downloads/polyphyletic_clades_counts.txt', header=T, sep='\t')

# Apply the italicizing function to the labels
#results_df$italicized_sister = sapply(results_df$sisterSpecies, italicize_taxon))
italicized_sister_list <- lapply(results_df$sisterSpecies, italicize_taxon)


fixed_label_position <- 300  # Adjust this value to the fixed distance you want
results_df <- results_df %>%
  mutate(label_position = n + fixed_label_position)  # Add fixed value above the bars

results_df$num_lines <- sapply(results_df$sisterSpecies, function(x) length(strsplit(x, "\n")[[1]]))
results_df$max_line_length <- sapply(results_df$sisterSpecies, function(x) max(nchar(strsplit(x, "\n")[[1]])))

# Normalize hjust based on the max_line_length (e.g., divide by max line length)
# You may adjust the normalization factor depending on the font and size.
results_df$hjust <- 1- 1 / (2*results_df$num_lines)


pdf('~/transfer/sister_clade_barplot.pdf',8,15)
ggplot(results_df, aes(x = reorder_within(sister, -n, clade), y = n, fill = color)) +
  geom_bar(stat = "identity") +         # Create bar plot
  facet_wrap(~ clade, scales = "free", ncol = 3) + # Create facets by clade
  scale_fill_identity() +                # Use manual fill colors
  scale_x_reordered() +                  # Reorder x-axis labels within facets
  labs(x = "Sister", y = "Count", title = "Sister Count per Clade") +
  theme_minimal() +                      # Use minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels
        legend.position = "none")        # Hide legend
ggplot(results_df, aes(x = reorder_within(sisterSpecies, -n, clade), y = n, fill = color)) +
  geom_bar(stat = "identity") +         # Create bar plot
  geom_text(aes(y = label_position, label = italicized_sister_list),        # Add labels on top of bars
 #           position = position_stack(vjust = 1.05), # Position text above the bars
            size = 2, angle = 90, ,vjust = results_df$num_lines/2, lineheight=0.7, parse=T,position = position_nudge(x = 0)) + 
  facet_wrap(~ clade, scales = "free_x", ncol = 3) + # Create facets by clade
  scale_fill_identity() +                # Use manual fill colors
  scale_x_reordered() +                  # Reorder x-axis labels within facets
  labs(x = "Sister", y = "Count", title = "Sister Count per Clade") +
  theme_minimal() +                      # Use minimal theme
  theme(axis.text.x = element_blank(), #element_text(angle = 90, hjust = 0.5),  # Rotate x-axis labels
        legend.position = "none")        # Hide legend


ggplot(results_df[results_df$sister!='Monophyletic',], aes(x = reorder_within(sister, -n, clade), y = n, fill = color)) +
  geom_bar(stat = "identity") +         # Create bar plot
  facet_wrap(~ clade, scales = "free", ncol = 3) + # Create facets by clade
  scale_fill_identity() +                # Use manual fill colors
  scale_x_reordered() +                  # Reorder x-axis labels within facets
  labs(x = "Sister", y = "Count", title = "Sister Count per Clade") +
  theme_minimal() +                      # Use minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels
        legend.position = "none")        # Hide legend


dev.off()


results_df$label_position_adj=results_df$label_position
results_df$label_position=2000

achinebp=ggplot(results_df[results_df$clade=='achineClade',], aes(x = reorder_within(sisterSpecies, -n, clade), y = n, fill = color)) +
  geom_bar(stat = "identity") +         # Create bar plot
  geom_text(aes(y = label_position, label = italicized_sister_list[results_df$clade=='achineClade']),        # Add labels on top of bars
            #           position = position_stack(vjust = 1.05), # Position text above the bars
            size = 2, angle = 90, vjust = results_df$num_lines[results_df$clade=='achineClade']/2, lineheight=0.7, parse=T,position = position_nudge(x = 0)) + 
#  facet_wrap(~ clade, scales = "free_x", ncol = 3) + # Create facets by clade
  scale_fill_identity() +                # Use manual fill colors
  scale_x_reordered() +                  # Reorder x-axis labels within facets
  labs(x = "A. chinensis Duplicates Clade", y = "Count") +
  theme(axis.text.x = element_blank(), #element_text(angle = 90, hjust = 0.5),  # Rotate x-axis labels
        legend.position = "none") +ylim(0,6200)       # Hide legend
achinebp

tripbp=ggplot(results_df[results_df$clade=='tdacs1Clade',], aes(x = reorder_within(sisterSpecies, -n, clade), y = n, fill = color)) +
  geom_bar(stat = "identity") +         # Create bar plot
  geom_text(aes(y = label_position, label = italicized_sister_list[results_df$clade=='tdacs1Clade']),        # Add labels on top of bars
            #           position = position_stack(vjust = 1.05), # Position text above the bars
            size = 2, angle = 90, vjust = results_df$num_lines[results_df$clade=='tdacs1Clade']/2, lineheight=0.7, parse=T,position = position_nudge(x = 0)) + 
  #  facet_wrap(~ clade, scales = "free_x", ncol = 3) + # Create facets by clade
  scale_fill_identity() +                # Use manual fill colors
  scale_x_reordered() +                  # Reorder x-axis labels within facets
  labs(x = "T. dactyloides Duplicates Clade", y = "Count") +
  theme(axis.text.x = element_blank(), #element_text(angle = 90, hjust = 0.5),  # Rotate x-axis labels
        legend.position = "none")  +ylim(0,6200)       # Hide legend
tripbp

sscopabp=ggplot(results_df[results_df$clade=='sscopaClade',], aes(x = reorder_within(sisterSpecies, -n, clade), y = n, fill = color)) +
  geom_bar(stat = "identity") +         # Create bar plot
  geom_text(aes(y = label_position, label = italicized_sister_list[results_df$clade=='sscopaClade']),        # Add labels on top of bars
            #           position = position_stack(vjust = 1.05), # Position text above the bars
            size = 2, angle = 90, vjust = results_df$num_lines[results_df$clade=='sscopaClade']/2, lineheight=0.7, parse=T,position = position_nudge(x = 0)) + 
  #  facet_wrap(~ clade, scales = "free_x", ncol = 3) + # Create facets by clade
  scale_fill_identity() +                # Use manual fill colors
  scale_x_reordered() +                  # Reorder x-axis labels within facets
  labs(x = "S. scoparium Duplicates Clade", y = "Count") +
  theme(axis.text.x = element_blank(), #element_text(angle = 90, hjust = 0.5),  # Rotate x-axis labels
        legend.position = "none")  +ylim(0,6200)       # Hide legend
sscopabp

agerarbp=ggplot(results_df[results_df$clade=='agerarClade',], aes(x = reorder_within(sisterSpecies, -n, clade), y = n, fill = color)) +
  geom_bar(stat = "identity") +         # Create bar plot
  geom_text(aes(y = label_position, label = italicized_sister_list[results_df$clade=='agerarClade']),        # Add labels on top of bars
            #           position = position_stack(vjust = 1.05), # Position text above the bars
            size = 2, angle = 90, vjust = results_df$num_lines[results_df$clade=='agerarClade']/2, lineheight=0.7, parse=T,position = position_nudge(x = 0)) + 
  #  facet_wrap(~ clade, scales = "free_x", ncol = 3) + # Create facets by clade
  scale_fill_identity() +                # Use manual fill colors
  scale_x_reordered() +                  # Reorder x-axis labels within facets
  labs(x = "A. gerardi Duplicates Clade", y = "Count") +
  theme(axis.text.x = element_blank(), #element_text(angle = 90, hjust = 0.5),  # Rotate x-axis labels
        legend.position = "none")  +ylim(0,6200)       # Hide legend
agerarbp

ccitrabp=ggplot(results_df[results_df$clade=='ccitraClade',], aes(x = reorder_within(sisterSpecies, -n, clade), y = n, fill = color)) +
  geom_bar(stat = "identity") +         # Create bar plot
  geom_text(aes(y = label_position, label = italicized_sister_list[results_df$clade=='ccitraClade']),        # Add labels on top of bars
            #           position = position_stack(vjust = 1.05), # Position text above the bars
            size = 2, angle = 90, vjust = results_df$num_lines[results_df$clade=='ccitraClade']/2, lineheight=0.7, parse=T,position = position_nudge(x = 0)) + 
  #  facet_wrap(~ clade, scales = "free_x", ncol = 3) + # Create facets by clade
  scale_fill_identity() +                # Use manual fill colors
  scale_x_reordered() +                  # Reorder x-axis labels within facets
  labs(x = "C. citratus Duplicates Clade", y = "Count") +
  theme(axis.text.x = element_blank(), #element_text(angle = 90, hjust = 0.5),  # Rotate x-axis labels
        legend.position = "none")  +ylim(0,6200)       # Hide legend
ccitrabp

blagurbp=ggplot(results_df[results_df$clade=='blagurClade',], aes(x = reorder_within(sisterSpecies, -n, clade), y = n, fill = color)) +
  geom_bar(stat = "identity") +         # Create bar plot
  geom_text(aes(y = label_position, label = italicized_sister_list[results_df$clade=='blagurClade']),        # Add labels on top of bars
            #           position = position_stack(vjust = 1.05), # Position text above the bars
            size = 2, angle = 90, vjust = results_df$num_lines[results_df$clade=='blagurClade']/2, lineheight=0.7, parse=T,position = position_nudge(x = 0)) + 
  #  facet_wrap(~ clade, scales = "free_x", ncol = 3) + # Create facets by clade
  scale_fill_identity() +                # Use manual fill colors
  scale_x_reordered() +                  # Reorder x-axis labels within facets
  labs(x = "B. laguroides Duplicates Clade", y = "Count") +
  theme(axis.text.x = element_blank(), #element_text(angle = 90, hjust = 0.5),  # Rotate x-axis labels
        legend.position = "none")  +ylim(0,6200)       # Hide legend
blagurbp


allsisters=plot_grid(tripbp, sscopabp, achinebp, agerarbp, ccitrabp, blagurbp, ncol=1, labels=c('D','E','F','G','H', 'I'), align='hv')
allsisters

plot_grid(ksplot, cpb ,  aapp + theme(axis.line.y=element_blank()) , allsisters, ncol=4,  rel_widths=c(1,0.2,0.2,0.5), align='h', labels='AUTO')

  
aapp=ggplot(tapp, aes(x=width/2, y=pct, fill=variable, width=width)) +
  geom_bar(stat='identity', position='fill') +
  coord_polar(theta='y') +
  facet_wrap(~speciesLabel, ncol=1, strip.position='left') +
  theme(strip.background = element_blank(),
        strip.text.y.left = element_blank(),
        panel.spacing = unit(3, "pt"),
        axis.text = element_text(size=9), legend.text = element_text(size = 8)) +
  guides(fill = guide_legend(title = "", nrow=2, override.aes = list(size = 0.5))) +
  theme(legend.position = c(0.2, 1.1), legend.justification = c(0,1),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank()) + #, axis.line.y=element_blank()) +
  scale_fill_manual(values = c('#5F4B8BFF', '#E69A8DFF')) +
  ylab('Proportion\nDuplicates\nMonophyletic') +
  xlab('') #+
  # geom_text(aes(label = label, x = width * 1.1), # Place text off to the side
  #           size = 5, color = 'black', fontface = 'bold', hjust = 0) +
  # geom_rect(data = tapp %>% filter(highlight), # Only draw boxes for highlighted pies
  #           aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 100), 
  #           fill = NA, color = 'black', size = 1) +  # Black box with no fill
#  geom_label(aes(x = x_label, y = y_label, label = highlight), size = 5, hjust = 0,
#             fill = 'white', color = 'black', label.padding = unit(0.25, "lines"))



plot_grid(ksplot, cpb ,  aapp + theme(axis.line.y=element_blank()) , allsisters, ncol=4,  rel_widths=c(1,0.2,0.2,0.5), align='h', axis = 'tb',labels=c('A', 'B', 'C', ''))

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
  geom_point(size=3, alpha=0.8)+
  scale_color_manual(values=ploidycolors)+theme(legend.position=c(0.75,0.9), legend.text = element_text(size=8)) +
  ylab('Ks between duplicates') + xlab('Heterozygosity') + labs(color='Ploidy', pch='') +guides(color='none')+
  theme() # Remove y-axis and other clutter

ae=read.csv('~/Downloads/ranges_panand_aubuchonelder2023.csv')
ae$het=het$het[match(ae$V2, het$genome)]
ae$ploidy=het$ploidy[match(ae$V2, het$genome)]
hetrange=ggplot(ae, aes(x=het, y=as.numeric(ConR.AOO.km2), color=ploidy, label=V2)) +
  geom_point(size=2, alpha=0.8)+ geom_text()+
  scale_color_manual(values=ploidycolors)+theme(legend.position='none') +
  ylab('Area of Occupancy (km2)') + xlab('Heterozygosity') + 
  theme()  # Remove y-axis and other clutter

allsisters=plot_grid(tripbp, sscopabp, achinebp, agerarbp, ccitrabp, blagurbp, ncol=1, labels=c('E','F','G','H','I','J'), align='hv')

#pdf('~/Downloads/polyploidy_prevalence.pdf',10,10)
plot_grid(plot_grid(heterozygosity, hethomeo, axis='tb', align='h', rel_widths=c(1,0.5), labels=c('A','B')), 
plot_grid(ksplot, 
          #cpb ,  
          aapp + theme(axis.line.y=element_blank()) , allsisters, ncol=3,  rel_widths=c(1,0.2,0.5), align='h', axis = 'tb',labels=c('C', 'D', '')),
ncol=1, rel_heights = c(0.35,1), align='v', axis='rl')
#dev.off()




mean(
  median(homeol$ks[homeol$ploidy=='Paleotetraploid'& substr(homeol$genome,1,1)=='z'], na.rm=T),
  median(homeol$ks[homeol$ploidy=='Paleotetraploid'& substr(homeol$genome,1,3)=='tda'], na.rm=T)
)
mean(
median(homeol$ks[homeol$ploidy=='Paleotetraploid'& substr(homeol$genome,1,1)=='z'], na.rm=T)/6.5e-9/2/1e6,
median(homeol$ks[homeol$ploidy=='Paleotetraploid'& substr(homeol$genome,1,3)=='tda'], na.rm=T)/6.5e-9/2/1e6
)



#### now prepare trees for just gerardi parentage
gerardiclade=c('avirgi', 'achine', 'agerar', 'sscopa', 'smicro', 'pvagin')
bas=sapply(1:length(filenames), function(i){ ## ep2 is not there
  
  awt=read.raxml(paste0('',filenames[i]))
  awt=as.phylo(awt)
  awt$tip.label=gsub('_R_', '', awt$tip.label)
  ## add to make sure outgroup is there, and is monophyletic - otherwise skip this tree
  if(any(substr(as.phylo(awt)$tip.label,1,6) %in% c('pvagin'))){
 #   if(is.monophyletic(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate', 'pvagin', 'Pavag0', 'Pavag1')])){
      awt=root(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('pvagin')], resolve.root=T)
      
      ## now keep only six digit code
      awt$tip.label=paste0(awt$tip.label, '_', substr(awt$tip.label,1,6))
      sixtips=substr(awt$tip.label,1,6)
      #awt$tip.label[grepl('Pavag', sixtips)]=paste0(awt$tip.label[grepl('Pavag', sixtips)], '_','paspal')

   #   print(i)
      if(any(gerardiclade%in%sixtips)){
      awt=keep.tip(awt, awt$tip.label[sixtips %in% gerardiclade])
      #awt$tip.label=paste0(awt$tip.label, 1:length(awt$tip.label),'_', awt$tip.label) ## grampa usese the species name after the _species, so just fake this
      return(as.phylo(awt))
      }
    }
  }
})

bas <- bas[!sapply(bas,is.null)]
d=do.call("c",bas)

write.tree(d, paste0('gerardiclade_anchors_grampa.', Sys.Date(), '.tre'))

rspsix=read.tree(text='((((((((((((zTIL11:1.0,zmB735:1.0):0.00399,zTIL01:1.0):0.123008,(zTIL25:1.0,zTIL18:1.0):0.051481):0.078959,zmhuet:1.0):0.315942,(zluxur:1.0,znicar:1.0):1.343):0.035313,(zdmomo:1.0,zdgigi:1.0):1.085725):4.966642,(((((tdacs1:1.0,tdacs2:1.0):0.104529,tdacn2:1.0):0.344325,tdacn1:1.0):0.085028,tdactm:1.0):3.804823,tzopol:1.0):0.163575):3.340689,((udigit:1.0,vcuspi:1.0):0.130867,rrottb:1.0):0.2067):0.532529,((rtuber:1.0,hcompr:1.0):1.047055,etrips:1.0):0.311621):0.382165,(((((((((((sscopa:1.0,smicro:1.0):0.375104,avirgi:1.0):0.464688,(achine:1.0,agerar:1.0):0.197343):2.589164,(crefra:1.0,ccitra:1.0):4.376438):0.274712,((hconto:1.0,ttrian:1.0):0.601483,blagur:1.0):0.896994):0.456651,ppanic:1.0):0.093559,sbicol:1.0):0.131326,irugos:1.0):0.006381,snutan:1.0):0.144602,atenui:1.0):0.524102,(telega:1.0,cserru:1.0):0.446963):0.337256):2.735692,((bdista:1.0,osativ:1.0):8.038372,svirid:1.0):0.322258):1.0,pvagin:1.0);')
rspger=keep.tip(rspsix, gerardiclade)
write.tree(rspger, paste0('gerardiclade_sptree_grampa.', Sys.Date(), '.tre'))





### for polyphest
#### now prepare trees for just gerardi parentage
gerardiclade=c('avirgi', 'achine', 'agerar', 'sscopa', 'smicro', 'pvagin')
bas=sapply(1:length(filenames), function(i){ ## ep2 is not there
  
  awt=read.raxml(paste0('',filenames[i]))
  awt=as.phylo(awt)
  awt$tip.label=gsub('_R_', '', awt$tip.label)
  ## add to make sure outgroup is there, and is monophyletic - otherwise skip this tree
  if(any(substr(as.phylo(awt)$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate', 'pvagin', 'Pavag0', 'Pavag1'))){
    if(is.monophyletic(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate', 'pvagin', 'Pavag0', 'Pavag1')])){
      awt=root(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate', 'pvagin', 'Pavag0', 'Pavag1')], resolve.root=T)
      

      sixtips=substr(awt$tip.label,1,6)
      #awt$tip.label[grepl('Pavag', sixtips)]=paste0(awt$tip.label[grepl('Pavag', sixtips)], '_','paspal')
      
      #   print(i)
 #     if(any(gerardiclade%in%sixtips)){
        awt=keep.tip(awt, awt$tip.label[sixtips %in% c(all$V2[!all$V2%in%c('tdacn2', 'tdacs2')], 'pvagin')])
        ## now keep only six digit code
        awt$tip.label=substr(awt$tip.label,1,6)
        #awt$tip.label=paste0(awt$tip.label, 1:length(awt$tip.label),'_', awt$tip.label) ## grampa usese the species name after the _species, so just fake this
        return(as.phylo(awt))
 #     }
    }
  }
})

bas <- bas[!sapply(bas,is.null)]
d=do.call("c",bas)

write.tree(d, paste0('syntenic_anchors_polyphest.', Sys.Date(), '.tre'))



### for polyphest - gerardi
#### now prepare trees for just gerardi parentage
gerardiclade=c('avirgi', 'achine', 'agerar', 'sscopa', 'smicro', 'pvagin')
bas=lapply(1:length(filenames), function(i){ ## ep2 is not there
  
  awt=read.raxml(paste0('',filenames[i]))
  awt=as.phylo(awt)
  awt$tip.label=gsub('_R_', '', awt$tip.label)
  ## add to make sure outgroup is there, and is monophyletic - otherwise skip this tree
  if(any(substr(as.phylo(awt)$tip.label,1,6) %in% c('pvagin'))){
 #   if(is.monophyletic(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate', 'pvagin', 'Pavag0', 'Pavag1')])){
 #     awt=root(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate', 'pvagin', 'Pavag0', 'Pavag1')], resolve.root=T)
      awt=root(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('pvagin')], resolve.root=T)
      
      
      sixtips=substr(awt$tip.label,1,6)
      #awt$tip.label[grepl('Pavag', sixtips)]=paste0(awt$tip.label[grepl('Pavag', sixtips)], '_','paspal')
      
      #   print(i)
      #     if(any(gerardiclade%in%sixtips)){
      awt=keep.tip(awt, awt$tip.label[sixtips %in% c(gerardiclade, 'pvagin')])
      ## now keep only six digit code
      awt$tip.label=substr(awt$tip.label,1,6)
      #awt$tip.label=paste0(awt$tip.label, 1:length(awt$tip.label),'_', awt$tip.label) ## grampa usese the species name after the _species, so just fake this
      return(as.phylo(awt))
      #     }
    }
#  }
})

bas <- bas[!sapply(bas,is.null)]
d=do.call("c",bas)

write.tree(d, paste0('gerardi_anchors_polyphest.', Sys.Date(), '.tre'))





#### udigit subgenomes - who is sister?

#udsg=read.table('~/Downloads/udigitk17_q50_f2.chrom-subgenome.tsv', header=F)
udsg=read.table('udigitk17_q50_f2.chrom-subgenome.tsv', header=F)

gerardiclade=c('avirgi', 'achine', 'agerar', 'sscopa', 'smicro', 'pvagin')



rotate_clades <- function(tree) {
  # Helper function to check if a node is a leaf
  is.leaf <- function(node) {
    return(!(node %in% tree$edge[, 1]))
  }
  
  # Helper function to get all tips of a clade
  get_tips <- function(tree, node) {
    if (is.leaf(node)) {
      return(tree$tip.label[node])
    } else {
      children <- tree$edge[tree$edge[, 1] == node, 2]
      return(unlist(sapply(children, function(x) get_tips(tree, x))))
    }
  }
  
  # Function to recursively sort the labels in a clade
  sort_clade <- function(node) {
    if (is.leaf(node)) {
      return(node)
    } else {
      children <- which(tree$edge[, 1] == node)
      child_nodes <- tree$edge[children, 2]
      
      # Recursively sort the child nodes
      sorted_children <- sapply(child_nodes, sort_clade)
      
      # Sort the child nodes alphabetically based on their tips
      tips <- sapply(sorted_children, function(x) min(get_tips(tree, x)))
      sorted_order <- order(tips)
      
      tree$edge[children, 2] <<- sorted_children[sorted_order]
      return(node)
    }
  }
  
  # Start sorting from the root
  root <- Ntip(tree) + 1
  sort_clade(root)
  return(tree)
}

udtree=lapply(1:length(filenames), function(i){ ## ep2 is not there
  
  awt=read.raxml(paste0('',filenames[i]))
  awt=as.phylo(awt)
  awt$tip.label=gsub('_R_', '', awt$tip.label)
  ## add to make sure outgroup is there, and is monophyletic - otherwise skip this tree
  if(any(substr(as.phylo(awt)$tip.label,1,6) %in% c('pvagin')) & any(substr(as.phylo(awt)$tip.label,1,6)=='udigit')){
    #   if(is.monophyletic(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate', 'pvagin', 'Pavag0', 'Pavag1')])){
    #     awt=root(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('osativ', 'bdista', 'pprate', 'pvagin', 'Pavag0', 'Pavag1')], resolve.root=T)
    awt=root(awt, as.phylo(awt)$tip.label[substr(as.phylo(awt)$tip.label,1,6) %in% c('pvagin')], resolve.root=T)
    
    
    sixtips=substr(awt$tip.label,1,6)
    #awt$tip.label[grepl('Pavag', sixtips)]=paste0(awt$tip.label[grepl('Pavag', sixtips)], '_','paspal')
    
    #   print(i)
    #     if(any(gerardiclade%in%sixtips)){
    awt=keep.tip(awt, awt$tip.label[sixtips %in% c('udigit', 'pvagin')])
    ut=awt
    ### get scaffold, so i can assign subgenome!!!
    ## udigit_Pavag01G000200_scaf_2_215007925-215011286
    ## now keep only six digit code
    awt$tip.label[substr(awt$tip.label,1,6)=='udigit']=udsg$V2[match(paste0('scaf_', str_split_fixed(awt$tip.label[substr(awt$tip.label,1,6)=='udigit'], '_',5)[,4]), udsg$V1)]
    awt$tip.label[substr(awt$tip.label,1,6)=='pvagin']='pvagin'
    #awt$tip.label=paste0(awt$tip.label, 1:length(awt$tip.label),'_', awt$tip.label) ## grampa usese the species name after the _species, so just fake this
    stripped_tree <- awt
    stripped_tree$edge.length <- NULL
    stripped_tree$node.label <- NULL
    
    sorted_tree <- rotate_clades(stripped_tree)
    topology <- write.tree(sorted_tree)
    
    ## now store 
    chr=paste0('scaf_', str_split_fixed(ut$tip.label[substr(ut$tip.label,1,6)=='udigit'], '_',5)[,4])
    startend=str_split_fixed(ut$tip.label[substr(ut$tip.label,1,6)=='udigit'], '_',5)[,5]
    start=str_split_fixed(startend, '-', 2)[,1]
    end=str_split_fixed(startend, '-', 2)[,2]
    pavag=str_split_fixed(ut$tip.label[substr(ut$tip.label,1,6)=='udigit'],'_',3)[,2]
    
     return(data.frame(tix=i, pavag=pavag, chr=chr, start=start, end=end, topology=topology))
    #     }
  }
  #  }
})

#sort(table(unlist(udtree)))

## so major topology is ((2,3),1) but ((1,2),3) is also really common, and other alt is also - where are these along chromsomes?


udig=do.call(rbind, udtree)


halloween_palette <- c("#FF7518",  # Pumpkin orange
                       "#4B0082",  # Dark purple (witch's robe)
                       "#000000",  # Black (night)
                       "#8A9A5B",  # Eerie green
                       "#FF6347",  # Blood red
                       "#FFD700")  # Gold (moonlight)


udig$pchr=substr(udig$pavag,6,7)
udig$subgenome=udsg$V2[match(udig$chr, udsg$V1)]

pdf('~/transfer/udig_topologies.pdf',20,20) 
ggplot(udig[udig$topology%in%names(tail(sort(table(udig$topology)))),], aes(x=as.numeric(start), fill=topology)) + geom_histogram(binwidth=10000000, position='fill')+ facet_wrap(pchr~paste(chr, subgenome), ncol=3) + scale_fill_manual(values=halloween_palette)
ggplot(udig[udig$topology%in%names(tail(sort(table(udig$topology)))),], aes(x=as.numeric(start), fill=topology)) + geom_histogram(binwidth=10000000)+ facet_wrap(pchr~paste(chr, subgenome), ncol=3) + scale_fill_manual(values=halloween_palette)
ggplot(udig[udig$topology%in%names(tail(sort(table(udig$topology)))),], aes(x=as.numeric(start), fill=topology)) + geom_histogram(binwidth=1000000, position='fill')+ facet_wrap(pchr~paste(chr, subgenome), ncol=3) + scale_fill_manual(values=halloween_palette)
dev.off()



pdf('~/transfer/udig_topologies.pdf',20,20) 
ggplot(udig[udig$topology%in%names(tail(sort(table(udig$topology)))),], aes(x=as.numeric(start), fill=topology)) + geom_histogram(binwidth=10000000, position='fill')+ facet_wrap(pchr~paste(chr, subgenome), ncol=3) + scale_fill_manual(values=halloween_palette)
ggplot(udig[udig$topology%in%names(tail(sort(table(udig$topology)))),], aes(x=as.numeric(start), fill=topology)) + geom_histogram(binwidth=10000000)+ facet_wrap(pchr~paste(chr, subgenome), ncol=3) + scale_fill_manual(values=halloween_palette)
ggplot(udig[udig$topology%in%names(tail(sort(table(udig$topology)))),], aes(x=as.numeric(start), fill=topology)) + geom_histogram(binwidth=1000000, position='fill')+ facet_wrap(pchr~paste(chr, subgenome), ncol=3) + scale_fill_manual(values=halloween_palette)
dev.off()




