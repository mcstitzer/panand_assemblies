

## XXXX NEVERMIND NOT WORKING on biohpc: singularity run --bind $PWD --pwd $PWD /programs/genespace-1.2.3/genespace.sif R
# conda activate genespace - to get the orthofinder install

### cbsuxm01 /workdir/mcs368/genespace/

library(GENESPACE)

genomeRepo <- "/workdir/mcs368/genespace/genomeRepo"
wd <- "/workdir/mcs368/genespace/"
path2mcscanx <- "/programs/MCScanX" ## directory that contains MCScanX_h


asize=read.table('../panand_sp_ploidy.txt')
#h=read.table('../general_summaries/panand_assembly_sizes.txt', header=T, sep='\t')
h=read.table('panand_assembly_sizes.txt', header=T, sep='\t')

#h$ploidy=asize$V3[match(h$V2, asize$V2)]*ifelse(h$haploid, 2,1) ## need to adjust for haploid assemblies!
asize$ploidy=asize$V3*ifelse(asize$V2%in%h$V2[h$haploid], 1,2)
asize$ploidy[asize$V2%in%c('tdacn1','tdacs1')]= 2

#### does NOT like helixer :(
#   parsedPaths <- parse_annotations(
#   rawGenomeRepo = genomeRepo, 
#   genomeDirs = c(asize$V2, 'pvagin'),
#   genomeIDs = c(asize$V2, 'pvagin'),
#   gffString = "gff",
#   faString = "fasta",
#   presets = c(rep("none", length(asize$V2)), 'phytozome'), 
#   genespaceWd = wd)


## instead make bed files on my own??


  
gpar <- init_genespace(
  wd = wd,
  path2mcscanx = path2mcscanx,
  genomeIDs=c(asize$V2, 'pvagin'),
    ploidy= c(asize$ploidy, 1),
    outgroup='pvagin',
    nCores=64)
  
  
  out <- run_genespace(gpar, overwrite = F)
  
  
  ## wait, maybe I don't want outgroup :(
  ## it keeps pvagin from being in the figures
  ## I can keep the already completed orthofinder results
  ## pvagin first makes plotting easier
    mkdir ../genespace_noog
  wd <- "/workdir/mcs368/genespace_noog" # no outgroup
  
  gpar <- init_genespace(
  wd = wd,
  path2mcscanx = path2mcscanx,
  genomeIDs=c('pvagin', asize$V2 ),
    ploidy= c(1, asize$ploidy),
#    outgroup='pvagin',
    nCores=64)
    
    out <- run_genespace(gpar, overwrite = F)
  
  ## dangit that stupid outgroup - this won't work


  
#     genomeDirs = c(asize$V2, 'pvagin'),
#   genomeIDs = c(asize$V2, 'pvagin'),
#   ploidy= c(asize$ploidy, 1),
#   outgroup='pvagin')
#   
#   
#   out <- run_genespace(gpar, overwrite = T)