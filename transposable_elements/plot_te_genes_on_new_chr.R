library(rtracklayer)
library(dplyr)
library(data.table)
library(stringr)
library(plyranges)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(ggtext)

#make new fasta
#source activate scaffolding
#agptools assemble ../genomes/Cr-AUB069-DRAFT-PanAnd-1.0.fasta crefra.agp > Cr-AUB069-DRAFT-PanAnd-1.0.chrSuperScaf.fa
## update position 
#python uplift_gff.py crefra.agp ../repeats/trash/Cr-AUB069-DRAFT-PanAnd-1.0_EDTAandTandemRepeat.gff3 Cr-AUB069-DRAFT-PanAnd-1.0.chrSuperScaf.EDTAandTandemRepeat.gff3 
#python uplift_gff.py crefra.agp ../../helixer_annotations/Cr-AUB069-DRAFT-PanAnd-1.0_helixer.gff Cr-AUB069-DRAFT-PanAnd-1.0.chrSuperScaf.helixer.gff 

## actually, incorporated into script, so need to use like this
## Rscript plot_te_genes_on_new_chr.R "$PREFIX" "$OUTPUT_TRASH_GFF" "$OUTPUT_HELIXER_GFF" "$AGP_FILE"
args <- commandArgs(trailingOnly = TRUE)
genotype <- args[1]       # The PREFIX value
tegff <- args[2]          # The OUTPUT_TRASH_GFF file
genegff <- args[3]        # The OUTPUT_HELIXER_GFF file
agpfile <- args[4]        # The agp file

# Print the Inputs for Debugging
cat("Genotype Prefix: ", genotype, "\n")
cat("TE GFF File: ", tegff, "\n")
cat("Gene GFF File: ", genegff, "\n")
cat("AGP File: ", agpfile, "\n")


## xm01 /workdir/mcs368/panand_assemblies/repeats/
all=read.table('../../panand_sp_ploidy.txt')
all=all[!all$V2 %in% c('pprate', 'tdactm', 'tzopol', 'osativ', 'bdista', 'agerjg', 'svirid', 'eophiu'),]

all=all[all$V2!='sbicol',] ## dummy me lost it

genomecountlist=vector(mode = "list", length = length(all$V2))
names(genomecountlist)=all$V2

repeatlengthslist=vector(mode = "list", length = length(all$V2))
names(repeatlengthslist)=all$V2

#genotype='crefraCHR'

## repeats
a=import.gff3(tegff) 
temp=data.frame(a) #%>% group_by(fam=gsub('_LTR', '', gsub('_INT', '', Name)), Classification) %>% dplyr::filter(!is.na(as.numeric(Identity))) 
if('Note' %in% colnames(temp)){temp$Note=''}
temp$genome=genotype

genomecountlist[[genotype]]=temp
repeatlengthslist[[genotype]]=sum(width(reduce(a, ignore.strand=T)))

genomecount=do.call(rbind, genomecountlist)
genomecount$Identity=as.numeric(genomecount$Identity)

#write.table(data.frame(genome=names(repeatlengthslist), repeatbp=unlist(repeatlengthslist)),'total_repeat_bp.txt', row.names=F, col.names=T, sep='\t', quote=F)

## add a superfamily based on matchign up to the classification field - note most of these NAs are relics from the B73 annotation being included!!!!!
genomecount$sup=c(NA, 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', 'DHH', NA,NA,NA,NA,NA,NA,'RLC', 'RLG', 'RLG', 'RLX', 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', NA,NA,NA)[match(genomecount$Classification, c("Cent/CentC", "DNA/DTA", "DNA/DTC", "DNA/DTH", "DNA/DTM", "DNA/DTT", 
"DNA/Helitron", "knob/knob180", "knob/TR-1", "LINE/L1", "LINE/RTE", 
"LINE/unknown", "Low_complexity", "LTR/Copia", "LTR/CRM", "LTR/Gypsy", 
"LTR/unknown", "MITE/DTA", "MITE/DTC", "MITE/DTH", "MITE/DTM", 
"MITE/DTT", "rDNA/spacer", "Simple_repeat", "subtelomere/4-12-1"))]

## tandem repeat maybe want to do by length class???
#genomecount$sup[genomecount$source=='TRASH']=genomecount$Name[genomecount$source=='TRASH']
genomecount$sup[genomecount$source=='TRASH']='TandemRepeat'

genomecount$sup[genomecount$source=='RepeatMasker']='TandemRepeat'
genomecount$Classification[genomecount$source=='RepeatMasker']='TandemRepeat'

## there's still knob in this one znicar family ...
genomecount=genomecount %>% group_by(genome, fam=gsub('_LTR', '', gsub('_INT', '', Name))) %>% filter(!(fam == "TE_00015576" & genome == "znicar")) %>% mutate(fam=fam)


## genes
a=import.gff(genegff, version='3')
#### read in genes!
genecountlist=vector(mode = "list", length = length(all$V2))
names(genecountlist)=all$V2

a$genome=genotype
a$Name=NULL ## remove Name from liftoff genes
genecountlist[[genotype]]=data.frame(a)

genes=do.call(rbind, genecountlist)
genes=genes[genes$type=='gene',]



### okay this will work!!! make a stacked bar plot of each of these genome-wide, using these colors
# Retrotransposon colors (pastel reds and pinks)
retro_colors <- c(RLG="#F7B7B7", RLC="#F28E8E", RLX="#E67373")
# DNA transposon colors (cohesive blue gradient)
dna_colors <- c(DHH="#377EB8", DTA="#4B8EBA", DTC="#5DA9BD", DTH="#70C2BF", 
                DTM="#83D6C1", DTT="#96E2C3", DTS="#A8EEC5")
# Tandem repeat color (pastel lavender)
tandem_color <- c(TandemRepeat="#D4B4F4", TR='#D4B4F4')
# Combine all colors
te_colors <- c(retro_colors, dna_colors, tandem_color)


# Short species label (replace as needed)
shortSpeciesLabel <- genotype

## get contig boundaris from agp file
crefra_agp <- read.table(agpfile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Ensure AGP file has appropriate column names (modify these if needed based on your AGP structure)
# Example column names: seqname, start, end, contig_name, ...
colnames(crefra_agp) <- c("seqname", "start", "end", "contig_name", "type", 'scaf', 'type2', 'length', 'strand')  # Adjust as necessary

intermediate_boundaries <- crefra_agp %>%
  group_by(seqname) %>%
  filter(start != min(start) & end != max(end)) %>%
  select(seqname, start, end)

# Combine the boundaries into a single column for plotting
intermediate_boundaries <- reshape2::melt(
  data = intermediate_boundaries,
  id.vars = "seqname",
  measure.vars = c("start", "end"),
  variable.name = "boundary_type",
  value.name = "boundary_position"
)





# Identify the primary chromosomes (chr1 to chr10) and additional scaffolds > 10Mb
chromosomes <- unique(genomecount$seqnames)
primary_chromosomes <- paste0("chr", 1:10)

# Identify additional scaffolds > 10Mb
additional_scaffolds <- unique(genomecount$seqnames[!genomecount$seqnames %in% primary_chromosomes & genomecount$width > 10e6])

# Combine primary chromosomes and additional scaffolds
all_seqs_to_plot <- c(primary_chromosomes, additional_scaffolds)

pdf(paste0('~/transfer/chromosomes_', genotype, '.pdf'), 20,10)
# Function to filter longest & genes data and create plots per chromosome
for (chromosome in all_seqs_to_plot) {
  # Filter for the current chromosome and sequences longer than 10Mb
  longest_filtered <- genomecount[genomecount$seqnames == chromosome, ]
  genes_filtered <- genes[genes$seqnames == chromosome , ]
  
  # Skip if no relevant data for the chromosome
  if (nrow(longest_filtered) == 0 | nrow(genes_filtered) == 0) {
    next
  }
  
    boundaries_filtered <- intermediate_boundaries[intermediate_boundaries$seqname == chromosome, ]


  # Plot using ggplot and plot_grid
  combined_plot <- plot_grid(
    ggplot(longest_filtered, aes(x=start, fill=factor(sup), weight=width)) +
      scale_fill_manual(values=te_colors) +
      geom_histogram(binwidth=1e6, position='stack') +
      ggtitle(paste0(shortSpeciesLabel, ' Longest Sequence: ', chromosome)) +
      xlab('Position on Sequence (1 Mb bins)') + ylab('Count') +
      theme(
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(1, "cm"),
        axis.text.y = element_text(size = 9),
        legend.justification = 'center',
        plot.title = element_text(hjust = 0.5)
      ) +
    # Add vertical lines for contig boundaries
    geom_vline(data = boundaries_filtered, aes(xintercept = boundary_position), 
               color = "black", linetype = "dashed", alpha = 0.7)+ 
      guides(fill = guide_legend(
        title.position = "top",
        label.position = "bottom",
        nrow = 1,
        byrow = TRUE,
        reverse = TRUE
      )),
    
    ggplot(genes_filtered, aes(x=start)) +
      geom_histogram(binwidth=1e6) +
      xlab('Position on Sequence (1 Mb bins)') +
      ylab('Count') +
      ggtitle('Helixer Genes') +
    # Add vertical lines for contig boundaries
    geom_vline(data = boundaries_filtered, aes(xintercept = boundary_position), 
               color = "black", linetype = "dashed", alpha = 0.7)+
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 9)
      ),
    
    align = 'hv',
    ncol = 1,
    rel_heights = c(1, 0.3)
  )
  
  # Print plot
  print(combined_plot)
}


dev.off()