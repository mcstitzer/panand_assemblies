library(stringr)
library(rtracklayer)
## reformat gff into a friendlier format
#ctg_1	RepeatMasker	similarity	1	150	 2.0	+	.	Target "Motif:achine_TandemRepeat1_183bp_repeat" 34 183
#ctg_1	EDTA	Gypsy_LTR_retrotransposon	1	15219	2950	+	.ID=TE_homo_0;Name=TE_00009095_LTR;Classification=LTR/Gypsy;Sequence_ontology=SO:0002265;Identity=0.967;Method=homology


all=read.table('../../panand_sp_ploidy.txt')
all=all[!all$V2 %in% c('pprate', 'tdactm', 'tzopol', 'osativ', 'bdista', 'agerjg', 'svirid', 'eophiu'),]


for(x in 1:nrow(all)) {
# a=read.csv(paste0(x, '/Summary.of.repetitive.regions.', all$V1[all$V2==x], '.fasta.csv'), header=T)
if(file.exists(paste0(all$V2[x], '_trrm/', all$V1[x], '.fasta.out.gff'))){
  a=read.table(paste0(all$V2[x], '_trrm/', all$V1[x], '.fasta.out.gff'), sep='\t')
 ## okay fine travis, I'll put something correct in field 3 :D http://www.sequenceontology.org/browser/current_release/term/SO:0000705
 gff=data.frame(seqname=a$V1, method='TRASH', kind='tandem_repeat', start=a$V4, end=a$V5, score=a$V6, strand=a$V7, noidea=a$V8, 
               field9=paste0('ID=TandemRepeat_', rownames(a), ';Name=', gsub('Motif:', '', str_split_fixed(a$V9, ' ', 3)[,2]), ';Classification=TandemRepeat'))

write.table(gff, paste0(all$V1[x], '_TandemRepeat.gff3'), row.names=F, col.names=F, sep='\t', quote=F)
  }
 print(x)
 }

### ughhhhh merge overlapping intervals taht come from the same repeat
### nevermind this is too hard, with all the score stuff - leaving here because it coul dbe fixed in the future
for(x in 1:nrow(all)){
  a=import.gff(paste0(all$V1[x], '_TandemRepeat.gff3'))
  b=unlist(reduce(split(a, ~Name))) ## merge bookended or overlapping TRs if they're the same repeat consensus
  
  }
