library(stringr)

## reformat gff into a friendlier format
#ctg_1	RepeatMasker	similarity	1	150	 2.0	+	.	Target "Motif:achine_TandemRepeat1_183bp_repeat" 34 183
#ctg_1	EDTA	Gypsy_LTR_retrotransposon	1	15219	2950	+	.ID=TE_homo_0;Name=TE_00009095_LTR;Classification=LTR/Gypsy;Sequence_ontology=SO:0002265;Identity=0.967;Method=homology

a=read.table('achine_trrm/Ac-Pasquet1232-DRAFT-PanAnd-1.0.fasta.out.gff', sep='\t')
## okay fine travis, I'll put something correct in field 3 :D http://www.sequenceontology.org/browser/current_release/term/SO:0000705
gff=data.frame(seqname=a$V1, method='TRASH', kind='tandem_repeat', start=a$V4, end=a$V5, score=a$V6, strand=a$V7, noidea=a$V8, 
               field9=paste0('ID=TandemRepeat_', rownames(a), ';Name=', gsub('Motif:', '', str_split_fixed(a$V9, ' ', 3)[,2]), ';Classification=TandemRepeat'))
