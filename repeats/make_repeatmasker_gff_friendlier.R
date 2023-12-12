## reformat gff into a friendlier format
#ctg_1	RepeatMasker	similarity	1	150	 2.0	+	.	Target "Motif:achine_TandemRepeat1_183bp_repeat" 34 183
#ctg_1	EDTA	Gypsy_LTR_retrotransposon	1	15219	2950	+	.ID=TE_homo_0;Name=TE_00009095_LTR;Classification=LTR/Gypsy;Sequence_ontology=SO:0002265;Identity=0.967;Method=homology

a=read.table('')
gff=data.frame(seqname=a$V1, method='TRASH', 'tandem_repeat'
