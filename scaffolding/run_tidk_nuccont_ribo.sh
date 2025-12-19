



six=crefra
FASTA=Cr-AUB069-DRAFT-PanAnd-1.0.chrSuperScaf.fa        
     
six=cserru
FASTA=Cs-KelloggPI219580-DRAFT-PanAnd-1.0.chrSuperScaf.fa  

six=hconto
FASTA=Hc-AUB53_1-DRAFT-PanAnd-1.0.chrSuperScaf.fa   

six=ppanic
FASTA=Pi-Clark-DRAFT-PanAnd-1.0.chrSuperScaf.fa    

## ribosomes
## where are the repeats on the genomes??
## ugh telomeres mess it up so make temp
conda activate panand_assemblies
GENOME=../updated_genomes/${FASTA}
sed 's/^[ACGT][ACGT][ACGT][ACGT]/GATC/' $GENOME > $GENOME.tmp
barrnap --quiet --kingdom euk $GENOME.tmp > ${FASTA%.fa}.rrna.gff3
rm $GENOME.tmp

### gc content
samtools faidx ../updated_genomes/${FASTA}
bedtools makewindows -g ../updated_genomes/${FASTA}.fai -w 1000000 > ${FASTA%.fa}.1mbwindows.bed
bedtools nuc -fi ../updated_genomes/${FASTA} -bed ${FASTA%.fa}.1mbwindows.bed > ${FASTA%.fa}.1mbnuccontent.bed

## run tidk
conda activate tidk
tidk find --clade Poales --output ${FASTA%.fa}.tidk --dir . ../updated_genomes/${FASTA}



## do later
## added 5 as a ploidy option!!
#cd ../../../
## in hic directory :(
six=hconto
Rscript generate_subphaser_input_cmdline.R hcontoCHR-Pv-4 4 ../scaffolding/Hc-AUB53_1-DRAFT-PanAnd-1.0.chrSuperScaf.fa.fai

six=tdacs1
Rscript generate_subphaser_input_cmdline.R tdacs1-Pv-2 2 ../genomes/Td-FL_9056069_6-REFERENCE-PanAnd-2.0a.fasta.fai

six=tdacn1
Rscript generate_subphaser_input_cmdline.R tdacn1-Pv-2 2 ../genomes/Td-KS_B6_1-REFERENCE-PanAnd-2.0a.fasta.fai



### also do on unchanged tripsacinae and avirgi

cd /workdir/mcs368/panand_assemblies/original_genomes

six=avirgi
FASTA=Av-Kellogg1287_8-REFERENCE-PanAnd-1.0.fasta

six=tdacn1
FASTA=Td-KS_B6_1-REFERENCE-PanAnd-2.0a.fasta

six=tdacn2
FASTA=Td-KS_B6_1-REFERENCE-PanAnd-2.0b.fasta

six=tdacs1
FASTA=Td-FL_9056069_6-REFERENCE-PanAnd-2.0a.fasta
six=tdacs2
FASTA=Td-FL_9056069_6-REFERENCE-PanAnd-2.0b.fasta

six=zdgigi
FASTA=Zd-Gigi-REFERENCE-PanAnd-1.0.fasta

six=zdmomo
FASTA=Zd-Momo-REFERENCE-PanAnd-1.0.fasta

six=zluxur
FASTA=Zl-RIL003-REFERENCE-PanAnd-1.0.fasta

six=zmhuet
FASTA=Zh-RIMHU001-REFERENCE-PanAnd-1.0.fasta

six=zTIL18
FASTA=Zx-TIL18-REFERENCE-PanAnd-1.0.fasta

six=zTIL25
FASTA=Zx-TIL25-REFERENCE-PanAnd-1.0.fasta

six=zTIL01
FASTA=Zv-TIL01-REFERENCE-PanAnd-1.0.fasta

six=zTIL11
FASTA=Zv-TIL11-REFERENCE-PanAnd-1.0.fasta

six=znicar
FASTA=Zn-PI615697-REFERENCE-PanAnd-1.0.fasta


six=sbicol
FASTA=Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fasta

six=zmB735
FASTA=Zm-B73-REFERENCE-NAM-5.0.fasta

conda activate panand_assemblies
GENOME=../genomes/${FASTA}
sed 's/^[ACGT][ACGT][ACGT][ACGT]/GATC/' $GENOME > $GENOME.tmp
barrnap --quiet --kingdom euk $GENOME.tmp > ${FASTA%.fasta}.rrna.gff3
rm $GENOME.tmp

### gc content
samtools faidx ../genomes/${FASTA}
bedtools makewindows -g ../genomes/${FASTA}.fai -w 1000000 > ${FASTA%.fasta}.1mbwindows.bed
bedtools nuc -fi ../genomes/${FASTA} -bed ${FASTA%.fasta}.1mbwindows.bed > ${FASTA%.fasta}.1mbnuccontent.bed

## run tidk
conda activate tidk
tidk find --clade Poales --output ${FASTA%.fasta}.tidk --dir . ../genomes/${FASTA}

