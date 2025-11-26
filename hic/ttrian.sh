## this should be relatively easy - 0.175% het, no Bionano
## 854 mb flow
conda activate panand_assemblies

hicR1=AV9240007_Themeda_triandra_HiC_pairsR1.fq.gz
hicR2=AV9240007_Themeda_triandra_HiC_pairsR2.fq.gz
six=ttrian

scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/082525_PanAnd/AV9240007_Themeda_triandra_HiC_pairsR1.fq.gz .
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/082525_PanAnd/AV9240007_Themeda_triandra_HiC_pairsR2.fq.gz .
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/082625_PanAnd/AV9590002_Themeda_triandra_HiFi_20250822_01_R84064-1_B01.filter.bam .
samtools bam2fq -@ 64 AV9590002_Themeda_triandra_HiFi_20250822_01_R84064-1_B01.filter.bam > ${six}_hifi.fastq
bgzip ${six}_hifi.fastq



## check kmers first!
mkdir tmp
/programs/kmc-3.2.4/kmc -k21 -t10 -m64 -ci1 -cs10000 ${six}_hifi.fastq.gz ${six}kmc tmp/
/programs/kmc-3.2.4/kmc_tools transform ${six}kmc histogram ${six}.histo -cx10000

../genomescope2.0/genomescope.R -i ${six}.histo -o ${six}_genomescope -k 21

export PATH=/programs/fastk-1.1/bin:$PATH
export PATH=/programs/smudgeplot-0.4.0/bin:$PATH
# run FastK to create a k-mer database
FastK -v -t4 -k31 -M16 -T4 ${six}_hifi.fastq.gz -N${six}_FastK_Table

# Find all k-mer pairs in the dataset using hetmer module
smudgeplot.py hetmers -L 12 -t 4 -o ${six}_kmerpairs --verbose ${six}_FastK_Table
# this now generated `data/Scer/kmerpairs_text.smu` file;
# it's a flat file with three columns; covB, covA and freq (the number of k-mer pairs with these respective coverages)

# use the .smu file to infer ploidy and create smudgeplot
smudgeplot.py all -o ${six}_smudgeplot ${six}_kmerpairs_text.smu


## dual-scaf, hic assembly
# hifiasm -o ${six}.asm -t64 --dual-scaf --h1 AV9240007_Themeda_triandra_HiC_pairsR1.fq.gz --h2 AV9240007_Themeda_triandra_HiC_pairsR2.fq.gz ${six}_hifi.fastq.gz
### this looks really weird - big assembly, also looks like a tetraploid? two paths in each haplotype visible in hic....
### hmm, maybe it's not getting homozygous coverage peak right?
# hifiasm -o ${six}.asm -t64 --hom-cov 22 --dual-scaf --h1 AV9240007_Themeda_triandra_HiC_pairsR1.fq.gz --h2 AV9240007_Themeda_triandra_HiC_pairsR2.fq.gz ${six}_hifi.fastq.gz
## back to normal because it was duplicating everyting???
#hifiasm -o ${six}.asm -t64  --dual-scaf --h1 AV9240007_Themeda_triandra_HiC_pairsR1.fq.gz --h2 AV9240007_Themeda_triandra_HiC_pairsR2.fq.gz ${six}_hifi.fastq.gz
## still big! try --hg-size to set haploid size?
#hifiasm -o ${six}.asm -t64  --dual-scaf --hg-size 850m --h1 AV9240007_Themeda_triandra_HiC_pairsR1.fq.gz --h2 AV9240007_Themeda_triandra_HiC_pairsR2.fq.gz ${six}_hifi.fastq.gz
## NOPE still big
##ummm chatgpt told me to use 55, so i'll try it...
# hifiasm -o ${six}.asm -t64 --hom-cov 48 --dual-scaf --h1 AV9240007_Themeda_triandra_HiC_pairsR1.fq.gz --h2 AV9240007_Themeda_triandra_HiC_pairsR2.fq.gz ${six}_hifi.fastq.gz

### this is beautiful!!!
 hifiasm -o ${six}.asm -t64 -s 0.2 --dual-scaf --h1 AV9240007_Themeda_triandra_HiC_pairsR1.fq.gz --h2 AV9240007_Themeda_triandra_HiC_pairsR2.fq.gz ${six}_hifi.fastq.gz
## 1.6 Gb a piece!


awk '/^S/{print ">"$2"\n"$3}' ${six}.asm.hic.hap1.p_ctg.gfa > ${six}.asm.hic.hap1.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' ${six}.asm.hic.hap2.p_ctg.gfa > ${six}.asm.hic.hap2.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' ${six}.asm.hic.p_ctg.gfa > ${six}.asm.hic.p_ctg.fa


## check assembly, how many contigs per chromosome??
/programs/bbmap-38.96/stats.sh -Xmx200g in="${six}.asm.hic.hap1.p_ctg.fa" out="${six}.asm.hic.hap1.p_ctg.stats.txt" overwrite='true'
/programs/bbmap-38.96/stats.sh -Xmx200g in="${six}.asm.hic.hap2.p_ctg.fa" out="${six}.asm.hic.hap2.p_ctg.stats.txt" overwrite='true'
/programs/bbmap-38.96/stats.sh -Xmx200g in="${six}.asm.hic.p_ctg.fa" out="${six}.asm.hic.p_ctg.stats.txt" overwrite='true'



### I DON"T GET IT... this is still huge, and what is happenoing??
# i'm going to build a tree using these assemblies and see if they're what we think they are?
## or like ahve they polyploidized??

## check out a dotplot
cd ..
./fa_to_dotplot.sh ${six}/${six}.asm.hic.hap1.p_ctg.fa ${six}Hap1 1
./fa_to_dotplot.sh ${six}/${six}.asm.hic.hap2.p_ctg.fa ${six}Hap2 1
cd ${six}


## combine haplotypes and haphic



cat ${six}.asm.hic.hap1.p_ctg.fa ${six}.asm.hic.hap2.p_ctg.fa >${six}.asm.hic.bothhaps.p_ctg.fa

bwa index ${six}.asm.hic.bothhaps.p_ctg.fa
bwa mem -5SP -t 128 ${six}.asm.hic.bothhaps.p_ctg.fa $hicR1 $hicR2 | samblaster | samtools view - -@ 64 -S -h -b -F 3340 -o ${six}BothHaps.bam
## duh it can't index if not sorted
##samtools index ${six}BothHaps.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
../HapHiC/utils/filter_bam ${six}BothHaps.bam 1 --nm 3 --threads 64 | samtools view - -b -@ 64 -o ${six}BothHaps.filtered.bam

## filtering reduces 27 Gb bam to 165 mb =--- seems dangerous! I think it's because this is so low heterozygosity.

## and can use the gfa to deal with depths and mis-pairing
../HapHiC/haphic pipeline ${six}.asm.hic.bothhaps.p_ctg.fa ${six}BothHaps.filtered.bam 40 --gfa "${six}.asm.hic.hap1.p_ctg.gfa,${six}.asm.hic.hap2.p_ctg.gfa" --outdir ${six}BothHaps

## skip filtering? potential rec here https://github.com/zengxiaofei/HapHiC/issues/86
## use gfa to be sure not to combine haps
../HapHiC/haphic pipeline ${six}.asm.hic.bothhaps.p_ctg.fa ${six}BothHaps.bam 40 --gfa "${six}.asm.hic.hap1.p_ctg.gfa,${six}.asm.hic.hap2.p_ctg.gfa" --outdir ${six}BothHapsUnfiltered



## to correct contigs with hic info
## --correct_nrounds 2 
## to remove hic links between alleles, so set to number of haplotypes present (polyploid)
## --remove_allelic_links 2

## try quickview 
#../HapHiC/haphic pipeline ${six}.asm.hic.bothhaps.p_ctg.fa ${six}BothHaps.filtered.bam 20 --quick_view --gfa "${six}.asm.hic.hap1.p_ctg.gfa,${six}.asm.hic.hap2.p_ctg.gfa" --outdir ${six}BothHapsQuickView
../HapHiC/haphic pipeline ${six}.asm.hic.bothhaps.p_ctg.fa ${six}BothHaps.bam 20 --quick_view  --gfa "${six}.asm.hic.hap1.p_ctg.gfa,${six}.asm.hic.hap2.p_ctg.gfa" --outdir ${six}BothHapsQuickViewUnfiltered


#########################
## scaffold each haplotype separately

conda activate haphic

bwa index ${six}.asm.hic.hap1.p_ctg.fa
bwa mem -5SP -t 128 ${six}.asm.hic.hap1.p_ctg.fa $hicR1 $hicR2 | samblaster | samtools view - -@ 64 -S -h -b -F 3340 -o ${six}Hap1.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
../HapHiC/utils/filter_bam ${six}Hap1.bam 1 --nm 3 --threads 64 | samtools view - -b -@ 64 -o ${six}Hap1.filtered.bam

## filtering reduces 27 Gb bam to 165 mb =--- seems dangerous! I think it's because this is so low heterozygosity.

## and can use the gfa to deal with depths and mis-pairing
../HapHiC/haphic pipeline ${six}.asm.hic.hap1.p_ctg.fa ${six}Hap1.filtered.bam 20 --outdir ${six}Hap1
cd ${six}/${six}Hap1/04.build/
bash juicebox.sh


## hap2
bwa index ${six}.asm.hic.hap2.p_ctg.fa
bwa mem -5SP -t 128 ${six}.asm.hic.hap2.p_ctg.fa $hicR1 $hicR2 | samblaster | samtools view - -@ 64 -S -h -b -F 3340 -o ${six}Hap2.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
../HapHiC/utils/filter_bam ${six}Hap2.bam 1 --nm 3 --threads 64 | samtools view - -b -@ 64 -o ${six}Hap2.filtered.bam

## filtering reduces 27 Gb bam to 165 mb =--- seems dangerous! I think it's because this is so low heterozygosity.

## and can use the gfa to deal with depths and mis-pairing
../HapHiC/haphic pipeline ${six}.asm.hic.hap2.p_ctg.fa ${six}Hap2.filtered.bam 20 --outdir ${six}Hap2
cd ${six}/${six}Hap2/04.build/
bash juicebox.sh


### then use corrected contact maps
# Generate the final FASTA file for the scaffolds
## from the 04 directory...
##python ../../../juicebox_scripts/juicebox_scripts/juicebox_assembly_converter.py -a out_JBAT.review.${six}Hap2.assembly -f scaffolds.fa -s


##### download out_JBAT.hic and out_JBAT.assembly, and correct in juicebox




###
juicer post -o ${six}Hap1_aggressivecorrection ${six}Hap1_aggressivecorrection.out_JBAT.review.assembly.assembly out_JBAT.liftover.agp ../../${six}.asm.hic.hap1.p_ctg.fa
cd ../../
./fa_to_dotplot.sh ${six}/${six}Hap1/04.build/${six}Hap1_aggressivecorrection.FINAL.fa ${six}Hap1aggressive 1

## ribosomes
## where are the repeats on the genomes??
## ugh telomeres mess it up so make temp
GENOME=${six}/${six}Hap1/04.build/${six}Hap1_aggressivecorrection.FINAL.fa 
sed 's/^[ACGT][ACGT][ACGT][ACGT]/GATC/' $GENOME > $GENOME.tmp
barrnap --quiet --kingdom euk $GENOME.tmp > ${GENOME%.fa}.rrna.gff3
rm $GENOME.tmp

## run trash2
conda activate TRASH_2
Rscript /workdir/mcs368/panand_assemblies/repeats/TRASH_2/src/TRASH.R -f ${six}Hap1_aggressivecorrection.FINAL.fa -o /workdir/mcs368/panand_assemblies/hic/${six}/${six}Hap1/04.build/${six}Hap1_aggressivecorrection_TRASH2 -p 24

## run tidk
conda activate tidk
tidk find --clade Poales --output ${six}Hap1_aggressivecorrection.tidk --dir . ${six}Hap1_aggressivecorrection.FINAL.fa 

## run helixer
ssh cbsugpu08
cd /workdir/mcs368/
scp cbsuxm01:/workdir/mcs368/panand_assemblies/hic/${six}/${six}Hap1/04.build/${six}Hap1_aggressivecorrection.FINAL.fa .
singularity run --nv --bind $PWD --pwd $PWD /programs/helixer-0.3.5/helixer.sif Helixer.py --fasta-path ${six}Hap1_aggressivecorrection.FINAL.fa --lineage land_plant --gff-output-path ${six}Hap1_aggressivecorrection.FINAL.helixer.gff3
scp ${six}Hap1_aggressivecorrection.FINAL.helixer.gff3 cbsuxm01:/workdir/mcs368/panand_assemblies/hic/${six}/${six}Hap1/04.build/${six}Hap1_aggressivecorrection.FINAL.helixer.gff3


### gc content
conda activate panand_assemblies
samtools faidx ${six}Hap1_aggressivecorrection.FINAL.fa
bedtools makewindows -g ${six}Hap1_aggressivecorrection.FINAL.fa.fai -w 1000000 > ${six}Hap1_aggressivecorrection.FINAL.1mbwindows.bed
bedtools nuc -fi ${six}Hap1_aggressivecorrection.FINAL.fa -bed ${six}Hap1_aggressivecorrection.FINAL.1mbwindows.bed > ${six}Hap1_aggressivecorrection.FINAL.1mbnuccontent.bed





############# NOW DO IT ALL WITH HAP2!!!!
###
juicer post -o ${six}Hap2_aggressivecorrection ${six}Hap2_aggressivecorrection.out_JBAT.review.assembly out_JBAT.liftover.agp ../../${six}.asm.hic.hap2.p_ctg.fa
cd ../../../
./fa_to_dotplot.sh ${six}/${six}Hap2/04.build/${six}Hap2_aggressivecorrection.FINAL.fa ${six}Hap2aggressive 1



## ribosomes
## where are the repeats on the genomes??
## ugh telomeres mess it up so make temp
GENOME=${six}/${six}Hap2/04.build/${six}Hap2_aggressivecorrection.FINAL.fa 
sed 's/^[ACGT][ACGT][ACGT][ACGT]/GATC/' $GENOME > $GENOME.tmp
barrnap --quiet --kingdom euk $GENOME.tmp > ${GENOME%.fa}.rrna.gff3
rm $GENOME.tmp

## run trash2
cd ${six}/${six}Hap2/04.build
mkdir -p ${six}Hap2_aggressivecorrection_TRASH2
Rscript /workdir/mcs368/panand_assemblies/repeats/TRASH_2/src/TRASH.R -f ${six}Hap2_aggressivecorrection.FINAL.fa -o /workdir/mcs368/panand_assemblies/hic/${six}/${six}Hap2/04.build/${six}Hap2_aggressivecorrection_TRASH2 -p 24

## run tidk
conda activate tidk
tidk find --clade Poales --output ${six}Hap2_aggressivecorrection.tidk --dir . ${six}Hap2_aggressivecorrection.FINAL.fa 

## run helixer
ssh cbsugpu08
cd /workdir/mcs368/
scp cbsuxm01:/workdir/mcs368/panand_assemblies/hic/${six}/${six}Hap2/04.build/${six}Hap2_aggressivecorrection.FINAL.fa .
singularity run --nv --bind $PWD --pwd $PWD /programs/helixer-0.3.5/helixer.sif Helixer.py --fasta-path ${six}Hap2_aggressivecorrection.FINAL.fa --lineage land_plant --gff-output-path ${six}Hap2_aggressivecorrection.FINAL.helixer.gff3
scp ${six}Hap2_aggressivecorrection.FINAL.helixer.gff3 cbsuxm01:/workdir/mcs368/panand_assemblies/hic/${six}/${six}Hap2/04.build/${six}Hap2_aggressivecorrection.FINAL.helixer.gff3


### gc content
samtools faidx ${six}Hap2_aggressivecorrection.FINAL.fa
bedtools makewindows -g ${six}Hap2_aggressivecorrection.FINAL.fa.fai -w 1000000 > ${six}Hap2_aggressivecorrection.FINAL.1mbwindows.bed
bedtools nuc -fi ${six}Hap2_aggressivecorrection.FINAL.fa -bed ${six}Hap2_aggressivecorrection.FINAL.1mbwindows.bed > ${six}Hap2_aggressivecorrection.FINAL.1mbnuccontent.bed





#####also map hic to existing referenceand try toscaffold??!?!


conda activate haphic

bwa index Tt-AUB21_1-DRAFT-PanAnd-1.0.fasta
bwa mem -5SP -t 96 Tt-AUB21_1-DRAFT-PanAnd-1.0.fasta $hicR1 $hicR2 | samblaster | samtools view - -@ 64 -S -h -b -F 3340 -o ${six}DIPLOID.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
../HapHiC/utils/filter_bam ${six}DIPLOID.bam 1 --nm 3 --threads 64 | samtools view - -b -@ 64 -o ${six}DIPLOID.filtered.bam

## filtering reduces 27 Gb bam to 165 mb =--- seems dangerous! I think it's because this is so low heterozygosity.

## well, it's allelic so i guess we say 40 chromsomes!!!
../HapHiC/haphic pipeline Tt-AUB21_1-DRAFT-PanAnd-1.0.fasta ${six}DIPLOID.filtered.bam 60 --outdir ${six}DIPLOID

### hmmm not the greatest - adding a lot of small (100kb) contigs to ends of each chromosome
### just be aggressive :)
juicer post -o ${six}DIPLOID_aggressivecorrection ${six}DIPLOID_aggressivecorrection.out_JBAT.review.assembly out_JBAT.liftover.agp ../../Tt-AUB21_1-DRAFT-PanAnd-1.0.fasta
cd ../../../
./fa_to_dotplot.sh ${six}/${six}/04.build/${six}DIPLOID_aggressivecorrection.FINAL.fa ${six}DIPLOIDaggressive 3











